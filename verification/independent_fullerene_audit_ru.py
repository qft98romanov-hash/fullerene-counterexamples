#!/usr/bin/env python3
"""
Независимый аудит сертификатов контрпримеров о фуллеренах.

В отличие от основных проверяющих программ, этот файл не импортирует их код:
- NetworkX используется для graph6, планарности и изоморфизма;
- SymPy — для точной целочисленной матричной алгебры;
- mpmath — только для недоказательной численной сверки;
- строгие знаки в C_34 получают точными Fraction-суммами и рациональной
  оценкой хвоста.

Запуск:
    python independent_fullerene_audit_ru.py /path/to/extracted/package

Зависимости:
    pip install networkx sympy mpmath
"""
from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import Counter
from fractions import Fraction
from math import factorial
from pathlib import Path

import mpmath as mp
import networkx as nx
import sympy as sp



# Читаем непустые graph6-строки из файла как ASCII-байты.
def read_graph6(path: Path) -> list[bytes]:
    return [
        line.strip().encode("ascii")
        for line in path.read_text("ascii").splitlines()
        if line.strip()
    ]



# Вычисляем SHA256 файла блоками, не загружая весь файл в память.
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()



# Проверяем планарность и перечисляем грани заданного вложения.
def planar_faces(graph: nx.Graph) -> tuple[nx.PlanarEmbedding, list[tuple[int, ...]]]:
    planar, embedding = nx.check_planarity(graph)
    if not planar:
        raise AssertionError("graph is not planar")
    marked: set[tuple[int, int]] = set()
    faces: list[tuple[int, ...]] = []
    for u, v in embedding.edges():
        if (u, v) not in marked:
            faces.append(tuple(embedding.traverse_face(u, v, marked)))
    return embedding, faces



# Строим геометрический двойственный граф по списку граней.
def geometric_dual(graph: nx.Graph) -> tuple[nx.Graph, list[tuple[int, ...]]]:
    _, faces = planar_faces(graph)
    dart_to_face: dict[tuple[int, int], int] = {}
    for index, face in enumerate(faces):
        for j, u in enumerate(face):
            v = face[(j + 1) % len(face)]
            if (u, v) in dart_to_face:
                raise AssertionError("directed edge belongs to two traversed faces")
            dart_to_face[(u, v)] = index

    dual = nx.Graph()
    dual.add_nodes_from(range(len(faces)))
    for u, v in graph.edges():
        left = dart_to_face[(u, v)]
        right = dart_to_face[(v, u)]
        if left == right:
            raise AssertionError("bridge in alleged triangulation")
        dual.add_edge(left, right)
    return dual, faces


def verify_dual_fullerene(raw: bytes, direct_order: int) -> tuple[nx.Graph, nx.Graph]:
    graph = nx.from_graph6_bytes(raw)
    if nx.number_of_selfloops(graph):
        raise AssertionError("loop in graph6 input")
    if not nx.is_connected(graph):
        raise AssertionError("dual graph is disconnected")

    dual_order = graph.number_of_nodes()
    dual_size = graph.number_of_edges()
    degrees = Counter(dict(graph.degree()).values())
    expected_degrees = Counter({5: 12, 6: dual_order - 12})
    if degrees != expected_degrees:
        raise AssertionError((degrees, expected_degrees))
    if direct_order != 2 * (dual_order - 2):
        raise AssertionError("wrong direct/dual order relation")

    _, triangular_faces = planar_faces(graph)
    if len(triangular_faces) != direct_order or any(len(face) != 3 for face in triangular_faces):
        raise AssertionError("not a spherical triangulation")
    if dual_order - dual_size + len(triangular_faces) != 2:
        raise AssertionError("Euler characteristic is not 2")

    direct, _ = geometric_dual(graph)
    if direct.number_of_nodes() != direct_order or not nx.is_connected(direct):
        raise AssertionError("bad geometric dual")
    if Counter(dict(direct.degree()).values()) != Counter({3: direct_order}):
        raise AssertionError("direct graph is not cubic")
    _, direct_faces = planar_faces(direct)
    face_profile = Counter(map(len, direct_faces))
    expected_profile = Counter({5: 12, 6: direct_order // 2 - 10})
    if face_profile != expected_profile:
        raise AssertionError((face_profile, expected_profile))
    return graph, direct


def sympy_adjacency(graph: nx.Graph) -> sp.Matrix:
    vertices = sorted(graph)
    position = {vertex: i for i, vertex in enumerate(vertices)}
    matrix = sp.zeros(len(vertices))
    for u, v in graph.edges():
        i, j = position[u], position[v]
        matrix[i, j] = matrix[j, i] = 1
    return matrix


def degree_matrix(graph: nx.Graph) -> sp.Matrix:
    return sp.diag(*[graph.degree(v) for v in sorted(graph)])


def enumerate_gsw_paths(graph: nx.Graph) -> list[tuple[int, ...]]:
    """Enumerate oriented gSW paths using the literal published Definition 1."""
    degree = dict(graph.degree())
    neighbors = {v: set(graph[v]) for v in graph}
    paths: set[tuple[int, ...]] = set()

    for p in graph:
        if degree[p] != 5:
            continue
        for h in neighbors[p]:
            if degree[h] != 6:
                continue
            for z in neighbors[p] & neighbors[h]:
                path = [p, h, z]
                used = set(path)
                while len(path) <= graph.number_of_nodes():
                    if (
                        len(path) >= 4
                        and len(path) % 2 == 0
                        and degree[path[-2]] == 6
                        and degree[path[-1]] == 5
                    ):
                        paths.add(tuple(path))

                    continuation = (
                        neighbors[path[-2]] & neighbors[path[-1]]
                    ) - {path[-3]}
                    if len(continuation) != 1:
                        raise AssertionError(
                            "zigzag continuation is not unique in the supplied triangulation"
                        )
                    next_vertex = next(iter(continuation))
                    if next_vertex in used:
                        break
                    path.append(next_vertex)
                    used.add(next_vertex)
    return sorted(paths)


def c32_certificate(package: Path) -> dict[str, object]:
    raw = read_graph6(package / "data/C32.g6")[0]
    graph, _ = verify_dual_fullerene(raw, 32)

    degree_six = [v for v, d in graph.degree() if d == 6]
    t6 = graph.subgraph(degree_six).copy()
    components = [t6.subgraph(c).copy() for c in nx.connected_components(t6)]
    if len(components) != 2 or not all(
        nx.is_isomorphic(component, nx.complete_graph(3)) for component in components
    ):
        raise AssertionError("T^6 is not K3 disjoint union K3")

    witness = (3, 0, 4, 5, 13, 14)  # zero-based form of (4,1,5,6,14,15)
    if len(set(witness)) != 6:
        raise AssertionError("repeated witness vertex")
    witness_degrees = [graph.degree(v) for v in witness]
    if witness_degrees != [5, 6, 5, 5, 6, 5]:
        raise AssertionError(witness_degrees)
    if not all(graph.has_edge(witness[i], witness[i + 1]) for i in range(5)):
        raise AssertionError("missing ordinary path edge")
    if not all(graph.has_edge(witness[i], witness[i + 2]) for i in range(4)):
        raise AssertionError("missing distance-two chord")

    all_paths = enumerate_gsw_paths(graph)
    if witness not in all_paths:
        raise AssertionError("explicit witness missing from independent enumeration")
    internal_six = [
        path
        for path in all_paths
        if all(graph.degree(v) == 6 for v in path[2:-2])
    ]

    return {
        "dual_order": graph.number_of_nodes(),
        "dual_size": graph.number_of_edges(),
        "T6_components_one_based": [
            [v + 1 for v in sorted(component.nodes())] for component in components
        ],
        "witness_one_based": [v + 1 for v in witness],
        "witness_degrees": witness_degrees,
        "oriented_gsw_paths_found": len(all_paths),
        "gsw_paths_with_all_internal_vertices_degree_6": len(internal_six),
    }


def c44_certificate(package: Path) -> dict[str, object]:
    raws = read_graph6(package / "data/C44_pair.g6")
    graphs = [verify_dual_fullerene(raw, 44)[0] for raw in raws]
    if nx.is_isomorphic(graphs[0], graphs[1]):
        raise AssertionError("C44 pair is isomorphic")

    adjacency = [sympy_adjacency(graph) for graph in graphs]
    trace_a6 = [int(sp.trace(matrix**6)) for matrix in adjacency]
    if trace_a6 != [40386, 40362]:
        raise AssertionError(trace_a6)

    q_matrices = [
        adjacency[i] + degree_matrix(graphs[i]) for i in range(2)
    ]
    charpolys = [matrix.charpoly().all_coeffs() for matrix in q_matrices]
    if charpolys[0] != charpolys[1]:
        raise AssertionError("charpoly(A+D) differs")

    variable = sp.symbols("x")
    moments = [
        sp.expand(
            sp.trace((adjacency[i] + variable * degree_matrix(graphs[i])) ** 6)
        )
        for i in range(2)
    ]
    mixed_difference = sp.expand(moments[0] - moments[1])
    if sp.expand(mixed_difference - 24 * (1 - variable**3)) != 0:
        raise AssertionError(mixed_difference)

    lam = sp.symbols("lambda")
    polynomial = sp.Poly(
        sum(c * lam ** (24 - i) for i, c in enumerate(charpolys[0])), lam
    )
    return {
        "trace_A6": trace_a6,
        "charpoly_equal": True,
        "charpoly_factorization": str(sp.factor(polynomial.as_expr())),
        "mixed_sixth_moment_difference": str(sp.factor(mixed_difference)),
    }


def exact_taylor_interval(
    t: Fraction,
    moment_differences: list[int],
    order: int = 19,
    maximum_row_sum: int = 18,
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    cutoff = len(moment_differences) - 1
    partial = sum(
        Fraction(moment_differences[k], factorial(k)) * (t / 2) ** k
        for k in range(cutoff + 1)
    )
    q = t * maximum_row_sum / 2
    if q >= cutoff + 2:
        raise AssertionError("Taylor cutoff too short")
    radius = (
        2
        * order
        * q ** (cutoff + 1)
        / factorial(cutoff + 1)
        / (1 - q / Fraction(cutoff + 2))
    )
    return partial - radius, partial + radius, partial, radius


def decimal(fraction: Fraction, digits: int = 80) -> str:
    with mp.workdps(digits + 10):
        return mp.nstr(mp.mpf(fraction.numerator) / fraction.denominator, digits)


def c34_certificate(package: Path) -> dict[str, object]:
    raws = read_graph6(package / "data/C34_pair.g6")
    graphs = [verify_dual_fullerene(raw, 34)[0] for raw in raws]
    if nx.is_isomorphic(graphs[0], graphs[1]):
        raise AssertionError("C34 pair is isomorphic")

    t6_triangle_counts: list[int] = []
    for graph in graphs:
        t6 = graph.subgraph([v for v, d in graph.degree() if d == 6])
        t6_triangle_counts.append(sum(nx.triangles(t6).values()) // 3)
    if t6_triangle_counts != [0, 1]:
        raise AssertionError(t6_triangle_counts)

    adjacency = [sympy_adjacency(graph) for graph in graphs]
    b_matrices = [
        2 * adjacency[i] + degree_matrix(graphs[i]) for i in range(2)
    ]

    powers = [sp.eye(19), sp.eye(19)]
    differences: list[int] = []
    for k in range(83):
        differences.append(int(sp.trace(powers[0]) - sp.trace(powers[1])))
        if k < 82:
            powers = [powers[i] * b_matrices[i] for i in range(2)]
    first_nonzero = next((k, value) for k, value in enumerate(differences) if value)
    if first_nonzero != (6, 960):
        raise AssertionError(first_nonzero)

    endpoint_data: dict[str, object] = {}
    for label, t, expected_sign in (
        ("half", Fraction(1, 2), 1),
        ("one", Fraction(1), -1),
    ):
        lower, upper, partial, radius = exact_taylor_interval(t, differences)
        sign = 1 if lower > 0 else -1 if upper < 0 else 0
        if sign != expected_sign:
            raise AssertionError((label, lower, upper))
        endpoint_data[label] = {
            "partial": decimal(partial, 65),
            "radius": decimal(radius, 20),
            "certified_sign": sign,
        }

    stored = json.loads((package / "data/certificates.json").read_text("utf-8"))
    bracket = stored["C34"]["alpha_root_bracket"]
    left = Fraction(bracket["left_fraction"])
    right = Fraction(bracket["right_fraction"])
    left_interval = exact_taylor_interval(left, differences)
    right_interval = exact_taylor_interval(right, differences)
    if not (Fraction(1, 2) < left < right < 1):
        raise AssertionError("root bracket outside Conjecture 3 domain")
    if left_interval[0] <= 0 or right_interval[1] >= 0:
        raise AssertionError("root bracket is not sign-changing")

    # Non-rigorous but independent numerical sanity check from eigenvalues.
    mp.mp.dps = 100
    eigenvalues: list[list[mp.mpf]] = []
    for matrix in b_matrices:
        mp_matrix = mp.matrix(
            [[mp.mpf(int(matrix[i, j])) for j in range(19)] for i in range(19)]
        )
        values, _ = mp.eigsy(mp_matrix)
        eigenvalues.append([values[i] for i in range(19)])

    def numerical_delta(t: Fraction) -> mp.mpf:
        tt = mp.mpf(t.numerator) / t.denominator
        return mp.fsum(mp.exp(tt * value / 2) for value in eigenvalues[0]) - mp.fsum(
            mp.exp(tt * value / 2) for value in eigenvalues[1]
        )

    return {
        "T6_triangle_counts": t6_triangle_counts,
        "first_nonzero_B_moment": first_nonzero,
        "endpoints": endpoint_data,
        "alpha_bracket": [decimal(left, 76), decimal(right, 76)],
        "beta_bracket": [decimal(left / 2, 76), decimal(right / 2, 76)],
        "bracket_width": decimal(right - left, 25),
        "numerical_delta_at_bracket": [
            mp.nstr(numerical_delta(left), 35),
            mp.nstr(numerical_delta(right), 35),
        ],
    }


def random_relabeling_checks(package: Path, trials: int = 100) -> None:
    rng = random.Random(0xC0FFEE)
    c32 = nx.from_graph6_bytes(read_graph6(package / "data/C32.g6")[0])
    c44 = [
        nx.from_graph6_bytes(raw) for raw in read_graph6(package / "data/C44_pair.g6")
    ]
    c34 = [
        nx.from_graph6_bytes(raw) for raw in read_graph6(package / "data/C34_pair.g6")
    ]
    target_c44 = [
        (sympy_adjacency(graph) + degree_matrix(graph)).charpoly().all_coeffs()
        for graph in c44
    ]

    for _ in range(trials):
        vertices = list(c32)
        shuffled = vertices[:]
        rng.shuffle(shuffled)
        relabeled = nx.relabel_nodes(c32, dict(zip(vertices, shuffled)), copy=True)
        degree_six = [v for v, d in relabeled.degree() if d == 6]
        components = [
            relabeled.subgraph(c)
            for c in nx.connected_components(relabeled.subgraph(degree_six))
        ]
        if not (
            len(components) == 2
            and all(nx.is_isomorphic(c, nx.complete_graph(3)) for c in components)
            and enumerate_gsw_paths(relabeled)
        ):
            raise AssertionError("C32 relabeling failure")

        relabeled_c44 = []
        for graph in c44:
            vertices = list(graph)
            shuffled = vertices[:]
            rng.shuffle(shuffled)
            relabeled_c44.append(
                nx.relabel_nodes(graph, dict(zip(vertices, shuffled)), copy=True)
            )
        observed = [
            (sympy_adjacency(graph) + degree_matrix(graph)).charpoly().all_coeffs()
            for graph in relabeled_c44
        ]
        if observed != target_c44 or nx.is_isomorphic(*relabeled_c44):
            raise AssertionError("C44 relabeling failure")

        counts = []
        for graph in c34:
            vertices = list(graph)
            shuffled = vertices[:]
            rng.shuffle(shuffled)
            relabeled = nx.relabel_nodes(graph, dict(zip(vertices, shuffled)), copy=True)
            t6 = relabeled.subgraph([v for v, d in relabeled.degree() if d == 6])
            counts.append(sum(nx.triangles(t6).values()) // 3)
        if counts != [0, 1]:
            raise AssertionError("C34 relabeling failure")



# Точка входа: читаем пакет, запускаем все независимые проверки и сохраняем JSON.
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("package", type=Path)
    parser.add_argument("--pdf", type=Path, default=None)
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()
    package_label = args.package.as_posix()
    package = args.package.resolve()

    result: dict[str, object] = {
        "package": package_label,
        "C32": c32_certificate(package),
        "C44": c44_certificate(package),
        "C34": c34_certificate(package),
    }
    random_relabeling_checks(package, 100)
    result["random_relabeling_trials"] = 100

    if args.pdf is not None:
        result["pdf_sha256"] = sha256(args.pdf)

    print(json.dumps(result, ensure_ascii=False, indent=2))
    print("INDEPENDENT AUDIT: PASS")
    if args.json is not None:
        args.json.write_text(
            json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )


if __name__ == "__main__":
    main()
