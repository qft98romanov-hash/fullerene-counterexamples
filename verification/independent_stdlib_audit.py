#!/usr/bin/env python3
"""
Независимая автономная проверка трёх основных сертификатов.

Файл не импортирует код из ``core`` и использует только
стандартную библиотеку Python. Все равенства целочисленные,
а знаки для C34 сертифицируются рациональными интервалами.

Запуск из корня пакета:

    python independent_stdlib_audit.py

Программа работает и при ``python -O``: доказательные
условия выполняются через явные исключения.
"""

from __future__ import annotations

from collections import Counter, deque
from decimal import Decimal, getcontext
from fractions import Fraction
from itertools import combinations
from math import factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

EXPECTED_C32 = ("Q|eMID@WH?e@E@B?_wGBB?MC?NW",)
EXPECTED_C44 = (
    "W|eMID@WH_a@E@B?__GM@?OK@_G@_G?wC?C@?@wG?@P_?@|",
    "W|eMID@WH_b@A@B?__GM@?OK@_G@oG?WC?E@??wG?@W_??~",
)
EXPECTED_C34 = (
    "R|eMID`GH_b@A@B?_wGAF?[C?QW?{G",
    "R|eMID`GH_a@E@B?_wGB`?FC?FG?Bw",
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_nonempty_ascii_lines(path: Path) -> tuple[str, ...]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return tuple(line.strip() for line in path.read_text("ascii").splitlines() if line.strip())


def decode_short_graph6(text: str) -> list[set[int]]:
    """Строгий декодер канонической короткой graph6-формы, n<=62."""
    if text.startswith(">>graph6<<"):
        text = text[len(">>graph6<<") :]
    if not text:
        raise ValueError("empty graph6 string")
    values = [ord(char) - 63 for char in text]
    if any(value < 0 or value > 63 for value in values):
        raise ValueError("invalid graph6 character")
    n = values[0]
    if n == 63:
        raise ValueError("only the canonical short graph6 header is accepted")
    bit_count = n * (n - 1) // 2
    payload_count = (bit_count + 5) // 6
    if len(values) != 1 + payload_count:
        raise ValueError("wrong graph6 payload length")
    bits: list[int] = []
    for value in values[1:]:
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    if any(bits[bit_count:]):
        raise ValueError("nonzero graph6 padding")
    adjacency = [set() for _ in range(n)]
    position = 0
    for right in range(1, n):
        for left in range(right):
            if bits[position]:
                adjacency[left].add(right)
                adjacency[right].add(left)
            position += 1
    return adjacency


def edges(adjacency: list[set[int]]) -> list[tuple[int, int]]:
    return [(u, v) for u, row in enumerate(adjacency) for v in row if u < v]


def triangles(adjacency: list[set[int]]) -> list[tuple[int, int, int]]:
    return [
        (u, v, w)
        for u in range(len(adjacency))
        for v in sorted(x for x in adjacency[u] if u < x)
        for w in sorted(x for x in adjacency[u] & adjacency[v] if v < x)
    ]


def connected_after_removing(adjacency: list[set[int]], removed: set[int]) -> bool:
    remaining = set(range(len(adjacency))) - removed
    if len(remaining) <= 1:
        return True
    start = min(remaining)
    seen = {start}
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in adjacency[u] & remaining - seen:
            seen.add(v)
            queue.append(v)
    return seen == remaining


def at_least_three_connected(adjacency: list[set[int]]) -> bool:
    vertices = range(len(adjacency))
    return all(
        connected_after_removing(adjacency, set(removed))
        for size in (0, 1, 2)
        for removed in combinations(vertices, size)
    )


def components(adjacency: list[set[int]], selected: set[int]) -> list[set[int]]:
    todo = set(selected)
    answer: list[set[int]] = []
    while todo:
        start = min(todo)
        todo.remove(start)
        component = {start}
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in sorted(adjacency[u] & todo):
                todo.remove(v)
                component.add(v)
                queue.append(v)
        answer.append(component)
    return answer


def sphere_and_dual_certificate(adjacency: list[set[int]], direct_order: int) -> dict[str, object]:
    n = len(adjacency)
    edge_list = edges(adjacency)
    face_list = triangles(adjacency)
    require(connected_after_removing(adjacency, set()), "dual graph is disconnected")
    require(n - len(edge_list) + len(face_list) == 2, "Euler characteristic is not 2")

    edge_faces: dict[tuple[int, int], list[int]] = {edge: [] for edge in edge_list}
    for face_index, (a, b, c) in enumerate(face_list):
        for edge in ((a, b), (a, c), (b, c)):
            edge_faces[tuple(sorted(edge))].append(face_index)
    require(all(len(face_ids) == 2 for face_ids in edge_faces.values()),
            "an edge is not contained in exactly two triangular faces")

    for vertex in range(n):
        neighbours = set(adjacency[vertex])
        link = {u: adjacency[u] & neighbours for u in neighbours}
        require(neighbours and all(len(link[u]) == 2 for u in neighbours),
                "a vertex link is not 2-regular")
        link_seen = {min(neighbours)}
        queue = deque(link_seen)
        while queue:
            u = queue.popleft()
            for v in link[u] - link_seen:
                link_seen.add(v)
                queue.append(v)
        require(link_seen == neighbours, "a vertex link is disconnected")

    direct = [set() for _ in face_list]
    for face_ids in edge_faces.values():
        left, right = face_ids
        require(left != right, "loop in geometric dual")
        require(right not in direct[left], "parallel edge in geometric dual")
        direct[left].add(right)
        direct[right].add(left)

    require(len(face_list) == direct_order, "wrong direct fullerene order")
    require(all(len(row) == 3 for row in direct), "geometric dual is not cubic")
    require(at_least_three_connected(adjacency), "dual graph is not 3-connected")
    require(at_least_three_connected(direct), "direct graph is not 3-connected")
    profile = Counter(map(len, adjacency))
    expected = Counter({5: 12, 6: n - 12})
    require(profile == expected, f"wrong dual degree profile: {profile}")
    require(Counter(map(len, adjacency)) == Counter({5: 12, 6: direct_order // 2 - 10}),
            "wrong 5/6 face profile of geometric dual")
    return {
        "dual_vertices": n,
        "dual_edges": len(edge_list),
        "triangular_faces": len(face_list),
        "direct_vertices": len(direct),
        "direct_edges": sum(map(len, direct)) // 2,
        "degree_profile": dict(sorted(profile.items())),
    }


def literal_gsw_paths(adjacency: list[set[int]]) -> list[tuple[int, ...]]:
    degrees = list(map(len, adjacency))
    answer: set[tuple[int, ...]] = set()

    def extend(path: list[int]) -> None:
        if len(path) >= 4 and len(path) % 2 == 0:
            if degrees[path[-2]] == 6 and degrees[path[-1]] == 5:
                answer.add(tuple(path))
        for candidate in sorted(adjacency[path[-1]] & adjacency[path[-2]]):
            if candidate not in path:
                path.append(candidate)
                extend(path)
                path.pop()

    for first in range(len(adjacency)):
        if degrees[first] == 5:
            for second in sorted(adjacency[first]):
                if degrees[second] == 6:
                    extend([first, second])
    return sorted(answer)


def induced_is_path_square(adjacency: list[set[int]], path: tuple[int, ...]) -> bool:
    """Дополнительный фильтр: порождённый фрагмент ровно равен P_{2w}^2."""
    expected = {
        tuple(sorted((path[i], path[j])))
        for i in range(len(path))
        for j in range(i + 1, len(path))
        if j - i in (1, 2)
    }
    observed = {
        tuple(sorted((u, v)))
        for position, u in enumerate(path)
        for v in path[position + 1 :]
        if v in adjacency[u]
    }
    return observed == expected


def matrix(adjacency: list[set[int]]) -> list[list[int]]:
    return [[int(j in adjacency[i]) for j in range(len(adjacency))] for i in range(len(adjacency))]


def identity(n: int) -> list[list[int]]:
    return [[int(i == j) for j in range(n)] for i in range(n)]


def multiply(left, right):
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("incompatible matrix dimensions")
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right)))
         for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def power(value, exponent: int):
    if not isinstance(exponent, int):
        raise TypeError("matrix exponent must be an integer")
    if exponent < 0:
        raise ValueError("negative matrix exponent is not supported")
    result = identity(len(value))
    base = value
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        exponent //= 2
    return result


def trace(value) -> int:
    return sum(value[i][i] for i in range(len(value)))


def a_plus_xd(adjacency_matrix, x=1, adjacency_coefficient=1):
    result = [[adjacency_coefficient * entry for entry in row] for row in adjacency_matrix]
    for i, row in enumerate(adjacency_matrix):
        result[i][i] += x * sum(row)
    return result


def characteristic_polynomial(value) -> list[int]:
    n = len(value)
    unit = identity(n)
    auxiliary = identity(n)
    coefficients = [1]
    for k in range(1, n + 1):
        product = multiply(value, auxiliary)
        coefficient, remainder = divmod(-trace(product), k)
        require(remainder == 0, "inexact Faddeev-LeVerrier division")
        coefficients.append(coefficient)
        auxiliary = [
            [product[i][j] + coefficient * unit[i][j] for j in range(n)]
            for i in range(n)
        ]
    require(not any(entry for row in auxiliary for entry in row),
            "Cayley-Hamilton check failed")
    return coefficients


def moment_differences(first, second, cutoff: int) -> list[int]:
    first_power = identity(len(first))
    second_power = identity(len(second))
    answer = []
    for k in range(cutoff + 1):
        answer.append(trace(first_power) - trace(second_power))
        if k < cutoff:
            first_power = multiply(first_power, first)
            second_power = multiply(second_power, second)
    return answer


def taylor_interval(t: Fraction, moments: list[int], first, second):
    cutoff = len(moments) - 1
    partial = sum(
        Fraction(moments[k], factorial(k)) * (t / 2) ** k
        for k in range(cutoff + 1)
    )
    max_row_sum = max(
        max(sum(abs(entry) for entry in row) for row in first),
        max(sum(abs(entry) for entry in row) for row in second),
    )
    q = abs(t) * max_row_sum / 2
    require(q < cutoff + 2, "Taylor cutoff is too short")
    radius = (
        (len(first) + len(second))
        * q ** (cutoff + 1)
        / factorial(cutoff + 1)
        / (1 - q / (cutoff + 2))
    )
    return partial - radius, partial + radius


def certified_sign(t: Fraction, moments, first, second) -> int:
    lower, upper = taylor_interval(t, moments, first, second)
    return 1 if lower > 0 else -1 if upper < 0 else 0


def verify_c32(text: str) -> dict[str, object]:
    adjacency = decode_short_graph6(text)
    topology = sphere_and_dual_certificate(adjacency, 32)
    degree_six = {v for v, row in enumerate(adjacency) if len(row) == 6}
    pieces = components(adjacency, degree_six)
    require(len(pieces) == 2, "C32: T^6 does not have two components")
    require(all(len(piece) == 3 and all(len(adjacency[v] & piece) == 2 for v in piece)
                for piece in pieces), "C32: T^6 is not K3 disjoint union K3")
    paths = literal_gsw_paths(adjacency)
    witness = (3, 0, 4, 5, 13, 14)
    require(witness in paths, "C32: explicit gSW witness is absent")
    require([len(adjacency[v]) for v in witness] == [5, 6, 5, 5, 6, 5],
            "C32: wrong witness degrees")
    unoriented = {min(path, tuple(reversed(path))) for path in paths}
    induced_path_squares = [path for path in paths if induced_is_path_square(adjacency, path)]
    internal_degree_six = [
        path for path in paths
        if all(len(adjacency[v]) == 6 for v in path[2:-2])
    ]
    require(len(paths) == 72, "C32: expected 72 oriented literal gSW paths")
    require(len(unoriented) == 36, "C32: expected 36 literal paths modulo reversal")
    require(len(induced_path_squares) == 60,
            "C32: expected 60 oriented fragments induced exactly as path squares")
    require(not internal_degree_six, "C32: unexpected all-internal-degree-6 path")
    return {
        **topology,
        "T6_components_one_based": [[v + 1 for v in sorted(piece)] for piece in pieces],
        "literal_oriented_gsw": len(paths),
        "literal_mod_reversal": len(unoriented),
        "induced_equals_path_square_oriented": len(induced_path_squares),
        "all_internal_degree6": len(internal_degree_six),
    }


def verify_c44(texts: tuple[str, str]) -> dict[str, object]:
    graphs = [decode_short_graph6(text) for text in texts]
    topology = [sphere_and_dual_certificate(graph, 44) for graph in graphs]
    matrices = list(map(matrix, graphs))
    sixth = [trace(power(value, 6)) for value in matrices]
    require(sixth == [40386, 40362], "C44: wrong sixth adjacency moments")
    q_matrices = [a_plus_xd(value) for value in matrices]
    charpolys = list(map(characteristic_polynomial, q_matrices))
    require(charpolys[0] == charpolys[1], "C44: charpoly(A+D) differ")
    for x in range(7):
        difference = trace(power(a_plus_xd(matrices[0], x), 6))
        difference -= trace(power(a_plus_xd(matrices[1], x), 6))
        require(difference == 24 * (1 - x**3),
                f"C44: mixed-moment polynomial failed at x={x}")
    return {
        "topology": topology,
        "trace_A6": sixth,
        "charpoly_A_plus_D_equal": True,
        "charpoly_degree": len(charpolys[0]) - 1,
        "mixed_identity": "tr((A1+xD1)^6)-tr((A2+xD2)^6)=24(1-x^3)",
    }


def verify_c34(texts: tuple[str, str]) -> dict[str, object]:
    graphs = [decode_short_graph6(text) for text in texts]
    topology = [sphere_and_dual_certificate(graph, 34) for graph in graphs]
    triangle_counts = []
    for graph in graphs:
        selected = {v for v, row in enumerate(graph) if len(row) == 6}
        triangle_counts.append(sum(set(face) <= selected for face in triangles(graph)))
    require(triangle_counts == [0, 1], "C34: T^6 triangle counts are not [0,1]")
    matrices = list(map(matrix, graphs))
    b_matrices = [a_plus_xd(value, 1, 2) for value in matrices]
    moments = moment_differences(b_matrices[0], b_matrices[1], 82)
    first_nonzero = next((k, value) for k, value in enumerate(moments) if value)
    require(first_nonzero == (6, 960), "C34: first nonzero B-moment is not (6,960)")
    require(certified_sign(Fraction(1, 2), moments, *b_matrices) == 1,
            "C34: Delta(1/2) is not certified positive")
    require(certified_sign(Fraction(1), moments, *b_matrices) == -1,
            "C34: Delta(1) is not certified negative")
    left, right = Fraction(1, 2), Fraction(1)
    for _ in range(100):
        midpoint = (left + right) / 2
        sign = certified_sign(midpoint, moments, *b_matrices)
        require(sign != 0, "C34: undecided bisection sign")
        if sign > 0:
            left = midpoint
        else:
            right = midpoint
    expected_left = Fraction(1492341137137844774488888961329, 2535301200456458802993406410752)
    expected_right = Fraction(746170568568922387244444480665, 1267650600228229401496703205376)
    require((left, right) == (expected_left, expected_right),
            "C34: unexpected exact root bracket")
    getcontext().prec = 78
    decimal = lambda value: str(Decimal(value.numerator) / Decimal(value.denominator))
    return {
        "topology": topology,
        "T6_triangle_counts": triangle_counts,
        "first_nonzero_B_moment": first_nonzero,
        "alpha_bracket_exact": [str(left), str(right)],
        "alpha_bracket_decimal": [decimal(left), decimal(right)],
    }


def main() -> None:
    c32_texts = read_nonempty_ascii_lines(DATA / "C32.g6")
    c44_texts = read_nonempty_ascii_lines(DATA / "C44_pair.g6")
    c34_texts = read_nonempty_ascii_lines(DATA / "C34_pair.g6")
    require(c32_texts == EXPECTED_C32, "C32 machine identifier differs")
    require(c44_texts == EXPECTED_C44, "C44 machine identifiers differ")
    require(c34_texts == EXPECTED_C34, "C34 machine identifiers differ")

    print("=== INDEPENDENT STDLIB AUDIT ===")
    print("C32:", verify_c32(c32_texts[0]))
    print("C44:", verify_c44(c44_texts))
    print("C34:", verify_c34(c34_texts))
    print("INDEPENDENT STDLIB AUDIT: PASS")


if __name__ == "__main__":
    main()
