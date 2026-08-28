#!/usr/bin/env python3
"""
Детерминированные property/regression-тесты без внешних зависимостей.

Проверяются два независимых graph6-декодера, матричные степени,
характеристические многочлены, инвариантность относительно
перенумерации, gSW-перебор и рациональные интервалы C34.
"""

from __future__ import annotations

import importlib.util
import random
import sys
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


core_c32 = load("property_core_c32", ROOT / "core/audit_c32_op3_ru.py")
core_op4 = load("property_core_op4", ROOT / "core/verify_op4_ru.py")
independent = load("property_independent", ROOT / "independent_stdlib_audit.py")


def encode_short_graph6(adjacency: list[set[int]]) -> str:
    n = len(adjacency)
    if n > 62:
        raise ValueError("short graph6 only")
    bits = [
        int(left in adjacency[right])
        for right in range(1, n)
        for left in range(right)
    ]
    bits.extend([0] * ((-len(bits)) % 6))
    payload = []
    for offset in range(0, len(bits), 6):
        value = 0
        for bit in bits[offset : offset + 6]:
            value = (value << 1) | bit
        payload.append(chr(value + 63))
    return chr(n + 63) + "".join(payload)


def random_graph(n: int, rng: random.Random) -> list[set[int]]:
    adjacency = [set() for _ in range(n)]
    threshold = 20 if n > 12 else 40
    for right in range(1, n):
        for left in range(right):
            if rng.randrange(100) < threshold:
                adjacency[left].add(right)
                adjacency[right].add(left)
    return adjacency


def as_matrix(adjacency: list[set[int]]) -> list[list[int]]:
    return [[int(j in adjacency[i]) for j in range(len(adjacency))] for i in range(len(adjacency))]


def naive_power(value, exponent: int):
    result = core_op4.identity(len(value))
    for _ in range(exponent):
        result = core_op4.matmul(result, value)
    return result


def permute_adjacency(adjacency: list[set[int]], permutation: list[int]) -> list[set[int]]:
    result = [set() for _ in adjacency]
    for old_u, row in enumerate(adjacency):
        for old_v in row:
            result[permutation[old_u]].add(permutation[old_v])
    return result


def expect_failure(function, *args) -> None:
    try:
        function(*args)
    except (ValueError, TypeError, RuntimeError, ArithmeticError):
        return
    raise RuntimeError(f"expected failure from {function.__name__}")


def main() -> None:
    rng = random.Random(20260806)

    decoder_cases = 0
    for n in (0, 1, 2, 5, 18, 24, 62):
        for _ in range(8):
            adjacency = random_graph(n, rng)
            encoded = encode_short_graph6(adjacency)
            decoded_a = core_c32.decode_graph6(encoded)
            decoded_b = independent.decode_short_graph6(encoded)
            decoded_c = core_op4.decode_graph6(encoded)
            require(decoded_a == adjacency, f"C32 decoder mismatch at n={n}")
            require(decoded_b == adjacency, f"independent decoder mismatch at n={n}")
            require(decoded_c == as_matrix(adjacency), f"OP4 decoder mismatch at n={n}")
            decoder_cases += 1

    for malformed in ("", "~", ">>graph6<<~", "A@", "B~", "C??"):
        expect_failure(independent.decode_short_graph6, malformed)

    matrix_cases = 0
    charpoly_cases = 0
    for n in range(1, 7):
        for exponent in range(10):
            value = [[rng.randrange(-2, 3) for _ in range(n)] for _ in range(n)]
            require(core_op4.matrix_power(value, exponent) == naive_power(value, exponent),
                    "binary matrix power differs from naive power")
            require(independent.power(value, exponent) == naive_power(value, exponent),
                    "independent matrix power differs from naive power")
            matrix_cases += 1
        for _ in range(12):
            value = [[rng.randrange(-3, 4) for _ in range(n)] for _ in range(n)]
            require(
                core_op4.characteristic_polynomial(value)
                == independent.characteristic_polynomial(value),
                "two characteristic-polynomial implementations disagree",
            )
            charpoly_cases += 1

    expect_failure(core_op4.matrix_power, [[1]], -1)
    expect_failure(independent.power, [[1]], -1)
    expect_failure(core_op4.matrix_power, [[1]], Fraction(1, 2))
    expect_failure(independent.power, [[1]], Fraction(1, 2))

    c32 = independent.decode_short_graph6(independent.EXPECTED_C32[0])
    original_paths = independent.literal_gsw_paths(c32)
    original_unoriented = {min(path, tuple(reversed(path))) for path in original_paths}
    require((len(original_paths), len(original_unoriented)) == (72, 36),
            "unexpected base gSW counts")
    relabeling_cases = 0
    for _ in range(30):
        permutation = list(range(len(c32)))
        rng.shuffle(permutation)
        relabeled = permute_adjacency(c32, permutation)
        paths = independent.literal_gsw_paths(relabeled)
        unoriented = {min(path, tuple(reversed(path))) for path in paths}
        require((len(paths), len(unoriented)) == (72, 36),
                "gSW count changed after relabeling")
        require(independent.sphere_and_dual_certificate(relabeled, 32)["dual_edges"] == 48,
                "topological certificate changed after relabeling")
        relabeling_cases += 1

    damaged = [set(row) for row in c32]
    u, v = independent.edges(damaged)[0]
    damaged[u].remove(v)
    damaged[v].remove(u)
    expect_failure(independent.sphere_and_dual_certificate, damaged, 32)

    c44 = [independent.decode_short_graph6(text) for text in independent.EXPECTED_C44]
    target_polynomials = [
        independent.characteristic_polynomial(
            independent.a_plus_xd(independent.matrix(graph))
        )
        for graph in c44
    ]
    for graph, target in zip(c44, target_polynomials):
        permutation = list(range(len(graph)))
        rng.shuffle(permutation)
        relabeled = permute_adjacency(graph, permutation)
        observed = independent.characteristic_polynomial(
            independent.a_plus_xd(independent.matrix(relabeled))
        )
        require(observed == target, "C44 charpoly changed after relabeling")

    c34 = [independent.decode_short_graph6(text) for text in independent.EXPECTED_C34]
    b_independent = [
        independent.a_plus_xd(independent.matrix(graph), 1, 2)
        for graph in c34
    ]
    b_core = [
        core_op4.adjacency_plus_degree(independent.matrix(graph), adjacency_coefficient=2)
        for graph in c34
    ]
    moments_independent = independent.moment_differences(*b_independent, 82)
    moments_core = core_op4.moment_differences(*b_core, 82)
    require(moments_independent == moments_core, "C34 moment implementations disagree")
    interval_cases = 0
    for t in (Fraction(0), Fraction(1, 100), Fraction(1, 2), Fraction(3, 5), Fraction(9, 10), Fraction(1)):
        interval_a = independent.taylor_interval(t, moments_independent, *b_independent)
        interval_b = core_op4.taylor_interval(t, moments_core, *b_core)[:2]
        require(interval_a == interval_b, f"C34 Taylor intervals disagree at t={t}")
        interval_cases += 1

    print("PROPERTY TESTS STDLIB: PASS")
    print("graph6 decoder cases:", decoder_cases)
    print("matrix power cases:", matrix_cases)
    print("characteristic polynomial comparisons:", charpoly_cases)
    print("C32 random relabelings:", relabeling_cases)
    print("C34 interval comparisons:", interval_cases)


if __name__ == "__main__":
    main()
