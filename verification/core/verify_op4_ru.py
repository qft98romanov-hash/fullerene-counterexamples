#!/usr/bin/env python3
"""
Точная проверка двух контрпримеров к Conjecture 3 из статьи
A. Bille, V. Buchstaber, E. Spodarev,
"Some Open Mathematical Problems on Fullerenes".

Проверка C_44 полностью целочисленная:
двойственные графы неизоморфны, но матрицы A+D имеют одинаковый
характеристический многочлен. Поэтому их характеры совпадают на
диагонали alpha=beta.

Проверка C_34 использует точные рациональные частичные суммы ряда
Тейлора и строгую оценку хвоста. Это позволяет сертифицировать
смену знака функции разности характеров на луче beta=alpha/2.

Используется только стандартная библиотека Python. Все доказательные
проверки явные и сохраняются в режиме ``python -O``.
"""

from decimal import Decimal, getcontext
from fractions import Fraction
from math import factorial


def require(condition, message):
    """Обязательная проверка, не отключаемая ключом ``python -O``."""
    if not condition:
        raise RuntimeError(message)


C44_G6 = (
    "W|eMID@WH_a@E@B?__GM@?OK@_G@_G?wC?C@?@wG?@P_?@|",
    "W|eMID@WH_b@A@B?__GM@?OK@_G@oG?WC?E@??wG?@W_??~",
)

C34_G6 = (
    "R|eMID`GH_b@A@B?_wGAF?[C?QW?{G",
    "R|eMID`GH_a@E@B?_wGB`?FC?FG?Bw",
)



# Декодируем короткий формат graph6 в точную матрицу смежности.
# Для наших графов число вершин меньше 63, поэтому короткого заголовка достаточно.
def decode_graph6(s):
    """Decode the short graph6 format (0 <= number of vertices <= 62)."""
    prefix = ">>graph6<<"
    if s.startswith(prefix):
        s = s[len(prefix) :]
    if not s:
        raise ValueError("empty graph6 string")

    n = ord(s[0]) - 63
    if not 0 <= n <= 62:
        raise ValueError("реализованы только короткие заголовки graph6")

    bits = []
    for ch in s[1:]:
        value = ord(ch) - 63
        if not 0 <= value <= 63:
            raise ValueError("invalid graph6 character")
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))

    needed = n * (n - 1) // 2
    padded_length = 6 * ((needed + 5) // 6)
    if len(bits) != padded_length or any(bits[needed:]):
        raise ValueError("bad graph6 length or nonzero padding")

    A = [[0] * n for _ in range(n)]
    pos = 0
    # graph6 order: (0,1), (0,2),(1,2), (0,3),(1,3),(2,3), ...
    for j in range(1, n):
        for i in range(j):
            A[i][j] = A[j][i] = bits[pos]
            pos += 1
    return A



# Единичная матрица порядка n.
def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]



# Точное матричное умножение по стандартной формуле (AB)_{ij}=sum_k A_{ik}B_{kj}.
def matmul(A, B):
    rows, middle, cols = len(A), len(B), len(B[0])
    if len(A[0]) != middle:
        raise ValueError("incompatible matrix dimensions")
    return [
        [sum(A[i][k] * B[k][j] for k in range(middle)) for j in range(cols)]
        for i in range(rows)
    ]



# Бинарное возведение матрицы в неотрицательную целую степень.
def matrix_power(A, exponent):
    if not isinstance(exponent, int):
        raise TypeError("matrix exponent must be an integer")
    if exponent < 0:
        raise ValueError("negative matrix exponent is not supported")
    result = identity(len(A))
    base = A
    while exponent:
        if exponent & 1:
            result = matmul(result, base)
        base = matmul(base, base)
        exponent //= 2
    return result



# След квадратной матрицы.
def trace(A):
    return sum(A[i][i] for i in range(len(A)))



# Степени вершин: сумма элементов соответствующей строки матрицы смежности.
def degrees(A):
    return [sum(row) for row in A]



# Строим cA+D, где D — диагональная матрица степеней.
def adjacency_plus_degree(A, adjacency_coefficient=1):
    """Return adjacency_coefficient*A + D, with D the degree matrix."""
    M = [[adjacency_coefficient * x for x in row] for row in A]
    for i, degree in enumerate(degrees(A)):
        M[i][i] += degree
    return M



# Строим A+xD; параметр x может быть точным рациональным числом.
def adjacency_plus_x_degree(A, x):
    """Return A+xD (x may be any exact scalar)."""
    M = [row[:] for row in A]
    for i, degree in enumerate(degrees(A)):
        M[i][i] += x * degree
    return M



# Точный характеристический многочлен алгоритмом Фаддеева—Леверье.
# Финальная нулевая матрица служит встроенной проверкой теоремы Кэли—Гамильтона.
def characteristic_polynomial(M):
    """Exact Faddeev-LeVerrier coefficients, in descending order.

    The return value [1,c1,...,cn] represents
    x^n+c1*x^(n-1)+...+cn.
    """
    n = len(M)
    I = identity(n)
    B = I
    coefficients = [1]
    for k in range(1, n + 1):
        MB = matmul(M, B)
        quotient, remainder = divmod(-trace(MB), k)
        if remainder:
            raise ArithmeticError("Faddeev-LeVerrier division was not exact")
        coefficient = quotient
        coefficients.append(coefficient)
        B = [
            [MB[i][j] + coefficient * I[i][j] for j in range(n)]
            for i in range(n)
        ]
    if any(entry for row in B for entry in row):
        raise ArithmeticError("Cayley-Hamilton check failed")
    return coefficients



# Точное перемножение многочленов, заданных списками коэффициентов.
def polynomial_product(p, q):
    """Multiply coefficient lists in descending order."""
    result = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            result[i + j] += a * b
    return result



# Выделяем вершины степени 6 и считаем треугольники в индуцированном подграфе T^6.
# Для C_34 различие этих чисел является простым сертификатом неизоморфности.
def induced_degree_six_triangles(A):
    vertices = [i for i, degree in enumerate(degrees(A)) if degree == 6]
    triangles = []
    for a, i in enumerate(vertices):
        for b in range(a + 1, len(vertices)):
            j = vertices[b]
            for k in vertices[b + 1 :]:
                if A[i][j] and A[i][k] and A[j][k]:
                    triangles.append((i, j, k))
    return vertices, triangles



# Рёбра индуцированного подграфа на вершинах степени 6.
def induced_degree_six_edges(A):
    vertices = [i for i, degree in enumerate(degrees(A)) if degree == 6]
    edges = [
        (i, j)
        for position, i in enumerate(vertices)
        for j in vertices[position + 1 :]
        if A[i][j]
    ]
    return vertices, edges



# Комбинаторный сертификат триангуляции двумерной сферы:
# связность, два треугольника на каждом ребре, циклические линки и chi=2.
def sphere_triangulation_certificate(A):
    """Check an explicit combinatorial 2-sphere certificate.

    We take every 3-clique as a triangular face, verify that every edge is in
    exactly two faces and every vertex link is one cycle, and check Euler
    characteristic 2.  Thus the clique complex is a connected closed
    triangulated surface with Euler characteristic 2, hence a 2-sphere.
    """
    n = len(A)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if A[i][j]]
    faces = [
        (i, j, k)
        for i in range(n)
        for j in range(i + 1, n)
        if A[i][j]
        for k in range(j + 1, n)
        if A[i][k] and A[j][k]
    ]
    face_sets = [set(face) for face in faces]

    graph_seen = set()
    graph_stack = [0]
    while graph_stack:
        u = graph_stack.pop()
        if u not in graph_seen:
            graph_seen.add(u)
            graph_stack.extend(v for v in range(n) if A[u][v] and v not in graph_seen)
    if len(graph_seen) != n:
        return False
    if n - len(edges) + len(faces) != 2:
        return False
    if any(sum({i, j} <= face for face in face_sets) != 2 for i, j in edges):
        return False

    for vertex in range(n):
        neighbours = [i for i in range(n) if A[vertex][i]]
        link_edges = {
            tuple(sorted((next(iter(face - {vertex, other})), other)))
            for face in face_sets
            if vertex in face
            for other in face - {vertex}
        }
        # The comprehension above produces each link edge twice, but the set
        # removes duplicates.  A finite 2-regular connected graph is a cycle.
        link_adjacency = {u: set() for u in neighbours}
        for u, v in link_edges:
            link_adjacency[u].add(v)
            link_adjacency[v].add(u)
        if any(len(link_adjacency[u]) != 2 for u in neighbours):
            return False
        seen = set()
        stack = [neighbours[0]]
        while stack:
            u = stack.pop()
            if u not in seen:
                seen.add(u)
                stack.extend(link_adjacency[u] - seen)
        if seen != set(neighbours):
            return False
    return True



# Точные разности моментов tr(B1^k)-tr(B2^k) до заданного порядка.
def moment_differences(B1, B2, cutoff):
    """Return tr(B1^k)-tr(B2^k), k=0,...,cutoff, exactly."""
    powers = [identity(len(B1)), identity(len(B2))]
    answer = []
    for k in range(cutoff + 1):
        answer.append(trace(powers[0]) - trace(powers[1]))
        if k != cutoff:
            powers = [matmul(powers[0], B1), matmul(powers[1], B2)]
    return answer



# Строгий рациональный интервал для разности следов матричных экспонент.
# Частичная сумма точна; хвост оценивается через норму матриц и геометрическую мажоранту.
def taylor_interval(t, moments, B1, B2):
    """Certified interval for

       Delta(t) = tr(exp(t*B1/2)) - tr(exp(t*B2/2)).

    The Taylor partial sum is a Fraction.  Since B1,B2 are symmetric,
    |tr((t*Bi/2)^k)| <= ni*q^k, where
    q=(|t|/2)*max_i ||Bi||_infinity.  Bounding the remaining exponential
    series by a geometric majorant gives the exact rational radius below.
    """
    if not isinstance(t, Fraction):
        t = Fraction(t)
    N = len(moments) - 1
    partial = sum(
        Fraction(moments[k], factorial(k)) * (t / 2) ** k
        for k in range(N + 1)
    )

    max_row_sum = max(
        max(sum(abs(x) for x in row) for row in B1),
        max(sum(abs(x) for x in row) for row in B2),
    )
    # Абсолютная величина делает оценку корректной и для отрицательных t,
    # хотя ниже используются только положительные значения параметра.
    q = abs(t) * max_row_sum / 2
    if q >= N + 2:
        raise ValueError("Taylor cutoff too small for the geometric tail bound")
    radius = (
        (len(B1) + len(B2))
        * q ** (N + 1)
        / factorial(N + 1)
        / (1 - q / (N + 2))
    )
    return partial - radius, partial + radius, partial, radius



# Возвращаем строгий знак функции, если весь гарантированный интервал
# лежит по одну сторону от нуля.
def certified_sign(t, moments, B1, B2):
    lower, upper, _, _ = taylor_interval(t, moments, B1, B2)
    if lower > 0:
        return 1
    if upper < 0:
        return -1
    return 0



# Десятичная печать точного рационального числа; в доказательных сравнениях не используется.
def decimal_string(x, digits=60):
    getcontext().prec = digits
    return str(Decimal(x.numerator) / Decimal(x.denominator))



# Основной целочисленный сертификат для пары C_44:
# неизоморфность по tr(A^6), затем коспектральность A+D.
def verify_c44():
    print("=== C44: exact integral certificate ===")
    graphs = [decode_graph6(s) for s in C44_G6]
    polynomials = []
    sixth_traces = []
    for index, A in enumerate(graphs, 1):
        deg = degrees(A)
        edge_number = sum(deg) // 2
        sphere_certificate = sphere_triangulation_certificate(A)
        sixth_trace = trace(matrix_power(A, 6))
        require(len(A) == 24, "C44: expected 24 dual vertices")
        require(edge_number == 66, "C44: expected 66 dual edges")
        require(deg.count(5) == 12 and deg.count(6) == 12,
                "C44: wrong degree multiplicities")
        require(set(deg) == {5, 6}, "C44: an unexpected vertex degree was found")
        require(sphere_certificate, "C44: triangulated-sphere certificate failed")
        print(
            "graph", index,
            "vertices", len(A),
            "edges", edge_number,
            "degree multiset", {d: deg.count(d) for d in sorted(set(deg))},
        )
        print("graph", index, "triangulated 2-sphere certificate", sphere_certificate)
        print("graph", index, "tr(A^6)", sixth_trace)
        sixth_traces.append(sixth_trace)
        polynomials.append(characteristic_polynomial(adjacency_plus_degree(A)))

    require(sixth_traces == [40386, 40362], "C44: unexpected tr(A^6) values")
    require(polynomials[0] == polynomials[1], "C44: charpoly(A+D) differ")
    print("charpoly(A+D) equal:", True)
    print("common coefficients:")
    print(polynomials[0])

    factors = (
        ([1, -4], 1),
        ([1, -18, 98, -164], 1),
        ([1, -15, 71, -107], 2),
        ([1, -26, 236, -892, 1192], 1),
        ([1, -27, 277, -1353, 3158, -2830], 2),
    )
    product = [1]
    for factor, multiplicity in factors:
        for _ in range(multiplicity):
            product = polynomial_product(product, factor)
    require(product == polynomials[0], "C44: displayed factorization is incorrect")
    print("displayed factorization verified:", True)

    # Both sides below are polynomials in x of degree at most six.  Checking
    # seven distinct integer x-values therefore proves the identity exactly.
    for x in range(7):
        difference = trace(matrix_power(adjacency_plus_x_degree(graphs[0], x), 6))
        difference -= trace(matrix_power(adjacency_plus_x_degree(graphs[1], x), 6))
        require(difference == 24 * (1 - x**3),
                f"C44: mixed sixth-moment identity failed at x={x}")
    print("tr((A1+xD1)^6)-tr((A2+xD2)^6) = 24*(1-x^3):", True)

    # Дополнительный сертификат согласует программу с замечанием в тексте:
    # именно эта C44-пара не сталкивается в рекомендованной точке
    # (alpha,beta)=(1/2,1/4). Здесь снова B=2A+D и t=1/2.
    B1, B2 = [
        adjacency_plus_degree(A, adjacency_coefficient=2)
        for A in graphs
    ]
    recommended_moments = moment_differences(B1, B2, 100)
    lower, upper, partial, radius = taylor_interval(
        Fraction(1, 2), recommended_moments, B1, B2
    )
    displayed_lower = Fraction(235725175552915, 10**17)
    displayed_upper = Fraction(235725175552916, 10**17)
    require(lower > displayed_lower, "C44: lower bound at (1/2,1/4) is too small")
    require(upper < displayed_upper, "C44: upper bound at (1/2,1/4) is too large")
    require(radius < Fraction(6, 10**93), "C44: Taylor radius is too large")
    print("recommended point (1/2,1/4) Taylor partial =", decimal_string(partial))
    print("recommended point rigorous radius =", decimal_string(radius))
    print("recommended point displayed interval verified:", True)
    print()



# Основной интервальный сертификат для пары C_34:
# смена знака Delta(t) и точная рациональная бисекция.
def verify_c34():
    print("=== C34: certified IVT collision on beta=alpha/2 ===")
    graphs = [decode_graph6(s) for s in C34_G6]
    expected_vertices = (
        [1, 3, 5, 8, 10, 13, 16],
        [1, 3, 4, 14, 15, 16, 17],
    )
    expected_edges = (
        [(1, 3), (1, 5), (3, 10), (5, 13), (8, 16), (13, 16)],
        [(1, 3), (1, 4), (3, 4), (14, 15), (15, 16), (16, 17)],
    )
    degree_six_triangle_counts = []
    for index, A in enumerate(graphs, 1):
        deg = degrees(A)
        edge_number = sum(deg) // 2
        sphere_certificate = sphere_triangulation_certificate(A)
        vertices, triangles = induced_degree_six_triangles(A)
        edge_vertices, induced_edges = induced_degree_six_edges(A)
        one_based_vertices = [v + 1 for v in vertices]
        one_based_edges = [(u + 1, v + 1) for u, v in induced_edges]
        require(len(A) == 19, "C34: expected 19 dual vertices")
        require(edge_number == 51, "C34: expected 51 dual edges")
        require(deg.count(5) == 12 and deg.count(6) == 7,
                "C34: wrong degree multiplicities")
        require(set(deg) == {5, 6}, "C34: an unexpected vertex degree was found")
        require(sphere_certificate, "C34: triangulated-sphere certificate failed")
        require(edge_vertices == vertices, "C34: inconsistent degree-six vertex sets")
        require(one_based_vertices == expected_vertices[index - 1],
                "C34: degree-six vertices differ from the signed certificate")
        require(one_based_edges == expected_edges[index - 1],
                "C34: T^6 edges differ from the signed certificate")
        print(
            "graph", index,
            "vertices", len(A),
            "edges", edge_number,
            "degree multiset", {d: deg.count(d) for d in sorted(set(deg))},
        )
        print("triangulated 2-sphere certificate:", sphere_certificate)
        print("degree-6 vertices (one-based):", one_based_vertices)
        print("edges in induced degree-6 graph (one-based):", one_based_edges)
        print(
            "triangles in induced degree-6 graph (one-based):",
            [tuple(v + 1 for v in triangle) for triangle in triangles],
        )
        degree_six_triangle_counts.append(len(triangles))

    require(degree_six_triangle_counts == [0, 1],
            "C34: expected 0 and 1 triangles in the two T^6 graphs")

    # t*A+(t/2)*D = (t/2)*(2*A+D).
    B1, B2 = [adjacency_plus_degree(A, adjacency_coefficient=2) for A in graphs]
    cutoff = 82
    moments = moment_differences(B1, B2, cutoff)
    first_nonzero_moment = next((k, x) for k, x in enumerate(moments) if x)
    require(first_nonzero_moment == (6, 960),
            "C34: first nonzero B-moment is not (6,960)")
    print("first nonzero moment:", first_nonzero_moment)

    for t in (Fraction(1, 2), Fraction(1)):
        lower, upper, partial, radius = taylor_interval(t, moments, B1, B2)
        print("t =", t)
        print("  Taylor partial =", decimal_string(partial))
        print("  rigorous radius =", decimal_string(radius))
        print("  certified sign =", "+" if lower > 0 else "-")
        require(lower > 0 if t == Fraction(1, 2) else upper < 0,
                f"C34: required strict sign was not certified at t={t}")

    # Exact-rational bisection.  Every sign decision uses a rigorous interval.
    left, right = Fraction(1, 2), Fraction(1)
    require(certified_sign(left, moments, B1, B2) == 1,
            "C34: Delta(1/2) was not certified positive")
    require(certified_sign(right, moments, B1, B2) == -1,
            "C34: Delta(1) was not certified negative")
    for _ in range(100):
        midpoint = (left + right) / 2
        sign = certified_sign(midpoint, moments, B1, B2)
        if sign == 1:
            left = midpoint
        elif sign == -1:
            right = midpoint
        else:
            raise ArithmeticError("increase Taylor cutoff")

    expected_left = Fraction(
        1492341137137844774488888961329,
        2535301200456458802993406410752,
    )
    expected_right = Fraction(
        746170568568922387244444480665,
        1267650600228229401496703205376,
    )
    require(left == expected_left, "C34: unexpected left bisection endpoint")
    require(right == expected_right, "C34: unexpected right bisection endpoint")
    require(certified_sign(left, moments, B1, B2) == 1,
            "C34: left endpoint is not certified positive")
    require(certified_sign(right, moments, B1, B2) == -1,
            "C34: right endpoint is not certified negative")

    print("certified sign-changing alpha bracket:")
    print(" exact left =", left)
    print(" exact right =", right)
    print(" ", decimal_string(left, 75))
    print(" ", decimal_string(right, 75))
    print("certified beta=alpha/2 bracket:")
    print(" ", decimal_string(left / 2, 75))
    print(" ", decimal_string(right / 2, 75))


if __name__ == "__main__":
    verify_c44()
    verify_c34()
