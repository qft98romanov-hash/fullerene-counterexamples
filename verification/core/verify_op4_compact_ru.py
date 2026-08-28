#!/usr/bin/env python3
"""
Компактная точная проверка двух контрпримеров к Conjecture 3.

Это сокращённая версия полного файла verify_op4_ru.py. Она проверяет:

1. Пару C_44:
   - разные значения tr(A^6), следовательно, графы неизоморфны;
   - одинаковый характеристический многочлен матриц A+D;
   - отсюда совпадение характеров ch_{t,t} для всех t.

2. Пару C_34:
   - разное число треугольников в индуцированных подграфах T^6;
   - строгую смену знака функции
         Delta(t)=ch_{t,t/2}(T_1)-ch_{t,t/2}(T_2);
   - сертифицированный интервал для корня t_0 in (1/2,1).

Используется только стандартная библиотека Python.

ОГРАНИЧЕНИЕ:
эта компактная версия не повторяет полный топологический сертификат
триангуляции сферы. Для основной проверки используйте verify_op4_ru.py.
"""

from decimal import Decimal, getcontext
from fractions import Fraction as Q
from math import factorial


def require(condition: bool, message: str) -> None:
    """Обязательная проверка, работающая и при ``python -O``."""
    if not condition:
        raise RuntimeError(message)


# Точные graph6-идентификаторы двух C_44- и двух C_34-графов.
G44 = (
    "W|eMID@WH_a@E@B?__GM@?OK@_G@_G?wC?C@?@wG?@P_?@|",
    "W|eMID@WH_b@A@B?__GM@?OK@_G@oG?WC?E@??wG?@W_??~",
)

G34 = (
    "R|eMID`GH_b@A@B?_wGAF?[C?QW?{G",
    "R|eMID`GH_a@E@B?_wGB`?FC?FG?Bw",
)


def graph6(s: str) -> list[list[int]]:
    """Декодировать короткую graph6-строку в матрицу смежности."""
    if not s:
        raise ValueError("Пустая graph6-строка.")

    n = ord(s[0]) - 63
    if not 0 <= n <= 62:
        raise ValueError("Поддерживается только короткий graph6-заголовок.")

    bits: list[int] = []
    for char in s[1:]:
        value = ord(char) - 63
        if not 0 <= value <= 63:
            raise ValueError(f"Недопустимый graph6-символ: {char!r}")
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))

    needed = n * (n - 1) // 2
    padded = 6 * ((needed + 5) // 6)
    require(len(bits) == padded, "Неверная длина graph6-данных.")
    require(not any(bits[needed:]), "Ненулевые биты заполнения graph6.")

    A = [[0] * n for _ in range(n)]
    position = 0
    for j in range(1, n):
        for i in range(j):
            A[i][j] = A[j][i] = bits[position]
            position += 1
    return A


def eye(n: int) -> list[list[int]]:
    """Единичная матрица порядка n."""
    return [[int(i == j) for j in range(n)] for i in range(n)]


def tr(A) -> int:
    """След квадратной матрицы."""
    return sum(A[i][i] for i in range(len(A)))


def deg(A) -> list[int]:
    """Степени вершин по строкам матрицы смежности."""
    return list(map(sum, A))


def mm(A, B):
    """Точное матричное умножение."""
    if len(A[0]) != len(B):
        raise ValueError("Несовместимые размеры матриц.")
    return [
        [
            sum(A[i][k] * B[k][j] for k in range(len(B)))
            for j in range(len(B[0]))
        ]
        for i in range(len(A))
    ]


def mpow(A, exponent: int):
    """Бинарное возведение матрицы в неотрицательную степень."""
    if exponent < 0:
        raise ValueError("Отрицательная степень не поддерживается.")
    result = eye(len(A))
    base = A
    while exponent:
        if exponent & 1:
            result = mm(result, base)
        base = mm(base, base)
        exponent //= 2
    return result


def cA_plus_D(A, coefficient=1):
    """Построить coefficient*A + D, где D — матрица степеней."""
    M = [[coefficient * x for x in row] for row in A]
    for i, degree in enumerate(deg(A)):
        M[i][i] += degree
    return M


def A_plus_xD(A, x):
    """Построить A+xD в точной арифметике."""
    M = [row[:] for row in A]
    for i, degree in enumerate(deg(A)):
        M[i][i] += x * degree
    return M


def charpoly(M) -> list[int]:
    """
    Точный характеристический многочлен алгоритмом Фаддеева–Леверье.

    [1,c1,...,cn] обозначает
        lambda^n+c1*lambda^(n-1)+...+cn.
    """
    n = len(M)
    I = eye(n)
    B = eye(n)
    coefficients = [1]

    for k in range(1, n + 1):
        MB = mm(M, B)
        coefficient, remainder = divmod(-tr(MB), k)
        require(remainder == 0, "Неточное деление в алгоритме Фаддеева–Леверье.")
        coefficients.append(coefficient)
        B = [
            [MB[i][j] + coefficient * I[i][j] for j in range(n)]
            for i in range(n)
        ]

    # Проверка Кэли–Гамильтона.
    require(not any(value for row in B for value in row),
            "Проверка Кэли–Гамильтона не пройдена.")
    return coefficients


def pmul(a, b):
    """Перемножить многочлены, заданные коэффициентами сверху вниз."""
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


# ============================================================
# C_44: полностью целочисленный сертификат
# ============================================================

A44 = list(map(graph6, G44))
polynomials = []
sixth_moments = []

for index, A in enumerate(A44, 1):
    degrees = deg(A)
    moment6 = tr(mpow(A, 6))
    sixth_moments.append(moment6)

    print(
        "C44 граф", index,
        "V,E,степени,tr(A^6) =",
        len(A),
        sum(degrees) // 2,
        {d: degrees.count(d) for d in sorted(set(degrees))},
        moment6,
    )

    require(len(A) == 24, "C44: ожидалось 24 вершины.")
    require(sum(degrees) // 2 == 66, "C44: ожидалось 66 рёбер.")
    require(degrees.count(5) == 12, "C44: ожидалось 12 вершин степени 5.")
    require(degrees.count(6) == 12, "C44: ожидалось 12 вершин степени 6.")

    polynomials.append(charpoly(cA_plus_D(A)))

# Различие шестых моментов доказывает неизоморфность.
require(sixth_moments == [40386, 40362], "C44: неверные tr(A^6).")
require(sixth_moments[0] != sixth_moments[1], "C44: не подтверждена неизоморфность.")

# Совпадение charpoly(A+D) доказывает коспектральность A+D.
require(polynomials[0] == polynomials[1], "C44: charpoly(A+D) различаются.")

factors = [
    ([1, -4], 1),
    ([1, -18, 98, -164], 1),
    ([1, -15, 71, -107], 2),
    ([1, -26, 236, -892, 1192], 1),
    ([1, -27, 277, -1353, 3158, -2830], 2),
]

factor_product = [1]
for factor, multiplicity in factors:
    for _ in range(multiplicity):
        factor_product = pmul(factor_product, factor)

require(factor_product == polynomials[0], "C44: факторизация неверна.")
print("C44: charpoly(A+D) совпадают, факторизация проверена.")

# Разность ниже — многочлен по x степени <=6.
# Семь точных значений x=0,...,6 доказывают тождество.
for x in range(7):
    delta6 = tr(mpow(A_plus_xD(A44[0], x), 6))
    delta6 -= tr(mpow(A_plus_xD(A44[1], x), 6))
    require(delta6 == 24 * (1 - x**3),
            f"C44: тождество шестого момента не выполнено при x={x}.")

print("C44: tr((A1+xD1)^6)-tr((A2+xD2)^6) = 24(1-x^3).")


# ============================================================
# C_34: строгая смена знака и теорема о промежуточном значении
# ============================================================

A34 = list(map(graph6, G34))
triangle_counts = []

for index, A in enumerate(A34, 1):
    degrees = deg(A)
    V6 = [v for v, degree in enumerate(degrees) if degree == 6]

    triangles = [
        (u + 1, v + 1, w + 1)
        for a, u in enumerate(V6)
        for b, v in enumerate(V6[a + 1:], a + 1)
        for w in V6[b + 1:]
        if A[u][v] and A[u][w] and A[v][w]
    ]
    triangle_counts.append(len(triangles))

    print(
        "C34 граф", index,
        "V,E,степени,треугольники T6 =",
        len(A),
        sum(degrees) // 2,
        {d: degrees.count(d) for d in sorted(set(degrees))},
        triangles,
    )

    require(len(A) == 19, "C34: ожидалось 19 вершин.")
    require(sum(degrees) // 2 == 51, "C34: ожидалось 51 ребро.")
    require(degrees.count(5) == 12, "C34: ожидалось 12 вершин степени 5.")
    require(degrees.count(6) == 7, "C34: ожидалось 7 вершин степени 6.")

# Разное число треугольников в T^6 — сертификат неизоморфности.
require(triangle_counts == [0, 1], "C34: ожидались 0 и 1 треугольников в T^6.")

# B_i=2A_i+D_i, поэтому
# ch_{t,t/2}(T_i)=tr exp((t/2)B_i).
B = [cA_plus_D(A, 2) for A in A34]
N = 82
powers = [eye(19), eye(19)]
moments = []

for k in range(N + 1):
    moments.append(tr(powers[0]) - tr(powers[1]))
    if k < N:
        powers = [mm(powers[i], B[i]) for i in range(2)]

first_nonzero = next((k, value) for k, value in enumerate(moments) if value)
require(first_nonzero == (6, 960), "C34: первый ненулевой B-момент не равен (6,960).")
print("C34: первый ненулевой B-момент =", first_nonzero)


def interval(t: Q):
    """
    Строгий интервал для
        Delta(t)=tr exp((t/2)B1)-tr exp((t/2)B2).
    """
    partial = sum(
        Q(moments[k], factorial(k)) * (t / 2) ** k
        for k in range(N + 1)
    )

    # ||B_i||_2 <= ||B_i||_infinity <= 18,
    # размеры обеих матриц равны 19.
    q = 9 * t
    radius = (
        38
        * q ** (N + 1)
        / factorial(N + 1)
        / (1 - q / Q(N + 2))
    )
    return partial - radius, partial + radius, partial, radius


def sign(t: Q) -> int:
    """Вернуть строгий знак Delta(t), либо 0 при недостаточной точности."""
    lower, upper, _, _ = interval(t)
    return 1 if lower > 0 else -1 if upper < 0 else 0


getcontext().prec = 70


def dec(value: Q) -> Decimal:
    """Десятичная печать точного рационального числа."""
    return Decimal(value.numerator) / Decimal(value.denominator)


for t in (Q(1, 2), Q(1)):
    lower, upper, partial, radius = interval(t)
    print(
        "t =", t,
        "partial =", dec(partial),
        "tail <=", dec(radius),
        "sign =", sign(t),
    )
    require(lower > 0 if t == Q(1, 2) else upper < 0,
            f"C34: не сертифицирован нужный знак при t={t}.")

# Точная рациональная бисекция.
left, right = Q(1, 2), Q(1)
require(sign(left) == 1 and sign(right) == -1, "C34: исходные концы не имеют разных знаков.")

for _ in range(100):
    midpoint = (left + right) / 2
    current_sign = sign(midpoint)
    require(current_sign != 0, "C34: не хватае точности для шага бисекции.")
    if current_sign > 0:
        left = midpoint
    else:
        right = midpoint

require(sign(left) == 1, "C34: левый конец не сертифицирован положительным.")
require(sign(right) == -1, "C34: правый конец не сертифицирован отрицательным.")

print("C34: сертифицированный интервал alpha:", dec(left), dec(right))
print("C34: сертифицированный интервал beta:", dec(left / 2), dec(right / 2))
print("КОМПАКТНЫЕ СПЕКТРАЛЬНЫЕ СЕРТИФИКАТЫ: PASS")
