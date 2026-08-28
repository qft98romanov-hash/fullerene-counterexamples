#!/usr/bin/env python3
"""
Проверка смысловой согласованности файлов проверочного пакета.

Эта программа НЕ доказывает гипотезы заново. Её задача иная:
убедиться, что одни и те же сертификаты одинаково записаны в:

- graph6-файлах каталога data/;
- константах основных проверяющих программ;
- файле data/certificates.json;
- значениях, повторно вычисленных основными модулями.

Это именно внутрипакетная сверка, а не независимое математическое доказательство.
Независимый аудит запускается отдельно.

Ожидаемая структура каталогов:

package_root/
    core/
        audit_c32_op3_ru.py
        verify_op4_ru.py
        check_package_consistency_ru.py   # этот файл можно положить сюда
    data/
        C32.g6
        C34_pair.g6
        C44_pair.g6
        certificates.json

Запуск из корня пакета:

    py -3 .\\core\\check_package_consistency_ru.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path


# Не создаём каталоги __pycache__ при динамической загрузке модулей.
sys.dont_write_bytecode = True


# Этот файл предполагается расположенным в package_root/core/.
ROOT = Path(__file__).resolve().parent.parent
CORE = ROOT / "core"
DATA = ROOT / "data"


def require(condition: bool, message: str) -> None:
    """Обязательная проверка с явным диагностическим сообщением."""
    if not condition:
        raise AssertionError(message)


def load_module(name: str, path: Path):
    """
    Динамически загрузить Python-файл как модуль.

    Это позволяет использовать функции основных проверяющих программ,
    не копируя их реализацию повторно.
    """
    if not path.is_file():
        raise FileNotFoundError(path)

    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Не удалось подготовить загрузку модуля: {path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def g6_lines(name: str) -> list[str]:
    """Прочитать непустые graph6-строки из data/<name>."""
    path = DATA / name
    return [
        line.strip()
        for line in path.read_text("ascii").splitlines()
        if line.strip()
    ]


def edge_count(A) -> int:
    """Число рёбер по симметричной матрице смежности."""
    return sum(map(sum, A)) // 2


def degree_counts(A) -> dict[str, int]:
    """Кратности степеней с текстовыми ключами как в JSON."""
    degree_list = [sum(row) for row in A]
    return {
        str(degree): degree_list.count(degree)
        for degree in sorted(set(degree_list))
    }


def main() -> None:
    print("Корень пакета: .")

    # Загружаем исходные основные проверяющие программы пакета.
    c32 = load_module(
        "audit_c32_data_check",
        CORE / "audit_c32_op3_ru.py",
    )
    op4 = load_module(
        "verify_op4_data_check",
        CORE / "verify_op4_ru.py",
    )

    certificate_path = DATA / "certificates.json"
    certificate = json.loads(certificate_path.read_text("utf-8"))

    # ========================================================
    # 1. Согласованность первичных graph6-идентификаторов
    # ========================================================

    require(
        g6_lines("C32.g6") == [c32.G6],
        "C32: graph6 в data/C32.g6 не совпадает с константой программы.",
    )
    require(
        g6_lines("C34_pair.g6") == list(op4.C34_G6),
        "C34: graph6-файл не совпадает с константами программы.",
    )
    require(
        g6_lines("C44_pair.g6") == list(op4.C44_G6),
        "C44: graph6-файл не совпадает с константами программы.",
    )

    require(
        certificate["C32"]["graph6"] == c32.G6,
        "C32: graph6 в certificates.json не совпадает с программой.",
    )
    require(
        certificate["C34"]["graph6"] == list(op4.C34_G6),
        "C34: graph6 в certificates.json не совпадает с программой.",
    )
    require(
        certificate["C44"]["graph6"] == list(op4.C44_G6),
        "C44: graph6 в certificates.json не совпадает с программой.",
    )

    print("graph6-идентификаторы: PASS")

    # ========================================================
    # 2. Повторная проверка C32-данных
    # ========================================================

    A32_sets = c32.decode_graph6(c32.G6)

    require(
        len(A32_sets) == certificate["C32"]["dual_vertices"],
        "C32: неверное число вершин.",
    )
    require(
        sum(map(len, A32_sets)) // 2 == certificate["C32"]["dual_edges"],
        "C32: неверное число рёбер.",
    )

    counts32 = {
        str(degree): sum(len(row) == degree for row in A32_sets)
        for degree in (5, 6)
    }
    require(
        counts32 == certificate["C32"]["degree_multiplicities"],
        "C32: не совпадают кратности степеней.",
    )
    require(
        len(c32.all_triangles(A32_sets))
        == certificate["C32"]["triangular_faces"],
        "C32: не совпадает число треугольников.",
    )

    print("C32 JSON/код: PASS")

    # ========================================================
    # 3. Повторная проверка C44-данных
    # ========================================================

    graphs44 = [op4.decode_graph6(s) for s in op4.C44_G6]

    require(
        all(len(A) == certificate["C44"]["dual_vertices"] for A in graphs44),
        "C44: неверный порядок графов.",
    )
    require(
        all(edge_count(A) == certificate["C44"]["dual_edges"] for A in graphs44),
        "C44: неверный размер графов.",
    )
    require(
        all(
            degree_counts(A) == certificate["C44"]["degree_multiplicities"]
            for A in graphs44
        ),
        "C44: не совпадают кратности степеней.",
    )

    traces = [
        op4.trace(op4.matrix_power(A, 6))
        for A in graphs44
    ]
    require(
        traces == certificate["C44"]["trace_A6"],
        "C44: не совпадают значения tr(A^6).",
    )

    charpolys = [
        op4.characteristic_polynomial(op4.adjacency_plus_degree(A))
        for A in graphs44
    ]
    require(
        charpolys[0] == charpolys[1],
        "C44: характеристические многочлены A+D различаются.",
    )
    require(
        charpolys[0]
        == certificate["C44"]["charpoly_A_plus_D_coefficients_descending"],
        "C44: характеристический многочлен не совпадает с JSON.",
    )

    print("C44 JSON/код: PASS")

    # ========================================================
    # 4. Повторная проверка C34-данных
    # ========================================================

    graphs34 = [op4.decode_graph6(s) for s in op4.C34_G6]

    require(
        all(len(A) == certificate["C34"]["dual_vertices"] for A in graphs34),
        "C34: неверный порядок графов.",
    )
    require(
        all(edge_count(A) == certificate["C34"]["dual_edges"] for A in graphs34),
        "C34: неверный размер графов.",
    )
    require(
        all(
            degree_counts(A) == certificate["C34"]["degree_multiplicities"]
            for A in graphs34
        ),
        "C34: не совпадают кратности степеней.",
    )

    triangle_counts = [
        len(op4.induced_degree_six_triangles(A)[1])
        for A in graphs34
    ]
    require(
        triangle_counts == certificate["C34"]["T6_triangle_counts"],
        "C34: не совпадают числа треугольников в T^6.",
    )

    B1, B2 = [
        op4.adjacency_plus_degree(A, adjacency_coefficient=2)
        for A in graphs34
    ]
    moments = op4.moment_differences(B1, B2, 82)

    first = next(
        (power, value)
        for power, value in enumerate(moments)
        if value
    )
    require(
        first
        == (
            certificate["C34"]["first_nonzero_B_moment"]["power"],
            certificate["C34"]["first_nonzero_B_moment"]["difference"],
        ),
        "C34: не совпадает первый ненулевой B-момент.",
    )

    # Сопоставляем точные частичные суммы, радиусы и знаки
    # в точках t=1/2 и t=1.
    for t, key in (
        (Fraction(1, 2), "taylor_at_half"),
        (Fraction(1), "taylor_at_one"),
    ):
        lower, upper, partial, radius = op4.taylor_interval(
            t,
            moments,
            B1,
            B2,
        )
        record = certificate["C34"][key]

        require(
            partial == Fraction(record["partial_fraction"]),
            f"{key}: не совпадает частичная сумма.",
        )
        require(
            radius == Fraction(record["radius_fraction"]),
            f"{key}: не совпадает радиус остатка.",
        )

        sign = 1 if lower > 0 else -1 if upper < 0 else 0
        require(
            sign == record["certified_sign"],
            f"{key}: не совпадает сертифицированный знак.",
        )

    bracket = certificate["C34"]["alpha_root_bracket"]
    left = Fraction(bracket["left_fraction"])
    right = Fraction(bracket["right_fraction"])

    require(
        Fraction(1, 2) < left < right < 1,
        "C34: интервал корня не лежит внутри (1/2,1).",
    )
    require(
        op4.certified_sign(left, moments, B1, B2) == 1,
        "C34: на левой границе интервала знак не положителен.",
    )
    require(
        op4.certified_sign(right, moments, B1, B2) == -1,
        "C34: на правой границе интервала знак не отрицателен.",
    )

    print("C34 JSON/код: PASS")
    print("СМЫСЛОВАЯ СОГЛАСОВАННОСТЬ ПАКЕТА: PASS")


if __name__ == "__main__":
    main()
