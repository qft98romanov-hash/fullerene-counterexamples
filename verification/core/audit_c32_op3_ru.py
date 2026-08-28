"""
Аннотированный точный сертификат для C32-контрпримера
к Conjectures 1 и 2 из статьи:

A. Bille, V. Buchstaber, E. Spodarev,
"Some Open Mathematical Problems on Fullerenes",
J. Chem. Inf. Model. 65 (2025), 2911–2923.
DOI: 10.1021/acs.jcim.4c01997

Навигация по статье:
- p. 2913: двойственный фуллерен T_n, подграфы T_n^5 и T_n^6;
- p. 2915: Definition 1 (gSW path);
- pp. 2916–2917: Definition 2 (t-triangle);
- pp. 2917–2918: Construction 1 и Definition 4 (cut-partition);
- p. 2918: Conjectures 1 и 2.

Скрипт использует только стандартную библиотеку Python.
"""

from collections import Counter, deque
from typing import Iterable


def require(condition: bool, message: str) -> None:
    """Обязательная доказательная проверка.

    Она не исчезает при ``python -O``.
    """
    if not condition:
        raise RuntimeError(message)


# Строка graph6 является первичным и однозначным машинным
# идентификатором рассматриваемого конечного графа.
#
# После декодирования получается граф на 18 вершинах.
# Он будет интерпретирован как двойственный граф T фуллерена C_32.
G6 = "Q|eMID@WH?e@E@B?_wGBB?MC?NW"


def decode_graph6(s: str) -> list[set[int]]:
    """
    Декодировать строку формата graph6 в список множеств соседей.

    Возвращаемый объект adj устроен так:
        adj[v] = множество всех вершин, смежных с v.

    Формат graph6 хранит верхний треугольник матрицы смежности
    в порядке
        (0,1), (0,2), (1,2), (0,3), (1,3), (2,3), ...

    Эта функция поддерживает все три стандартных варианта заголовка
    порядка графа и дополнительно проверяет длину полезной нагрузки.
    """
    s = s.strip()

    # Иногда строка явно снабжается префиксом формата.
    if s.startswith(">>graph6<<"):
        s = s[len(">>graph6<<") :]

    if not s:
        raise ValueError("Пустая строка graph6.")

    # Каждый символ graph6 кодирует шестибитное число:
    # из ASCII-кода вычитается 63.
    vals = [ord(c) - 63 for c in s]
    if any(not 0 <= x <= 63 for x in vals):
        raise ValueError("В строке присутствует недопустимый graph6-символ.")

    # Читаем число вершин n.
    if vals[0] != 63:
        # Короткий заголовок: n <= 62.
        n, pos = vals[0], 1
    elif len(vals) >= 4 and vals[1] != 63:
        # Средний заголовок: 18 бит на n.
        n = (vals[1] << 12) | (vals[2] << 6) | vals[3]
        pos = 4
    elif len(vals) >= 8:
        # Длинный заголовок: 36 бит на n.
        n = 0
        for x in vals[2:8]:
            n = (n << 6) | x
        pos = 8
    else:
        raise ValueError("Повреждённый заголовок graph6.")

    # Для простого неориентированного графа надо закодировать
    # n(n-1)/2 потенциальных рёбер.
    need_bits = n * (n - 1) // 2
    need_chars = (need_bits + 5) // 6

    payload = vals[pos:]
    if len(payload) != need_chars:
        raise ValueError(
            f"Неверная длина кодовой части строки graph6: "
            f"ожидалось {need_chars} символов, получено {len(payload)}."
        )

    bits: list[int] = []
    for x in payload:
        bits.extend((x >> shift) & 1 for shift in range(5, -1, -1))

    # Последний символ при необходимости дополняется нулевыми битами.
    if any(bits[need_bits:]):
        raise ValueError("Ненулевые биты заполнения в конце graph6-строки.")

    # Строим симметричный список смежности.
    adj = [set() for _ in range(n)]
    k = 0
    for v in range(1, n):
        for u in range(v):
            if bits[k]:
                adj[u].add(v)
                adj[v].add(u)
            k += 1

    return adj


def connected_components(
    adj: list[set[int]],
    vertices: Iterable[int],
) -> list[set[int]]:
    """
    Найти компоненты связности подграфа, индуцированного vertices.

    Используется обычный обход в ширину. Мы не создаём отдельную
    матрицу смежности подграфа: при обходе просто пересекаем
    adj[u] с множеством ещё не посещённых выбранных вершин.
    """
    todo = set(vertices)
    result: list[set[int]] = []

    while todo:
        root = min(todo)
        component = {root}
        todo.remove(root)

        queue = deque([root])
        while queue:
            u = queue.popleft()
            for v in sorted(adj[u] & todo):
                todo.remove(v)
                component.add(v)
                queue.append(v)

        result.append(component)

    return result


def all_triangles(adj: list[set[int]]) -> list[tuple[int, int, int]]:
    """
    Перечислить все 3-клики графа, то есть все тройки попарно
    смежных вершин.

    Для данного конкретного графа дальнейшие проверки показывают,
    что эти 3-клики образуют двумерные симплексы триангуляции сферы:
    каждое ребро лежит ровно в двух треугольниках, линк каждой
    вершины является циклом, а характеристика Эйлера равна 2.
    """
    triangles: list[tuple[int, int, int]] = []

    for a in range(len(adj)):
        for b in sorted(v for v in adj[a] if a < v):
            for c in sorted(v for v in (adj[a] & adj[b]) if b < v):
                triangles.append((a, b, c))

    return triangles


def induced_degree(
    adj: list[set[int]],
    vertex: int,
    vertices: set[int],
) -> int:
    """Степень vertex в подграфе, индуцированном множеством vertices."""
    return len(adj[vertex] & vertices)


def is_k3_component(adj: list[set[int]], component: set[int]) -> bool:
    """
    Проверить, что индуцированная компонента является K_3.

    Для трёх вершин это эквивалентно тому, что каждая из них
    имеет внутреннюю степень 2.
    """
    return (
        len(component) == 3
        and all(induced_degree(adj, u, component) == 2 for u in component)
    )


def is_gsw_path(adj: list[set[int]], path: list[int]) -> bool:
    """
    Проверить Definition 1(a) из статьи буквально.

    Для path = (v_1, ..., v_{2w}) требуется:
    1) w >= 2;
    2) v_1 и v_{2w} имеют степень 5;
    3) v_2 и v_{2w-1} имеют степень 6;
    4) все вершины различны;
    5) последовательные пары являются рёбрами, поскольку это путь;
    6) (v_i, v_{i+2}) является ребром для i=1,...,2w-2.

    В опубликованном Definition 1 нет условия, что промежуточные
    вершины v_3,...,v_{2w-2} обязаны иметь степень 6.
    """
    if len(path) % 2 != 0:
        return False

    w = len(path) // 2
    if w < 2:
        return False

    if len(set(path)) != len(path):
        return False

    degree = [len(neighbors) for neighbors in adj]

    if degree[path[0]] != 5 or degree[path[-1]] != 5:
        return False

    if degree[path[1]] != 6 or degree[path[-2]] != 6:
        return False

    # Обычные рёбра пути (v_i, v_{i+1}).
    if any(path[i + 1] not in adj[path[i]] for i in range(len(path) - 1)):
        return False

    # Дополнительные рёбра (v_i, v_{i+2}) из Definition 1.
    if any(path[i + 2] not in adj[path[i]] for i in range(len(path) - 2)):
        return False

    return True


def sphere_triangulation_certificate(
    adj: list[set[int]],
) -> dict[str, object]:
    """
    Построить комбинаторный сертификат триангуляции S^2.

    Математическая мотивация:
    статья рассматривает двойственный фуллерен T_n как
    триангуляцию сферы. Чтобы не принимать это на веру для
    входной graph6-строки, мы независимо проверяем:

    - граф связен;
    - все двумерные клетки берутся как 3-клики;
    - каждое ребро входит ровно в два треугольника;
    - линк каждой вершины является одним циклом;
    - V - E + F = 2.

    Первые четыре свойства дают связное замкнутое
    комбинаторное 2-многообразие, а характеристика Эйлера 2
    идентифицирует его со сферой.
    """
    n = len(adj)
    edges = sum(len(neighbors) for neighbors in adj) // 2
    triangles = all_triangles(adj)

    components = connected_components(adj, range(n))
    connected = len(components) == 1

    # Сколько треугольников инцидентно каждому ребру.
    edge_triangle_count: Counter[tuple[int, int]] = Counter()
    for a, b, c in triangles:
        for u, v in ((a, b), (a, c), (b, c)):
            edge_triangle_count[tuple(sorted((u, v)))] += 1

    each_edge_in_two_triangles = (
        len(edge_triangle_count) == edges
        and set(edge_triangle_count.values()) == {2}
    )

    # Линк вершины v: граф на соседях v, где два соседа
    # соединены, если они смежны в исходном графе.
    # Для триангулированного 2-многообразия линк должен быть циклом.
    vertex_links_are_cycles = True

    for v in range(n):
        neighbors = set(adj[v])

        if not neighbors:
            vertex_links_are_cycles = False
            break

        link_adj = {
            u: adj[u] & neighbors
            for u in neighbors
        }

        # В цикле каждая вершина имеет степень 2.
        if any(len(link_adj[u]) != 2 for u in neighbors):
            vertex_links_are_cycles = False
            break

        # И весь линк должен быть связен.
        seen = set()
        queue = deque([min(neighbors)])
        seen.add(queue[0])

        while queue:
            u = queue.popleft()
            for z in link_adj[u]:
                if z not in seen:
                    seen.add(z)
                    queue.append(z)

        if seen != neighbors:
            vertex_links_are_cycles = False
            break

    return {
        "vertices": n,
        "edges": edges,
        "triangles": len(triangles),
        "connected": connected,
        "each_edge_in_two_triangles": each_edge_in_two_triangles,
        "vertex_links_are_cycles": vertex_links_are_cycles,
        "euler_characteristic": n - edges + len(triangles),
    }


def main() -> None:
    # ------------------------------------------------------------
    # 1. Декодируем точный граф и вычисляем базовые данные.
    # ------------------------------------------------------------
    adj = decode_graph6(G6)
    n = len(adj)
    edges = sum(len(neighbors) for neighbors in adj) // 2
    degree = [len(neighbors) for neighbors in adj]

    # В двойственном фуллерене:
    # - вершины степени 5 соответствуют пятиугольникам;
    # - вершины степени 6 соответствуют шестиугольникам.
    #
    # Именно эти индуцированные подграфы статья обозначает T_n^5 и T_n^6.
    v5 = {v for v, d in enumerate(degree) if d == 5}
    v6 = {v for v, d in enumerate(degree) if d == 6}

    # ------------------------------------------------------------
    # 2. Независимо подтверждаем, что входной граф является
    #    триангуляцией сферы нужного степенного типа.
    # ------------------------------------------------------------
    sphere = sphere_triangulation_certificate(adj)

    # Для C_32 двойственный граф должен иметь:
    # |V(T)| = 32/2 + 2 = 18,
    # |E(T)| = 3*32/2 = 48,
    # 32 треугольные грани,
    # 12 вершин степени 5 и 6 вершин степени 6.
    require(n == 18, f"Ожидалось 18 вершин, получено {n}.")
    require(edges == 48, f"Ожидалось 48 рёбер, получено {edges}.")
    require(degree.count(5) == 12, "Должно быть ровно 12 вершин степени 5.")
    require(degree.count(6) == 6, "Должно быть ровно 6 вершин степени 6.")
    require(set(degree) == {5, 6}, "Обнаружены недопустимые степени вершин.")

    require(sphere["connected"] is True, "Граф должен быть связен.")
    require(sphere["triangles"] == 32, "Для C_32 ожидаются 32 треугольные грани.")
    require(sphere["each_edge_in_two_triangles"] is True,
            "Каждое ребро должно входить ровно в два треугольника.")
    require(sphere["vertex_links_are_cycles"] is True,
            "Линк каждой вершины должен быть циклом.")
    require(sphere["euler_characteristic"] == 2,
            "Характеристика Эйлера должна быть равна 2.")

    # ------------------------------------------------------------
    # 3. Строим T^6 — подграф на вершинах степени 6.
    #
    #    Это именно тот объект, к которому статья применяет
    #    Definition 2, Construction 1 и Definition 4.
    # ------------------------------------------------------------
    t6_components = connected_components(adj, v6)

    each_component_is_k3 = (
        len(t6_components) == 2
        and all(is_k3_component(adj, comp) for comp in t6_components)
    )

    require(each_component_is_k3, "Ожидалось T^6 ≅ K_3 ⊔ K_3.")

    t6_edges = {
        tuple(sorted((u, v)))
        for u in v6
        for v in (adj[u] & v6)
        if u < v
    }
    require(len(t6_edges) == 6, "В T^6 должно быть ровно 6 рёбер.")

    # ------------------------------------------------------------
    # 4. Почему каждая компонента K_3 является 1-треугольником?
    #
    #    Definition 2 при t=1 требует:
    #    - t^2 = 1 внутреннюю треугольную грань;
    #    - (t+1)(t+2)/2 = 3 вершины;
    #    - ровно 3 вершины степени 2;
    #    - 3(t-1) = 0 вершин степени 4;
    #    - остальных вершин степени 6 нет.
    #
    #    Компонента K_3 удовлетворяет этим условиям буквально.
    # ------------------------------------------------------------
    each_component_is_t1_triangle = all(
        len(comp) == 3
        and all(induced_degree(adj, u, comp) == 2 for u in comp)
        for comp in t6_components
    )
    require(each_component_is_t1_triangle,
            "Обе компоненты T^6 должны быть 1-треугольниками.")

    # ------------------------------------------------------------
    # 5. Почему Construction 1 ничего не делает?
    #
    #    Это СПЕЦИАЛЬНЫЙ аргумент для T^6 = K_3 ⊔ K_3,
    #    а не общий программный алгоритм Construction 1.
    #
    #    Стадия 1:
    #    две непересекающиеся треугольные компоненты на сфере имеют
    #    общую внешнюю область с двумя компонентами границы; каждая
    #    вершина принадлежит границе не более чем одной грани,
    #    большей треугольника. Поэтому 2-facet- и 3-facet-вершин нет.
    #
    #    Стадия 2:
    #    каждая вершина имеет степень 2 внутри T^6, а разрезаются
    #    вершины степени 5. Следовательно, стадия 2 также пуста.
    # ------------------------------------------------------------
    construction_stage_1_empty_by_special_case = each_component_is_k3

    construction_stage_2_empty = all(
        induced_degree(adj, u, comp) != 5
        for comp in t6_components
        for u in comp
    )

    require(construction_stage_1_empty_by_special_case,
            "Не удалось сертифицировать пустоту первой стадии Construction 1.")
    require(construction_stage_2_empty,
            "В T^6 обнаружена вершина степени 5.")

    # Поэтому cut-partition совпадает с исходными компонентами:
    # Cut(T^6) = {K_3, K_3}, то есть состоит из двух 1-треугольников.

    # ------------------------------------------------------------
    # 6. Проверяем явный gSW-путь.
    #
    #    В статье вершины обычно нумеруются с 1.
    #    В Python мы используем нумерацию с нуля, поэтому
    #    (4,1,5,6,14,15) превращается в:
    #    (3,0,4,5,13,14).
    # ------------------------------------------------------------
    path = [3, 0, 4, 5, 13, 14]
    require(is_gsw_path(adj, path), "Явный путь не удовлетворяет Definition 1.")

    # ------------------------------------------------------------
    # 7. Печатаем читаемый сертификат.
    # ------------------------------------------------------------
    print("=== C32: точный комбинаторный сертификат ===")
    print("graph6:", G6)
    print("Порядок |V|:", n)
    print("Размер |E|:", edges)
    print("Распределение степеней:", dict(sorted(Counter(degree).items())))
    print("Сферический сертификат:", sphere)

    print(
        "Вершины степени 5 (нумерация с 1):",
        [v + 1 for v in sorted(v5)],
    )
    print(
        "Вершины степени 6 (нумерация с 1):",
        [v + 1 for v in sorted(v6)],
    )
    print(
        "Компоненты T^6 (нумерация с 1):",
        [[v + 1 for v in sorted(comp)] for comp in t6_components],
    )
    print(
        "Рёбра T^6 (нумерация с 1):",
        [(u + 1, v + 1) for u, v in sorted(t6_edges)],
    )
    print("T^6 ≅ K3 ⊔ K3:", each_component_is_k3)
    print("Обе компоненты являются 1-треугольниками:", each_component_is_t1_triangle)
    print("Стадия 1 Construction 1 пуста:", construction_stage_1_empty_by_special_case)
    print("Стадия 2 Construction 1 пуста:", construction_stage_2_empty)

    print("gSW-путь (нумерация с 1):", [v + 1 for v in path])
    print("Степени вершин пути:", [degree[v] for v in path])
    print("Definition 1 выполнено:", is_gsw_path(adj, path))

    print()
    print("ИТОГ:")
    print("- cut-partition состоит из двух 1-треугольников;")
    print("- граф при этом содержит gSW-путь;")
    print("- все машинно-проверяемые сертификаты: PASS.")


if __name__ == "__main__":
    main()
