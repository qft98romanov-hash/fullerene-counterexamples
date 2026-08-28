#!/usr/bin/env python3
"""
Единый офлайн-запуск всех стандартных сертификатов.

Каждый файл запускается в обычном и оптимизированном режимах.
Выводы должны совпасть побайтово: это регрессионная проверка
отсутствия доказательных условий, исчезающих при ``python -O``.

Все изменяемые результаты пишутся только в ``reproduced/`` в корне
пакета. Файлы ``data/`` и ``manifest.sha256`` не изменяются.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
OUTPUTS = ROOT / "reproduced"

SCRIPTS = (
    "core/audit_c32_op3_ru.py",
    "core/verify_op4_ru.py",
    "core/verify_op4_compact_ru.py",
    "core/check_package_consistency_ru.py",
    "independent_stdlib_audit.py",
    "property_tests_stdlib.py",
)


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def run_script(relative: str, optimized: bool) -> subprocess.CompletedProcess[str]:
    path = ROOT / relative
    if not path.is_file():
        raise FileNotFoundError(path)
    command = [sys.executable]
    if optimized:
        command.append("-O")
    command.append(str(path))
    environment = dict(os.environ)
    environment.update({"PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0"})
    return subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def main() -> int:
    OUTPUTS.mkdir(exist_ok=True)
    summary: dict[str, object] = {
        "python": sys.version,
        "executable": Path(sys.executable).name,
        "package_root": ".",
        "scripts": [],
    }
    failed = False

    for relative in SCRIPTS:
        normal = run_script(relative, optimized=False)
        optimized = run_script(relative, optimized=True)
        stem = relative.replace("/", "__").removesuffix(".py")
        normal_path = OUTPUTS / f"{stem}.txt"
        optimized_path = OUTPUTS / f"{stem}.optimized.txt"
        normal_path.write_text(normal.stdout, encoding="utf-8")
        optimized_path.write_text(optimized.stdout, encoding="utf-8")
        outputs_equal = normal.stdout == optimized.stdout
        passed = normal.returncode == 0 and optimized.returncode == 0 and outputs_equal
        record = {
            "script": relative,
            "normal_returncode": normal.returncode,
            "optimized_returncode": optimized.returncode,
            "normal_output": str(normal_path.relative_to(ROOT)),
            "optimized_output": str(optimized_path.relative_to(ROOT)),
            "outputs_byte_identical": outputs_equal,
            "normal_output_sha256": digest(normal_path),
        }
        summary["scripts"].append(record)
        print(f"[{'PASS' if passed else 'FAIL'}] {relative} (normal/-O)")
        if not passed:
            failed = True
            if normal.returncode:
                print(normal.stdout)
            if optimized.returncode:
                print(optimized.stdout)
            if not outputs_equal:
                print("Выводы обычного и -O режимов различаются.")

    summary_path = OUTPUTS / "run_summary_ru.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if failed:
        print("ЕДИНЫЙ ОФЛАЙН-ЗАПУСК: FAIL")
        return 1
    print("ЕДИНЫЙ ОФЛАЙН-ЗАПУСК: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
