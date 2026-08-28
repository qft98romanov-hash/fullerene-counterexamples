#!/usr/bin/env python3
r"""
Проверить SHA256-манифест пакета без его перезаписи.

Формат каждой строки manifest.sha256:

    <64 hex-символа><два пробела><относительный путь>

Скрипт следует положить в корень распакованного пакета рядом с
manifest.sha256 и запустить:

    py -3 .\verify_manifest_ru.py
"""

from __future__ import annotations

import hashlib
import string
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "manifest.sha256"


def digest(path: Path) -> str:
    """Вычислить SHA256 файла блоками по 1 МиБ."""
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def main() -> int:
    if not MANIFEST.is_file():
        print("Файл manifest.sha256 не найден:", MANIFEST)
        return 1

    failed = False
    checked = 0
    listed: set[Path] = set()

    for line_number, raw_line in enumerate(
        MANIFEST.read_text("utf-8").splitlines(),
        start=1,
    ):
        line = raw_line.strip("\n")
        if not line.strip():
            continue

        try:
            expected, relative = line.split("  ", 1)
        except ValueError:
            print(f"Строка {line_number}: неверный формат")
            failed = True
            continue

        relative_path = Path(relative)
        if (
            len(expected) != 64
            or any(character not in string.hexdigits for character in expected)
        ):
            print(f"Строка {line_number}: неверный SHA-256")
            failed = True
            continue
        if relative_path.is_absolute() or ".." in relative_path.parts:
            print(f"Строка {line_number}: небезопасный путь {relative!r}")
            failed = True
            continue
        normalized = Path(*[part for part in relative_path.parts if part != "."])
        if normalized in listed:
            print(f"Строка {line_number}: повторный путь {relative!r}")
            failed = True
            continue
        listed.add(normalized)

        path = ROOT / normalized
        if not path.is_file():
            print(f"MISSING: {relative}")
            failed = True
            continue

        actual = digest(path)
        checked += 1

        if actual.lower() != expected.lower():
            print(f"FAILED:  {relative}")
            print(f"  ожидается: {expected}")
            print(f"  получено:  {actual}")
            failed = True
        else:
            print(f"OK:      {relative}")

    actual_files = {
        path.relative_to(ROOT)
        for path in ROOT.rglob("*")
        if path.is_file()
        and path != MANIFEST
        and "reproduced" not in path.relative_to(ROOT).parts
        and "__pycache__" not in path.relative_to(ROOT).parts
        and path.suffix not in {".pyc", ".pyo"}
    }
    for unexpected in sorted(actual_files - listed, key=lambda path: path.as_posix()):
        print(f"UNLISTED: {unexpected}")
        failed = True
    for stale in sorted(listed - actual_files, key=lambda path: path.as_posix()):
        print(f"STALE:    {stale}")
        failed = True

    print()
    print("Проверено файлов:", checked)

    if failed:
        print("ПРОВЕРКА МАНИФЕСТА: FAILED")
        return 1

    print("ПРОВЕРКА МАНИФЕСТА: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
