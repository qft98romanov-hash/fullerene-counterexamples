# Certified Counterexamples to Three Conjectures on Dual Fullerene Graphs

This repository contains the reproducibility materials for the manuscript
**“Certified Counterexamples to Three Conjectures on Dual Fullerene Graphs.”**

The package verifies three results:

1. a dual `C32` fullerene graph whose hexagonal subgraph is `K3 ⊔ K3`, whose cut-partition consists of two non-degenerate `1`-triangles, and which nevertheless contains a generalized Stone–Wales path under the published definition;
2. a pair of nonisomorphic dual `C44` fullerene graphs for which the signless Laplacian matrices `A + D` are exactly cospectral;
3. a pair of nonisomorphic dual `C34` fullerene graphs for which exact rational interval arithmetic certifies a collision of the `(α,β)`-character on the ray `β = α/2`.

All proof-bearing comparisons in the core package use integer or rational arithmetic. Floating-point values are used only for diagnostic output.

## Repository layout

```text
.
├── README.md
├── CITATION.cff
├── .gitignore
├── LICENSE_TODO.md
├── verification/
│   ├── data/
│   ├── core/
│   ├── validation_logs/
│   ├── verify_manifest_ru.py
│   ├── independent_stdlib_audit.py
│   ├── independent_fullerene_audit_ru.py
│   ├── property_tests_stdlib.py
│   └── manifest.sha256
└── formal/
    └── README.md
```

The directory `verification/` is an unchanged copy of the audited verification package. Its internal `manifest.sha256` therefore continues to identify the original proof-bearing files exactly.

## Requirements

The autonomous core checks require:

- Python 3.12;
- no network access;
- no third-party Python packages.

The optional independent audit additionally uses the versions listed in
`verification/requirements-independent.txt`.

## Quick verification

From the repository root:

```bash
cd verification
python verify_manifest_ru.py
python core/run_all_local_ru.py
```

The expected final lines are:

```text
ПРОВЕРКА МАНИФЕСТА: PASS
ЕДИНЫЙ ОФЛАЙН-ЗАПУСК: PASS
```

The second command runs each proof-oriented program in ordinary and optimized (`python -O`) modes and requires the corresponding outputs to agree byte for byte. Mutable outputs are written only to `verification/reproduced/`.

## Optional independent audit

After installing the pinned optional dependencies:

```bash
python -m pip install -r verification/requirements-independent.txt
python verification/independent_fullerene_audit_ru.py verification \
  --json verification/reproduced/independent_fullerene_audit_results.json
```

### Windows PowerShell

On Windows, enable UTF-8 mode before running the verification scripts:

```powershell
$env:PYTHONUTF8 = "1"
$env:PYTHONIOENCODING = "utf-8"

py -3.12 .\verify_manifest_ru.py
py -3.12 .\core\run_all_local_ru.py


## Exact machine-readable graph descriptions

The five graphs underlying the counterexamples are specified by the exact
ASCII `graph6` encodings stored in:

- `verification/data/C32.g6`;
- `verification/data/C44_pair.g6`;
- `verification/data/C34_pair.g6`.

Each nonempty line contains one complete encoding and decodes unambiguously to
an adjacency matrix with the vertex ordering used by the verification code.
The encodings are exact but not canonical under relabeling: another ordering
of the same abstract graph may produce a different `graph6` string. No
catalogue numbering is required to reconstruct or verify any graph.

## Scope and limitations

The package verifies the claims made in the manuscript but does not claim:

- a general implementation or choice-independence theorem for Construction 1;
- uniqueness or simplicity of the certified `C34` root;
- a character collision at the separately recommended point `(1/2,1/4)`;
- minimality or uniqueness of the displayed counterexamples, unless such claims are separately added and supported by reproducible census code.

## Formal verification

Independent Lean and Rocq source files were supplied during an external review. They are not included in this public-upload draft until permission, attribution, and licensing have been confirmed. See `formal/README.md`.

## Citation

Citation metadata are provided in `CITATION.cff`. The repository URL, release date, journal citation, and DOI should be added after the repository and article records have been created.

## License

No public license has yet been selected for this draft. Resolve `LICENSE_TODO.md` before making the repository public.
