# Rank-88 extremal even unimodular lattices: computations

This repository distributes the calculations supporting the constructions of two rank-88 lattices, `L` and `L0`. It contains executable source, exact inputs, verification programs, and small certificates. The large exhaustive certificates are packaged separately and indexed by `full_certificate_index.jsonl`.

## Quick verification

With Python 3.9 or later on Linux or macOS, run:

```sh
python3 scripts/quick_verify.py
```

This checks distribution hashes, input identities, exact formula calculations, positive/even/unimodular Gram matrices, the two supplied norm-8 witness sets, and their independent perfection rank certificates.

## Full reproduction

Install PARI/GP, a C++17 compiler, GMP, and GAP with SmallGrp and GAPDoc. Then run:

```sh
python3 scripts/reproduce_all.py --output /absolute/path/to/new-rank88-run
```

The output directory must not exist. See [the reproducibility guide](docs/reproducibility.md) for stage commands, dependencies, and the difference between exhaustive searches and witness checks. Allow several gigabytes of free disk space. Recorded timings in `environment.json` are examples, not guarantees.

## Contents

- `data/`: exact Gram matrices, order-115 actions, and factor inputs.
- `certificates/`: independently checkable perfection certificates for both lattices.
- `compute/`: source and inputs for the `L`, `J0`, orbit-table, class-number, unit-group, and factor-group computations. Legacy computational directory names are retained so the drivers continue to work.
- `scripts/`: input reconstruction, quick checks, and full reproduction drivers.
- `verification/`: compact recorded verification results.
- `theorem_code_map.*` and `computation_dependency_map.*`: computational claims, source, input, evidence, expected output, and prerequisites.
- `full_certificate_index.jsonl`: every member of the separate certificate archive with size and SHA-256.
- `zenodo-large-file-manifest.json`: archive identity and checksum.

A fresh full rerun does not require the large archive. All code and data needed for the source-only rerun are included here.

Code licensing is in `LICENSE`. Rights in third-party references are described in `NOTICE.md`. [日本語の説明](README_JA.md).
