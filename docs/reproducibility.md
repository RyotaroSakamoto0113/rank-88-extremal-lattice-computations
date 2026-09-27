# Reproducibility guide

## Requirements and installation

Tested: macOS/arm64, Python 3.9.6, PARI/GP 2.17.4, Apple clang 21, GMP, GAP 4.14.0. Exact versions and stage timings are in `environment.json`. Python scripts use the standard library; no SageMath or Magma license is needed.

On Ubuntu install `python3 pari-gp g++ libgmp-dev gap gap-smallgrp gap-gapdoc`. On macOS install the command-line compiler tools, PARI/GP and GMP, and GAP with its packages. Set `CXX=g++` if a suitable `clang++` is unavailable. Pass `--gap /absolute/path/to/gap` when GAP is not on PATH. The drivers use POSIX process accounting, so use Linux/macOS or Linux under WSL.

Do not run Python with `-O` or `PYTHONOPTIMIZE`: the historical programs use assertions as mathematical checks. Never reuse a completed output directory. Every command below writes to a new directory outside the source tree.

## Levels of verification

1. `python3 scripts/quick_verify.py`: distribution hashes plus independent explicit certificates. No PARI/GAP/compiler. It verifies the two Gram matrices and rank certificates in exact integer/finite-field arithmetic.
2. `python3 scripts/reproduce_all.py --output /absolute/new/run`: all exhaustive minimum computations, field arithmetic, factor groups, L0 reconstruction, and certificate checks.

The full driver runs `scripts/verify_factor_claims.py` after the factor-group stage. It checks shell counts and stabilizers, Sylow-23 actions, natural conjugations, norm-seven exclusion, and a Steinitz ideal identity.

## Independent stage commands

```sh
python3 compute/L/REPRODUCE_FULL.py --output /absolute/new/L
python3 compute/j0/reproduce.py --run-dir /absolute/new/J0
python3 compute/groups/reproduce.py --output /absolute/new/groups --gap /path/to/gap
python3 scripts/verify_factor_claims.py --repository . --groups /absolute/new/groups --gap /path/to/gap --output /absolute/new/factor_claims
python3 scripts/build_L0.py --output /absolute/new/L0_input
python3 scripts/verify_split_certificate.py --input data/L/lattice.json --certificate certificates/L_perfection.json --output /absolute/new/check_L
python3 scripts/verify_split_certificate.py --input data/L0/lattice.json --certificate certificates/L0_perfection.json --output /absolute/new/check_L0
python3 scripts/verify_formulas.py
```

`compute/L/REPRODUCE_FULL.py --smoke` checks prerequisites, h=3, perfection, and group actions. The full L driver also regenerates the 78 branches with final cutoff 47 in each ideal class, and rechecks 324468 exact quotient forms.

J0's exclusion is a separate calculation with 1556 binary survivors and 453 source representatives. These counts are different from J's 1376 and 362. J0's positive-control norm-four target has eight integral embeddings. The positive control and separate rational affine enumeration guard against an incorrectly empty search.

## Perfection witness generation

Perfection is an existence assertion, so verifying the supplied 44 seeds suffices. The seed finder is optional and randomized: different runs may produce different vectors and may need more trials. The source of the finder and the exact split-rank producer are included as `scripts/find_seeds.gp` and `scripts/split_rank.cpp`; the preserved L source includes `workspace/output/perfect_reaudit/search_run.py` and `rank_run.py`. A complete portable wrapper is `scripts/find_perfection.py`:

```sh
python3 scripts/find_perfection.py --lattice L0 --trials 160 --seed 20260926 --output /absolute/new/L0_search
```

The resulting certificate is then checked with `scripts/verify_split_certificate.py`. The independent verifier imports no production search or rank code. It checks positive leading Gram minors by Bareiss elimination, order and isometry equations over the integers, the left eigenbasis equation, the pair partition, exact norm 8 for every seed, and 115 nonzero determinants modulo 461. Full modular rank implies full rational rank because an integer minor nonzero modulo a prime is nonzero over Q. Perfection additionally uses minimum 8 proved by the independent exclusion pipeline.

## Certificate archive

The large archive named in `zenodo-large-file-manifest.json` is distributed separately. Check its archive SHA-256, extract only into an empty directory, and check its individual entries against `full_certificate_index.jsonl`. `scripts/verify_release_archive.py` verifies the compressed archive and all member hashes without extracting. `external:` paths in the claim map refer to archive members.

Some saved successful outputs contain original absolute paths; rerun drivers rebuild working paths dynamically. A complete source-only rerun requires no historical project, `/mnt/data`, Desktop folder, or cloud account. Raw outputs may differ in timings, PIDs, temporary names, GAP generator ordering, and equivalent basis choices. Compare mathematical invariants and exact lattice bindings rather than raw output hashes across fresh runs.

The exact witnesses for the projective minima of `E` are checked by `compute/L/src/verify_additional.py`.

The MIT license covers original code. Third-party references are described in `NOTICE.md`.
