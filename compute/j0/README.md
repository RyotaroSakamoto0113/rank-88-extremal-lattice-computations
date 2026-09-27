# J₀⊗E: exact nonexistence of rank-4 tensors of h=3

Result: no such tensor exists. The complete exclusion is
8159 → 1556 → 453 → 1 → 0. The last affine completion is 47 → 77 → 25 → 0.

See `report_ja.md` for the proof and its dependence on the previously proved
projective minima (6,16,27,48). Combined with the previous norm-4 witness,
the trace lattice of E⊗J₀ is even unimodular of rank 88 and minimum 8.

The algebraic identification of the projective minima of J and J₀ is proved
in Proposition 5.4 of the accompanying paper. It includes nonfree sublattices
and is independent of the computations establishing the values for J.

## Fresh reproduction

Requirements: Python 3.9+ standard library, PARI/GP (`gp` in PATH),
Clang C++17 (`clang++` in PATH), and GMP headers/libraries.
The supplied compiler configuration includes `/opt/homebrew/include` and
`/opt/homebrew/lib`; standard compiler paths are also searched. Adjust those
two include/library arguments in `run_initial.py` if GMP is elsewhere.

From this directory:

```sh
python3 reproduce.py
```

Optionally specify a new, nonexistent directory:

```sh
python3 reproduce.py --run-dir /absolute/path/to/new-audit
```

Every run copies only the source files and the original lattice definitions,
then regenerates the shell, charts, candidates, obstructions, embeddings and
affine spheres. Saved candidate lists, shells and exclusion results are not
used as computational input. No proof assertions are disabled.

The last message must be `PROOF_COMPLETE`. Check `verified: true` and
`rank4_h3_exists: false` in the fresh `FULL_RUN_SUMMARY.json`.

## Evidence included

`run/input/` identifies the original J and the normalized ideal twist J₀.
`run/output/rank4_audit/` contains source, exact intermediate data, logs,
independent verifications and complete candidate coverage.
`FRESH_RUN_SUMMARY.json` records the successful fresh reproduction.
`MANIFEST.json` lists SHA-256 hashes for the packaged files.
The fresh run completed in 43.60 seconds on the original machine, including
compilation and independent verifications, excluding reproof of the previously
established projective minima. All 23 compared evidence files agreed.

`norm4_certificate/verify.py` (standard-library Python) and
`norm4_certificate/verify.gp` independently verify the previous explicit h=4
witness from its cyclotomic polynomial coordinates. Run the GP verifier from
inside that directory.

For a saved-evidence coverage check only, run `python3 finalize.py`.
This does not rerun the computations; use `reproduce.py` for that.

The proof uses F-linear ambient isometries only. No unverified full
automorphism-group identification or anti-linear orbit identification is used.
