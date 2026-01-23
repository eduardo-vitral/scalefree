# scalefree tests (updated for intrinsic VP)

This test suite is designed to regression-test the structured output parsing for the Fortran backend, including the new intrinsic VP diagnostics.

## What is covered

For each algorithm in `(1, 2, 3)` we run:

- `kinematics="intrinsic"` with `average=False`  → `intrinsic_point` (+ `vp_intrinsic`, `vp_table_intrinsic`)
- `kinematics="intrinsic"` with `average=True`   → `intrinsic_shell_average` (+ beta, + `vp_intrinsic`, `vp_table_intrinsic`)
- `kinematics="projected"` with `average=False`  → `projected_point` (+ `vp`, `vp_table`)
- `kinematics="projected"` with `average=True`   → `projected_circle_average` (+ `vp`, `vp_table`)

Optionally (off by default) we also run `kinematics="both"` for `average=False/True`.

## Golden reference files

Reference outputs live in `tests/data/*_ref.txt`.

Generate/refresh them locally (after changing Fortran or the parser) with:

```bash
python tests/make_vprofile_refs.py
```

To also generate the optional `both` references:

```bash
python tests/make_vprofile_refs.py --include-both
```

## CI behaviour

- If reference files are missing **locally**, tests will skip and print the generation command.
- If reference files are missing in **CI** (`CI` env var set), tests will **fail**.
