# Troubleshooting

## “gfortran not found” / backend cannot be built

`vprofile()` (ScaleFree) and `hermite.hermite()` (fitvp) can auto-compile backends if needed,
but require `gfortran`.

Install it, or point to a precompiled executable:

```bash
export SCALEFREE_EXE=/path/to/scalefree.e
export SCALEFREE_FITVP_EXE=/path/to/fitvp.e
```

## “Fortran returned success but no structured output was detected”

This typically means:
- the backend executable is not the expected modified version that prints structured `# kind=...` blocks, or
- VP output was not requested (try `usevp=True`), or
- you are using a custom build that writes to an output file only.

Fixes:
- Use the packaged source compilation route, or
- Pass `output_path="out.txt"` to force file parsing.

## Unexpected NaNs / infinities in Gauss–Hermite moments

This can occur for extreme parameter combinations, or if VP reconstruction becomes numerically unstable.

Mitigations:
- Try a different `algorithm` (1 vs 3).
- Reduce extremes (e.g., move `beta` closer to 0, increase `q`).
- Increase numerical robustness (change `integration`, adjust `ngl_or_eps`).
