# HypoShrink 1.0.

## Companion construction

- Replaced the compact-root construction based on `crossprod(H)` and LDL with a direct SVD of `H`.
- Added stable lower- and upper-trapezoidal companion construction using orthogonal Givens rotations.
- Added `trapez = "lower"`, `"upper"`, or `"none"`; retained `utrapez` for backward compatibility.
- Transform non-zero right-hand sides directly through the compact SVD and the same orthogonal rotations.
- Added an explicit check that `H %*% theta = y` has a non-empty solution set.
- Removed the `fastmatrix` dependency.

## New features

- Added `checkCompanion()` to assess the numerical validity and minimality of companion representations.

## Corrections

- Corrected the `ATS_f` formula to match the documented adjusted ATS definition.
- Corrected the scaled-equivalence logic in `CompareHypothesis()`, including the right-hand-side norm condition and numerical comparisons of small values.
- Corrected lower/upper trapezoidal terminology for `centeringCompanion()`.
- Harmonized function names and examples in README/manuscript files.
- Harmonized the documented software license with `DESCRIPTION` (GPL-3).
- Strengthened input validation and tests.
