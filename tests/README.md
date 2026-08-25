# Tests

Regression tests for the XML readers, using real anonymised device files as fixtures
(see `fixtures/README.md`).

## Running

From the repo root in MATLAB:

```matlab
results = runtests('tests');
table(results)
```

## What's covered

`testReadBPplusReaders.m` exercises `read_BPplus.m`'s dispatch to
`read_BPplusBPplus.m` / `read_BPplusCardioScope.m` against all six fixture
files, and specifically regression-tests the beta7 fixes:

- the duplicate-block bug that overwrote `ao.averagePulsePointsIndexes`
  with the wrong source field
- the NIBP mode guard that previously never set `metadata.mode`
- `ao.ed` resolving to `NaN` (not `-1` or `[]`) when unavailable
- the `patienet_id` -> `patient_id` typo
- `ao.lag` resolving correctly via the firmware-build lookup table

## What's not covered yet

`ai_v2.m`, `BPfitres_v1.m`, and `kreservoir_v15.m` have no tests yet -
these are the better candidates for synthetic-input unit tests, since
they're pure numerical functions with no file I/O. `bpp_Res2.m` itself
isn't tested (it's a script with GUI dialogs and side effects, not a
function) - splitting it into functions would make that practical.

## A note on these specific tests

These assertions were written by working through the raw XML files by
hand (see the beta7 fixes and the conversation that produced them) and
have not been executed against MATLAB. Run them and report back anything
that fails, rather than assuming the test file itself is correct.
