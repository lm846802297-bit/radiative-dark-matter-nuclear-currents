# Failed-trial and out-of-scope archive exclusions

## Superseded C1 `nw8-fail` trial

The source archive records this run as **FAIL_CONVERGENCE** for the superseded production quadrature `(8,6,4)`. It is not a production input and is not used by C2. The accepted C1 archive is `output/stage-c1-full-lorentz-grid/`, with production order `(14,6,4)` and comparison order `(16,14,12)`; see `C1_PASS_REPORT.md`.

The following five provenance records are retained in the target archive:

- `output/stage-c1-full-lorentz-grid-nw8-fail/C1_FAILURE_REPORT.md`
- `output/stage-c1-full-lorentz-grid-nw8-fail/SHA256SUMS.txt`
- `output/stage-c1-full-lorentz-grid-nw8-fail/c1_run.log`
- `output/stage-c1-full-lorentz-grid-nw8-fail/convergence_report.json`
- `output/stage-c1-full-lorentz-grid-nw8-fail/run_manifest.json`

The following three failed-trial data payloads are deliberately excluded. Their exact source identities and SHA-256 values remain recorded in the retained source `SHA256SUMS.txt` and in this note:

| Excluded source file | Size (bytes) | Source Git blob SHA | SHA-256 |
|---|---:|---|---|
| `output/stage-c1-full-lorentz-grid-nw8-fail/full_lorentz_tensor_grid.npz` | 16,067,801 | `c813dd4cf222530ed2adf336212f82fba6735556` | `6f87674f449680ed7170f3c2b602f43902135a6d0bd19cfa6895b11c56ec2121` |
| `output/stage-c1-full-lorentz-grid-nw8-fail/point_coverage.csv` | 83,981 | `2a63032d4be06d230085b81c7554ea258d5916aa` | `cce0dd1630596185d3ddb2262350d582e141b199499092dbf230f4d616227638` |
| `output/stage-c1-full-lorentz-grid-nw8-fail/tensor_index.csv` | 7,417,828 | `1bad061eaa07bc804baa330a81fe7d374f01bc5d` | `081874056b005a30af501a4d01c6073fe3685f845227764bb5e502bc49747cdb` |

The exclusion is deliberate: failed tensor/index payloads are not needed to reproduce or consume the accepted C1/C2 production chain. The failure report and provenance are retained so the failed convergence attempt remains auditable.

## Separate `pilot_248` package

The source subtree `pilot_248/` contains 58 files for a separate 248-keV pilot/re-coupling study. It is not part of the selected Stage A/B/C1/C2 production chain and was therefore excluded as a separate package, not classified as a failed production run. All 58 files remain at the pinned source commit.

## Counts

- Source P4 subtree: 139 files.
- Imported production package: 78 files.
- Excluded separate `pilot_248` package: 58 files.
- Excluded failed C1 data payloads: 3 files.
- Failed-trial provenance records retained: 5 files.
- Nuclear calculations rerun during migration: **0**.
