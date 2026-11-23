# Raw Native AFH zero-recoil normalization audit

This checks immutable AFH v6.0 native CSVs already imported into the VV
baseline. It does not alter the old data or claim SN100PN physics from it.

The `W_M00_AFH` CSV column is the direct `ResponseNuclear[...,0,0]`
output; it is not the reduced spin-normalized quantity. The inherited
[conversion](https://github.com/lm846802297-bit/radiative-dark-matter-nuclear-currents/blob/5204f1a78b32b6ff5b423cc5d278d781e88add5e/nuclear/afh/data/native_v6_2026-11-22/source/README.md)
is `R_raw = W_native * 4*pi/(2*J+1)`. Only **then** does the equal
physical proton/neutron convention multiply by four: `R_equal = 4*R_raw`.
These factors are distinct. The ideal zero-momentum count limits are
`R_raw=A²/4` and `R_equal=A²`, not `W_native=A²/4`.

The independent audit hashes all seven native CSVs against the pinned
import manifest, rejects ambiguous zero rows, and measures their actual
normalizations without an additional count rescaling.

| Isotope | W_native(0) | R_equal(0) | ideal A² |
|---|---:|---:|---:|
| Xe-129 | 662.05131246 | 16639.16431618 | 16641 |
| Xe-131 | 1365.51629340 | 17159.58382286 | 17161 |

Across seven isotopes the largest fractional difference from the count
limit is 2.068e-4 (Xe-134). The replay threshold 2.5e-4 is a regression
check for this pinned dataset, not a physical uncertainty or permission
to rescale the data. The audit does not attribute these small discrepancies
to a particular cause. A W(0)-normalized shape is a separate observable
and must be labeled as such. See the accompanying JSON for all seven
isotopes and hashes.

Replay from a clone of the research branch:
```sh
audit=nuclear/sn100pn/xe129-xe131-rank0/native_afh_audit
python "$audit/audit_native_afh_charge_normalization.py" \
  nuclear/afh/data/native_v6_2026-11-22 /tmp/native-afh-zero-q-replay.json
cmp /tmp/native-afh-zero-q-replay.json "$audit/native-afh-zero-q-audit.json"
```
