# Radiative dark-matter nuclear-currents storage audit — 2026-11-24

**Scope:** read-only current-tree audit and adoption of central PR #42's three-tier scientific data storage policy.  
**Repository:** `lm846802297-bit/radiative-dark-matter-nuclear-currents`  
**Audited branch:** `main`

## 1. Repository baseline

- GitHub repository metadata size: **3,836 KiB**.
- Complete recursive main tree: **102 blobs**, approximately **3.556 MiB** summed uncompressed file bytes.
- Recursive tree response was complete (`truncated=false`).
- No file exceeds 1 MiB.
- No root `.gitattributes` or `.lfsconfig` was found.
- The tree contains multiple SHA-256 manifests and a native AFH run manifest.

GitHub metadata size and current-tree uncompressed bytes are different measures; the latter is not a historical Git-object-size breakdown.

## 2. Native scientific data already archived

The current repository contains native AFH xenon response data, including:

- Eight-response native CSVs for Xe-128, Xe-129, Xe-130, Xe-131, Xe-132, Xe-134, and Xe-136.
- Native Xe-129/Xe-131 spatial multipole CSV.
- A source manifest, data SHA-256 manifests, and a native AFH run manifest.

Largest current files:

| Path | Size |
|---|---:|
| `nuclear/afh/data/native_v6_2026-11-22/native/spatial_multipoles/afh_v6_xe129_xe131_spatial_multipoles_native.csv` | 680,119 bytes |
| `nuclear/afh/data/native_v6_2026-11-22/native/afh_v6_xe131_eight_response_native.csv` | 538,596 bytes |
| `nuclear/afh/data/native_v6_2026-11-22/native/afh_v6_xe129_eight_response_native.csv` | 493,415 bytes |
| `schemas/research-result-record-v1.0.1.schema.json` | 30,073 bytes |

These are already-native scientific files with provenance and checksums. Keep them unchanged in ordinary Git.

## 3. Storage-tier decision

### Ordinary Git — retain current native dataset

Ordinary Git is appropriate for the current AFH native CSVs, run manifest, source manifests, SHA-256 lists, theory notes, schemas, scripts, and validation reports. The largest file is approximately 0.65 MiB.

### Git LFS — not enabled

No current file approaches GitHub's ordinary-Git per-file warning threshold. No LFS configuration was found. Do not enable LFS for the current AFH dataset.

### External object storage — not configured

No current dataset requires external object storage based on the audited snapshot. For future bulk detector outputs or large raw scans, evaluate durable versioned object storage using immutable object/version IDs and repository-side manifests with exact URI, size, source/run ID, SHA-256, retrieval instructions, and retention responsibility.

## 4. Integrity and migration boundary

- Existing AFH native files and historical manifests remain in place.
- No files were modified, deleted, recompressed, or moved.
- No Git history was rewritten.
- No nuclear calculation was run or rerun.
- Existing checksum manifests were inventoried but not all native-file SHA-256 values were independently recomputed in this read-only audit.
- No LFS or external storage was configured.

## 5. Standards adoption

The canonical lock is `docs/standards-lock.yaml`. It is updated to central commit `14155487afc075d301fa7f6cf452c7bff137e869`, which includes central PR #42's three-tier storage policy. Project-specific nuclear-current channels, AFH response conventions, chiral-EFT extensions, and physical validation gates remain local.

## 6. Conclusion

The repository already stores its AFH native nuclear data in ordinary Git with run and SHA-256 manifests. Its current size is small. **Keep the existing native data unchanged; no LFS or external-storage migration is warranted.** Apply the three-tier policy prospectively when a concrete larger dataset is produced.
