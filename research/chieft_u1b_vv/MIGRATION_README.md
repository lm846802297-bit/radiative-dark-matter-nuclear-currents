# U(1)_B VV chiEFT audit migration

This is a source-preserving migration of the literature/source-matching work that was developed as an exploratory extension in `LZ-3-body-for-hdm`.

## Origin

- Repository: https://github.com/lm846802297-bit/LZ-3-body-for-hdm
- Source branch: `research/chiral-vv-source-matching-20261124`
- Source PR: https://github.com/lm846802297-bit/LZ-3-body-for-hdm/pull/133
- Source PR head: `4f63f37d7ca339c505977367911beb64248c2bfb`
- Scope: literature and source-matching audit only; no numerical chiEFT current, NN+3N production, xenon many-body calculation, or full-current completion is claimed.

## Migration policy

The source documents and original task/validation records are preserved without changing their scientific claims. The provenance JSON files intentionally retain their original source-repository identities and paths as historical provenance. This README records their new archival location.

## Files

| Original source file | Original Git blob SHA | Destination |
|---|---|---|
| `research/three_body_full_current_nucleon_match/CHI_EFT_VV_SOURCE_MATCHING_AUDIT_2026-11-24.md` | `b06ae5923ca6f97e2ec8316ec788a01b71df19d8` | `research/chieft_u1b_vv/CHI_EFT_VV_SOURCE_MATCHING_AUDIT_2026-11-24.md` |
| `research/three_body_full_current_nucleon_match/FULL_CURRENT_CHI_EFT_INPUT_AUDIT.md` | `34f739178fa5832af5f4db546dab47b334107e4b` | `research/chieft_u1b_vv/FULL_CURRENT_CHI_EFT_INPUT_AUDIT.md` |
| `research/three_body_full_current_nucleon_match/FULL_CURRENT_MODEL_CLOSURE_AUDIT.md` | `4a34d5c20f714be206f1c923eb5607ae87b85cc1` | `research/chieft_u1b_vv/FULL_CURRENT_MODEL_CLOSURE_AUDIT.md` |
| `research/three_body_full_current_nucleon_match/task_briefs/20261124-chiral-vv-source-matching.md` | `a48bd6b8f1d672d45e0b98891dea38979351e6da` | `research/chieft_u1b_vv/task_brief.md` |
| `research/standards-provenance/tasks/20261124-chiral-vv-source-matching.json` | `67819e7b4d891808fb288f62f84784d2a46e81dd` | `research/chieft_u1b_vv/provenance/task_provenance.json` |
| `research/standards-provenance/reports/20261124-chiral-vv-source-matching.json` | `b69324e430b28e183c1a4866fb29e01fa830c00e` | `research/chieft_u1b_vv/provenance/validation_report.json` |

## Relationship to the LZ calculation

The validated one-body impulse baseline and the original-paper three-body recoil predictions remain in `LZ-3-body-for-hdm`. This archive does not change those results or their completion status. The chiEFT program is optional follow-up research in this general nuclear-current project.
