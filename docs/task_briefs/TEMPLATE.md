# Research Task Brief — [short task title]

简体中文版：[research-task-brief.zh-CN.md](research-task-brief.zh-CN.md)。

> Complete every field before substantive work. If a field does not apply, write N/A and explain why. Do not leave ambiguous blanks.

**Task ID:** YYYYMMDD-short-kebab-case  
**Status:** draft | assigned | in-progress | blocked | completed  
**Prepared on:** YYYY-MM-DD  
**Prepared by:** [person/model/agent]  
**Assigned to:** [person/model/agent/team]  
**Research topic:** [topic or cross-topic]  
**Baseline repository:** [owner/repository]  
**Baseline branch and commit:** [branch] @ [full SHA]  
**Related issue / PR / task:** [URL or N/A]

---

## 1. Research question

**Primary question**

[One precise, answerable question. Define the physical/computational object and requested parameter or kinematic range.]

**Scope**
- Included:
- Excluded:
- Intended claim this task may support:

## 2. Primary sources, files, and repositories

List primary sources first. Include exact versions and locations within each source.

| Type | Exact source and version | Relevant pages / equations / paths / symbols | Provenance, URL, or SHA-256 |
|---|---|---|---|
| Primary paper/source | | | |
| Secondary/reference source | | | |
| Repository and code | | | |
| Supplied file / dataset | | | |
| Other | | | |

Unavailable sources and their impact:

## 3. Git baseline and versions

- Repository:
- Baseline branch:
- Baseline commit SHA:
- Relevant code path(s) and commit(s):
- Input / dataset / release version:
- Submodule or external package versions:
- Working-tree state (clean or list local changes):
- Exact revision against which results must be compared:

## 3A. Central standards provenance

- Central repository: lm846802297-bit/high-energy-theory-research
- Exact central source commit (40-character SHA):
- Standards lock path:
- Standards lock SHA-256 (computed from actual file bytes):
- Applied profile IDs (select only task-relevant profiles):
- Project-local extension paths and scope:
- Task brief SHA-256 (computed after the brief is finalized):
- Machine-readable task provenance record path:

Do not infer applied profiles from the lock alone. Record only profiles actually consulted and used for this task. The task provenance record must bind the exact lock and brief bytes by SHA-256.

## 4. Conventions

State definitions explicitly. Do not rely on labels such as “standard sign” or “usual normalization” without defining them.

- Metric signature and four-vector ordering:
- Units and conversion constants:
- Momentum/coordinate definitions and directions:
- Spatial magnitudes versus Lorentz invariants:
- Fourier / phase conventions:
- Algebra, matrix, spinor, polarization, or trace conventions:
- Interaction, propagator, boundary-condition, and regulator signs:
- State/current/observable normalization:
- Coupling and internal-symmetry normalization:
- Averages, sums, and symmetry factors:
- Reference frames, axes, angles, and energy definitions:
- Other conventions:
- Unresolved conventions:

## 5. Assumptions and approximations

| Assumption / approximation | Reason and source | Validity range | Possible impact / omitted contribution |
|---|---|---|---|
| | | | |

Also state the model/content, degrees of freedom, benchmark parameters, truncation order, and deliberately omitted effects.

## 6. Known validated results

Only list results whose evidence and scope can be identified. Unsupported prior conclusions remain hypotheses.

| Result / claim | Exact scope and benchmark | Evidence, code/input revision | Status and tolerance |
|---|---|---|---|
| | | | |

Relevant historical, provisional, diagnostic, or superseded results:

## 7. Required independent derivation

Specify the exact derivations required and their starting sources.

1. [Equation/result to derive from the primary source or model]
2. [Intermediate identities, signs, factors, or approximations to expose]
3. [Independent comparison target; do not use it as a substitute for derivation]

The derivation must show intermediate steps, assumptions, conventions, normalization factors, and where each approximation is applied.

## 8. Required numerical and technical checks

Set methods and acceptance criteria before interpreting outputs. Mark irrelevant checks N/A with a reason.

| Check | Method / independent implementation | Acceptance criterion / tolerance | Result: PASS / FAIL / BLOCKED / NOT RUN / N/A |
|---|---|---|---|
| Dimensional and unit consistency | | | |
| Analytic limits / thresholds / endpoints | | | |
| Conservation and exact identities | | | |
| Signs, normalization, averages, symmetry factors | | | |
| Independent code or primary-source comparison | | | |
| Precision / quadrature / grid / interpolation convergence | | | |
| Regression against frozen benchmark | | | |
| Uncertainty and assumption sensitivity | | | |
| Output integrity and checksums | | | |

Record actual residuals, tolerances, commands, logs, and exact code/input revisions. A check not run must remain NOT RUN.

## 9. Expected outputs

| Deliverable | Required path / format | Required contents and metadata | Acceptance condition |
|---|---|---|---|
| Independent derivation | | | |
| Code / implementation | | | |
| Tests / validation | | | |
| Machine-readable results | | | |
| Figures / report | | | |
| Run manifest / provenance | | | |
| PR / archive record | | | |

**Reproduction commands and environment:**

**Required conclusion level:** exploratory / diagnostic / reproducible / validated / other.

## 10. Unresolved questions

| Open question or discrepancy | Competing sources / conventions | Expected impact | Blocking? | Next discriminating check |
|---|---|---|---|---|
| | | | | |

---

## Execution instructions

- Derive independently from the supplied primary sources; do not assume earlier model conclusions are correct.
- Do not silently change metric, coordinate, sign, unit, or normalization conventions.
- If sources disagree, preserve the alternatives and quantify their consequences before adopting one.
- Separate mathematical derivation, code correctness, numerical convergence, and physical validity.
- State what was not run, what failed, and what remains blocked.
- CI is optional unless the task explicitly requests a particular workflow.

## Return record

- Answer and scope actually resolved:
- Independent derivation / key equations:
- Conventions confirmed or changed (with reasons):
- Assumptions and validity range:
- Numerical/technical checks, tolerances, and residuals:
- Unresolved discrepancies and next checks:
- Output paths, checksums, reproduction commands:
- Exact completion commit SHA / PR:
- Final status and the precise claim it supports:

## Change log

| Date | Revision | Change and reason | Affected conclusions / outputs |
|---|---|---|---|
| | | | |
