# Validation and test strategy

Tests are layered:

1. Algebra: gamma-matrix identities, spin sums, metric signs, exact kinematic identities.
2. Amplitudes: both scalar-emission diagrams, coherent interference, direct spinor-amplitude evaluation.
3. Tensors: direct amplitude square versus integrated three-body tensor; Hermiticity and positivity where expected.
4. Matching: relativistic current versus declared NR expansion; vector Ward or axial PCAC checks.
5. Nuclear: AFH multipole reconstruction, normalization, signed interference, positivity, Cauchy-Schwarz bounds.
6. Numerics: quadrature, phase-space boundary, recoil-grid, precision, and interpolation refinements.
7. Independent validation: distinct implementation or independent source data at selected points.
8. Regression: frozen VV benchmark remains unchanged when VA/AA/SS modules are added.

Each report states tested scope, source/input hashes, tolerances, and actual result.
