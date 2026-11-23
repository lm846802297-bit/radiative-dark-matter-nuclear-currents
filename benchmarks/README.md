# Benchmarks

Each benchmark is a frozen, machine-readable definition of a physical point or scan. Specify channel, dark-sector masses/couplings, scalar mass/emission vertex, mediator and width prescription, isotope/current model, \(E_R\), speed/halo model, integration bounds, \(q_A,q_N,t_A,t_N,y_A\), form-factor/response normalization, precision, quadrature, interpolation, and tolerances.

The first benchmark recovers the original VV model. Add VA only after the VV recovery gate passes. Keep AA and SS separate rather than reusing the VV kernel.

## Primary paper recovery target

The original paper's Fig. 3 natural-xenon VV benchmark is frozen in
\`vv_recovery/paper_fig3_vector_natural_xenon.json\`. The imported source
snapshot and Native AFH data are inputs; the independent Fig. 3 reproduction
is still pending. See \`docs/VV_BASELINE_MIGRATION.md\`.
