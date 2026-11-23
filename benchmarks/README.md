# Benchmarks

Each benchmark is a frozen, machine-readable definition of a physical point
or scan. Specify channel, dark-sector masses/couplings, scalar mass/emission
vertex, mediator and width prescription, isotope/current model, \(E_R\),
speed/halo model, integration bounds, \(q_A,q_N,t_A,t_N,y_A\),
form-factor/response normalization, precision, quadrature, interpolation,
and tolerances.

## Primary recovery target

\`vv_recovery/paper_fig3_vector_natural_xenon.json\` freezes the original
paper's Fig. 3 natural-xenon VV benchmark. It includes the stated dark-sector
parameters, coupling convention, 2.84 tonne-year exposure, and the distinct
normalization prescriptions for the comparison curves.

The first numerical task is an independent reproduction of this benchmark.
The imported source and Native AFH data are inputs, not a completed Fig. 3
reproduction. Add VA only after the VV recovery gate passes. Keep AA and SS
separate rather than reusing the VV kernel by changing a response label.
