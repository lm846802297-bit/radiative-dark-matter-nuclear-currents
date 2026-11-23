# Original-paper VV recovery benchmark

The machine-readable definition in
\`paper_fig3_vector_natural_xenon.json\` freezes the natural-xenon VV
benchmark used in Fig. 3 of the original \`3body_for_HDM.pdf\`.

## Fixed parameters

\[
m_{\chi_1}=1\ {\rm TeV},\quad m_{A'}=5\ {\rm GeV},\quad
m_s=300\ {\rm keV},\quad \Delta m=10\ {\rm MeV},
\]
\[
g_1=0.1,\quad g_2=1.725,\quad g_{A'N}=g_{VN}=0.04857.
\]

The coupling mapping uses \(g_{VN}=3g_{Vq}\), with \(g_{Vq}=0.01619\).
Under this convention,
\((g_1g_2g_{VN})^2=7.0196329805625\times10^{-5}\).

The plotted observable is \(dN/dE_R\) per keV for a 2.84 tonne-year
natural-xenon exposure. A faithful event-spectrum reproduction must include
the paper's LZ nuclear-recoil efficiency; exposure alone does not determine
that efficiency.

The endothermic \(2\to2\) comparison uses \(\Delta_{\rm iDM}=300\) keV
and is normalized to the same integrated rate over the plotted interval.
The loop-induced elastic curve retains its predicted relative normalization
through the loop function \(F_3\).

## Scope discipline

This is the primary recovery target. The later P4 natural-xenon spectrum
in the source repository used a different secondary parameter point and is
a leading-static one-body candidate without the paper's LZ efficiency or
2.84 tonne-year exposure. It is not a reproduction of Fig. 3.

The benchmark definition is frozen. Independent numerical reproduction
in this new repository remains pending; this migration claims no new
event rate.
