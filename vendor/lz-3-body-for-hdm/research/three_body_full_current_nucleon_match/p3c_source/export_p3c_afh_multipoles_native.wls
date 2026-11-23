(* Native AFH v6.0 rank-resolved amplitudes at the four P3-C recoil points. *)
scriptDir = DirectoryName[$InputFileName];
repoRoot = DirectoryName[DirectoryName[DirectoryName[DirectoryName[scriptDir]]]];
envPath = Environment["DMFormFactorPath"];
packagePath = If[StringQ[envPath] && envPath =!= "", envPath,
  FileNameJoin[{repoRoot, "v6", "dmformfactor-V6.m"}]];
If[!FileExistsQ[packagePath], Print["ERROR: AFH source not found: ", packagePath]; Exit[2]];
Get[packagePath];

outDir = Environment["P3C_NATIVE_OUTDIR"];
If[!StringQ[outDir] || outDir === "", outDir = FileNameJoin[{scriptDir, "output"}]];
If[!DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

isotopes = {129, 131};
recoils = {0.0002, 0.000267};
uGeV = 0.93149410242;
rows = {};
termRows = {};
maxImaginary = 0.;
rankMap = {};

Do[
  A = isotope; Z = 54;
  SetIsotope[Z, A, "default", "default"];
  bFM = Sqrt[41.467/(45.0 A^(-1.0/3.0) - 25.0 A^(-2.0/3.0))];
  bGeVInv = bFM * 5.0677;
  isoFactor = Sqrt[2]/Sqrt[2 dmformfactor`Private`TIso + 1];
  chargeRanks = Range[0, Round[2 dmformfactor`Private`JIso]];
  (* Avoid relying on the package's response-sum loop for the rank map. *)
  chargeRanks = Sort@DeleteDuplicates@Cases[
    Table[If[dmformfactor`Private`DensityMatrix[i][[3]] == 0,
      dmformfactor`Private`DensityMatrix[i][[2]], -1],
      {i, dmformfactor`Private`DMmatrixLength}], n_Integer /; 0 <= n <= 2 dmformfactor`Private`JIso];
  spatialRanks = Sort@DeleteDuplicates@Cases[
    Table[If[dmformfactor`Private`DensityMatrix[i][[3]] == 0 &&
        dmformfactor`Private`DensityMatrix[i][[2]] >= 1,
      dmformfactor`Private`DensityMatrix[i][[2]], 0],
      {i, dmformfactor`Private`DMmatrixLength}], n_Integer /; n > 0];
  If[chargeRanks =!= Range[0, Round[2 dmformfactor`Private`JIso]],
    Print["ERROR: incomplete charge rank map for Xe-", A, ": ", chargeRanks]; Exit[3]];
  AppendTo[rankMap, ToString[A] <> ":charge=" <> StringRiffle[ToString /@ chargeRanks, ","] <>
    ";spatial=" <> StringRiffle[ToString /@ spatialRanks, ","]];

  Do[
    ER = recoil; Q = Sqrt[2.0 A uGeV ER + ER^2]; yValue = (Q bGeVInv/2.0)^2;
    Do[
      Do[
        Module[{channel, op, value},
          {channel, op} = pair;
          value = Sum[
            If[dmformfactor`Private`DensityMatrix[i][[3]] == 0 &&
                dmformfactor`Private`DensityMatrix[i][[2]] == J,
              dmformfactor`Private`DensityMatrix[i][[1]] * isoFactor *
                op[yValue,
                  {dmformfactor`Private`DensityMatrix[i][[4]], dmformfactor`Private`DensityMatrix[i][[5]]/2},
                  {dmformfactor`Private`DensityMatrix[i][[6]], dmformfactor`Private`DensityMatrix[i][[7]]/2}, J], 0],
            {i, dmformfactor`Private`DMmatrixLength}];
          value = N[value, 30];
          If[!NumericQ[value] || !FreeQ[value, Indeterminate | ComplexInfinity | DirectedInfinity[___]],
            Print["ERROR: invalid ", channel, " amplitude for Xe-", A, " ER=", ER, " J=", J]; Exit[4]];
          maxImaginary = Max[maxImaginary, Abs[Im[value]]];
          AppendTo[rows, {A, ER, Q, yValue, J, channel, Re[value], Im[value]}]
        ],
        {pair, {{"M", dmformfactor`Private`MJ}, {"PhiPP", dmformfactor`Private`PhiPPJ}}}],
      {J, chargeRanks}];
    Do[
      Do[
        Module[{channel, op, value},
          {channel, op} = pair;
          value = Sum[
            If[dmformfactor`Private`DensityMatrix[i][[3]] == 0 &&
                dmformfactor`Private`DensityMatrix[i][[2]] == J,
              dmformfactor`Private`DensityMatrix[i][[1]] * isoFactor *
                op[yValue,
                  {dmformfactor`Private`DensityMatrix[i][[4]], dmformfactor`Private`DensityMatrix[i][[5]]/2},
                  {dmformfactor`Private`DensityMatrix[i][[6]], dmformfactor`Private`DensityMatrix[i][[7]]/2}, J], 0],
            {i, dmformfactor`Private`DMmatrixLength}];
          value = N[value, 30];
          If[!NumericQ[value] || !FreeQ[value, Indeterminate | ComplexInfinity | DirectedInfinity[___]],
            Print["ERROR: invalid ", channel, " amplitude for Xe-", A, " ER=", ER, " J=", J]; Exit[5]];
          maxImaginary = Max[maxImaginary, Abs[Im[value]]];
          AppendTo[rows, {A, ER, Q, yValue, J, channel, Re[value], Im[value]}]
        ],
        {pair, {{"Delta", dmformfactor`Private`DeltaJ}, {"SigmaPrime", dmformfactor`Private`SigmaPJ}}}],
      {J, spatialRanks}];
    Do[
      Do[
        If[dmformfactor`Private`DensityMatrix[i][[3]] == 0 &&
            dmformfactor`Private`DensityMatrix[i][[2]] == J,
          Module[{np, jp, n, j, coeff, mAmp, phiAmp, deltaAmp, sigmaAmp, sigmaMinus, sigmaPlus, sigmaRebuilt},
            np = dmformfactor`Private`DensityMatrix[i][[4]];
            jp = dmformfactor`Private`DensityMatrix[i][[5]]/2;
            n = dmformfactor`Private`DensityMatrix[i][[6]];
            j = dmformfactor`Private`DensityMatrix[i][[7]]/2;
            coeff = dmformfactor`Private`DensityMatrix[i][[1]];
            mAmp = dmformfactor`Private`MJ[yValue, {np, jp}, {n, j}, J];
            phiAmp = dmformfactor`Private`PhiPPJ[yValue, {np, jp}, {n, j}, J];
            deltaAmp = If[J > 0, dmformfactor`Private`DeltaJ[yValue, {np, jp}, {n, j}, J], 0];
            sigmaAmp = If[J > 0, dmformfactor`Private`SigmaPJ[yValue, {np, jp}, {n, j}, J], 0];
            sigmaMinus = If[J > 0,
              dmformfactor`Private`MJLSigma[yValue,
                {dmformfactor`Private`Nodal[np, jp], dmformfactor`Private`Lnumber[np, jp], jp},
                {dmformfactor`Private`Nodal[n, j], dmformfactor`Private`Lnumber[n, j], j}, J, J - 1], 0];
            sigmaPlus = If[J > 0,
              dmformfactor`Private`MJLSigma[yValue,
                {dmformfactor`Private`Nodal[np, jp], dmformfactor`Private`Lnumber[np, jp], jp},
                {dmformfactor`Private`Nodal[n, j], dmformfactor`Private`Lnumber[n, j], j}, J, J + 1], 0];
            sigmaRebuilt = If[J > 0, -Sqrt[J]/Sqrt[2 J + 1] sigmaPlus + Sqrt[J + 1]/Sqrt[2 J + 1] sigmaMinus, 0];
            AppendTo[termRows, {A, ER, Q, yValue, J, i, coeff, np, 2 jp, n, 2 j,
              Re[mAmp], Re[phiAmp], Re[deltaAmp], Re[sigmaAmp], Re[sigmaMinus], Re[sigmaPlus], Re[sigmaRebuilt]}]
          ]
        ],
        {i, dmformfactor`Private`DMmatrixLength}],
      {J, chargeRanks}];
  , {recoil, recoils}];
, {isotope, isotopes}];

If[maxImaginary > 2*10^-12, Print["ERROR: imaginary residual exceeds tolerance: ", maxImaginary]; Exit[6]];
csvPath = FileNameJoin[{outDir, "afh_v6_p3c_xe129_xe131_native_multipoles.csv"}];
Export[csvPath, Prepend[rows, {"A", "ER_GeV", "Q_A_GeV", "y", "J", "channel", "amplitude_real", "amplitude_imag"}], "CSV"];
termCsvPath = FileNameJoin[{outDir, "afh_v6_p3c_xe129_xe131_native_obdme_terms.csv"}];
Export[termCsvPath, Prepend[termRows, {"A", "ER_GeV", "Q_A_GeV", "y", "J", "density_index", "obdme_coefficient",
  "Nprime", "two_jprime", "N", "two_j", "M_element", "PhiPP_element", "Delta_element", "SigmaPrime_element",
  "SigmaPrime_LJminus1_element", "SigmaPrime_LJplus1_element", "SigmaPrime_rebuilt_element"}], "CSV"];
metadata = StringRiffle[{
  "source=AFH DMFormFactor v6.0 native Wolfram execution",
  "source_path=" <> packagePath,
  "source_sha256=" <> IntegerString[FileHash[packagePath, "SHA256"], 16, 64],
  "wolfram_version=" <> $Version,
  "isotopes=129,131",
  "recoils_GeV=0.0002,0.000267",
  "mass_scheme=historical benchmark mA=A*0.93149410242 GeV; Q_A=sqrt(2*mA*ER+ER^2)",
  "channels=raw native M, PhiPP, Delta, SigmaPrime reduced multipole OBDME sums",
  "isospin=density-matrix isospin-rank index T=0 rows only; AFH target TIso=|A/2-Z|; apply AFH factor Sqrt[2/(2 TIso+1)]",
  "normalization=raw AFH package convention; factor-four and C0 not applied",
  "rank_map=" <> StringRiffle[rankMap, ";"],
  "max_imaginary_residual=" <> ToString[N[maxImaginary, 16], InputForm],
  "rows=" <> ToString[Length[rows]]
  , "obdme_term_rows=" <> ToString[Length[termRows]]
}, "\n"];
Export[FileNameJoin[{outDir, "native_metadata.txt"}], metadata, "Text"];
Print["P3-C native AFH multipoles complete: rows=", Length[rows], "; max imaginary=", maxImaginary];
