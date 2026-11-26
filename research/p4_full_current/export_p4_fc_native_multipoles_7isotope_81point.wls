(* Native AFH v6.0 rank-resolved charge/spatial multipoles for P4 full-current work. *)
fail[msg_, code_:1] := (Print["ERROR: ", msg]; Exit[code]);
env[name_] := Module[{s=Environment[name]}, If[StringQ[s] && s=!="",s,""]];
scriptDir=DirectoryName[$InputFileName];
repoRoot=ExpandFileName[FileNameJoin[{scriptDir,"..","..",".."}]];
packagePath=env["DMFormFactorPath"];
If[packagePath==="",packagePath=FileNameJoin[{repoRoot,"v6","dmformfactor-V6.m"}]];
massInputPath=env["P4_XENON_MASS_INPUT"];
If[massInputPath==="",massInputPath=FileNameJoin[{repoRoot,"research","afh_default_density_natural_xenon_halo","input","natural_xenon_isotope_weights.csv"}]];
outDir=env["P4_FC_NATIVE_OUTDIR"];
If[outDir==="",outDir=FileNameJoin[{scriptDir,"output","native-seven-isotope-81-point"}]];
expectedAFHSHA="dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372";
expectedMassSHA="fd1bb996e5f582c74c3c86d129f1d2dded510c65237d0a7b2ba08861f1d8c991";
If[!FileExistsQ[packagePath],fail["AFH source not found: "<>packagePath,2]];
If[!FileExistsQ[massInputPath],fail["mass input not found: "<>massInputPath,2]];
afhSHA=FileHash[packagePath,"SHA256","HexString"];
massSHA=FileHash[massInputPath,"SHA256","HexString"];
If[afhSHA=!=expectedAFHSHA,fail["AFH source SHA-256 mismatch: "<>afhSHA,3]];
If[massSHA=!=expectedMassSHA,fail["P4 mass input SHA-256 mismatch: "<>massSHA,3]];
If[!DirectoryQ[outDir],CreateDirectory[outDir,CreateIntermediateDirectories->True]];
Get[packagePath];

massData=Import[massInputPath,"CSV"];
massHeader=First[massData];
col[name_]:=Module[{p=FirstPosition[massHeader,name,Missing["NotFound"]]},
  If[MissingQ[p],fail["mass input missing column "<>name,4],First[p]]];
colA=col["A"]; colAtomic=col["atomic_mass_u"]; colAvailable=col["AFH_default_density_available"];
asNumber[x_]:=If[NumberQ[x],N[x,30],N[ToExpression[StringTrim[ToString[x]]],30]];
isTrue[x_]:=ToLowerCase[StringTrim[ToString[x]]]==="true";
uGeV=0.93149410242; electronMassU=0.000548579909065; Z=54;
massRows=Rest[massData];
massTable=Association[];
Do[
  If[isTrue[row[[colAvailable]]],
    a=Round[asNumber[row[[colA]]]];
    massTable[a]=(asNumber[row[[colAtomic]]]-Z*electronMassU)*uGeV
  ],
  {row,massRows}
];
isotopes={128,129,130,131,132,134,136};
If[Sort[Keys[massTable]]=!=isotopes,fail["unexpected AFH isotope mass coverage",5]];
recoilKeV=N[40.+9.5 Range[0,80],18];
If[Length[recoilKeV]!=81 || Last[recoilKeV]!=800.,fail["invalid 81-point recoil grid",5]];

responseNames={"M","SigmaDoublePrime","SigmaPrime","PhiPP","PhiTildePrime","Delta","MPhiPP","SigmaPrimeDelta"};
responseIndices=Range[8];
multipoleRows={}; responseRows={}; termRows={};
maxImaginary=0.; maxResponseImaginary=0.;
rankMap={}; c0Rows={}; totalPoints=0;

Do[
  A=isotope; mA=massTable[A];
  Print["Native AFH P4-FC export: Xe-",A];
  SetIsotope[Z,A,"default","default"];
  JA=N[dmformfactor`Private`JIso,20];
  Tiso=N[dmformfactor`Private`TIso,20];
  bFM=Sqrt[41.467/(45.0 A^(-1.0/3.0)-25.0 A^(-2.0/3.0))];
  bGeVInv=bFM*5.0677;
  isoFactor=Sqrt[2]/Sqrt[2 dmformfactor`Private`TIso+1];
  densityLength=dmformfactor`Private`DMmatrixLength;
  chargeRanks=Sort@DeleteDuplicates@Cases[
    Table[If[dmformfactor`Private`DensityMatrix[i][[3]]==0,
      dmformfactor`Private`DensityMatrix[i][[2]],-1],{i,densityLength}],
    n_Integer/;0<=n<=2 JA];
  spatialRanks=Sort@DeleteDuplicates@Cases[
    Table[If[dmformfactor`Private`DensityMatrix[i][[3]]==0 &&
      dmformfactor`Private`DensityMatrix[i][[2]]>=1,
      dmformfactor`Private`DensityMatrix[i][[2]],0],{i,densityLength}],
    n_Integer/;n>0];
  expectedChargeRanks=Range[0,Round[2 JA]];
  expectedSpatialRanks=Switch[A,129,{1},131,{1,2,3},_,{}];
  If[chargeRanks=!=expectedChargeRanks,fail["charge rank map mismatch Xe-"<>ToString[A]<>": "<>ToString[chargeRanks],6]];
  If[spatialRanks=!=expectedSpatialRanks,fail["spatial rank map mismatch Xe-"<>ToString[A]<>": "<>ToString[spatialRanks],6]];
  AppendTo[rankMap,ToString[A]<>":J_A="<>ToString[JA]<>":charge="<>StringRiffle[ToString/@chargeRanks,","]<>":spatial="<>StringRiffle[ToString/@spatialRanks,","]];
  m0=N[ResponseNuclear[y,1,0,0]/.y->0,40];
  c0=A^2/(4*(4 Pi/(2 JA+1))*m0);
  AppendTo[c0Rows,{A,JA,Tiso,mA,bGeVInv,c0}];

  Do[
    ER=energy*10^-6;
    Q=Sqrt[2*mA*ER+ER^2];
    yValue=(Q*bGeVInv/2)^2;
    responseValues=N[Table[ResponseNuclear[y,index,0,0],{index,responseIndices}]/.y->yValue,40];
    If[!VectorQ[responseValues,NumericQ] || !FreeQ[responseValues,Indeterminate|ComplexInfinity|DirectedInfinity[___]],
      fail["non-finite native response Xe-"<>ToString[A]<>" ER="<>ToString[energy],7]];
    maxResponseImaginary=Max[maxResponseImaginary,Max[Abs[Im[#]]&/@responseValues]];
    AppendTo[responseRows,Join[{A,energy,ER,-2*mA*ER,mA,Q,yValue,JA,c0},Re/@responseValues,Im/@responseValues]];

    Do[
      Module[{channel,op,value},
        {channel,op}=pair;
        value=Sum[
          If[dmformfactor`Private`DensityMatrix[i][[3]]==0 &&
            dmformfactor`Private`DensityMatrix[i][[2]]==J,
            dmformfactor`Private`DensityMatrix[i][[1]]*isoFactor*
              op[yValue,
                {dmformfactor`Private`DensityMatrix[i][[4]],dmformfactor`Private`DensityMatrix[i][[5]]/2},
                {dmformfactor`Private`DensityMatrix[i][[6]],dmformfactor`Private`DensityMatrix[i][[7]]/2},J],0],
          {i,densityLength}];
        value=N[value,40];
        If[!NumericQ[value] || !FreeQ[value,Indeterminate|ComplexInfinity|DirectedInfinity[___]],
          fail["invalid "<>channel<>" amplitude Xe-"<>ToString[A]<>" ER="<>ToString[energy]<>" J="<>ToString[J],8]];
        maxImaginary=Max[maxImaginary,Abs[Im[value]]];
        AppendTo[multipoleRows,{A,energy,ER,-2*mA*ER,mA,Q,yValue,JA,J,channel,Re[value],Im[value]}]
      ],
      {J,chargeRanks},{pair,{{"M",dmformfactor`Private`MJ},{"PhiPP",dmformfactor`Private`PhiPPJ}}}
    ];
    Do[
      Module[{channel,op,value},
        {channel,op}=pair;
        value=Sum[
          If[dmformfactor`Private`DensityMatrix[i][[3]]==0 &&
            dmformfactor`Private`DensityMatrix[i][[2]]==J,
            dmformfactor`Private`DensityMatrix[i][[1]]*isoFactor*
              op[yValue,
                {dmformfactor`Private`DensityMatrix[i][[4]],dmformfactor`Private`DensityMatrix[i][[5]]/2},
                {dmformfactor`Private`DensityMatrix[i][[6]],dmformfactor`Private`DensityMatrix[i][[7]]/2},J],0],
          {i,densityLength}];
        value=N[value,40];
        If[!NumericQ[value] || !FreeQ[value,Indeterminate|ComplexInfinity|DirectedInfinity[___]],
          fail["invalid "<>channel<>" amplitude Xe-"<>ToString[A]<>" ER="<>ToString[energy]<>" J="<>ToString[J],9]];
        maxImaginary=Max[maxImaginary,Abs[Im[value]]];
        AppendTo[multipoleRows,{A,energy,ER,-2*mA*ER,mA,Q,yValue,JA,J,channel,Re[value],Im[value]}]
      ],
      {J,spatialRanks},{pair,{{"Delta",dmformfactor`Private`DeltaJ},{"SigmaPrime",dmformfactor`Private`SigmaPJ}}}
    ];

    Do[
      If[dmformfactor`Private`DensityMatrix[i][[3]]==0 &&
        dmformfactor`Private`DensityMatrix[i][[2]]==J,
        Module[{np,jp,n,j,coefficient,mAmp,phiAmp,deltaAmp,sigmaAmp,sigmaMinus,sigmaPlus},
          np=dmformfactor`Private`DensityMatrix[i][[4]];
          jp=dmformfactor`Private`DensityMatrix[i][[5]]/2;
          n=dmformfactor`Private`DensityMatrix[i][[6]];
          j=dmformfactor`Private`DensityMatrix[i][[7]]/2;
          coefficient=dmformfactor`Private`DensityMatrix[i][[1]];
          mAmp=dmformfactor`Private`MJ[yValue,{np,jp},{n,j},J];
          phiAmp=dmformfactor`Private`PhiPPJ[yValue,{np,jp},{n,j},J];
          deltaAmp=If[J>0,dmformfactor`Private`DeltaJ[yValue,{np,jp},{n,j},J],0];
          sigmaAmp=If[J>0,dmformfactor`Private`SigmaPJ[yValue,{np,jp},{n,j},J],0];
          sigmaMinus=If[J>0,dmformfactor`Private`MJLSigma[yValue,
            {dmformfactor`Private`Nodal[np,jp],dmformfactor`Private`Lnumber[np,jp],jp},
            {dmformfactor`Private`Nodal[n,j],dmformfactor`Private`Lnumber[n,j],j},J,J-1],0];
          sigmaPlus=If[J>0,dmformfactor`Private`MJLSigma[yValue,
            {dmformfactor`Private`Nodal[np,jp],dmformfactor`Private`Lnumber[np,jp],jp},
            {dmformfactor`Private`Nodal[n,j],dmformfactor`Private`Lnumber[n,j],j},J,J+1],0];
          AppendTo[termRows,{A,energy,ER,-2*mA*ER,Q,yValue,J,i,coefficient,np,2 jp,n,2 j,
            Re[mAmp],Re[phiAmp],Re[deltaAmp],Re[sigmaAmp],Re[sigmaMinus],Re[sigmaPlus]}]
        ]],
      {J,chargeRanks},{i,densityLength}
    ];
    totalPoints++,
    {energy,recoilKeV}
  ],
  {isotope,isotopes}
];

If[maxImaginary>2*10^-12,fail["multipole imaginary residual exceeds tolerance: "<>ToString[maxImaginary],10]];
If[maxResponseImaginary>2*10^-12,fail["native response imaginary residual exceeds tolerance: "<>ToString[maxResponseImaginary],10]];

multipolePath=FileNameJoin[{outDir,"afh_v6_p4_fc_7isotope_81point_native_multipoles.csv"}];
responsePath=FileNameJoin[{outDir,"afh_v6_p4_fc_7isotope_81point_native_responses.csv"}];
termPath=FileNameJoin[{outDir,"afh_v6_p4_fc_7isotope_81point_native_obdme_terms.csv"}];
Export[multipolePath,Prepend[multipoleRows,{"A","ER_keV","ER_GeV","t_A_GeV2","mA_GeV","Q_A_GeV","y","J_A","J","channel","amplitude_real","amplitude_imag"}],"CSV"];
Export[responsePath,Prepend[responseRows,Join[
  {"A","ER_keV","ER_GeV","t_A_GeV2","mA_GeV","Q_A_GeV","y","J_A","C0"},
  ("W_"<>#<>"_raw")&/@responseNames,
  ("Im_W_"<>#<>"_raw")&/@responseNames
]],"CSV"];
Export[termPath,Prepend[termRows,{"A","ER_keV","ER_GeV","t_A_GeV2","Q_A_GeV","y","J","density_index","obdme_coefficient",
  "Nprime","two_jprime","N","two_j","M_element","PhiPP_element","Delta_element","SigmaPrime_element",
  "SigmaPrime_LJminus1_element","SigmaPrime_LJplus1_element"}],"CSV"];

sourceSHA=FileHash[packagePath,"SHA256","HexString"];
inputSHA=FileHash[massInputPath,"SHA256","HexString"];
scriptSHA=FileHash[$InputFileName,"SHA256","HexString"];
metadata={
  "status=NATIVE_AFH_P4_FC_SEVEN_ISOTOPE_EXPORT_COMPLETE",
  "wolfram_version="<>$Version,
  "system="<>$System,
  "processor_type="<>ToString[$ProcessorType],
  "processor_count="<>ToString[$ProcessorCount],
  "machine_name="<>$MachineName,
  "system_id="<>$SystemID,
  "afh_source_path="<>packagePath,
  "afh_source_sha256="<>sourceSHA,
  "mass_input_path="<>massInputPath,
  "mass_input_sha256="<>inputSHA,
  "mass_formula=(atomic_mass_u-54*0.000548579909065)*0.93149410242_GeV",
  "isotopes=128,129,130,131,132,134,136",
  "recoil_grid_keV=40 to 800 inclusive, 81 points, 9.5 keV step",
  "transfer=Q_A=sqrt(2*m_A*E_R+E_R^2); t_A=-2*m_A*E_R",
  "charge_multipoles=Native AFH MJ and PhiPPJ, all T=0 OBDME ranks",
  "spatial_multipoles=Native AFH DeltaJ and SigmaPJ, all nonzero T=0 OBDME ranks",
  "native_responses=ResponseNuclear[y,index,0,0], indices 1..8, raw package convention",
  "normalization=raw amplitudes/responses; no physical factor-four or C0 rescaling applied",
  "rank_map="<>StringRiffle[rankMap,";"],
  "c0_by_isotope="<>StringRiffle[(ToString[#[[1]]]<>"="<>ToString[InputForm[N[#[[6]],16]]])&/@c0Rows,";"],
  "recoil_isotope_points="<>ToString[totalPoints],
  "multipole_rows="<>ToString[Length[multipoleRows]],
  "native_response_rows="<>ToString[Length[responseRows]],
  "obdme_term_rows="<>ToString[Length[termRows]],
  "max_multipole_imaginary_residual="<>ToString[InputForm[N[maxImaginary,16]]],
  "max_response_imaginary_residual="<>ToString[InputForm[N[maxResponseImaginary,16]]],
  "multipole_csv_sha256="<>FileHash[multipolePath,"SHA256","HexString"],
  "response_csv_sha256="<>FileHash[responsePath,"SHA256","HexString"],
  "obdme_terms_csv_sha256="<>FileHash[termPath,"SHA256","HexString"],
  "script_sha256="<>scriptSHA
};
Export[FileNameJoin[{outDir,"native_export_metadata.txt"}],StringRiffle[metadata,"\n"]<>"\n","Text"];
Print[StringRiffle[metadata,"\n"]];
Exit[0];
