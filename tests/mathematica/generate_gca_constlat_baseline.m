#!/usr/bin/env math -script

defaultOutputPath = "tests/input/gca_constlat_cases_with_baseline.csv";
targetPerBin = 100;
seed = 20260410;

getCsvArgs[] := Module[
  {
    raw = DeleteDuplicates @ Join[
      If[ValueQ[$ScriptCommandLine], $ScriptCommandLine, {}],
      If[ValueQ[$CommandLine], $CommandLine, {}]
    ]
  },
  Select[raw, StringQ[#] && StringMatchQ[ToLowerCase[#], ___ ~~ ".csv"] &]
];

csvArgs = getCsvArgs[];
Which[
  Length[csvArgs] == 0,
    outputPath = defaultOutputPath;,
  Length[csvArgs] == 1,
    outputPath = First[csvArgs];,
  True,
    Print["Usage: math -script tests/mathematica/generate_gca_constlat_baseline.m [output_csv]"];
    Print["Default output_csv = ", defaultOutputPath];
    Print["Received CSV-like args: ", csvArgs];
    Quit[1];
];

wp = 100;
zeroTol = SetPrecision[10^-50, wp];
coplanarTol = SetPrecision[10^-35, wp];
arcTol = SetPrecision[10^-35, wp];
constLatTol = SetPrecision[10^-35, wp];
sphereTol = SetPrecision[10^-30, wp];
significandBits = 53;
maxAttemptsPerBin = 10000;

asInt[x_] := Module[{y = StringTrim[ToString[x]]},
  Quiet @ Check[ToExpression[y], $Failed]
];

fromSigExpPair[{sig_, exp_}] := SetPrecision[asInt[sig] * 2^asInt[exp], wp];

toMachineReal[x_] := N[x, MachinePrecision];

decomposeDouble[x_] := Module[{xm = toMachineReal[x], mant, exp2, sig},
  If[xm == 0., Return[{"0", "0"}]];
  {mant, exp2} = MantissaExponent[xm, 2];
  sig = Round[mant * 2^significandBits];
  {
    ToString[InputForm[sig]],
    ToString[InputForm[exp2 - significandBits]]
  }
];

serializeScalar[x_] := decomposeDouble[x];

roundTripScalar[x_] := fromSigExpPair[serializeScalar[x]];

serializeVec3[v_List] := Flatten[serializeScalar /@ v];

norm3[v_List] := Sqrt[v.v];

onMinorArc[q_List, a_List, b_List] := Module[
  {ab, abNorm, coplanar, abDot},
  If[norm3[a - b] <= zeroTol || norm3[a + b] <= zeroTol, Return[False]];
  ab = Cross[a, b];
  abNorm = norm3[ab];
  If[abNorm <= zeroTol, Return[False]];
  coplanar = Abs[q.ab];
  If[coplanar > coplanarTol, Return[False]];
  abDot = a.b;
  q.a >= abDot - arcTol && q.b >= abDot - arcTol
];

latLonToVec[latDeg_, lonDeg_] := Module[
  {lat = SetPrecision[latDeg Degree, wp], lon = SetPrecision[lonDeg Degree, wp]},
  {
    Cos[lat] Cos[lon],
    Cos[lat] Sin[lon],
    Sin[lat]
  }
];

computeBaseline[a0_List, a1_List, z0_] := Module[
  {
    normal, nx, ny, nz, denom, normalSq, planarSq, planar,
    candPos, candNeg, posValid, negValid, baseline
  },

  normal = Cross[a0, a1];
  {nx, ny, nz} = normal;
  denom = nx^2 + ny^2;
  normalSq = normal.normal;

  If[norm3[normal] <= zeroTol || denom <= zeroTol, Return[Missing["Degenerate"]]];

  planarSq = denom - z0^2 normalSq;
  If[planarSq <= zeroTol, Return[Missing["NoRealIntersection"]]];
  planar = Sqrt[planarSq];

  candPos = {
    -((nx nz z0 - ny planar) / denom),
    -((ny nz z0 + nx planar) / denom),
    z0
  };
  candNeg = {
    -((nx nz z0 + ny planar) / denom),
    -((ny nz z0 - nx planar) / denom),
    z0
  };

  posValid =
    Abs[candPos[[3]] - z0] <= constLatTol &&
    Abs[candPos.candPos - 1] <= sphereTol &&
    Abs[normal.candPos] <= coplanarTol &&
    onMinorArc[candPos, a0, a1];

  negValid =
    Abs[candNeg[[3]] - z0] <= constLatTol &&
    Abs[candNeg.candNeg - 1] <= sphereTol &&
    Abs[normal.candNeg] <= coplanarTol &&
    onMinorArc[candNeg, a0, a1];

  Which[
    posValid && !negValid, baseline = candPos,
    negValid && !posValid, baseline = candNeg,
    True, Return[Missing["NoUniqueMinorArcIntersection"]]
  ];

  baseline
];

generateOneCase[caseId_Integer, latMin_, latMax_] := Module[
  {
    lat0, lat1, lon0, lonDelta, lon1, mix, a0Exact, a1Exact, z0Exact,
    a0, a1, z0, baseline
  },

  lat0 = SetPrecision[RandomReal[{latMin, latMax}], wp];
  lat1 = SetPrecision[RandomReal[{latMin, latMax}], wp];
  lon0 = SetPrecision[RandomReal[{-180, 180}], wp];
  lonDelta = SetPrecision[RandomChoice[{-1, 1}] RandomReal[{1, 20}], wp];
  lon1 = lon0 + lonDelta;
  mix = SetPrecision[RandomReal[{0.01, 0.99}], wp];

  a0Exact = latLonToVec[lat0, lon0];
  a1Exact = latLonToVec[lat1, lon1];
  z0Exact = mix a0Exact[[3]] + (1 - mix) a1Exact[[3]];

  a0 = roundTripScalar /@ a0Exact;
  a1 = roundTripScalar /@ a1Exact;
  z0 = roundTripScalar[z0Exact];

  baseline = computeBaseline[a0, a1, z0];
  If[Head[baseline] === Missing, Return[Missing["InvalidCase"]]];

  Join[
    {ToString[caseId]},
    serializeVec3[a0],
    serializeVec3[a1],
    serializeScalar[z0],
    serializeScalar[baseline[[1]]],
    serializeScalar[baseline[[2]]]
  ]
];

generateBin[latMin_, latMax_, startCaseId_Integer] := Module[
  {rows = {}, attempts = 0, caseId = startCaseId, row},
  While[Length[rows] < targetPerBin && attempts < maxAttemptsPerBin,
    attempts++;
    row = generateOneCase[caseId, latMin, latMax];
    If[ListQ[row],
      AppendTo[rows, row];
      caseId++;
    ];
  ];
  If[Length[rows] < targetPerBin,
    Print["[ERROR] Failed to generate ", targetPerBin,
          " valid cases for latitude range [", latMin, ", ", latMax,
          "] after ", attempts, " attempts."];
    Quit[2];
  ];
  <|"rows" -> rows, "attempts" -> attempts, "nextCaseId" -> caseId|>
];

SeedRandom[seed];

equator = generateBin[0, 1, 1];
pole = generateBin[89, 90, equator["nextCaseId"]];

validRows = Join[equator["rows"], pole["rows"]];
totalAttempts = equator["attempts"] + pole["attempts"];
dropped = totalAttempts - Length[validRows];

outputHeader = {
  "case_id",
  "a0_x_sig", "a0_x_exp", "a0_y_sig", "a0_y_exp", "a0_z_sig", "a0_z_exp",
  "a1_x_sig", "a1_x_exp", "a1_y_sig", "a1_y_exp", "a1_z_sig", "a1_z_exp",
  "z0_sig", "z0_exp",
  "baseline_x_sig", "baseline_x_exp", "baseline_y_sig", "baseline_y_exp"
};

Print["[INFO] Total attempted: ", totalAttempts];
Print["[INFO] Valid kept: ", Length[validRows]];
Print["[INFO] Dropped: ", dropped];

Export[
  outputPath,
  Prepend[validRows, outputHeader],
  "CSV",
  "TextDelimiters" -> None
];

Print["[INFO] Wrote ", Length[validRows], " baseline rows to ", outputPath];
