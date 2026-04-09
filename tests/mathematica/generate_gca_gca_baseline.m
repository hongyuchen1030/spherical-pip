#!/usr/bin/env math -script

defaultInputPath =
  "/pscratch/sd/h/hyvchen/regrid_bench/GCAGCA/arcs/gca_gca_pairs_seed20251104_N100.csv";
defaultOutputPath =
  "tests/input/gca_gca_pairs_seed20251104_N100_with_baseline.csv";

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
    inputPath = defaultInputPath;
    outputPath = defaultOutputPath;,
  Length[csvArgs] == 2,
    {inputPath, outputPath} = csvArgs;,
  True,
    Print["Usage: math -script tests/mathematica/generate_gca_gca_baseline.wls [input_csv output_csv]"];
    Print["Defaults:"];
    Print["  input_csv  = ", defaultInputPath];
    Print["  output_csv = ", defaultOutputPath];
    Print["Received CSV-like args: ", csvArgs];
    Quit[1];
];

wp = 80;
zeroTol = SetPrecision[10^-40, wp];
coplanarTol = SetPrecision[10^-30, wp];
arcTol = SetPrecision[10^-30, wp];
significandBits = 53;

safeLines[path_] := Module[{lines},
  lines = Quiet @ Check[Import[path, "Lines"], $Failed];
  If[lines === $Failed || !ListQ[lines], {}, Select[lines, StringLength[StringTrim[#]] > 0 &]]
];

asInt[x_] := Module[{y = StringTrim[x]},
  Quiet @ Check[ToExpression[y], $Failed]
];

fromSigExp2[sig_, exp_] := SetPrecision[asInt[sig] * 2^asInt[exp], wp];

vec3FromRow[row_List, start_Integer] := {
  fromSigExp2[row[[start]], row[[start + 1]]],
  fromSigExp2[row[[start + 2]], row[[start + 3]]],
  fromSigExp2[row[[start + 4]], row[[start + 5]]]
};

norm3[v_List] := Sqrt[v.v];

normalize3[v_List] := Module[{n = norm3[v]},
  If[n <= zeroTol, Throw["normalize3: zero vector"]];
  v/n
];

onMinorArc[q_List, a_List, b_List] := Module[
  {ab, abNorm, coplanar, dir1, dir2},
  If[norm3[a - b] <= zeroTol || norm3[a + b] <= zeroTol, Return[False]];
  ab = Cross[a, b];
  abNorm = norm3[ab];
  If[abNorm <= zeroTol, Return[False]];
  coplanar = Abs[q.ab];
  If[coplanar > coplanarTol, Return[False]];
  dir1 = Cross[a, q].ab;
  dir2 = Cross[q, b].ab;
  dir1 >= -arcTol && dir2 >= -arcTol
];

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

appendBaseline[row_List] := Module[
  {a0, a1, b0, b1, w, wNorm, cand1, cand2, cand1On, cand2On, baseline},

  If[Length[row] != 26, Return[Missing["BadRow"]]];

  a0 = vec3FromRow[row, 3];
  a1 = vec3FromRow[row, 9];
  b0 = vec3FromRow[row, 15];
  b1 = vec3FromRow[row, 21];

  w = Cross[Cross[a0, a1], Cross[b0, b1]];
  wNorm = norm3[w];

  (* Degenerate great circles → skip *)
  If[wNorm <= zeroTol, Return[Missing["Degenerate"]]];

  cand1 = normalize3[w];
  cand2 = -cand1;

  cand1On = onMinorArc[cand1, a0, a1] && onMinorArc[cand1, b0, b1];
  cand2On = onMinorArc[cand2, a0, a1] && onMinorArc[cand2, b0, b1];

  Which[
    cand1On && !cand2On,
      baseline = cand1;,

    cand2On && !cand1On,
      baseline = cand2;,

    (* No valid arc intersection OR ambiguous *)
    True,
      Return[Missing["NoIntersection"]]
  ];

  Join[row, Sequence @@ (decomposeDouble /@ baseline)]
];

lines = safeLines[inputPath];
If[Length[lines] < 2,
  Print["[ERROR] Failed to read non-empty CSV: ", inputPath];
  Quit[5];
];

header = StringSplit[First[lines], ","];
If[Length[header] != 26,
  Print["[ERROR] Expected 26 input columns, got ", Length[header]];
  Quit[6];
];

rows = StringSplit[#, ","] & /@ Rest[lines];

processed = appendBaseline /@ rows;

validRows = Select[processed, ListQ];

dropped = Length[processed] - Length[validRows];

Print["[INFO] Total rows: ", Length[processed]];
Print["[INFO] Valid arc intersections: ", Length[validRows]];
Print["[INFO] Dropped rows: ", dropped];

outputHeader = Join[
  header,
  {
    "baseline_x_sig", "baseline_x_exp",
    "baseline_y_sig", "baseline_y_exp",
    "baseline_z_sig", "baseline_z_exp"
  }
];

If[Length[validRows] == 0,
  Print["[WARNING] No valid arc-arc intersections found. Output CSV not written."];
  Quit[0];
];

Export[
  outputPath,
  Prepend[validRows, outputHeader],
  "CSV",
  "TextDelimiters" -> None
];

Print["[INFO] Wrote ", Length[validRows], " baseline rows to ", outputPath];