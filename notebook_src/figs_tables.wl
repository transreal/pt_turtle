(* figs_tables.wl : tables (rendered as figures) and survey plots *)
sup[n_Integer] := StringJoin[Characters[ToString[n]] /. Thread[Characters["0123456789"] -> Characters["⁰¹²³⁴⁵⁶⁷⁸⁹"]]];
mirror[p_String] := StringReplace[p, {"R" -> "L", "L" -> "R"}];
rots[p_String] := Table[StringRotateLeft[p, i], {i, 0, StringLength[p] - 1}];
canon[p_String] := First[Sort[Join[rots[p], rots[mirror[p]], rots[StringReverse[p]], rots[mirror[StringReverse[p]]]]]];   (* rotation, mirror, reversal *)
primitive[p_String] := Module[{n = StringLength[p]}, SelectFirst[StringTake[p, #] & /@ Divisors[n], Function[q, StringJoin[ConstantArray[q, n/StringLength[q]]] === p]]];
canonP[p_String] := canon[primitive[p]];
lead[w_String] := StringLength[First[StringCases[w, StartOfString ~~ "S" ...]]];
rotS[p_String] := First[SortBy[rots[p], {-lead[#], #} &]];
(* net exit sequence of an R-only program (None if not of that form) *)
exitsOf[p0_String] := Module[{p = If[StringContainsQ[p0, "L"], mirror[p0], p0], w, k, seq = {}, i = 1, j, b, m, n},
   If[StringContainsQ[p, "L"] || StringContainsQ[p <> p, "RR"], Return[None]];
   If[! StringContainsQ[p <> p, "SS"], Return["(SR)^k"]];
   n = StringLength[p]; k = First[First[StringPosition[p <> p, "SSR"]]];
   w = StringTake[p <> p <> p, {k + 3, k + 2 + n}];
   While[i <= n,
    j = i; While[j + 1 <= n && StringTake[w, {j, j + 1}] === "SR", j += 2]; b = (j - i)/2;
    m = 0; While[j + m <= n && StringTake[w, {j + m}] === "S", m++];
    If[m < 2 || j + m > n || StringTake[w, {j + m}] =!= "R", Return[None]];
    seq = Join[seq, {3 + b}, ConstantArray[2, m - 2]]; i = j + m + 1];
   First[Sort[Table[RotateLeft[seq, t], {t, 0, Length[seq] - 1}]]]];
progFromExits[seq_List] := Module[{i0 = First[FirstPosition[seq, _?(# >= 3 &)]], s, blocks = {}},
   s = RotateLeft[seq, i0 - 1];
   Do[If[k >= 3, AppendTo[blocks, {k, 0}], blocks[[-1, 2]]++], {k, s}];
   StringJoin[Table[Which[b[[1]] == 3, "", b[[1]] == 4, "(SR)", True, "(SR)" <> sup[b[[1]] - 3]] <> "(" <> If[b[[2]] + 2 >= 5, "S" <> sup[b[[2]] + 2], StringJoin[ConstantArray["S", b[[2]] + 2]]] <> "R)", {b, blocks}]]];
blockForm[seq_List, turn_String] := Module[{b = StringReplace[progFromExits[seq], "R" -> turn]}, If[StringCount[b, "("] == 1, b <> "*", "[" <> b <> "]*"]];
showProg[p_String] := Module[{ex = exitsOf[p]},
   Which[StringContainsQ[p, "L"] && ! StringContainsQ[p, "R"] && ListQ[ex], blockForm[ex, "L"],
    ListQ[ex] && ! StringContainsQ[p, "L"], blockForm[ex, "R"],
    True, "(" <> rotS[p] <> ")*"]];
showExits[p_String] := With[{ex = exitsOf[p]}, If[ListQ[ex], "(" <> StringRiffle[ToString /@ ex, ", "] <> ")", "—"]];
(* the survey programs are L-forms; mirror(reverse(p)) has identical statistics on the same tiling (Remark 1) *)
dual[p_String] := rotS[mirror[StringReverse[p]]];
cellS[s_, opts___] := Style[s, FontFamily -> ptFont, FontSize -> 6.3, opts];
tableFig[name_, header_, rows_, widths_, align_ : Automatic] := Module[{g},
   g = Grid[Prepend[Map[cellS, rows, {2}], cellS[#, Bold] & /@ header], ItemSize -> {widths, Automatic}, Alignment -> {If[align === Automatic, Left, align], Center},
     Dividers -> {None, {1 -> Directive[Black, AbsoluteThickness[0.9]], 2 -> Directive[Black, AbsoluteThickness[0.5]], -1 -> Directive[Black, AbsoluteThickness[0.9]]}},
     Background -> {None, {None, {GrayLevel[0.955], None}}}, Spacings -> {0.9, 0.55}];
   ptSave[name, g]];
fnum[x_] := ToString[NumberForm[x, {5, 1}]];
int[n_] := ToString[n];

(* ---------- bounded programs (Theorem 4) ---------- *)
vc = Import[FileNameJoin[{$dataDir, "verify_classes.wxf"}]];
pm = Import[FileNameJoin[{$repoDir, "results", "proof_many.json"}], "RawJSON"];
bkeys = SortBy[Keys[pm], {StringLength[#], #} &];
Print["bounded table rows: ", Length[bkeys], "; agreement with certificate: ", And @@ Table[vc[p]["starts"] == pm[p]["starts"] && vc[p]["escape"] == 0 && Keys[vc[p]["spectrum"]] === Sort[ToExpression /@ Keys[pm[p]["spectrum"]]], {p, bkeys}]];
brows = Table[{showProg[p], showExits[p], int[vc[p]["n"]], int[vc[p]["starts"]], int[vc[p]["nclasses"]], StringRiffle[ToString /@ Keys[vc[p]["spectrum"]], ", "]}, {p, bkeys}];
tableFig["tab-bounded", {"命令列", "出口列", "n", "初期状態", "類数", "周期の集合"}, brows, {17, 7, 1.5, 4.5, 2.5, 32}];
Print["classes of the four main programs: ", {#, vc[#]["nclasses"], vc[#]["classes"]} & /@ {"SSR", "SSRSSSSR", "SRSRSRSSR", "SSRSRSRSRSSR"}];
Print["extra single-exit programs: ", Table[{k, With[{r = vc[StringJoin[ConstantArray["SR", k - 3]] <> "SSR"]}, {r["starts"], r["escape"], r["timeout"], r["nclasses"], Keys[r["spectrum"]]}]}, {k, {35, 40, 45, 50}}]];
boundedCanon = Union[canonP /@ Join[bkeys, Table[StringJoin[ConstantArray["SR", k - 3]] <> "SSR", {k, {35, 40, 45, 50}}]]];

(* ---------- unbounded programs (Theorem 6) ---------- *)
sm = Import[FileNameJoin[{$repoDir, "results", "unbounded_certs", "summary.json"}], "RawJSON"];
sm = SortBy[Select[sm, ! MemberQ[{"SSRSR", "SLSRSLSR", "SSSRSSSR"}, #["prog"]] &], {StringLength[#["prog"]], #["prog"]} &];
unboundedCanon = Union[canonP[#["prog"]] & /@ sm];
urows = Table[{showProg[r["prog"]], showExits[r["prog"]], int[r["B"]], int[r["npass"]], int[r["nports"]], int[r["identities"]],
    If[r["root_type"][[1]] == 1, "鈍角", "鋭角"] <> If[r["root_type"][[2]] > 0, " (+)", " (−)"] <> ", " <> int[r["m0"]], StringRiffle[ToString /@ r["N"], ", "] <> ", …",
    "N(m+2) = " <> int[r["a"]] <> "N(m) " <> If[r["b"] >= 0, "+ ", "− "] <> int[Abs[r["b"]]]}, {r, sm}];
tableFig["tab-unbounded", {"命令列", "出口列", "B", "パス族", "ポート族", "恒等式", "根 (型, m₀)", "長さ N(m₀), N(m₀+2), …", "漸化式"}, urows, {13, 6.5, 1.5, 3, 3.5, 3, 6.5, 15, 11}];
Print["unbounded rows: ", Length[urows]];

(* ---------- survey of the 43 programs ---------- *)
sv = Association[Table[n -> Import[FileNameJoin[{$dataDir, "survey_exact_all43_" <> ToString[n] <> ".wxf"}]], {n, {46, 160, 300}}]];
progs43 = Keys[sv[300]];
pct[r_] := 100. r["closed"]/r["starts"];
verdict[p_] := Module[{c = pct[sv[#][p]] & /@ {46, 160, 300}, mr = sv[#][p]["maxradius"] & /@ {46, 160, 300}},
   Which[mr[[3]] < 25 && Abs[mr[[3]] - mr[[2]]] < 0.1 && c[[3]] > 95, "一様有界", c[[3]] < 30 && Abs[c[[3]] - c[[2]]] < 12, "拡散型", mr[[3]] > 100, "階層型", True, "混合"]];
proof[p_] := Which[MemberQ[boundedCanon, canonP[p]], "有界 (定理 4)", MemberQ[unboundedCanon, canonP[p]], "非有界 (定理 6)", True, "—"];
srows = Table[{showProg[dual[p]], showExits[dual[p]], StringRiffle[fnum[pct[sv[#][p]]] & /@ {46, 160, 300}, " / "], StringRiffle[int[sv[#][p]["nclasses"]] & /@ {46, 160, 300}, " / "], int[sv[300][p]["maxperiod"]],
    StringRiffle[fnum[sv[#][p]["maxradius"]] & /@ {46, 160, 300}, " / "], verdict[p], proof[p]}, {p, progs43}];
tableFig["tab-survey43", {"命令列", "出口列", "閉軌道率 % (窓 46 / 160 / 300)", "類数 (46 / 160 / 300)", "最大周期 (300)", "最大半径 (46 / 160 / 300)", "様相", "証明"}, srows, {11, 6, 11.5, 8.5, 5, 11.5, 4.5, 7}];
Print["verdict tally: ", Tally[verdict /@ progs43], "  proof tally: ", Tally[proof /@ progs43]];
Print["survey rows:"]; Scan[Print, srows];
Export[FileNameJoin[{$dataDir, "survey_rows.wxf"}], <|"rows" -> srows, "verdict" -> AssociationMap[verdict, progs43], "proof" -> AssociationMap[proof, progs43]|>];

(* ---------- survey regimes plot ---------- *)
Module[{col = <|"一様有界" -> RGBColor[0.15, 0.45, 0.85], "階層型" -> RGBColor[0.9, 0.5, 0.1], "拡散型" -> RGBColor[0.8, 0.15, 0.2], "混合" -> GrayLevel[0.5]|>, g1, g2, order, leg},
  order = SortBy[progs43, Position[{"階層型", "拡散型", "一様有界", "混合"}, verdict[#]] &];
  g1 = ListLogLogPlot[Table[{#, Max[sv[#][p]["maxradius"], 0.5]} & /@ {46, 160, 300}, {p, order}], Joined -> True, PlotStyle -> (Directive[col[verdict[#]], AbsoluteThickness[0.8], Opacity[0.75]] & /@ order),
    Frame -> True, FrameLabel -> {ptText["窓の一辺", 9], ptText["閉軌道の最大半径", 9]}, FrameTicks -> {{Automatic, None}, {{46, 160, 300}, None}}, PlotRange -> {{40, 340}, {0.7, 600}}, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8, Black}, AspectRatio -> 0.85];
  g2 = ListLogLinearPlot[Table[{#, pct[sv[#][p]]} & /@ {46, 160, 300}, {p, order}], Joined -> True, PlotStyle -> (Directive[col[verdict[#]], AbsoluteThickness[0.8], Opacity[0.75]] & /@ order),
    Frame -> True, FrameLabel -> {ptText["窓の一辺", 9], ptText["閉軌道率 (%)", 9]}, FrameTicks -> {{Automatic, None}, {{46, 160, 300}, None}}, PlotRange -> {{40, 340}, {0, 102}}, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8, Black}, AspectRatio -> 0.85];
  leg = LineLegend[Values[KeyTake[col, {"一様有界", "階層型", "拡散型"}]], {"一様有界", "階層型", "拡散型"}, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8.5}, LegendLayout -> "Row"];
  ptSave["survey-regimes", Column[{GraphicsRow[{g1, g2}, Spacings -> 10, ImageSize -> 600], leg}, Alignment -> Center]]];

(* ---------- single-exit family and pairs (if available) ---------- *)
fS = FileNameJoin[{$dataDir, "survey_exact_single_300.wxf"}];
If[FileExistsQ[fS],
  Module[{ss = Import[fS], ks = Join[Range[3, 30], {35, 40, 45, 50}], rec, stat, rows, col},
   rec[k_] := ss[StringJoin[ConstantArray["SR", k - 3]] <> "SSR"];
   ks = Select[ks, KeyExistsQ[ss, StringJoin[ConstantArray["SR", # - 3]] <> "SSR"] &];
   stat[k_] := With[{p = StringJoin[ConstantArray["SR", k - 3]] <> "SSR"}, Which[MemberQ[boundedCanon, canonP[p]], "有界", MemberQ[unboundedCanon, canonP[p]], "非有界", rec[k]["maxradius"] > 100, "階層型", True, "?"]];
   rows = Table[{int[k], fnum[pct[rec[k]]], int[rec[k]["nclasses"]], int[rec[k]["maxperiod"]], fnum[rec[k]["maxradius"]], stat[k]}, {k, ks}];
   Print["single-exit rows:"]; Scan[Print, rows];
   tableFig["tab-single", {"k", "閉軌道率 %", "類数", "最大周期", "最大半径", "判定"}, rows, {3, 7, 5, 7, 7, 7}, Right];
   col = <|"有界" -> RGBColor[0.15, 0.45, 0.85], "非有界" -> RGBColor[0.8, 0.15, 0.2], "階層型" -> RGBColor[0.9, 0.5, 0.1], "?" -> GrayLevel[0.5]|>;
   ptSave["single-exit", Column[{
      Graphics[{Table[With[{k = ks[[i]], h = Log10[Max[rec[ks[[i]]]["maxradius"], 1.01]]}, {col[stat[k]], Rectangle[{i - 0.4, 0}, {i + 0.4, h}], Black, Text[ptText[ToString[k], 6.5], {i, -0.12}]}], {i, Length[ks]}],
        {GrayLevel[0.3], AbsoluteThickness[0.5], Line[{{0.3, 0}, {Length[ks] + 0.7, 0}}], Table[{GrayLevel[0.8], Line[{{0.3, y}, {Length[ks] + 0.7, y}}], Black, Text[ptText[ToString[10^y], 7], {0.1, y}, {1, 0}]}, {y, {1, 2}}]},
        Text[ptText["閉軌道の最大半径 (窓 300)", 8], {Length[ks]/2, 2.95}], Text[ptText["k", 8, Italic], {Length[ks] + 1.2, -0.12}]}, PlotRange -> {{-2.2, Length[ks] + 2}, {-0.3, 3.1}}, AspectRatio -> 0.36, ImageSize -> 600],
      LineLegend[Values[KeyTake[col, {"有界", "非有界", "階層型"}]], {"有界 (証明済)", "非有界 (証明済)", "階層型 (走査による分類)"}, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8.5}, LegendLayout -> "Row"]}, Alignment -> Center]];
   Print["single closed counts: ", Table[{k, rec[k]["closed"], rec[k]["nclasses"], rec[k]["maxperiod"], Round[rec[k]["maxradius"], 0.1], rec[k]["timeout"]}, {k, ks}]];
   Print["identical statistics: ", Select[Subsets[ks, {2}], Function[pr, rec[pr[[1]]]["closed"] == rec[pr[[2]]]["closed"] && rec[pr[[1]]]["classes"] === rec[pr[[2]]]["classes"]]]];
   Print["coincidences: ", Table[{pr, rec[pr[[1]]]["nclasses"] == rec[pr[[2]]]["nclasses"], rec[pr[[1]]]["maxperiod"] == rec[pr[[2]]]["maxperiod"], Abs[rec[pr[[1]]]["maxradius"] - rec[pr[[2]]]["maxradius"]] < 0.05}, {pr, Select[{{6, 18}, {12, 24}, {10, 50}, {20, 40}, {15, 45}, {14, 16}}, SubsetQ[ks, #] &]}]]]];
fP = FileNameJoin[{$dataDir, "survey_exact_pairs_300.wxf"}];
If[FileExistsQ[fP],
  Module[{pp = Import[fP], rows},
   rows = Table[With[{r = pp[p]}, {showExits[p], showProg[p], fnum[pct[r]], int[r["nclasses"]], int[r["maxperiod"]], fnum[r["maxradius"]], Which[MemberQ[boundedCanon, canonP[p]], "有界 (定理 4)", MemberQ[unboundedCanon, canonP[p]], "非有界 (定理 6)", r["maxradius"] > 100, "階層型", True, "?"]}], {p, Keys[pp]}];
   Print["pairs rows:"]; Scan[Print, rows];
   tableFig["tab-pairs", {"出口列", "命令列", "閉軌道率 %", "類数", "最大周期", "最大半径", "判定"}, rows, {5, 16, 6, 4, 6, 6, 9}]]];
