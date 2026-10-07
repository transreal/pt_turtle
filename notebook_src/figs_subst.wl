(* figs_subst.wl : Robinson triangles, substitution, supertiles, corona atlas, closure, verification patch, ownership *)
triCol[col_, sg_] := Which[col == 0 && sg > 0, RGBColor[0.96, 0.60, 0.33], col == 0, RGBColor[0.995, 0.85, 0.70], col == 1 && sg > 0, RGBColor[0.40, 0.60, 0.86], True, RGBColor[0.78, 0.87, 0.97]];
tx[s_, pos_, size_ : 9, opts___] := Text[ptText[s, size, opts], pos];
mx[s_, pos_, size_ : 11, opts___] := Text[ptMath[s, size, opts], pos];
title[s_] := ptText[s, 9.5, Black];
drawTris[t_, edge_ : Directive[GrayLevel[0.3], AbsoluteThickness[0.4]]] := Module[{xy = ptTriXY[t], sg = ptSign[t], col = t["col"]},
   {EdgeForm[edge], Table[With[{idx = Flatten[Position[Transpose[{col, sg}], {c, s}]]}, If[idx === {}, {}, {triCol[c, s], Polygon[xy[[idx]]]}]], {c, {0, 1}}, {s, {1, -1}}]}];
legend[] := Graphics[{Table[{EdgeForm[GrayLevel[0.3]], triCol[it[[1]], it[[2]]], Rectangle[{it[[3]], 0}, {it[[3]] + 0.35, 0.22}], Black, Text[ptText[it[[4]], 8], {it[[3]] + 0.42, 0.11}, {-1, 0}]},
     {it, {{0, 1, 0., "鋭角 (+)"}, {0, -1, 1.5, "鋭角 (−)"}, {1, 1, 3.0, "鈍角 (+)"}, {1, -1, 4.5, "鈍角 (−)"}}}]}, PlotRange -> {{-0.1, 6.}, {-0.1, 0.32}}, ImageSize -> 520];

(* ---------- Robinson triangles in the two rhombs ---------- *)
Module[{rhombPanel, g},
  rhombPanel[col_, ttl_] := Module[{t = ptStdTriangle[col, 1], A, B, C, As, cen, pA, pB, pC, pAs},
    {A, B, C} = Rest[t]; As = B + C - A; {pA, pB, pC, pAs} = ptEmbed /@ {A, B, C, As}; cen = Mean[{pA, pB, pC, pAs}];
    Graphics[{EdgeForm[{GrayLevel[0.2], AbsoluteThickness[0.8]}], {triCol[col, 1], Polygon[{pA, pB, pC}]}, {triCol[col, -1], Polygon[{pAs, pB, pC}]},
      {Black, Dashed, AbsoluteThickness[1.], Line[{pB, pC}]},
      mx["A", pA + 0.13 Normalize[pA - cen], 12], mx["B", pB + 0.13 Normalize[pB - cen], 12], mx["C", pC + 0.13 Normalize[pC - cen], 12], mx["A*", pAs + 0.14 Normalize[pAs - cen], 12],
      tx["(+)", Mean[{pA, pB, pC}], 9], tx["(−)", Mean[{pAs, pB, pC}], 9],
      {GrayLevel[0.1], AbsoluteThickness[0.8], Arrowheads[0.04], Arrow[Table[Mean[{pA, pB, pC}] + 0.085 {Cos[a], Sin[a]}, {a, 0.3, 5.2, 0.2}]], Arrow[Table[Mean[{pAs, pB, pC}] + 0.085 {Cos[a], Sin[a]}, {a, 5.2, 0.3, -0.2}]]}},
     PlotLabel -> title[ttl], PlotRange -> {cen[[1]] + {-1.2, 1.2}, cen[[2]] + {-0.72, 0.72}}]];
  ptSave["robinson-triangles", GraphicsRow[{rhombPanel[0, "細い菱形 = 鋭角三角形 (36°, 72°, 72°) 2 枚"], rhombPanel[1, "太い菱形 = 鈍角三角形 (108°, 36°, 36°) 2 枚"]}, Spacings -> 10, ImageSize -> 560]]];

(* ---------- substitution rule ---------- *)
Module[{panel},
  panel[col_, sign_, ttl_] := Module[{t0 = ptTri[{ptStdTriangle[col, sign]}], t1, xy0, xy1, sg1, cen0, cen1, names, pts, shift, g0, g1, lab0},
    t1 = ptSubdivide[t0]; xy0 = First[ptTriXY[t0]]; xy1 = ptTriXY[t1]; sg1 = ptSign[t1];
    cen0 = Mean[xy0]; cen1 = Mean[ptPhi xy0];
    shift = If[col == 0, {1.55, 0}, {2.05, 0}];
    lab0 = MapThread[mx[#1, #2 + 0.11 Normalize[#2 - cen0], 11] &, {{"A", "B", "C"}, xy0}];
    (* named points of the inflated triangle *)
    names = If[col == 0,
      {{"A'", ptPhi xy0[[1]]}, {"B'", ptPhi xy0[[2]]}, {"C'", ptPhi xy0[[3]]}, {"P", ptPhi xy0[[1]] + (xy0[[2]] - xy0[[1]])}},
      {{"A'", ptPhi xy0[[1]]}, {"B'", ptPhi xy0[[2]]}, {"C'", ptPhi xy0[[3]]}, {"Q", ptPhi xy0[[2]] + (xy0[[1]] - xy0[[2]])}, {"R", ptPhi xy0[[2]] + (xy0[[3]] - xy0[[2]])}}];
    Graphics[{EdgeForm[{GrayLevel[0.2], AbsoluteThickness[0.8]}],
      {triCol[col, sign], Polygon[xy0]}, lab0,
      {GrayLevel[0.3], AbsoluteThickness[1.2], Arrowheads[0.035], Arrow[{{1.12, 0.32}, {1.42, 0.32}}], mx["σ", {1.27, 0.4}, 12, FontSlant -> Plain]},
      Table[{triCol[t1["col"][[i]], sg1[[i]]], Polygon[(# + shift) & /@ xy1[[i]]]}, {i, Length[xy1]}],
      Table[With[{c = Mean[xy1[[i]]]}, MapThread[Text[Style[#1, FontFamily -> "Times New Roman", Italic, FontSize -> 7.5, GrayLevel[0.05]], shift + #2 + 0.13 Normalize[c - #2]] &, {{"A", "B", "C"}, xy1[[i]]}]], {i, Length[xy1]}],
      Table[tx["子 " <> ToString[t1["slot"][[i]]], shift + Mean[xy1[[i]]] + {0, -0.0}, 7.5, GrayLevel[0.1]], {i, Length[xy1]}],
      Table[{Black, PointSize[0.012], Point[shift + nm[[2]]], mx[nm[[1]], shift + nm[[2]] + 0.12 Normalize[nm[[2]] - cen1 + {0.001, 0.}], 11]}, {nm, names}]},
     PlotLabel -> title[ttl], PlotRange -> If[col == 0, {{-0.55, 4.0}, {-0.22, 1.12}}, {{-0.55, 4.0}, {-0.22, 1.72}}]]];
  ptSave["substitution-rule", GraphicsColumn[{panel[0, 1, "鋭角三角形 (+): P = φA + (B − A)"], panel[1, 1, "鈍角三角形 (+): Q = φB + (A − B),  R = φB + (C − B)"], legend[]}, Spacings -> 4, ImageSize -> 560]]];

(* ---------- supertiles ---------- *)
Module[{panel},
  panel[m_] := Module[{t = ptSubstitute[ptTri[{ptStdTriangle[1, 1]}], m]},
    Graphics[{drawTris[t, Directive[GrayLevel[0.25], AbsoluteThickness[If[m > 4, 0.15, 0.3]]]], {Black, AbsoluteThickness[1.1], Line[ptPhi^m Append[#, First[#]] &[First[ptTriXY[ptTri[{ptStdTriangle[1, 1]}]]]]]}},
     PlotLabel -> title["レベル " <> ToString[m] <> " (" <> ToString[ptNTri[t]] <> " 枚)"], PlotRangePadding -> Scaled[0.03]]];
  ptSave["supertiles", GraphicsGrid[{{panel[1], panel[2], panel[3]}, {panel[4], panel[5], panel[6]}}, Spacings -> {4, 6}, ImageSize -> 600]]];

(* ---------- corona atlas ---------- *)
keys = Import[FileNameJoin[{$dataDir, "atlas.wxf"}]];
ctype[k_] := With[{t = ptCanonToTri[k], c = ptCanonCentre[k]}, {t["col"][[c]], -ptSign[t][[c]]}];
keysS = SortBy[keys, {ctype[#], Length[#], #} &];
Module[{panel, allxy, xr, yr},
  allxy = Flatten[ptTriXY[ptCanonToTri[#]] & /@ keysS, 2];
  xr = MinMax[allxy[[All, 1]]] + {-0.08, 0.08}; yr = MinMax[allxy[[All, 2]]] + {-0.08, 0.08};
  panel[k_, i_] := Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], xy},
    xy = ptTriXY[t];
    Graphics[{drawTris[t, Directive[GrayLevel[0.3], AbsoluteThickness[0.25]]], {EdgeForm[{Black, AbsoluteThickness[1.3]}], FaceForm[None], Polygon[xy[[c]]]}, {Black, PointSize[0.025], Point[Mean[xy[[c]]]]},
      Text[ptText[ToString[i] <> " (" <> ToString[Length[k]] <> ")", 7.5], {xr[[1]] + 0.05, yr[[2]] - 0.05}, {-1, 1}]},
     PlotRange -> {xr, yr}, PlotRangePadding -> 0]];
  ptSave["corona-atlas", Column[{GraphicsGrid[Partition[MapIndexed[panel[#1, #2[[1]]] &, keysS], 6], Spacings -> {1, 1}, ImageSize -> 600], legend[]}, Alignment -> Center]]];
Print["atlas order types: ", ctype /@ keysS];

(* ---------- closure illustration ---------- *)
Module[{k = SelectFirst[keysS, ctype[#] == {1, -1} && Length[#] == 17 &], t, c, st, inc, kids, cors, xy0, xy1, g0, g1, g2, fade, j = 1},
  t = ptCanonToTri[k]; c = ptCanonCentre[k]; st = ptSubstitute[t, 1]; inc = ptTriIncidence[st];
  kids = Flatten[Position[st["anc"], c]]; cors = ptCorona[inc, #] & /@ kids;
  xy0 = ptTriXY[t]; xy1 = ptTriXY[st];
  g0 = Graphics[{drawTris[t], {EdgeForm[{Black, AbsoluteThickness[1.8]}], FaceForm[None], Polygon[xy0[[c]]]}}, PlotLabel -> title["(a) コロナ C と中心 t"], PlotRangePadding -> Scaled[0.04]];
  g1 = Graphics[{drawTris[st, Directive[GrayLevel[0.4], AbsoluteThickness[0.25]]],
     {GrayLevel[0.15], AbsoluteThickness[0.8], Line[Append[#, First[#]]] & /@ (ptPhi xy0)},
     {EdgeForm[{Black, AbsoluteThickness[1.8]}], FaceForm[None], Polygon[xy1[[kids]]]}},
    PlotLabel -> title["(b) σ(C) と σ(t) の子 3 枚"], PlotRangePadding -> Scaled[0.04]];
  fade = Complement[Range[ptNTri[st]], cors[[j]]];
  g2 = Graphics[{{EdgeForm[{GrayLevel[0.75], AbsoluteThickness[0.25]}], GrayLevel[0.94], Polygon[xy1[[fade]]]},
     drawTris[KeyTake[Map[#[[cors[[j]]]] &, KeyTake[st, {"col", "A", "B", "C"}]], {"col", "A", "B", "C"}], Directive[GrayLevel[0.3], AbsoluteThickness[0.3]]],
     {EdgeForm[{Black, AbsoluteThickness[1.8]}], FaceForm[None], Polygon[xy1[[kids[[j]]]]]}},
    PlotLabel -> title["(c) 子 t' のコロナ (σ(C) の内部で完全)"], PlotRangePadding -> Scaled[0.04]];
  Print["closure example: corona size ", Length[k], ", sigma(C) has ", ptNTri[st], " triangles; child corona sizes ", Length /@ cors];
  ptSave["closure", GraphicsRow[{g0, g1, g2}, Spacings -> 4, ImageSize -> 600]]];

(* ---------- verification patch ---------- *)
Module[{k = SelectFirst[keysS, ctype[#] == {1, -1} && Length[#] == 17 &], t, c, n = 7, st, rh, nav, xy, cent, mask, inside, starts, res, reps, visited, tri0, outline, orbs, cols, prog = ptProg["SSR"], g},
  t = ptCanonToTri[k]; c = ptCanonCentre[k]; st = ptSubstitute[t, n]; rh = ptRhombs[st]; nav = ptNav[rh]; xy = ptRhombXY[rh]; cent = Mean /@ xy;
  mask = Unitize[st["anc"] - c];
  inside = Flatten[Position[mask[[rh["tri"][[All, 1]]]]*mask[[rh["tri"][[All, 2]]]], 0]];
  starts = Flatten[Table[4 (r - 1) + e + 1, {r, inside}, {e, 0, 3}]];
  res = ptScanC[nav["sNext"], prog, starts, 100000];
  reps = DeleteDuplicates[res[[All, 3]]];
  orbs = Association[Table[s0 -> ptRun[nav["sNext"], prog, s0], {s0, reps}]];
  visited = Union @@ (DeleteDuplicates[#["rhombs"]] & /@ Values[orbs]);
  outline = ptPhi^n Append[#, First[#]] &[ptTriXY[t][[c]]];
  cols = <|100 -> RGBColor[0.85, 0.15, 0.15], 80 -> RGBColor[0.95, 0.55, 0.1], 50 -> RGBColor[0.1, 0.6, 0.3], 30 -> RGBColor[0.55, 0.2, 0.75], 28 -> RGBColor[0.1, 0.35, 0.8]|>;
  g = Graphics[{EdgeForm[{GrayLevel[0.7], AbsoluteThickness[0.12]}],
     {GrayLevel[0.96], Polygon[xy]},
     {RGBColor[0.82, 0.90, 1.], Polygon[xy[[visited]]]},
     {RGBColor[0.55, 0.72, 0.95], Polygon[xy[[inside]]]},
     MapIndexed[Function[{pc, ii}, With[{tgt = Most[outline][[Mod[ii[[1]], 3] + 1]] + ({0.25, 0.5, 0.72, 0.4, 0.6}[[ii[[1]]]]) (Mean[Most[outline]] - Most[outline][[Mod[ii[[1]], 3] + 1]])},
         With[{o = First[SortBy[Select[Values[orbs], #["cycles"] == pc[[1]] &], Norm[Mean[cent[[#["rhombs"]]]] - tgt] &]]}, {pc[[2]], EdgeForm[{GrayLevel[0.2], AbsoluteThickness[0.15]}], Polygon[xy[[DeleteDuplicates[o["rhombs"]]]]]}]]], Normal[cols] /. Rule -> List],
     {Black, AbsoluteThickness[1.3], Line[outline]}}, PlotRangePadding -> Scaled[0.01], ImageSize -> 600];
  Print["verification patch: rhombs ", Length[xy], " starts ", Length[starts], " flags ", Tally[res[[All, 1]]], " distinct orbits ", Length[reps], " visited ", Length[visited], " periods ", Sort[Tally[#["cycles"] & /@ Values[orbs]]]];
  ptSave["verification-patch", g]];

(* ---------- ownership ---------- *)
Module[{m = 4, t0 = ptTri[{ptStdTriangle[1, 1]}], t, rhG, rhP, xyT, sg, posIdx, ghost, xyR, outline, neg, negUnowned, paired},
  t = ptSubstitute[t0, m]; xyT = ptTriXY[t]; sg = ptSign[t];
  rhG = ptRhombs[t, "Ghost" -> True]; xyR = ptRhombXY[rhG];
  ghost = Flatten[Position[rhG["tri"][[All, 2]], 0]];
  paired = Complement[Range[Length[xyR]], ghost];
  negUnowned = Complement[Flatten[Position[sg, -1]], rhG["tri"][[paired, 2]]];
  outline = ptPhi^m Append[#, First[#]] &[First[ptTriXY[t0]]];
  Print["ownership: level ", m, " triangles ", ptNTri[t], " positive ", Count[sg, 1], " owned rhombs ", Length[xyR], " (ghost-completed ", Length[ghost], "), negative triangles of non-owned rhombs ", Length[negUnowned]];
  ptSave["ownership", Graphics[{
      (* ghost halves outside X *)
      {EdgeForm[{GrayLevel[0.35], AbsoluteThickness[0.5], Dashing[{0.006, 0.006}]}], RGBColor[1, 0.95, 0.75], Polygon[Table[{xyR[[r, 2]], xyR[[r, 3]], xyR[[r, 4]]}, {r, ghost}]]},
      (* triangles of X *)
      {EdgeForm[{GrayLevel[0.45], AbsoluteThickness[0.3]}], {RGBColor[0.55, 0.72, 0.93], Polygon[xyT[[Flatten[Position[sg, 1]]]]]}, {RGBColor[0.84, 0.91, 0.99], Polygon[xyT[[rhG["tri"][[paired, 2]]]]]}, {GrayLevel[0.88], Polygon[xyT[[negUnowned]]]}},
      (* owned rhombs outlines *)
      {EdgeForm[{RGBColor[0.1, 0.2, 0.55], AbsoluteThickness[0.9]}], FaceForm[None], Polygon[xyR]},
      {Black, AbsoluteThickness[1.8], Line[outline]},
      MapIndexed[{EdgeForm[GrayLevel[0.3]], #1[[1]], Rectangle[{3.6, 6.7 - 0.42 #2[[1]]}, {4.05, 6.95 - 0.42 #2[[1]]}], Black, Text[ptText[#1[[2]], 8.5], {4.2, 6.82 - 0.42 #2[[1]]}, {-1, 0}]} &,
       {{RGBColor[0.55, 0.72, 0.93], "X の正の三角形 (所有する菱形の正の半分)"}, {RGBColor[0.84, 0.91, 0.99], "その鏡像の相手 (X の内部)"}, {RGBColor[1, 0.95, 0.75], "その鏡像の相手 (X の外部, 幽霊半分)"}, {GrayLevel[0.88], "X の負の三角形で相手が X の外のもの (所有しない)"}}]}, PlotRangePadding -> Scaled[0.03], ImageSize -> 560]]];
