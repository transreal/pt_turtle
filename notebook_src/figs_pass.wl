(* figs_pass.wl : pass-substitution figures for (SSSR)* from the certificate cert2_SSSR_1p_1.json *)
tx[s_, pos_, size_ : 9, opts___] := Text[ptText[s, size, opts], pos];
mx[s_, pos_, size_ : 11, opts___] := Text[ptMath[s, size, opts], pos];
title[s_] := ptText[s, 9.5, Black];
cert = Import[FileNameJoin[{$repoDir, "results", "unbounded_certs", "cert2_SSSR_1p_1.json"}], "RawJSON"];
prog = ptProg[cert["prog"]]; nP = Length[prog]; B0 = cert["B"]; root = First[cert["roots"]];
c54[v_] := Most[v] - Last[v];
phi2[v_] := v . ptPhiM . ptPhiM;
portAt[pid_, m_] := Module[{f = cert["ports"][ToString[pid]], mm, t, h, dt, dh},
   mm = If[Mod[B0, 2] == f["parity"], B0, B0 + 1];
   t = c54[f["ref"][[1]]]; h = c54[f["ref"][[2]]]; dt = c54[f["delta"][[1]]]; dh = c54[f["delta"][[2]]];
   While[mm < m, t = phi2[t] + dt; h = phi2[h] + dh; mm += 2];
   {t, h, f["ref"][[3]]}];
ptype[pid_] := cert["ports"][ToString[pid]]["typ"];
pass[pid_] := cert["passes"][ToString[pid]];
(* standard supertile with owned rhombs *)
super[typ_, m_] := super[typ, m] = Module[{t, rh, nav, q, dmap},
    t = ptSubstituteSlots[ptTri[{ptStdTriangle[typ[[1]], typ[[2]]]}], m];
    rh = ptRhombs[t, "Ghost" -> True]; nav = ptNav[rh]; q = rh["quad"];
    dmap = Association[Flatten[Table[{ptKey1[q[[r, e + 1]]], ptKey1[q[[r, Mod[e + 1, 4] + 1]]]} -> 4 (r - 1) + e + 1, {r, Length[q]}, {e, 0, 3}]]];
    <|"tri" -> t, "rh" -> rh, "nav" -> nav, "xy" -> ptRhombXY[rh], "dmap" -> dmap, "m" -> m, "typ" -> typ|>];
(* run a pass from an entry crossing {tail, head, i}; returns <|"rhombs", "N", "exit" (state, phase)|> *)
runPass[st_, cross_] := Module[{sNext = st["nav"]["sNext"], s, ph, rh, n = 0, t, c, done = False},
   s = st["dmap"][{ptKey1[cross[[2]]], ptKey1[cross[[1]]]}];
   s = 4 Quotient[s - 1, 4] + Mod[Mod[s - 1, 4] + 2, 4] + 1; ph = cross[[3]] + 1;
   rh = Internal`Bag[{Quotient[s - 1, 4] + 1}];
   While[! done,
    c = prog[[Mod[ph, nP] + 1]];
    If[c == 0, n++; t = sNext[[s]]; If[t == 0, done = True, s = t; Internal`StuffBag[rh, Quotient[s - 1, 4] + 1]],
     s = 4 Quotient[s - 1, 4] + Mod[Mod[s - 1, 4] + If[c == 1, 3, 1], 4] + 1];
    If[! done, ph++]];
   <|"rhombs" -> Internal`BagPart[rh, All], "N" -> n, "exit" -> {st["rh"]["quad"][[Quotient[s - 1, 4] + 1, Mod[s - 1, 4] + 1]], st["rh"]["quad"][[Quotient[s - 1, 4] + 1, Mod[Mod[s - 1, 4] + 1, 4] + 1]], Mod[ph, nP]}|>];
rtyp = ptype[root];
Print["root ", root, " type ", rtyp, " children ", pass[root]["children"]];
outline[typ_, m_] := ptPhi^m Append[#, First[#]] &[First[ptTriXY[ptTri[{ptStdTriangle[typ[[1]], typ[[2]]]}]]]];

(* ---------- root pass at levels 7, 9, 11 ---------- *)
runs = Association[Table[m -> Module[{st = super[rtyp, m], r}, r = runPass[st, portAt[root, m]];
      Print["level ", m, ": owned rhombs ", Length[st["xy"]], "  S-moves ", r["N"], "  distinct rhombs ", Length[Union[r["rhombs"]]], "  exit ok ", r["exit"] === portAt[pass[root]["out"], m]]; r], {m, {7, 9, 11}}]];
Module[{panel},
  panel[m_] := Module[{st = super[rtyp, m], r = runs[m], xy, cen},
    xy = st["xy"]; cen = Mean /@ xy;
    Graphics[{{If[m > 9, EdgeForm[None], EdgeForm[{GrayLevel[0.75], AbsoluteThickness[0.15]}]], GrayLevel[0.95], Polygon[xy]},
      {EdgeForm[None], RGBColor[0.35, 0.58, 0.9], Polygon[xy[[Union[r["rhombs"]]]]]},
      {ptRed, AbsoluteThickness[If[m > 9, 0.25, 0.5]], Line[cen[[r["rhombs"]]]]},
      {Black, AbsoluteThickness[1.], Line[outline[rtyp, m]]},
      {ptGreen, PointSize[0.03], Point[cen[[First[r["rhombs"]]]]]}, {Black, EdgeForm[None], Rectangle[cen[[Last[r["rhombs"]]]] - 0.012 ptPhi^m {1, 1}, cen[[Last[r["rhombs"]]]] + 0.012 ptPhi^m {1, 1}]}},
     PlotLabel -> title["レベル " <> ToString[m] <> ": " <> ToString[r["N"]] <> " 回"], PlotRangePadding -> Scaled[0.03]]];
  ptSave["root-pass", GraphicsRow[{panel[7], panel[9], panel[11]}, Spacings -> 4, ImageSize -> 600]]];

(* ---------- decomposition ---------- *)
Module[{segs, st7 = super[rtyp, 7], st9 = super[rtyp, 9], r7 = runs[7], r9 = runs[9], slot1, slot12, seg7, seg9, gkids, pal, g7, g9, childOutline, isRoot},
  (* slot of an owned rhomb = slot chain of its positive triangle *)
  slot1[st_] := st["tri"]["slots"][[1]][[st["rh"]["tri"][[All, 1]]]];
  slot12[st_] := Transpose[{st["tri"]["slots"][[1]][[st["rh"]["tri"][[All, 1]]]], st["tri"]["slots"][[2]][[st["rh"]["tri"][[All, 1]]]]}];
  seg7 = Split[r7["rhombs"], slot1[st7][[#1]] == slot1[st7][[#2]] &];
  seg9 = Split[r9["rhombs"], slot12[st9][[#1]] === slot12[st9][[#2]] &];
  gkids = Flatten[Table[pass[ch[[2]]]["children"][[All, 2]], {ch, pass[root]["children"]}]];
  isRoot = (# == root) & /@ gkids;
  Print["level 7: child segments ", Length[seg7], " lengths ", Length /@ seg7, " slots ", slot1[st7][[First[#]]] & /@ seg7];
  Print["level 9: grandchild segments ", Length[seg9], " lengths ", Length /@ seg9, "  expected ", Length[gkids], " root positions ", Flatten[Position[isRoot, True]]];
  pal = {RGBColor[0.30, 0.55, 0.88], RGBColor[0.2, 0.65, 0.35], RGBColor[0.9, 0.6, 0.1], RGBColor[0.6, 0.35, 0.8], RGBColor[0.1, 0.7, 0.75], RGBColor[0.85, 0.35, 0.55], RGBColor[0.5, 0.5, 0.2], RGBColor[0.2, 0.3, 0.7], RGBColor[0.75, 0.45, 0.25], RGBColor[0.4, 0.7, 0.2], RGBColor[0.55, 0.55, 0.6]};
  childOutline[m_, depth_] := Module[{t = ptSubstitute[ptTri[{ptStdTriangle[rtyp[[1]], rtyp[[2]]]}], depth]}, {Black, AbsoluteThickness[0.6], Line[ptPhi^(m - depth) Append[#, First[#]]] & /@ ptTriXY[t]}];
  g7 = Module[{xy = st7["xy"], cen}, cen = Mean /@ xy;
    Graphics[{{EdgeForm[{GrayLevel[0.75], AbsoluteThickness[0.15]}], GrayLevel[0.96], Polygon[xy]},
      MapIndexed[{EdgeForm[{GrayLevel[0.3], AbsoluteThickness[0.15]}], pal[[Mod[#2[[1]] - 1, Length[pal]] + 1]], Polygon[xy[[Union[#1]]]]} &, seg7],
      {GrayLevel[0.1], AbsoluteThickness[0.5], Line[cen[[r7["rhombs"]]]]},
      MapIndexed[If[Length[#1] >= 20, {Black, Text[ptText[ToString[#2[[1]]], 9, Bold], Mean[cen[[#1]]]]}, {}] &, seg7],
      childOutline[7, 1], {Black, AbsoluteThickness[1.2], Line[outline[rtyp, 7]]}},
     PlotLabel -> title["レベル 7: 子パス " <> ToString[Length[seg7]] <> " 本"], PlotRangePadding -> Scaled[0.03]]];
  g9 = Module[{xy = st9["xy"], cen, k = 0}, cen = Mean /@ xy;
    Graphics[{{EdgeForm[{GrayLevel[0.8], AbsoluteThickness[0.08]}], GrayLevel[0.96], Polygon[xy]},
      MapIndexed[{EdgeForm[None], If[TrueQ[isRoot[[#2[[1]]]]], If[(k++) == 0, RGBColor[0.85, 0.15, 0.15], RGBColor[0.98, 0.6, 0.1]], Blend[{RGBColor[0.55, 0.7, 0.92], RGBColor[0.25, 0.4, 0.75]}, Mod[#2[[1]], 3]/2.]], Polygon[xy[[Union[#1]]]]} &, seg9],
      {GrayLevel[0.1], AbsoluteThickness[0.3], Line[cen[[r9["rhombs"]]]]},
      childOutline[9, 2], {Black, AbsoluteThickness[1.2], Line[outline[rtyp, 9]]}},
     PlotLabel -> title["レベル 9: 孫パス " <> ToString[Length[seg9]] <> " 本"], PlotRangePadding -> Scaled[0.03]]];
  ptSave["pass-decomposition", GraphicsRow[{g7, g9}, Spacings -> 6, ImageSize -> 600]]];

(* ---------- ports: affine law ---------- *)
Module[{ms = Range[7, 21, 2], tri = First[ptTriXY[ptTri[{ptStdTriangle[rtyp[[1]], rtyp[[2]]]}]]], pin, pout, lim, cOf, g1, g2},
  cOf[pid_] := -ptEmbed[c54[cert["ports"][ToString[pid]]["delta"][[1]]]]/ptPhi;
  lim[pid_] := (ptEmbed[portAt[pid, 7][[1]]] - cOf[pid])/ptPhi^7;
  pin = Table[ptEmbed[portAt[root, m][[1]]]/ptPhi^m, {m, ms}];
  pout = Table[ptEmbed[portAt[pass[root]["out"], m][[1]]]/ptPhi^m, {m, ms}];
  Print["limit points: entry ", lim[root], " exit ", lim[pass[root]["out"]], "  |c| entry ", Norm[cOf[root]], " exit ", Norm[cOf[pass[root]["out"]]]];
  Print["deviation check (should equal |c|/phi^m): ", Table[{Norm[pin[[i]] - lim[root]] ptPhi^ms[[i]], Norm[pout[[i]] - lim[pass[root]["out"]]] ptPhi^ms[[i]]}, {i, 3}]];
  g1 = Graphics[{{EdgeForm[{Black, AbsoluteThickness[1.2]}], RGBColor[0.93, 0.96, 1.], Polygon[tri]},
     {ptGreen, PointSize[0.035], Point[lim[root]]}, {Black, PointSize[0.035], Point[lim[pass[root]["out"]]]},
     {ptGreen, AbsoluteThickness[0.8], Line[{pin[[1]], lim[root]}], PointSize[0.018], Point[Take[pin, 2]]}, {Black, AbsoluteThickness[0.8], Line[{pout[[1]], lim[pass[root]["out"]]}], PointSize[0.018], Point[Take[pout, 2]]},
     tx["入口の極限点", lim[root] + {0.2, 0.06}, 8.5, ptGreen], tx["出口の極限点", lim[pass[root]["out"]] + {-0.2, 0.07}, 8.5],
     tx["m = 7", pin[[1]] + {-0.03, -0.06}, 7.5], tx["m = 7", pout[[1]] + {0.12, -0.03}, 7.5]}, PlotRange -> {{-0.55, 1.2}, {-0.1, 1.05}}];
  g2 = ListLogPlot[{Transpose[{ms, Norm /@ ((# - lim[root]) & /@ pin)}], Transpose[{ms, Norm /@ ((# - lim[pass[root]["out"]]) & /@ pout)}]}, Joined -> True, PlotMarkers -> {Automatic, 5}, PlotStyle -> {ptGreen, Black}, Frame -> True,
    FrameLabel -> {ptText["レベル m", 9], ptText["極限点からの距離", 9]}, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8, Black}, GridLines -> Automatic, GridLinesStyle -> GrayLevel[0.9], PlotRange -> All, AspectRatio -> 0.75];
  ptSave["ports", GraphicsRow[{g1, g2}, Spacings -> 10, ImageSize -> 600]]];

(* ---------- growth of pass lengths (all certificates) ---------- *)
Module[{sm = Import[FileNameJoin[{$repoDir, "results", "unbounded_certs", "summary.json"}], "RawJSON"], data, skip = {"SSRSR", "SLSRSLSR", "SSSRSSSR"}},
  sm = Select[sm, ! MemberQ[skip, #["prog"]] &];
  data = Table[With[{n0 = First[r["N"]]}, Table[{r["m0"] + 2 k, Nest[r["a"] # + r["b"] &, n0, k]}, {k, 0, 5}]], {r, sm}];
  ptSave["growth", ListLogPlot[data, Joined -> True, PlotMarkers -> {Automatic, 4}, PlotStyle -> AbsoluteThickness[0.7], Frame -> True, FrameLabel -> {ptText["レベル m", 9], ptText["パスの長さ N(m)", 9]},
    Epilog -> {GrayLevel[0.3], Dashed, AbsoluteThickness[1], Line[{{7, Log[40.]}, {19, Log[40. 2^12]}}], Text[ptText["傾き: 1 レベルあたり 2 倍", 8], {15.5, Log[1500.]}]},
    PlotRange -> All, ImageSize -> 430, LabelStyle -> {FontFamily -> ptFont, FontSize -> 8, Black}, GridLines -> Automatic, GridLinesStyle -> GrayLevel[0.9]]]];

(* ---------- exit lemma schematic ---------- *)
Module[{t0 = ptStdTriangle[1, 1], Z, kids, ch, U, W, y, d, r, rp, c},
  Z = 6 ptPhi (ptEmbed /@ Rest[t0]);
  kids = 6 ptTriXY[ptSubdivide[ptTri[{t0}]]];
  ch = kids[[2]];                         (* child 1 = obtuse (Q, R, B'): its base R B' lies on the base of the parent *)
  U = ch[[2]]; W = ch[[3]]; y = U + 0.46 (W - U); d = Min[Norm[y - U], Norm[y - W]]; r = d Sin[36 Degree];
  c = y + 0.5 Normalize[W - U];
  rp = With[{p = c, u = Normalize[W - U], v = RotationMatrix[72 Degree] . Normalize[W - U]}, {p, p + u, p + u + v, p + v}];
  ptSave["exit-lemma", Graphics[{
     {EdgeForm[{GrayLevel[0.5], AbsoluteThickness[0.5]}], GrayLevel[0.95], Polygon[kids]},
     {EdgeForm[{GrayLevel[0.3], AbsoluteThickness[0.6]}], RGBColor[0.78, 0.87, 0.98], Polygon[ch]},
     {EdgeForm[{Black, AbsoluteThickness[1.6]}], FaceForm[None], Polygon[Z]},
     {EdgeForm[{ptGreen, AbsoluteThickness[1.], Dashed}], RGBColor[0.3, 0.75, 0.4, 0.12], Disk[y, r]},
     {EdgeForm[{ptRed, AbsoluteThickness[1.]}], RGBColor[0.95, 0.5, 0.45, 0.6], Polygon[rp]},
     {Black, PointSize[0.014], Point[{y, U, W, c}]},
     mx["Z", Z[[3]] + {-0.45, -0.2}, 13], tx["子の超三角形", Mean[ch] + {0.9, -0.15}, 9], mx["E", (y + W)/2 + 1.6 Normalize[W - U] + {0.25, 0.28}, 12],
     mx[Subscript["y", "n"], y + {-0.38, -0.3}, 12], mx[Subscript["a", "n"], c + {-0.1, -0.42}, 12], mx["r'", Mean[rp] + {0.85, 0.1}, 12, ptRed], mx["D", y + r {Cos[1.9], Sin[1.9]} + {0.05, 0.35}, 12, ptGreen],
     mx["U", U + {0.25, 0.3}, 11], mx["W", W + {0.35, 0.0}, 11],
     {GrayLevel[0.2], AbsoluteThickness[0.7], Arrowheads[{-0.02, 0.02}], Arrow[{y + 0.45 {0.588, 0.809}, U + 0.45 {0.588, 0.809}}], mx["≥ d", (y + U)/2 + 0.85 {0.588, 0.809}, 11]}},
    PlotRange -> {{Min[Z[[All, 1]]] - 0.6, Max[Z[[All, 1]]] + 1.6}, {-2.2, Max[Z[[All, 2]]] + 0.5}}, ImageSize -> 480]]];
