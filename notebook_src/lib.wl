(* lib.wl -- exact Penrose rhomb machinery in Z[zeta], zeta = Exp[2 Pi I/5].
   A point is an integer 4-vector {a,b,c,d} = a + b zeta + c zeta^2 + d zeta^3 (row-vector convention). *)

ptZ = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {-1, -1, -1, -1}};   (* v.ptZ = zeta v *)
ptPhiM = -(ptZ . ptZ) - (ptZ . ptZ . ptZ);                              (* phi = -zeta^2 - zeta^3 *)
ptRhoM = -(ptZ . ptZ . ptZ);                                            (* rho = -zeta^3 = Exp[I Pi/5] *)
ptRhoPow = Table[MatrixPower[ptRhoM, k], {k, 0, 9}];
ptPhiPow[n_] := MatrixPower[ptPhiM, n];
ptE4 = N[Table[{Cos[2 Pi j/5], Sin[2 Pi j/5]}, {j, 0, 3}]];
ptEmbed[v_] := v . ptE4;
ptUnit[k_] := {1, 0, 0, 0} . ptRhoPow[[Mod[k, 10] + 1]];
ptUnits = Table[ptUnit[k], {k, 0, 9}];
ptPhi = N[GoldenRatio];

(* ---- vertex keys ---- *)
ptKey[v_] := Module[{w = v + 16384}, ((w[[All, 1]]*32768 + w[[All, 2]])*32768 + w[[All, 3]])*32768 + w[[All, 4]]];
ptKey1[v_] := Module[{w = v + 16384}, ((w[[1]]*32768 + w[[2]])*32768 + w[[3]])*32768 + w[[4]]];

(* ---- triangles: <|"col","A","B","C"|> ; col 0 acute, 1 obtuse ---- *)
ptTri[list_] := <|"col" -> list[[All, 1]], "A" -> list[[All, 2]], "B" -> list[[All, 3]], "C" -> list[[All, 4]]|>;
ptStdTriangle[col_, sign_] := Module[{b = ptUnit[0], c = If[col == 0, ptUnit[1], ptUnit[3]]},
   If[sign > 0, {col, {0, 0, 0, 0}, b, c}, {col, {0, 0, 0, 0}, c, b}]];
ptWheel[] := ptTri[Table[
    With[{b = ptUnit[i], c = ptUnit[i + 1]}, If[EvenQ[i], {0, {0, 0, 0, 0}, c, b}, {0, {0, 0, 0, 0}, b, c}]], {i, 0, 9}]];
ptNTri[t_] := Length[t["col"]];
ptSign[t_] := Module[{a = ptEmbed[t["A"]], b = ptEmbed[t["B"]], c = ptEmbed[t["C"]]},
   Sign[(b[[All, 1]] - a[[All, 1]]) (c[[All, 2]] - a[[All, 2]]) - (b[[All, 2]] - a[[All, 2]]) (c[[All, 1]] - a[[All, 1]])]];

(* one substitution step; returns triangles with "par" (parent index) and "slot" *)
ptSubdivide[t_] := Module[{col = t["col"], ia, io, A, B, C, A2, B2, C2, P, Q, R, res, e4 = {{0, 0, 0, 0}}},
   ia = Flatten[Position[col, 0]]; io = Flatten[Position[col, 1]];
   res = {};
   If[ia =!= {},
    A = t["A"][[ia]]; B = t["B"][[ia]]; C = t["C"][[ia]];
    A2 = A . ptPhiM; B2 = B . ptPhiM; C2 = C . ptPhiM; P = A2 + (B - A);
    AppendTo[res, {ConstantArray[0, Length[ia]], C2, P, B2, ia, ConstantArray[0, Length[ia]]}];
    AppendTo[res, {ConstantArray[1, Length[ia]], P, C2, A2, ia, ConstantArray[1, Length[ia]]}]];
   If[io =!= {},
    A = t["A"][[io]]; B = t["B"][[io]]; C = t["C"][[io]];
    A2 = A . ptPhiM; B2 = B . ptPhiM; C2 = C . ptPhiM; Q = B2 + (A - B); R = B2 + (C - B);
    AppendTo[res, {ConstantArray[1, Length[io]], R, C2, A2, io, ConstantArray[0, Length[io]]}];
    AppendTo[res, {ConstantArray[1, Length[io]], Q, R, B2, io, ConstantArray[1, Length[io]]}];
    AppendTo[res, {ConstantArray[0, Length[io]], R, Q, A2, io, ConstantArray[2, Length[io]]}]];
   <|"col" -> Developer`ToPackedArray[Join @@ res[[All, 1]]], "A" -> Developer`ToPackedArray[Join @@ res[[All, 2]]],
    "B" -> Developer`ToPackedArray[Join @@ res[[All, 3]]], "C" -> Developer`ToPackedArray[Join @@ res[[All, 4]]],
    "par" -> Developer`ToPackedArray[Join @@ res[[All, 5]]], "slot" -> Developer`ToPackedArray[Join @@ res[[All, 6]]]|>];

(* n steps, tracking the ancestor index (into the original list) *)
ptSubstitute[t_, n_] := Module[{cur = t, anc = Range[ptNTri[t]], nx},
   Do[nx = ptSubdivide[cur]; anc = anc[[nx["par"]]]; cur = nx, {n}];
   Append[KeyDrop[cur, {"par", "slot"}], "anc" -> anc]];
(* n >= 1 steps, tracking for every fine triangle the chain of slots (list of n slot arrays, coarsest first) *)
ptSubstituteSlots[t_, n_] := Module[{cur = t, anc = Range[ptNTri[t]], nx, slots = {}},
   Do[nx = ptSubdivide[cur]; anc = anc[[nx["par"]]]; slots = Append[(#[[nx["par"]]] &) /@ slots, nx["slot"]]; cur = nx, {n}];
   Join[KeyDrop[cur, {"par", "slot"}], <|"anc" -> anc, "slots" -> slots|>]];

(* ---- rhombs ---- *)
(* returns <|"quad" (n x 4 x 4 int, CCW), "fat", "tri" (n x 2: positive, negative triangle index; 0 = ghost), "vid" ...|>
   ghost -> True : unpaired POSITIVE triangles are completed by their mirror image (ownership convention) *)
ptRhombs[t_, OptionsPattern[{"Ghost" -> False}]] := Module[{n = ptNTri[t], kB, kC, lo, hi, sg, ord, same, f, p1, p2, pos, neg, bad, quad, tri, single, Ap, Bp, Cp},
   kB = ptKey[t["B"]]; kC = ptKey[t["C"]];
   lo = MapThread[Min, {kB, kC}]; hi = MapThread[Max, {kB, kC}];
   sg = ptSign[t];
   ord = Ordering[Transpose[{lo, hi}]];
   same = With[{los = lo[[ord]], his = hi[[ord]]},
     Flatten[Position[Unitize[Most[los] - Rest[los]] + Unitize[Most[his] - Rest[his]], 0]]];
   p1 = ord[[same]]; p2 = ord[[same + 1]];
   bad = Count[sg[[p1]] + sg[[p2]], Except[0]];
   If[bad > 0, Print["ptRhombs: ", bad, " base pairs are not mirror images!"]];
   pos = MapThread[If[sg[[#1]] > 0, #1, #2] &, {p1, p2}]; neg = p1 + p2 - pos;
   quad = Transpose[{t["A"][[pos]], t["B"][[pos]], t["A"][[neg]], t["C"][[pos]]}];
   tri = Transpose[{pos, neg}];
   If[TrueQ[OptionValue["Ghost"]],
    single = Complement[Flatten[Position[sg, 1]], pos];
    If[single =!= {},
     Ap = t["A"][[single]]; Bp = t["B"][[single]]; Cp = t["C"][[single]];
     quad = Join[quad, Transpose[{Ap, Bp, Bp + Cp - Ap, Cp}]];
     tri = Join[tri, Transpose[{single, ConstantArray[0, Length[single]]}]]]];
   <|"quad" -> Developer`ToPackedArray[quad], "fat" -> t["col"][[tri[[All, 1]]]], "tri" -> tri|>];
ptRhombXY[rh_] := Map[ptEmbed, rh["quad"]];              (* n x 4 x 2 *)
ptRhombCent[rh_] := Mean /@ ptRhombXY[rh];

(* ---- navigation: state s = 4 (r-1) + e + 1, e in 0..3 ; sNext[s] = state after S (0 = no neighbour) ---- *)
ptNav[rh_] := Module[{q = rh["quad"], n, vk, u, vid, V, ek, rk, assoc, nb, r, e},
   n = Length[q];
   vk = ptKey[Flatten[q, 1]]; u = Union[vk];
   vid = Partition[Lookup[AssociationThread[u -> Range[Length[u]]], vk], 4]; V = Length[u] + 1;
   ek = Flatten[vid*V + RotateLeft[vid, {0, 1}]];      (* directed edge q_e -> q_{e+1}, in state order *)
   rk = Flatten[RotateLeft[vid, {0, 1}]*V + vid];
   If[Length[Union[ek]] != Length[ek], Print["ptNav: duplicated directed edges (overlapping rhombs)!"]];
   assoc = AssociationThread[ek -> Range[4 n]];
   nb = Lookup[assoc, rk, 0];
   (* neighbour state index nb = 4 (r'-1) + e' + 1  -> new heading e'+2 *)
   <|"sNext" -> Developer`ToPackedArray[If[# == 0, 0, 4 Quotient[# - 1, 4] + Mod[Mod[# - 1, 4] + 2, 4] + 1] & /@ nb],
     "nbr" -> Developer`ToPackedArray[nb], "vid" -> vid, "nv" -> Length[u]|>];

ptProg[s_String] := Characters[s] /. {"S" -> 0, "R" -> 1, "L" -> 2};

(* single run: {flag, cycles}; flag 0 closed, 1 escaped, 2 timeout *)
ptRunC = Compile[{{sNext, _Integer, 1}, {prog, _Integer, 1}, {s0, _Integer}, {maxCyc, _Integer}},
   Module[{s = s0, cyc = 0, np = Length[prog], flag = 2, j = 1, c = 0, t = 0},
    While[cyc < maxCyc && flag == 2,
     j = 1;
     While[j <= np,
      c = prog[[j]];
      If[c == 0,
       t = sNext[[s]]; If[t == 0, flag = 1; j = np + 1, s = t],
       s = 4 Quotient[s - 1, 4] + Mod[Mod[s - 1, 4] + If[c == 1, 3, 1], 4] + 1];
      j++];
     cyc++;
     If[flag == 2 && s == s0, flag = 0]];
    {flag, cyc}], RuntimeOptions -> "Speed"];

(* scan of all given start states with orbit sharing.  returns n x 3 array {flag, period, orbitId} for the starts *)
ptScanC = Compile[{{sNext, _Integer, 1}, {prog, _Integer, 1}, {starts, _Integer, 1}, {maxCyc, _Integer}},
   Module[{ns = Length[sNext], np = Length[prog], fl, per, oid, buf, nb = 0, s = 0, s0 = 0, cyc = 0, flag = 2, j = 1, c = 0, t = 0, k = 0, hit = 0, out},
    fl = Table[-1, {ns}]; per = Table[0, {ns}]; oid = Table[0, {ns}];
    buf = Table[0, {maxCyc + 1}];
    Do[
     s0 = starts[[i]];
     If[fl[[s0]] < 0,
      s = s0; cyc = 0; flag = 2; nb = 0; hit = 0;
      While[cyc < maxCyc && flag == 2,
       nb++; buf[[nb]] = s;
       j = 1;
       While[j <= np,
        c = prog[[j]];
        If[c == 0,
         t = sNext[[s]]; If[t == 0, flag = 1; j = np + 1, s = t],
         s = 4 Quotient[s - 1, 4] + Mod[Mod[s - 1, 4] + If[c == 1, 3, 1], 4] + 1];
        j++];
       cyc++;
       If[flag == 2,
        If[s == s0, flag = 0,
         If[fl[[s]] >= 0, hit = s; flag = 3]]]];
      Which[
       flag == 0, Do[fl[[buf[[k]]]] = 0; per[[buf[[k]]]] = cyc; oid[[buf[[k]]]] = s0, {k, 1, nb}],
       flag == 3, Do[fl[[buf[[k]]]] = fl[[hit]]; per[[buf[[k]]]] = per[[hit]]; oid[[buf[[k]]]] = oid[[hit]], {k, 1, nb}],
       True, Do[fl[[buf[[k]]]] = flag; per[[buf[[k]]]] = 0; oid[[buf[[k]]]] = s0, {k, 1, nb}]]],
     {i, Length[starts]}];
    out = Table[{fl[[starts[[i]]]], per[[starts[[i]]]], oid[[starts[[i]]]]}, {i, Length[starts]}];
    out], RuntimeOptions -> "Speed"];

(* recorded run (uncompiled): returns <|"flag","cycles","states" (state before every command), "rhombs" (sequence of rhombs entered, starting with the first)|> *)
ptRun[sNext_, prog_, s0_, maxCyc_ : 100000] := Module[{s = s0, cyc = 0, flag = 2, np = Length[prog], states, rh, t, c, j},
   states = Internal`Bag[]; rh = Internal`Bag[{Quotient[s0 - 1, 4] + 1}];
   While[cyc < maxCyc && flag == 2,
    j = 1;
    While[j <= np,
     c = prog[[j]]; Internal`StuffBag[states, s];
     If[c == 0,
      t = sNext[[s]]; If[t == 0, flag = 1; j = np + 1, s = t; Internal`StuffBag[rh, Quotient[s - 1, 4] + 1]],
      s = 4 Quotient[s - 1, 4] + Mod[Mod[s - 1, 4] + If[c == 1, 3, 1], 4] + 1];
     j++];
    cyc++;
    If[flag == 2 && s == s0, flag = 0]];
   <|"flag" -> flag, "cycles" -> cyc, "states" -> Internal`BagPart[states, All], "rhombs" -> Internal`BagPart[rh, All], "last" -> s|>];

(* isometry-invariant signature of a tile set *)
ptSignature[cent_, cells_] := Module[{cs = cent[[cells]], m}, m = Mean[cs]; Sort[Round[100 (Norm /@ ((# - m) & /@ cs))]]];

(* ---- vertex -> triangles incidence and coronas (marked triangle level) ---- *)
ptAngUnits = {{1, 2, 2}, {3, 1, 1}};   (* units of 36 deg at A,B,C for acute / obtuse *)
ptTriIncidence[t_] := Module[{n = ptNTri[t], k, u, vid, inc, ang},
   k = Join[ptKey[t["A"]], ptKey[t["B"]], ptKey[t["C"]]];
   u = Union[k]; vid = Lookup[AssociationThread[u -> Range[Length[u]]], k];
   ang = Join[ptAngUnits[[t["col"] + 1, 1]], ptAngUnits[[t["col"] + 1, 2]], ptAngUnits[[t["col"] + 1, 3]]];
   inc = GroupBy[Transpose[{vid, Join[Range[n], Range[n], Range[n]], ang}], First -> Rest];
   <|"vid" -> Transpose[Partition[vid, n]], "inc" -> inc,
     "full" -> Association[KeyValueMap[#1 -> (Total[#2[[All, 2]]] == 10) &, inc]]|>];
ptCorona[inc_, i_] := Module[{vs = inc["vid"][[i]]},
   If[And @@ (inc["full"] /@ vs), Union @@ (inc["inc"][#][[All, 1]] & /@ vs), None]];
(* canonical form in the frame of the centre triangle: A -> 0, B - A -> 1 *)
ptCoronaCanon[t_, idx_, c_] := Module[{A0 = t["A"][[c]], d, k, M},
   d = t["B"][[c]] - A0; k = First[FirstPosition[ptUnits, d]] - 1; M = ptRhoPow[[Mod[-k, 10] + 1]];
   Sort[Table[Flatten[{t["col"][[j]], (t["A"][[j]] - A0) . M, (t["B"][[j]] - A0) . M, (t["C"][[j]] - A0) . M, Boole[j == c]}], {j, idx}]]];
ptCanonToTri[canon_] := ptTri[{#[[1]], #[[2 ;; 5]], #[[6 ;; 9]], #[[10 ;; 13]]} & /@ canon];
ptCanonCentre[canon_] := First[FirstPosition[canon[[All, 14]], 1]];
ptAtlas[t_] := Module[{inc = ptTriIncidence[t], atlas = <||>, cor, key},
   Do[cor = ptCorona[inc, i];
    If[cor =!= None, key = ptCoronaCanon[t, cor, i]; atlas[key] = Lookup[atlas, Key[key], 0] + 1], {i, ptNTri[t]}];
   atlas];

(* ---- drawing helpers ---- *)
ptColFat = RGBColor[0.985, 0.975, 0.955]; ptColThin = RGBColor[0.86, 0.87, 0.89];
ptEdgeCol = GrayLevel[0.55];
ptColAcute = RGBColor[0.98, 0.80, 0.62]; ptColObtuse = RGBColor[0.70, 0.82, 0.95];
ptBlue = RGBColor[0.30, 0.55, 0.88]; ptRed = RGBColor[0.82, 0.20, 0.17]; ptGreen = RGBColor[0.15, 0.6, 0.3];
ptDrawRhombs[rh_, sel_ : All, opts___] := Module[{xy = ptRhombXY[rh], fat = rh["fat"], idx},
   idx = If[sel === All, Range[Length[xy]], sel];
   {EdgeForm[{ptEdgeCol, AbsoluteThickness[0.5]}], opts,
    {ptColFat, Polygon[xy[[Select[idx, fat[[#]] == 1 &]]]]}, {ptColThin, Polygon[xy[[Select[idx, fat[[#]] == 0 &]]]]}}];
ptTriXY[t_] := Transpose[{ptEmbed[t["A"]], ptEmbed[t["B"]], ptEmbed[t["C"]]}];
ptSave[name_, g_, w_ : 1800] := Module[{img = Rasterize[g, RasterSize -> w, Background -> White]},
   Export[FileNameJoin[{$figDir, name <> ".png"}], img];
   Print["saved ", name, " ", ImageDimensions[img]]; img];
ptFont = "Yu Gothic UI";
ptLabel[s_, pos_, size_ : 12, opts___] := Text[Style[s, FontFamily -> ptFont, FontSize -> size, opts], pos];

(* ---- pentagrid (de Bruijn) generator; same conventions as GeneratePenroseRhombs ---- *)
ptEV5 = N[Table[{Cos[2 Pi j/5], Sin[2 Pi j/5]}, {j, 0, 4}]];
ptOffsets = {0, 0.1, -0.1, 0.2, -0.2};
ptPentagrid[{xmin_, ymin_}, {w_, h_}, off_ : ptOffsets] := Module[{xmax = xmin + w, ymax = ymin + h, corners, proj, lo, hi, parts, K, r, s, id5 = IdentityMatrix[5], er, es, k5, quad, lines},
   corners = {{xmin, ymin}, {xmax, ymin}, {xmin, ymax}, {xmax, ymax}};
   proj = corners . Transpose[ptEV5] - ConstantArray[off, 4];
   lo = Floor[Min /@ Transpose[proj]] - 1; hi = Ceiling[Max /@ Transpose[proj]] + 1;
   parts = Flatten[Table[
      Module[{M, Nn, Ainv, P, in, KK, rr, ss},
       {M, Nn} = Transpose[Tuples[{Range[lo[[i]], hi[[i]]], Range[lo[[j]], hi[[j]]]}]];
       Ainv = Inverse[{ptEV5[[i]], ptEV5[[j]]}];
       P = Transpose[Ainv . {M + off[[i]], Nn + off[[j]]}];
       in = Flatten[Position[UnitStep[P[[All, 1]] - xmin] UnitStep[xmax - P[[All, 1]]] UnitStep[P[[All, 2]] - ymin] UnitStep[ymax - P[[All, 2]]], 1]];
       P = P[[in]]; M = M[[in]]; Nn = Nn[[in]];
       KK = Ceiling[P . Transpose[ptEV5] - ConstantArray[off, Length[P]]];
       KK[[All, i]] = M; KK[[All, j]] = Nn;
       {rr, ss} = If[MemberQ[{3, 4}, j - i], {j, i}, {i, j}];
       {KK, ConstantArray[rr, Length[KK]], ConstantArray[ss, Length[KK]], Transpose[{ConstantArray[i, Length[KK]], M, ConstantArray[j, Length[KK]], Nn}]}],
      {i, 1, 4}, {j, i + 1, 5}], 1];
   K = Join @@ parts[[All, 1]]; r = Join @@ parts[[All, 2]]; s = Join @@ parts[[All, 3]]; lines = Join @@ parts[[All, 4]];
   er = id5[[r]]; es = id5[[s]];
   k5 = Transpose[{K, K + er, K + er + es, K + es}];            (* n x 4 x 5 *)
   quad = k5[[All, All, 1 ;; 4]] - k5[[All, All, {5, 5, 5, 5}]];
   <|"quad" -> Developer`ToPackedArray[quad], "fat" -> 1 - Unitize[Mod[s - r, 5] - 1], "lines" -> lines, "rs" -> Transpose[{r, s}]|>];

(* ---- turtle state geometry ---- *)
ptStateRhomb[s_] := Quotient[s - 1, 4] + 1;
ptStateEdge[s_] := Mod[s - 1, 4];
(* directed edge (tail, head) of a state as 4-vectors / as xy *)
ptStateDE[rh_, s_] := With[{q = rh["quad"][[ptStateRhomb[s]]], e = ptStateEdge[s]}, {q[[e + 1]], q[[Mod[e + 1, 4] + 1]]}];
ptStateDExy[rh_, s_] := ptEmbed /@ ptStateDE[rh, s];
ptApply[nav_, s_, c_] := Which[c == 0, nav["sNext"][[s]], c == 1, 4 Quotient[s - 1, 4] + Mod[Mod[s - 1, 4] + 3, 4] + 1, True, 4 Quotient[s - 1, 4] + Mod[Mod[s - 1, 4] + 1, 4] + 1];

(* states at given cycle checkpoints (sorted ascending); 0 after escape; also returns {flag, cycles} in the last two slots *)
ptRunCheckC = Compile[{{sNext, _Integer, 1}, {prog, _Integer, 1}, {s0, _Integer}, {chk, _Integer, 1}},
   Module[{s = s0, cyc = 0, np = Length[prog], flag = 2, j = 1, c = 0, t = 0, nc = Length[chk], out, ic = 1, maxCyc = 0},
    out = Table[0, {nc + 2}]; maxCyc = chk[[nc]];
    While[cyc < maxCyc && flag == 2,
     j = 1;
     While[j <= np,
      c = prog[[j]];
      If[c == 0,
       t = sNext[[s]]; If[t == 0, flag = 1; j = np + 1, s = t],
       s = 4 Quotient[s - 1, 4] + Mod[Mod[s - 1, 4] + If[c == 1, 3, 1], 4] + 1];
      j++];
     cyc++;
     If[flag == 2,
      If[ic <= nc && cyc == chk[[ic]], out[[ic]] = s; ic++];
      If[s == s0, flag = 0]]];
    out[[nc + 1]] = flag; out[[nc + 2]] = cyc;
    out], RuntimeOptions -> "Speed"];

(* classes of closed orbits in a scan: list of <|period, ncells, count(starts), rep state, radius, cells|> *)
ptClasses[nav_, cent_, prog_, starts_, res_] := Module[{ok, reps, cls = <||>},
   ok = Flatten[Position[res[[All, 1]], 0]];
   reps = DeleteDuplicates[res[[ok, 3]]];
   Do[Module[{run = ptRun[nav["sNext"], prog, s0], cells, sig, key, m, rad},
     cells = DeleteDuplicates[run["rhombs"]];
     sig = ptSignature[cent, cells]; key = {run["cycles"], Length[cells], sig};
     m = Mean[cent[[cells]]]; rad = Max[Norm /@ ((# - m) & /@ cent[[cells]])];
     If[KeyExistsQ[cls, key], cls[key] = MapAt[# + 1 &, cls[key], Key["norbits"]],
      cls[key] = <|"period" -> run["cycles"], "ncells" -> Length[cells], "norbits" -> 1, "rep" -> s0, "radius" -> rad, "cells" -> cells, "centre" -> m|>]], {s0, reps}];
   SortBy[Values[cls], {#period, #ncells} &]];

(* draw an orbit panel: background tiles in the bounding box, swept tiles highlighted.
   Figures are composed at ImageSize ~600 (fonts 8-11 pt) and rasterised to 1800 px. *)
ptPanel[rh_, xy_, cent_, cells_, OptionsPattern[{"Margin" -> 1.3, "Path" -> None, "Label" -> None, "Color" -> Automatic, "Extra" -> {}, "Square" -> True, "Edges" -> True, "LabelSize" -> 9}]] := Module[{pts = Flatten[xy[[cells]], 1], x0, x1, y0, y1, bg, m = OptionValue["Margin"], col, hw},
   {x0, x1} = MinMax[pts[[All, 1]]] + {-m, m}; {y0, y1} = MinMax[pts[[All, 2]]] + {-m, m};
   If[TrueQ[OptionValue["Square"]], hw = Max[x1 - x0, y1 - y0]/2; {x0, x1} = Mean[{x0, x1}] + {-hw, hw}; {y0, y1} = Mean[{y0, y1}] + {-hw, hw}];
   bg = Flatten[Position[UnitStep[cent[[All, 1]] - x0 + 1.5] UnitStep[x1 + 1.5 - cent[[All, 1]]] UnitStep[cent[[All, 2]] - y0 + 1.5] UnitStep[y1 + 1.5 - cent[[All, 2]]], 1]];
   col = If[OptionValue["Color"] === Automatic, ptBlue, OptionValue["Color"]];
   Graphics[{
     {EdgeForm[{GrayLevel[0.72], AbsoluteThickness[0.3]}], {ptColFat, Polygon[xy[[Select[bg, rh["fat"][[#]] == 1 &]]]]}, {ptColThin, Polygon[xy[[Select[bg, rh["fat"][[#]] == 0 &]]]]}},
     {If[TrueQ[OptionValue["Edges"]], EdgeForm[{GrayLevel[0.15], AbsoluteThickness[0.4]}], EdgeForm[None]], col, Opacity[0.85], Polygon[xy[[cells]]]},
     If[OptionValue["Path"] =!= None, {ptRed, AbsoluteThickness[0.7], Line[OptionValue["Path"]]}, {}],
     OptionValue["Extra"]},
    PlotRange -> {{x0, x1}, {y0, y1}}, PlotRangePadding -> 0, PlotRangeClipping -> True,
    PlotLabel -> If[OptionValue["Label"] === None, None, Style[OptionValue["Label"], FontFamily -> ptFont, FontSize -> OptionValue["LabelSize"], Black]]]];
ptMath[s_, size_ : 11, opts___] := Style[s, FontFamily -> "Times New Roman", Italic, FontSize -> size, opts];
ptText[s_, size_ : 9, opts___] := Style[s, FontFamily -> ptFont, FontSize -> size, opts];

(* rhombs entered along a closed orbit of known period (sequence, one entry per S move) *)
ptOrbitRhombsC = Compile[{{sNext, _Integer, 1}, {prog, _Integer, 1}, {s0, _Integer}, {cycles, _Integer}},
   Module[{s = s0, np = Length[prog], nS = 0, out, k = 0, c = 0, t = 0},
    Do[If[prog[[j]] == 0, nS++], {j, np}];
    out = Table[0, {nS cycles}];
    Do[
     Do[c = prog[[j]];
      If[c == 0, t = sNext[[s]]; If[t > 0, s = t]; k++; out[[k]] = Quotient[s - 1, 4] + 1,
       s = 4 Quotient[s - 1, 4] + Mod[Mod[s - 1, 4] + If[c == 1, 3, 1], 4] + 1], {j, np}], {cyc, cycles}];
    out], RuntimeOptions -> "Speed"];
(* per-orbit data for all closed orbits of a scan: list of <|period, ncells, rep, centre, radius, sig|> *)
ptOrbitData[sNext_, cent_, prog_, res_, minPeriod_ : 1] := Module[{ok, reps},
   ok = Select[res, #[[1]] == 0 && #[[2]] >= minPeriod &];
   reps = DeleteDuplicatesBy[ok, #[[3]] &];
   Table[Module[{cells = DeleteDuplicates[ptOrbitRhombsC[sNext, prog, r[[3]], r[[2]]]], cs, m},
     cs = cent[[cells]]; m = Mean[cs];
     <|"period" -> r[[2]], "ncells" -> Length[cells], "rep" -> r[[3]], "centre" -> m, "radius" -> Max[Sqrt[Total[(Transpose[cs] - m)^2]]],
       "sig" -> Hash[Sort[Round[100 Sqrt[Total[(Transpose[cs] - m)^2]]]]]|>], {r, reps}]];
ptClassesOf[orb_] := SortBy[Values[GroupBy[orb, {#period, #ncells, #sig} &, Append[First[#], "norbits" -> Length[#]] &]], {#period, #ncells} &];

(* ---- exact orbit classes: congruence of swept tile sets under translations, rotations by 36 deg and reflections ---- *)
ptConjM = {{1, 0, 0, 0}, {-1, -1, -1, -1}, {0, 0, 0, 1}, {0, 0, 1, 0}};      (* complex conjugation on Z[zeta] *)
ptSym20 = Join[ptRhoPow, (ptConjM . #) & /@ ptRhoPow];
ptTransCanon[q_] := Module[{mn = First[Sort[Flatten[q, 1]]]}, Sort[Sort /@ (q - ConstantArray[mn, {Length[q], 4}])]];
ptD10Canon[q_] := First[Sort[Table[ptTransCanon[q . M], {M, ptSym20}]]];
(* exact classification of the closed orbits of a scan.  returns list of <|period, ncells, norbits, ntrans, rep, radius, centre|> *)
ptExactClasses[quad_, sNext_, cent_, prog_, res_] := Module[{ok, reps, tr = <||>, cls = <||>},
   ok = Select[res, #[[1]] == 0 &];
   reps = DeleteDuplicatesBy[ok, #[[3]] &];
   (* translation classes *)
   Do[Module[{cells = DeleteDuplicates[ptOrbitRhombsC[sNext, prog, r[[3]], r[[2]]]], key},
     key = {r[[2]], ptTransCanon[quad[[cells]]]};   (* the canonical form itself is the key (no hashing) *)
     If[KeyExistsQ[tr, key], tr[key] = {tr[key][[1]] + 1, tr[key][[2]]}, tr[key] = {1, r[[3]]}]], {r, reps}];
   (* D10 classes of the translation classes *)
   KeyValueMap[Function[{key, val}, Module[{cells = DeleteDuplicates[ptOrbitRhombsC[sNext, prog, val[[2]], key[[1]]]], k2, cs, m},
       k2 = {key[[1]], ptD10Canon[quad[[cells]]]};
       cs = cent[[cells]]; m = Mean[cs];
       If[KeyExistsQ[cls, k2], cls[k2] = Join[cls[k2], <|"norbits" -> cls[k2]["norbits"] + val[[1]], "ntrans" -> cls[k2]["ntrans"] + 1|>],
        cls[k2] = <|"period" -> key[[1]], "ncells" -> Length[cells], "norbits" -> val[[1]], "ntrans" -> 1, "rep" -> val[[2]], "centre" -> m, "radius" -> Max[Sqrt[Total[(Transpose[cs] - m)^2]]]|>]]], tr];
   SortBy[Values[cls], {#period, #ncells} &]];
