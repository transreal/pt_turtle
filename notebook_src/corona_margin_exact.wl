(* corona_margin_exact.wl : exact value of delta_C = distance from the centre triangle t_C to the boundary of the
   union of the corona C, for the 30 coronas of the atlas, in exact arithmetic without any numerical approximation:
   vertices are elements of Z[zeta] embedded with the exact Cos/Sin of multiples of 72 degrees, boundary edges are
   found by exact matching of vertex keys, and every quantity that is compared (squared distances, dot products,
   segment parameters) lies in Q(sqrt 5); it is written as a + b sqrt 5 with rational a, b (q5) and its sign is
   decided with rational arithmetic only (q5Sign).  For each corona the script checks exactly that every candidate
   squared distance is >= (5 - sqrt 5)/8 = sin^2 36 and that at least one candidate is equal to it.
   usage: run.wl corona_margin_exact.wl *)
keys = Import[FileNameJoin[{$dataDir, "atlas.wxf"}]];
e4x = Table[{Cos[2 Pi k/5], Sin[2 Pi k/5]}, {k, 0, 3}];                       (* exact basis 1, zeta, zeta^2, zeta^3 *)
embX[v_] := v . e4x;
q5[x_] := Module[{r = RootReduce[x], e, a, b},
   e = Expand[r]; b = Coefficient[e, Sqrt[5]]; a = e /. Sqrt[5] -> 0;
   If[! (Element[a, Rationals] && Element[b, Rationals] && RootReduce[r - a - b Sqrt[5]] === 0),
      Print["not an element of Q(sqrt 5): ", x]; Exit[1]];
   {a, b}];
(* sign of a + b sqrt 5 (a, b rational): if a and b have different signs, compare a^2 with 5 b^2 *)
q5Sign[{a_, b_}] := Which[b == 0, Sign[a], a == 0, Sign[b], Sign[a] == Sign[b], Sign[a], True, Sign[a] Sign[a^2 - 5 b^2]];
sgn[x_] := q5Sign[q5[x]];
exactMin[list_] := Fold[If[sgn[#2 - #1] < 0, #2, #1] &, list];
segDist2[p_, {a_, b_}] := Module[{d = b - a, t},
   t = ((p - a) . d)/(d . d);
   Which[sgn[t] <= 0, (p - a) . (p - a), sgn[t - 1] >= 0, (p - b) . (p - b), True, (p - a - t d) . (p - a - t d)]];
delta2 = (5 - Sqrt[5])/8;                                                      (* sin^2 36 *)
res = Table[Module[{t = ptCanonToTri[k], c = ptCanonCentre[k], tris, edges, bnd, cv, cands, signs},
    tris = Table[{t["A"][[j]], t["B"][[j]], t["C"][[j]]}, {j, ptNTri[t]}];                    (* exact 4-vectors *)
    edges = Flatten[Table[Sort /@ Partition[Append[tris[[j]], tris[[j, 1]]], 2, 1], {j, Length[tris]}], 1];
    bnd = Select[Tally[edges], #[[2]] == 1 &][[All, 1]];                                       (* boundary edges, exact *)
    cv = embX /@ tris[[c]];
    cands = Join[Flatten[Table[segDist2[p, embX /@ e], {p, cv}, {e, bnd}]],
       Flatten[Table[segDist2[embX[q], {cv[[i]], cv[[Mod[i, 3] + 1]]}], {q, Union[Flatten[bnd, 1]]}, {i, 3}]]];
    signs = sgn[# - delta2] & /@ cands;
    <|"ncands" -> Length[cands], "allGE" -> AllTrue[signs, # >= 0 &], "someEQ" -> MemberQ[signs, 0],
      "min2" -> RootReduce[exactMin[cands]]|>], {k, keys}];
Print["coronas: ", Length[res], "; candidate distances per corona: ", MinMax[res[[All, "ncands"]]]];
Print["every candidate squared distance >= (5 - sqrt 5)/8, for all coronas (exact): ", AllTrue[res[[All, "allGE"]], TrueQ]];
Print["some candidate equal to (5 - sqrt 5)/8, for all coronas (exact): ", AllTrue[res[[All, "someEQ"]], TrueQ]];
Print["delta_C^2 (exact minimum of the candidates) for the 30 coronas: ", Union[res[[All, "min2"]]]];
Print["all delta_C^2 equal to sin^2 36 = (5 - sqrt 5)/8 (exact): ", AllTrue[res[[All, "min2"]], RootReduce[# - delta2] === 0 &]];
Print["delta = sqrt(min) >= 0.58 (exact, i.e. min >= (29/50)^2): ", sgn[delta2 - (29/50)^2] >= 0, "   numerically delta = ", N[Sqrt[delta2], 12], " (for information only)"];
