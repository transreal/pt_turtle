(* wl_verify.wl -- independent Wolfram Language re-check of the proof patches.
   Reads patches_<PROG>_n<n>.json (exported from Python): each patch is a list of rhombs
   given by 4 integer Z^5 vertex keys (CCW), and a list of start tiles (1-based).
   Uses the user's own TurtleTiling.wl (BuildTurtleNav / TurtleRun) on a graph built from
   the keys, and counts closed / boundary-escape / non-converged orbits from every start tile
   and every heading.  Usage: wolframscript -file wl_verify.wl patches_SSR_n7.json *)
(* directory holding the author's CellularAutomata.wl and TurtleTiling.wl (not part of this repository):
   environment variable MYPACKAGES, or the directory of this script *)
$packageDirectory = With[{e = Environment["MYPACKAGES"]}, If[StringQ[e], e, DirectoryName[$InputFileName]]];
AppendTo[$Path, $packageDirectory];
Block[{$CharacterEncoding = "UTF-8"},
  Needs["CellularAutomata`", FileNameJoin[{$packageDirectory, "CellularAutomata.wl"}]];
  Needs["TurtleTiling`", FileNameJoin[{$packageDirectory, "TurtleTiling.wl"}]]];
file = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "patches_SSR_n7.json"];
data = Import[file, "RawJSON"];
prog = data["prog"]; n = data["n"];
Print["program (", prog, ")*  n=", n, "  patches: ", Length[data["patches"]]];
eVec = Table[N@{Cos[2 Pi i/5], Sin[2 Pi i/5]}, {i, 0, 4}];
t0 = AbsoluteTime[];
totals = <|"starts" -> 0, "closed" -> 0, "boundary" -> 0, "open" -> 0|>;
spectrum = <||>;
Do[
  keys = patch["keys"]; starts = patch["starts"];
  tiles = Table[<|"label" -> Join[keys[[i, 1]], {0, 0}],
     "type" -> "unknown",
     "vertices" -> (keys[[i]] . eVec),
     "vertexKeys" -> keys[[i]]|>, {i, Length[keys]}];
  graph = <|"tiles" -> tiles, "neumannNeighbors" -> Null|>;
  (* neumann neighbours via shared edge keys *)
  edgeMap = <||>;
  Do[Do[Module[{k = Sort[{keys[[i, e]], keys[[i, Mod[e, 4] + 1]]}]},
      edgeMap[k] = Append[Lookup[edgeMap, Key[k], {}], i]], {e, 4}], {i, Length[keys]}];
  nb = Table[{}, Length[keys]];
  Do[If[Length[v] == 2, nb[[v[[1]]]] = Union[nb[[v[[1]]]], {v[[2]]}];
       nb[[v[[2]]]] = Union[nb[[v[[2]]]], {v[[1]]}]], {v, Values[edgeMap]}];
  graph["neumannNeighbors"] = Association[Table[i -> nb[[i]], {i, Length[keys]}]];
  nav = BuildTurtleNav[graph];
  If[! AllTrue[Values[nav["windings"]], # == 1 &], Print["  WARNING: non-CCW tile in patch ", patch["corona"]]];
  res = Flatten[Table[TurtleRun[nav, prog, {c, e}, 200000], {c, starts}, {e, 1, 4}], 1];
  cl = Count[res, r_ /; r["closed"]]; bd = Count[res, r_ /; r["boundary"]]; op = Length[res] - cl - bd;
  totals["starts"] += Length[res]; totals["closed"] += cl; totals["boundary"] += bd; totals["open"] += op;
  Do[If[r["closed"], spectrum[r["period"]] = Lookup[spectrum, Key[r["period"]], 0] + 1], {r, res}];
  Print["  corona ", patch["corona"], ": rhombs ", Length[keys], "  starts ", Length[res],
    "  closed ", cl, "  boundary ", bd, "  open ", op, "   [", Round[AbsoluteTime[] - t0], " s]"];
  , {patch, data["patches"]}];
Print["TOTAL: ", totals];
Print["period spectrum: ", KeySort[spectrum]];
Print[If[totals["boundary"] == 0 && totals["open"] == 0,
   "ALL CLOSED: independent WL verification passed", "FAILED"]];
