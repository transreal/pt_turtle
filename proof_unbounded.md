# Unboundedness of (SSSR)* on Penrose rhomb tilings — statement, proof and certificate format

## Statement

**Theorem U.** Let T be any Penrose rhomb tiling and P = SSSR.
1. For every M there is a turtle state of T whose orbit under (SSSR)* is either infinite or has
   period > M.  More precisely, for every m ≥ 7 odd and every positively oriented obtuse level-m
   supertriangle X of T there is a state in a rhomb owned by X whose forward trajectory performs
   N_m S-moves, the first N_m − 1 between rhombs owned by X and the last one leaving R⁺(X), all states
   being pairwise distinct, where
   N_7 = 202 and N_{m+2} = 4 N_m + 66 (N_9 = 874, N_11 = 3562, N_13 = 14314, ...).
   Hence the orbit of that state has period ≥ N_m / 3 (cycles) or is infinite.
2. There exists a Penrose rhomb tiling and a state whose (SSSR)* orbit is infinite, hence unbounded.

So (SSSR)* is *not* bounded: no Penrose tiling admits a uniform bound on orbit sizes, and some
Penrose tiling carries a genuinely unbounded orbit.  (Whether a given regular pentagrid tiling
carries an infinite orbit is not decided here.)

## Setting and conventions

* Marked Robinson triangles: acute (col 0; apex A, angles 36°,72°,72°; legs 1, base 1/φ) and obtuse
  (col 1; apex A, 108°,36°,36°; legs 1, base φ); a triangle is a tuple (col, A, B, C), the order of
  B, C carrying the marking.  Points are elements of Z[ζ], ζ = e^{2πi/5}, stored as integer 5-vectors
  modulo (1,1,1,1,1).  Multiplication by φ = 1+ζ+ζ^4 and rotation by 36° (= multiplication by −ζ^3)
  are integer linear maps; all coordinates in the proof are exact.
* Substitution σ (inflate by φ, then subdivide), with X′ = φX:
  acute (A,B,C): P = A′+(B−A); children, in this order (slots 0,1): acute(C′,P,B′), obtuse(P,C′,A′).
  obtuse (A,B,C): Q = B′+(A−B), R = B′+(C−B); slots 0,1,2: obtuse(R,C′,A′), obtuse(Q,R,B′), acute(R,Q,A′).
* A type is τ = (col, s), s = orientation of (A,B,C).  The standard triangle of type (col, +1) has A = 0,
  B = 1, C = ρ^θ (θ = 1 for acute, 3 for obtuse; ρ = rotation by 36°); the standard triangle of type
  (col, −1) is the same point set with B and C exchanged (A = 0, B = ρ^θ, C = 1).  The standard level-m supertile X_m^τ is σ^m of it (as a set D of fine triangles);
  coordinates scale exactly by φ per level.  A child in slot j of a level-m supertile of type τ is
  G_{j,m}(X_{m−1}^{τ_j}) with G_{j,m}(x) = ρ^{k_j} x + φ^{m−1} t_j (ρ = rotation by 36°), where
  (τ_j, k_j, t_j) is read off at level 1.
* Classical facts used: every Penrose rhomb tiling is obtained from a marked Robinson-triangle tiling
  by merging the two triangles sharing each base (they are mirror images of each other); the triangle
  tilings form one local-isomorphism class, closed under composition; hence T is, for every m, a
  union of level-m supertiles G(X_m^τ) (G an orientation-preserving isometry with rotation a multiple
  of 36°), and all four types occur.  The hull (all Penrose rhomb tilings) is compact in the local topology.
* Turtle: state (rhomb, edge e ∈ Z_4 in counter-clockwise order, command index).  S crosses edge e
  into the adjacent rhomb and takes the opposite edge of that rhomb; R: e → e−1; L: e → e+1.
  Every rhomb edge is a leg of a fine triangle.
* **Ownership.** A rhomb is owned by the supertile that contains its positively oriented half.  For a
  supertile X, R⁺(X) is the set of rhombs owned by X.  The geometry of every rhomb in R⁺(X) is
  determined by X alone (the other half is the mirror image across the base), and two rhombs of
  R⁺(X) are adjacent in T iff they share a leg.  Ownership at level m−1 refines ownership at level m.
* **Crossing.** (tail, head, i): the turtle executes command number i (an S) and crosses the leg
  tail→head, written in the counter-clockwise order of the rhomb it leaves.
* A turtle move from a rhomb of R⁺(X) to a rhomb not in R⁺(X) is an *exit* of X.  The motion of the
  turtle inside R⁺(X) until its exit depends only on X (on D and the ghost halves), therefore it can
  be computed in the standard supertile.

## Certificate

* Base level B (here B = 6); R(q) = the smallest level ≥ B of parity q.
* Port families: id ↦ (type τ, parity q, reference crossing χ at level R(q), δ = (δ_tail, δ_head)).
  The port at level m = R(q)+2s is defined by χ_{m+2} = (φ² tail + δ_tail, φ² head + δ_head, i).
* Pass families: entry port id ↦ (exit port id, children [(slot j, child entry port, child exit port)]).
  A pass has the type and parity of its ports; children have the opposite parity and the type of the slot.

**H(m)** (for m ≥ B and every pass family π of parity m mod 2, X = X_m^{τ(π)}): the crossing
port(π.in, m) leads into a rhomb owned by X; starting there, the turtle performs exactly N(π,m)
S-moves; the first N−1 of them are between rhombs owned by X and the N-th is the crossing
port(π.out, m), which leads to a rhomb not owned by X.  Here N(π,B) is obtained by simulation (families
of parity B mod 2) and N(π,m) = Σ_children N(child, m−1) for every m > B.  (Owned rhombs protrude beyond
the triangle X by their ghost halves; "inside X" always means "in rhombs of R⁺(X)".)

## Obligations checked by `verify_cert.py` (all exact except the distance margins)

* **O1 (base).** For every pass family of parity B mod 2: H(B) by direct simulation in X_B^τ.
* **O3 (closure).** Every child entry port is itself a pass family of the table, with the stated exit;
  types and parities are consistent; consecutive children occupy different slots.
* **O2 (port identities).** For every pass family π of parity q with children (j_1,a_1,b_1),…,(j_r,a_r,b_r),
  and for the two lowest levels m ≡ q (mod 2) above B:
  (a) port(π.in, m) = G_{j_1,m}(port(a_1, m−1));
  (b) G_{j_i,m}(port(b_i, m−1)) = G_{j_{i+1},m}(port(a_{i+1}, m−1)) for 1 ≤ i < r;
  (c) port(π.out, m) = G_{j_r,m}(port(b_r, m−1)).
  *Lemma 1.* Both sides of each identity are sequences u_m with u_{m+2} = φ² u_m + const (left side by
  definition; right side because G_{j,m+2}(φ² x + δ) = φ² G_{j,m}(x) + ρ^{k_j} δ).  Two such sequences
  that agree at two consecutive levels of the parity class have the same constant and therefore agree
  at all levels.
* **O5 (exit property above the base).** For every pass family π (parity q), let (j_r, ·, b) be its last
  child, n_0 the lowest child level used in an inductive step, p_n the tail point of port(b, n),
  c = −δ_tail/φ (so p_n − c = φ^{n−n_0}(p_{n_0} − c)), y = p_{n_0} − c.  The checker verifies that there is
  an edge UW of the standard child triangle such that (i) y lies on the open segment φ^{n_0}U—φ^{n_0}W
  (exact collinearity in Z[ζ]); (ii) the placed edge G_j(U)G_j(W) lies on a side of the parent triangle
  (exact); (iii) d := min(|y−φ^{n_0}U|, |y−φ^{n_0}W|) ≥ (1.903 + |c|)/sin 36°.
  *Lemma 2.* Under (i)–(iii), for every child level n ≥ n_0 of that parity the rhomb r′ entered by
  the crossing port(b, n) (in the parent frame) is not owned by the parent X_{n+1}.
  Proof. y_n = φ^{n−n_0} y lies on the child edge E′ at distance ≥ d from its endpoints, so the ball
  B(y_n, d sin 36°) meets the parent triangle only inside the child (all triangle angles are ≥ 36° and
  E′ lies on a side of the parent).  r′ contains the leg tail→head, has diameter ≤ 2cos18° < 1.903 and
  p_n = y_n + c, hence r′ ⊂ B(y_n, 1.903 + |c|) ⊂ B(y_n, d sin 36°).  By H(n) for the child, r′ is not
  owned by the child, i.e. its positive half t⁺ is not a triangle of the child; t⁺ is a fine triangle
  inside that ball, so it cannot lie in the parent triangle.  ∎
* **O4 (growth).** Some pass family occurs μ ≥ 1 times in its own two-level expansion, which has more
  than μ elements; hence N(π, m+2) ≥ μ N(π, m) + 1 and N(π, m) → ∞.  The root of the certificate must be a
  pass family of the table.
* **O6 (exact lengths).** The lengths of the families of one parity satisfy N(m+2) = M N(m) with a fixed
  non-negative integer matrix M (two-level expansion).  Hence y_s = N(π, m_0+2s+2) − a N(π, m_0+2s) − b
  obeys a linear recurrence of order ≤ dim M + 1 (Cayley–Hamilton), and if it vanishes for dim M + 2
  consecutive s it vanishes identically.  The checker verifies this for dim M + 3 steps, which proves the
  stated closed recurrences (e.g. N(m+2) = 4 N(m) + 66 for the obtuse (+) root at odd levels).

## Proof of H(m) for all m ≥ B

Base: O1.  Step m−1 → m for a family π of parity m mod 2: by (a) the entry crossing of π at level m
is the entry crossing of the first child pass, so by H(m−1) (applied in the frame of child j_1, which
is an orientation-preserving isometric copy of the standard supertile and whose owned rhombs are owned
by X) the turtle traverses child j_1 and leaves it by G(port(b_1, m−1)); by (b) this crossing is the
entry crossing of the second child pass, which leads into a rhomb owned by child j_2; and so on.  The
last child pass ends with the crossing G(port(b_r, m−1)) = port(π.out, m) by (c), and by Lemma 2 this
crossing leads to a rhomb not owned by X.  All earlier S-moves are between rhombs owned by children of
X, hence owned by X.  The number of S-moves is the sum of the children's.  ∎

## Proof of Theorem U

(1) Let π be a self-reproducing pass family (O4) and X ⊂ T a level-m supertile of type τ(π), m of the
right parity.  By H(m) the state right after the entry crossing has a forward trajectory with N(π,m)
S-moves inside R⁺(X) followed by an exit.  If a state occurred twice before the exit, determinism
would force the trajectory to stay in R⁺(X) forever; hence the N(π,m) states are pairwise distinct
and the orbit is infinite or has at least N(π,m) S-moves per period, i.e. period ≥ N(π,m)/3 cycles.
N(π,m) → ∞ by O4.
(2) Let A_M be the set of pairs (tiling in the hull with a marked state at the origin) whose orbit does
not close within M cycles.  "Closing within M cycles" depends on a bounded patch, so A_M is closed;
the A_M are nested and non-empty by (1); the pointed hull is compact, so ∩_M A_M ≠ ∅.  An element has
an orbit that never closes; since the step map is a bijection on states, a non-closing orbit visits
infinitely many states, hence infinitely many rhombs, and is unbounded.  ∎

## Files

`lsys.py` (scattering through standard supertiles), `extract2.py` (constructs the certificate),
`verify_cert.py` (independent checker; shares no code with the extractor),
`cert2_SSSR_*.json` (certificates; e.g. `cert2_SSSR_1p_1.json`: 60 pass families, root = obtuse (+), odd levels).

## Review (2026-10-01)

An independent checker written only from this specification (own Z[zeta] arithmetic, substitution,
ownership rule and turtle; `indep_check.py`) accepts the certificate: O1 30/30, O2 320 identities,
O3, O5 for all 60 families (worst margin ratio 1.21), O4; direct simulations of the root pass at levels
9, 11 and 13 give exactly 874, 3562 and 14314 S-moves and the predicted exit ports; 180 further direct
simulations of all families at levels 6–11 agree with the recursion.  Four adversarial reviews (induction,
exit lemma, deduction of the theorem, checker code) found no soundness gap.  Their remarks were
incorporated: the convention for negatively oriented standard triangles, the definition of N for all
m > B, the wording of Theorem U(1), and the new obligation O6 for the exact recurrences; the checker now
also validates the program string and the root.
