(* Pythagorean Null Vector Lemma Verification *)
(* Claim: v = (E, px, py, pz) with the given parametrization satisfies E^2 - px^2 - py^2 - pz^2 = 0 *)

ClearAll[m, n, p, q];

E0 = m^2 + n^2 + p^2 + q^2;
px = m^2 + n^2 - p^2 - q^2;
py = 2*(m*q + n*p);
pz = 2*(n*q - m*p);

(* Minkowski norm: E^2 - px^2 - py^2 - pz^2 *)
norm = E0^2 - px^2 - py^2 - pz^2;

(* Expand and simplify *)
result = Simplify[Expand[norm]];

Print["========================================"];
Print["PYTHAGOREAN NULL VECTOR LEMMA"];
Print["========================================"];
Print[""];
Print["Parametrization:"];
Print["  E  = m^2 + n^2 + p^2 + q^2"];
Print["  px = m^2 + n^2 - p^2 - q^2"];
Print["  py = 2(mq + np)"];
Print["  pz = 2(nq - mp)"];
Print[""];
Print["Computing E^2 - px^2 - py^2 - pz^2..."];
Print[""];
Print["Result: ", result];
Print[""];
If[result === 0, 
  Print["ALGEBRAICALLY VERIFIED: This parametrization always produces null vectors."];
  Exit[0],
  Print["ERROR: Minkowski norm is not zero!"];
  Exit[1]
]
