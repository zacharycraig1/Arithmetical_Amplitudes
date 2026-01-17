(* 
REFEREE VERIFICATION: Gram-Levi-Civita Identity
================================================

This script algebraically verifies:
    D(v1, v2, v3, v4) = -epsilon(v1, v2, v3, v4)^2

where D is the Gram determinant and epsilon is the Levi-Civita contraction.

Usage:
    wolframscript -file verify_gram_levi.wl

Expected output:
    ALGEBRAICALLY VERIFIED: D = -epsilon^2
*)

Print["================================================"]
Print["REFEREE VERIFICATION: Gram-Levi-Civita Identity"]
Print["================================================"]
Print[""]

(* Define 4 generic 4-vectors *)
v1 = Array[Subscript[v, 1, #] &, 4];
v2 = Array[Subscript[v, 2, #] &, 4];
v3 = Array[Subscript[v, 3, #] &, 4];
v4 = Array[Subscript[v, 4, #] &, 4];

(* Minkowski metric: signature (+,-,-,-) *)
eta = DiagonalMatrix[{1, -1, -1, -1}];

(* Minkowski inner product *)
MinkowskiDot[a_, b_] := a . eta . b

(* Gram matrix *)
G = {{MinkowskiDot[v1, v1], MinkowskiDot[v1, v2], MinkowskiDot[v1, v3], MinkowskiDot[v1, v4]},
     {MinkowskiDot[v2, v1], MinkowskiDot[v2, v2], MinkowskiDot[v2, v3], MinkowskiDot[v2, v4]},
     {MinkowskiDot[v3, v1], MinkowskiDot[v3, v2], MinkowskiDot[v3, v3], MinkowskiDot[v3, v4]},
     {MinkowskiDot[v4, v1], MinkowskiDot[v4, v2], MinkowskiDot[v4, v3], MinkowskiDot[v4, v4]}};

Print["Computing Gram determinant D = det(G)..."]
D = Det[G] // Simplify;

(* 4x4 matrix of momentum components *)
M = {v1, v2, v3, v4};

Print["Computing Levi-Civita contraction epsilon = det(M)..."]
epsilon = Det[M];

(* The identity to verify: D = -epsilon^2 *)
Print["Verifying D + epsilon^2 = 0..."]
identity = Simplify[D + epsilon^2];

Print[""]
Print["Result: D + epsilon^2 = ", identity]
Print[""]

If[identity === 0,
   Print["================================================"]
   Print["ALGEBRAICALLY VERIFIED: D(v1,v2,v3,v4) = -epsilon(v1,v2,v3,v4)^2"]
   Print["================================================"]
   Print[""]
   Print["This identity is valid for ANY four 4-vectors in Minkowski space."]
   Print[""]
   Print["Consequence for CHY:"]
   Print["  Since D = -epsilon^2 and epsilon != 0 for generic kinematics,"]
   Print["  we have D < 0, hence sqrt(D) = i*|epsilon|, and the"]
   Print["  splitting field is Q(i)."]
   ,
   Print["VERIFICATION FAILED!"]
   Print["The identity D = -epsilon^2 does not hold."]
   Exit[1]
]
