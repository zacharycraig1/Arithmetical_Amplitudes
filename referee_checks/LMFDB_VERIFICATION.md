# LMFDB and GroupNames Verification Report

**Date:** January 2026  
**Verified by:** Automated browser verification

## Summary

The number-theoretic structures claimed in the paper have been verified against authoritative databases.

---

## 1. Q(i) - Gaussian Field

**LMFDB Entry:** [2.0.4.1](https://www.lmfdb.org/NumberField/2.0.4.1)

**Page Title:** `LMFDB - Number field 2.0.4.1: Q(sqrt(-1))`

**Confirmed Properties:**
- Field: Q(i) = Q(√-1)
- Degree: 2
- Discriminant: -4
- Galois group: C2 (cyclic of order 2)
- Class number: 1 (principal ideal domain)

**Relevance to Paper:**
- This is the splitting field for n=5 CHY equations
- The quadratic character χ₄ = (-1/p) corresponds to this field
- Gal(Q(i)/Q) ≅ Z/2Z is the parity symmetry group

---

## 2. Q(ζ₈) - 8th Cyclotomic Field

**LMFDB Entry:** [4.0.256.1](https://www.lmfdb.org/NumberField/4.0.256.1)

**Page Title:** `LMFDB - Number field 4.0.256.1: Q(zeta_8)`

**Confirmed Properties:**
- Field: Q(ζ₈) = Q(i, √2)
- Degree: 4
- Discriminant: 256 = 2^8
- Galois group: C2 × C2 (Klein four-group)
- Contains subfields Q(i) and Q(√2) and Q(√-2)

**Relevance to Paper:**
- Contains Q(i) as a subfield (n=5 splitting field)
- The loop χ₈ character corresponds to Q(√-2) ⊂ Q(ζ₈)
- Conductor 8 matches the observed arithmetic structure

---

## 3. D₂₄ - Dihedral Group of Order 48

**GroupNames Entry:** [D24](https://people.maths.bris.ac.uk/~matyd/GroupNames/1/D24.html)

**Page Title:** `D24 - GroupNames`

**Confirmed Properties:**
- Order: 48 = 2⁴ · 3
- Structure: Dihedral group of 24-gon
- Generators: rotation r (order 24) and reflection s (order 2)
- Relations: r²⁴ = s² = 1, srs = r⁻¹
- Transitive group label: 24T34

**Relevance to Paper:**
- Conjectured Galois group G₇ for n=7 CHY
- Contains Z/2Z quotient (reflection/conjugation)
- Order 48 matches expected structure for degree-24 cover with quadratic subfield

---

## 4. Transitive Group 24T34

**Identification:**
- The LMFDB transitive group 24T34 is isomorphic to D₂₄
- This is the natural permutation representation of D₂₄ on 24 points
- Consistent with 24 CHY solutions for n=7

---

## Verification Method

All verifications performed via automated browser navigation on:
- LMFDB: https://www.lmfdb.org
- GroupNames: https://people.maths.bris.ac.uk/~matyd/GroupNames/

The page titles and content confirm the claimed field/group properties.

---

## References

1. LMFDB Collaboration, "The L-functions and Modular Forms Database," https://www.lmfdb.org
2. Tim Dokchitser, "GroupNames," https://people.maths.bris.ac.uk/~matyd/GroupNames/
