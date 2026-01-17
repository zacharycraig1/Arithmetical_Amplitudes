# Candidate Galois Group: C2 × C2 × S4

## LMFDB Identification

Based on the hint "C×C×S4" and LMFDB analysis, the most likely candidate is:

**C2² × S4 = C2 × C2 × S4**

This appears in LMFDB as:
- **24T125**: C2² × S4 acting on 24 points (one representation)
- **24T126**: C2² × S4 acting on 24 points (another representation)
- **12T48**: C2² × S4 acting on 12 points
- **16T182**: C2² × S4 acting on 16 points

## Key Properties

| Property | Value |
|----------|-------|
| Abstract Group | C2² × S4 |
| Order | 96 = 2⁵ × 3 |
| Solvable | Yes |
| Degree (as 24T125) | 24 |
| Primitive | No (imprimitive) |

## Subfield Structure of 24T125

From LMFDB:
- **Degree 2**: C2 ← **This explains Q(i) obstruction!**
- **Degree 3**: S3
- **Degree 6**: S3, S4×C2 (×6)
- **Degree 12**: C2×S4 (×3)

## Why This Makes Sense for CHY n=7

1. **24 = 4! solutions**: S4 permutes the 4 unfixed particles naturally

2. **Q(i) subfield**: The C2 subfield at degree 2 explains why:
   - N_p = 0 at inert primes (Frobenius in non-identity C2 coset)
   - N_{p²} > 0 (Frobenius² in identity coset)

3. **C2 × C2 factor**: Could arise from:
   - Conjugation symmetry of solutions
   - Parity/helicity structure
   - Two independent quadratic extensions

4. **Imprimitive structure**: The 24 solutions may decompose into blocks related to particle subsets

## Relation to Original D24 Conjecture

D24 (dihedral of order 48) ≠ C2² × S4 (order 96):
- D24 has order 48, C2² × S4 has order 96
- D24 is dihedral, C2² × S4 is a direct product
- Both have index-2 subgroups giving Q(i) structure

The original conjecture may have confused:
- The action (dihedral-like behavior)
- With the actual group (which is larger and richer)

## Verification Path

To confirm this is the Galois group:

1. **Compute eliminant f(x) ∈ Q[x]** of degree 24
2. **Use Sage/Magma**: `f.galois_group()` 
3. **Check**: Should return 24T125 or 24T126

## LMFDB Links

- [24T125](https://www.lmfdb.org/GaloisGroup/24T125)
- [24T126](https://www.lmfdb.org/GaloisGroup/24T126)
- [12T48](https://www.lmfdb.org/GaloisGroup/12T48)

## Open Question

Our Frobenius data shows N_{p^{24}} < 24, which is puzzling if the Galois group acts transitively. Possible explanations:

1. **Bad reduction**: The eliminant has denominators divisible by test primes
2. **Factorization**: The eliminant factors, and each factor has smaller Galois group
3. **Computational issue**: The variety computation misses some solutions

The HPC eliminant computation will resolve this.
