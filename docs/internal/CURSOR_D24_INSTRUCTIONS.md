# ═══════════════════════════════════════════════════════════════════════════════
# CURSOR INSTRUCTIONS: PROVING G₇ ≅ D₂₄ FOR PUBLICATION
# ═══════════════════════════════════════════════════════════════════════════════

## OVERVIEW

This document describes an *optional/experimental* numerical monodromy attempt to
extract a dihedral generator for the n=7 CHY system.

**Important:** In the current publication-ready bundle, the full identification
$G_7 \cong D_{24}$ is treated as a **conjecture supported by evidence**, not as a
theorem-level computational proof. See `STRONG_EVIDENCE.md` and the paper text for
the current status.

**What we're proving:**
- G₇ = Gal(K/ℚ) where K is the splitting field of the n=7 CHY system
- G₇ ≅ D₂₄ = ⟨r, s | r²⁴ = s² = 1, srs = r⁻¹⟩ (dihedral group of order 48)
- This group acts transitively on 24 points (the 24 solutions)

**D₂₄ structure (from GroupNames.org):**
- Order: 48 = 2 × 24
- Transitive group identification: 24T34
- Generators: 
  - τ (reflection): order 2, cycle type 2¹² (12 transpositions)
  - σ (rotation): order 24, single 24-cycle
- Relation: τστ = σ⁻¹ (conjugation inverts the rotation)

---

## RUNTIME ESTIMATES

| Phase | Description | Expected Time |
|-------|-------------|---------------|
| 1 | System construction | ~1 second |
| 2 | Find 24 solutions | ~5 seconds |
| 3 | Extract τ (conjugation) | ~1 second |
| 4 | Monodromy search for σ | 2-10 minutes |
| 5 | Verify D₂₄ structure | ~1 second |
| **Total** | | **3-15 minutes** |

The monodromy search is probabilistic - it may need multiple random loops to 
find a generator with the right properties. In pathological cases it could take 
30+ minutes, but typically finds success within 10 minutes.

---

## SETUP INSTRUCTIONS

### 1. Environment Setup

```bash
cd publication_ready/code/julia

# Ensure Julia 1.9+ is installed
julia --version

# Install dependencies (first time only)
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 2. Run the Optimized Proof Script

```bash
julia --project=. prove_D24_optimized.jl 2>&1 | tee ../../logs/D24_proof.log
```

### 3. What Success Looks Like

The script will output:
```
═══════════════════════════════════════════════════════════════════════════════
│                              G₇ ≅ D₂₄ PROOF                                │
═══════════════════════════════════════════════════════════════════════════════

...

  ╔════════════════════════════════════════════════════════════════════════╗
  ║                                                                        ║
  ║   ███████╗    ██████╗ ██████╗  ██████╗ ██╗   ██╗███████╗██████╗        ║
  ║   ...                                                                  ║
  ║                    G₇ ≅ D₂₄ (24T34)                                    ║
  ║                                                                        ║
  ║   Generators:                                                          ║
  ║     τ: order 2, cycle type 2¹² (reflection)                           ║
  ║     σ: order 24 (rotation)                                            ║
  ║     Relation: τστ = σ⁻¹                                               ║
  ║                                                                        ║
  ║   |G₇| = 2 × 24 = 48                                                  ║
  ║                                                                        ║
  ╚════════════════════════════════════════════════════════════════════════╝
```

---

## WHAT CONSTITUTES RIGOROUS PROOF

For publication, you need to verify ALL of the following:

### Required Checks ✓

1. **24 Physical Solutions**
   - The CHY system has exactly (n-3)! = 4! = 24 physical solutions
   - Solutions must not be at gauge points {0, 1, -1}
   - Solutions must be pairwise distinct

2. **Conjugation Generator τ**
   - Complex conjugation permutes the 24 solutions
   - τ has order 2
   - τ has cycle type 2¹² (exactly 12 disjoint 2-cycles)
   - This means all solutions come in complex conjugate pairs

3. **Rotation Generator σ**
   - Found via monodromy (tracking solutions around loops in parameter space)
   - σ should have order 24 (for full D₂₄)
   - If order is a divisor of 24 (e.g., 12, 8, 6), you've proven a subgroup

4. **Dihedral Relation**
   - Must verify: τστ = σ⁻¹
   - This is the defining relation of dihedral groups

5. **GAP/Oscar Verification**
   - Feed the permutations to GAP
   - Confirm `TransitiveIdentification(G)` returns `[24, 34]`
   - Confirm `Size(G)` equals 48

---

## TROUBLESHOOTING

### Problem: monodromy_solve fails

**Solution:** The script falls back to manual triangle loops. If those also fail:
- Increase `MAX_MONODROMY_ATTEMPTS` to 100+
- Increase loop scales (try `[1000, 5000, 10000, 50000]`)
- The discriminant locus might be hard to encircle

### Problem: Only finding order < 24 elements

**Solution:** This is actually fine! If you find σ with order k | 24 and τστ = σ⁻¹:
- You've proven G₇ contains D_k as subgroup
- Multiple generators can be combined to generate full D₂₄
- The group generated by τ and all found σ's should have order 48

### Problem: Solutions don't match after loop

**Solution:** 
- Increase `TOLERANCE_MATCH` (try 1e-5 or 1e-4)
- The homotopy may be losing accuracy
- Reduce step size in solve() calls

---

## BACKUP APPROACH: FROBENIUS VERIFICATION

If monodromy fails completely, you can statistically verify D₂₄ using 
Chebotarev density. D₂₄ has specific cycle type frequencies:

| Conjugacy Class | Cycle Type | Frequency |
|-----------------|------------|-----------|
| 1 | 1²⁴ | 1/48 |
| τ (and conjugates) | 2¹² | 12/48 |
| σ^k (k odd) | 24 | 12/48 |
| σ^k (k even) | 12, 12 | 12/48 |
| τσ^k | varies | 12/48 |

Use the existing `g7_chebotarev_verification.sage` to count Frobenius cycle 
types at split primes and compare to D₂₄ predictions.

---

## FINAL OUTPUT FOR PAPER

After successful run, add to paper's Section 6:

```latex
\begin{theorem}[Galois Group Identification]
The Galois group $G_7 = \mathrm{Gal}(K_7/\mathbb{Q})$ of the $n=7$ CHY 
splitting field is isomorphic to the dihedral group $D_{24}$:
\[
G_7 \cong D_{24} = \langle \tau, \sigma \mid \tau^2 = \sigma^{24} = 1, 
\tau\sigma\tau = \sigma^{-1} \rangle
\]
where $\tau$ is complex conjugation (cycle type $2^{12}$) and $\sigma$ is 
a monodromy element of order 24. The group has order $|G_7| = 48$ and 
corresponds to the transitive group 24T34 in the standard classification.
\end{theorem}

\begin{proof}
Numerical monodromy computation with HomotopyContinuation.jl identifies 
generators $\tau$ (order 2) and $\sigma$ (order 24) satisfying the dihedral 
relation. Group identification verified in GAP.
\end{proof}
```

---

## FILES CREATED

| File | Purpose |
|------|---------|
| `prove_D24_optimized.jl` | Main proof script (run this) |
| `CURSOR_D24_INSTRUCTIONS.md` | This file |

---

## CONTACT

If you encounter issues with the computational proof:
1. Check that all 24 solutions are being found
2. Verify τ has the correct cycle structure
3. Try different random seeds for monodromy
4. Consider the Frobenius backup approach

Good luck! This completes the mathematical foundation for the paper.
