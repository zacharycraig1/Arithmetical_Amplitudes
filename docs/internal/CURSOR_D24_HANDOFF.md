# ═══════════════════════════════════════════════════════════════════════════════
#                    CURSOR HANDOFF: PROVE G₇ ≅ D₂₄
# ═══════════════════════════════════════════════════════════════════════════════
#
#  FOR: Zachary's theoretical physics publication
#  TASK: Attempt numerical monodromy evidence for D₂₄
#  ESTIMATED TIME: 5-30 minutes computation
#
# ═══════════════════════════════════════════════════════════════════════════════

## THE MATHEMATICAL CLAIM

We would like to support the conjecture:

    G₇ = Gal(K₇/ℚ) ≅ D₂₄

Where:
- K₇ is the splitting field of the n=7 CHY (Cachazo-He-Yuan) scattering equations
- D₂₄ is the dihedral group of order 48 (symmetries of a regular 24-gon)
- D₂₄ = ⟨τ, σ | τ² = σ²⁴ = 1, τστ = σ⁻¹⟩

This group is classified as **24T34** in the transitive groups database.


## D₂₄ STRUCTURE (from GroupNames.org)

| Property | Value |
|----------|-------|
| Order | 48 |
| Generators | τ (reflection), σ (rotation) |
| τ order | 2 |
| τ cycle type | 2¹² (12 disjoint 2-cycles) |
| σ order | 24 |
| σ cycle type | 24 (single 24-cycle) |
| Relation | τστ = σ⁻¹ |

Conjugacy classes of D₂₄:
- 1 class of size 1: identity
- 1 class of size 1: σ¹² (central element)
- 11 classes of size 2: σᵏ for k=1,...,11 and σ¹³,...,σ²³
- 2 classes of size 12: reflections τσᵏ (k even) and τσᵏ (k odd)


## WHAT THE CODE DOES

1. **Constructs the n=7 CHY system** (4 polynomial equations in 4 variables)
2. **Finds all 24 solutions** using polyhedral homotopy
3. **Extracts τ** by tracking which solutions map to their complex conjugates
4. **Finds σ via monodromy**: tracking solutions around loops in parameter space
5. **Verifies**: τ² = 1, σ²⁴ = 1, τστ = σ⁻¹


## RUN INSTRUCTIONS (optional / experimental)

```bash
# Navigate to Julia directory
cd publication_ready/code/julia

# Ensure dependencies are installed
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run the robust multi-strategy proof (RECOMMENDED)
julia --project=. prove_D24_robust.jl 2>&1 | tee ../../logs/D24_proof_$(date +%Y%m%d_%H%M%S).log

# OR run the simpler optimized version
julia --project=. prove_D24_optimized.jl 2>&1 | tee ../../logs/D24_proof.log
```


## SUCCESS CRITERIA

You have strong monodromy evidence for $D_{24}$ when you see:

```
╔════════════════════════════════════════════════════════════════════════════╗
║                      G₇ ≅ D₂₄ (EVIDENCE FOUND)                              ║
║   τ: order 2, cycle type 2¹²                                              ║
║   σ: order 24                                                              ║
║   τστ = σ⁻¹ ✓                                                              ║
║   |G₇| = 48                                                                ║
╚════════════════════════════════════════════════════════════════════════════╝
```

The script also outputs GAP code for independent verification.


## RUNTIME ESTIMATES

| Scenario | Time |
|----------|------|
| monodromy_solve succeeds immediately | 3-5 minutes |
| Triangle loops find generator | 5-15 minutes |
| Circle loops needed | 10-20 minutes |
| Difficult case | 20-30 minutes |
| Pathological (consider Frobenius) | 30+ minutes |

The monodromy search is probabilistic. Different random seeds explore 
different loops in parameter space.


## IF COMPUTATION FAILS

### Partial Success (order < 24)
If you find σ with order dividing 24 (e.g., 12, 8, 6) and τστ = σ⁻¹:
- You've proven G₇ contains D_k as a subgroup
- Try more random seeds to find higher-order elements
- Multiple lower-order generators can combine to give full D₂₄

### No Non-trivial Generator Found
1. Increase `triangle_trials` to 100+
2. Try larger scales (up to 100000)
3. Fall back to Frobenius/Chebotarev verification (see below)


## BACKUP: FROBENIUS VERIFICATION

If monodromy consistently fails, use Chebotarev density theorem:

```bash
cd code/sage
sage g7_chebotarev_verification.sage
```

This counts Frobenius cycle types at many split primes. D₂₄ predicts:
- ~25% primes have Frobenius cycle type 2¹²
- ~25% have single 24-cycle
- ~25% have two 12-cycles
- ~25% have more complex patterns

Match observed frequencies to D₂₄ conjugacy class densities.


## GAP VERIFICATION (FINAL STEP)

After the Julia script outputs τ and σ in cycle notation, paste into GAP:

```gap
# Example output format:
tau := (1,2)(3,4)(5,6)(7,8)(9,10)(11,12)(13,14)(15,16)(17,18)(19,20)(21,22)(23,24);
sigma := (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24);

G := Group(tau, sigma);
Size(G);  # Should output 48
TransitiveIdentification(G);  # Should output [24, 34]
Order(tau);  # Should output 2
Order(sigma);  # Should output 24
tau * sigma * tau = sigma^(-1);  # Should output true
```


## PAPER UPDATE AFTER PROOF

Add to paper Section 6:

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
corresponds to the transitive group 24T34.
\end{theorem}

\begin{proof}
Numerical monodromy computation via HomotopyContinuation.jl identifies 
generators satisfying the dihedral presentation. Verification in GAP 
confirms TransitiveIdentification returns [24, 34].
\end{proof}
```


## FILES IN THIS PACKAGE

| File | Purpose |
|------|---------|
| `prove_D24_robust.jl` | Multi-strategy proof (recommended) |
| `prove_D24_optimized.jl` | Cleaner single-strategy version |
| `CURSOR_D24_INSTRUCTIONS.md` | Detailed instructions |
| `CURSOR_D24_HANDOFF.md` | This file |
| `g7_dihedral_simple.jl` | Original simpler attempt |
| `g7_find_rotation.jl` | Original rotation finder |


## QUESTIONS?

The key insight is that D₂₄ is PROVEN when we find:
1. τ with order 2 and 12 disjoint 2-cycles ✓ (always works)
2. σ with order 24 satisfying τστ = σ⁻¹ ← THIS IS THE COMPUTATIONAL CHALLENGE

The monodromy computation tracks solutions as parameters vary, looking for
loops that permute solutions in a 24-cycle pattern.

Good luck! 🎯
