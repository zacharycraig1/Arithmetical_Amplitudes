# 8D/11D Generic Kinematics + χ₈ Pipeline

This directory contains new tools for testing the CHY scattering equation arithmetic structure using generic higher-dimensional kinematics.

## Motivation

The existing 4D kinematics generator produces points in a specific subspace of kinematic moduli space. The key question is:

> **Is the observed inert vanishing (N_p = 0 for p ≡ 3 mod 4) a genuine Galois/Frobenius phenomenon, or an artifact of 4D kinematics hitting a special locus?**

These tools provide alternative kinematics generators to test this hypothesis.

## New Files

### `generic_kinematics_8d11d.sage`

The main module providing:

1. **8D/11D Lightcone Null Kinematics** (`generate_null_kinematics_GFp`, `generate_null_kinematics_GFp2`)
   - Uses lightcone coordinates to avoid square roots over finite fields
   - Generates n massless momenta in D dimensions with momentum conservation
   - Works over GF(p) or GF(p²)

2. **Dimension-Free CHY Mandelstams** (`random_chy_mandelstams`, `random_chy_mandelstams_QQ`)
   - Generates random s_ij satisfying only the CHY row-sum constraint
   - Fastest approach for prime sweeps (no momentum-space machinery)

3. **χ₈ Character Utilities** (`chi8m2`, `chi8_2`, `chi4`, `classify_prime`)
   - χ₈⁻(p) = (-2/p) and χ₈⁺(p) = (2/p)
   - `chi8(p)` is an alias for χ₈⁻ to match the paper
   - χ₄(p) = (-1/p): distinguishes inert (p ≡ 3 mod 4) from split primes

4. **CHY Solution Counting** (`count_chy_solutions`)
   - Groebner saturation approach for counting solutions over finite fields

5. **Prime Sweep with χ₈ Correlation** (`chi8_prime_sweep`, `analyze_chi8_correlation`)
   - Systematic testing of χ₈ patterns across many primes

### `chi8_correlation_test.sage`

Command-line tool for testing χ₈⁻ correlation:

```bash
# Use 4D kinematics (existing approach)
sage chi8_correlation_test.sage --mode 4d --p_max 100

# Use 11D kinematics
sage chi8_correlation_test.sage --mode 8d11d --dim 11 --p_max 100

# Use dimension-free Mandelstams (fastest)
sage chi8_correlation_test.sage --mode mandelstam --p_max 200
```

### `frob_cycletypes_generic.sage`

Frobenius cycle type extraction with generic kinematics:

```bash
sage frob_cycletypes_generic.sage --kin_mode mandelstam --p_max 200
sage frob_cycletypes_generic.sage --kin_mode 8d11d --dim 11 --p_max 100
```

### `test_8d11d_kinematics.sage`

Smoke test for validating all generators work correctly:

```bash
sage test_8d11d_kinematics.sage
```

## Makefile Targets

```bash
make test-8d11d           # Smoke test
make chi8-mandelstam      # χ₈ test with dimension-free Mandelstams
make chi8-11d             # χ₈ test with 11D kinematics
make chi8-4d              # χ₈ test with 4D kinematics (comparison)
make chi8-full            # Full χ₈ sweep (all modes)
make frob-generic         # Frobenius cycle types with generic kinematics
make frob-11d             # Frobenius with 11D kinematics
make verify-generic-quick # Quick verification of new tools
```

## Key Findings from Initial Tests

1. **Inert Vanishing Confirmed**: ~86% of inert primes (p ≡ 3 mod 4) have N=0
2. **χ₈⁻ Correlation Not Significant**: In the tested range, no clear correlation with χ₈⁻(p)
3. **4D vs 11D**: Both approaches show similar inert vanishing behavior

## Important Note on Frobenius Analysis

For Frobenius/Chebotarev analysis, the correct workflow is:

1. Generate kinematics **over QQ** (a fixed point in moduli space)
2. Reduce mod p for each prime
3. Observe fiber structure (eliminant factor degrees = Frobenius cycle types)

The existing 4D-over-QQ approach is correct for this. The new 8D/11D-over-GF(p) generators test something different: they sample different kinematic points for each prime, which is useful for understanding "typical behavior" but not for Frobenius analysis of a fixed cover.

## Physics Story (Gut Instinct)

From the handoff document:

> My gut says you're not "seeing a random mod-p coincidence." You are watching **Frobenius split a branched cover**: the CHY scattering equation solutions form a cover of kinematic space, and when you reduce mod p, Frobenius chooses how the fiber splits/fuses.

The χ₄ (inert vanishing) and potential χ₈⁻ behavior is the signature of a **dyadic quadratic twist** living in the geometry. 8D/11D kinematics doesn't change CHY—it removes the "stage fog" so you can see the arithmetic structure cleanly.

## Mathematical Background

### Lightcone Coordinates

Represent momentum as p = (u, v, x₁, ..., x_{D-2}) with metric:
```
p² = -2uv + Σ xᵢ²
```

To generate a null vector (p² = 0) without square roots:
1. Set u = 1
2. Choose random x₁, ..., x_{D-2}
3. Set v = (Σ xᵢ²) / 2

### χ₈ Characters (Conductor 8)

We distinguish two quadratic characters:

```
χ₈⁻(p) = (-2/p) = { +1 if p ≡ 1, 3 (mod 8)
                  { -1 if p ≡ 5, 7 (mod 8)

χ₈⁺(p) = (2/p)  = { +1 if p ≡ 1, 7 (mod 8)
                  { -1 if p ≡ 3, 5 (mod 8)
```

Relation:
```
χ₈⁻(p) = χ₄(p) · χ₈⁺(p)
```
