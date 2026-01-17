# REFEREE VERIFICATION: n=7 Conjugate Pair Structure
# ====================================================
#
# This script verifies:
# 1. The n=7 CHY system has exactly 24 physical solutions
# 2. All 24 solutions come in conjugate pairs (12 pairs)
# 3. The conjugation permutation has cycle type 2^12
#
# Usage:
#     julia verify_conjugates.jl
#
# Requirements:
#     Julia 1.9+
#     HomotopyContinuation.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "code", "julia"))
Pkg.instantiate()

using HomotopyContinuation
using LinearAlgebra

println("="^70)
println("REFEREE VERIFICATION: n=7 Conjugate Pair Structure")
println("="^70)

# =============================================================================
# BUILD CHY SYSTEM
# =============================================================================

@var x4 x5 x6 x7

# Fixed rational Mandelstams (verified CHY-consistent)
s = Dict{Tuple{Int,Int}, Int}(
    (1,2) => -6728, (1,3) => -8080, (1,4) => -14144, (1,5) => -22672, 
    (1,6) => 17440, (1,7) => 34184,
    (2,3) => 2740, (2,4) => 15304, (2,5) => 17032, (2,6) => -9524, (2,7) => -18824,
    (3,4) => 4112, (3,5) => 6408, (3,6) => -1768, (3,7) => -3412,
    (4,5) => 6032, (4,6) => -3904, (4,7) => -7400,
    (5,6) => -2248, (5,7) => -4552,
    (6,7) => 4
)

# Gauge: x1=0, x2=1, x3=-1
sigmas = Dict(1=>0, 2=>1, 3=>-1, 4=>x4, 5=>x5, 6=>x6, 7=>x7)

function build_chy_poly(a, sigmas, s_dict)
    get_s(i,j) = s_dict[(min(i,j), max(i,j))]
    poly = 0
    for b in 1:7
        if b == a continue end
        term = get_s(a, b)
        for c in 1:7
            if c == a || c == b continue end
            term = term * (sigmas[a] - sigmas[c])
        end
        poly = poly + term
    end
    return poly
end

println("\nBuilding CHY system...")
f4 = build_chy_poly(4, sigmas, s)
f5 = build_chy_poly(5, sigmas, s)
f6 = build_chy_poly(6, sigmas, s)
f7 = build_chy_poly(7, sigmas, s)

F = System([f4, f5, f6, f7])

# =============================================================================
# SOLVE AND FILTER
# =============================================================================

println("Solving...")
result = solve(F)

all_sols = solutions(result)
println("Total solutions: ", length(all_sols))

function is_physical(sol; tol=1e-6)
    vals = [sol[1], sol[2], sol[3], sol[4]]
    for v in vals
        if abs(v) < tol || abs(v-1) < tol || abs(v+1) < tol
            return false
        end
    end
    for i in 1:4, j in (i+1):4
        if abs(vals[i] - vals[j]) < tol
            return false
        end
    end
    return true
end

physical_sols = filter(is_physical, all_sols)
println("Physical solutions: ", length(physical_sols))

# Sort for consistent indexing
sorted_sols = sort(physical_sols, by = sol -> (real(sol[1]), imag(sol[1]), real(sol[2])))

# =============================================================================
# VERIFY CONJUGATE PAIRS
# =============================================================================

println("\n" * "="^70)
println("CONJUGATE PAIR VERIFICATION")
println("="^70)

function find_conjugate_permutation(sols)
    n = length(sols)
    perm = zeros(Int, n)
    
    for i in 1:n
        conj_i = conj.(sols[i])
        for j in 1:n
            if norm(sols[j] .- conj_i) < 1e-8
                perm[i] = j
                break
            end
        end
    end
    return perm
end

τ = find_conjugate_permutation(sorted_sols)

# Check τ² = identity
τ_squared = [τ[τ[i]] for i in 1:length(τ)]
is_involution = all(τ_squared[i] == i for i in 1:length(τ))

# Count fixed points and 2-cycles
num_fixed = count(i -> τ[i] == i, 1:length(τ))
num_2cycles = (length(τ) - num_fixed) ÷ 2

println("Conjugation permutation τ:")
println("  Length: ", length(τ))
println("  τ² = identity: ", is_involution)
println("  Fixed points: ", num_fixed)
println("  2-cycles: ", num_2cycles)

# =============================================================================
# RESULTS
# =============================================================================

println("\n" * "="^70)
println("VERIFICATION RESULTS")
println("="^70)

tests_passed = true

# Test 1: Exactly 24 solutions
if length(physical_sols) == 24
    println("✓ Test 1: Exactly 24 physical solutions")
else
    println("✗ Test 1: Expected 24 solutions, got ", length(physical_sols))
    tests_passed = false
end

# Test 2: All in conjugate pairs
if num_fixed == 0 && num_2cycles == 12
    println("✓ Test 2: All 24 solutions in 12 conjugate pairs")
else
    println("✗ Test 2: Expected 12 pairs with 0 fixed points")
    tests_passed = false
end

# Test 3: τ is involution
if is_involution
    println("✓ Test 3: τ² = identity (τ is an involution)")
else
    println("✗ Test 3: τ is not an involution")
    tests_passed = false
end

# Test 4: Cycle type is 2^12
if num_2cycles == 12 && num_fixed == 0
    println("✓ Test 4: Cycle type is 2¹² (12 disjoint 2-cycles)")
else
    println("✗ Test 4: Incorrect cycle type")
    tests_passed = false
end

println("\n" * "="^70)
if tests_passed
    println("✓ ALL TESTS PASSED")
    println("  The n=7 CHY variety has the expected conjugate structure.")
else
    println("✗ SOME TESTS FAILED")
end
println("="^70)
