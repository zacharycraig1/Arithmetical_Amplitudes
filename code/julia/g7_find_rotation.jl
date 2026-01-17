# G₇ FIND ROTATION GENERATOR
# ===========================
# Try many random loops to find a non-trivial monodromy permutation

using HomotopyContinuation
using LinearAlgebra
using Random

println("="^70)
println("SEARCHING FOR ROTATION GENERATOR")
println("="^70)

# Setup system
@var x4 x5 x6 x7
@var p45 p46 p47

s = Dict{Tuple{Int,Int}, Any}(
    (1,2) => -6728, (1,3) => -8080, (1,4) => -14144, (1,5) => -22672, 
    (1,6) => 17440, (1,7) => 34184,
    (2,3) => 2740, (2,4) => 15304, (2,5) => 17032, (2,6) => -9524, (2,7) => -18824,
    (3,4) => 4112, (3,5) => 6408, (3,6) => -1768, (3,7) => -3412,
    (4,5) => p45, (4,6) => p46, (4,7) => p47,
    (5,6) => -2248, (5,7) => -4552,
    (6,7) => 4
)

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

g4 = build_chy_poly(4, sigmas, s)
g5 = build_chy_poly(5, sigmas, s)
g6 = build_chy_poly(6, sigmas, s)
g7 = build_chy_poly(7, sigmas, s)

G = System([g4, g5, g6, g7]; parameters=[p45, p46, p47])

# Base parameters
p0 = [6032.0 + 0.0im, -3904.0 + 0.0im, -7400.0 + 0.0im]

# First, get base solutions
println("\nGetting base solutions...")
@var x4_base x5_base x6_base x7_base
s_base = Dict{Tuple{Int,Int}, Int}(
    (1,2) => -6728, (1,3) => -8080, (1,4) => -14144, (1,5) => -22672, 
    (1,6) => 17440, (1,7) => 34184,
    (2,3) => 2740, (2,4) => 15304, (2,5) => 17032, (2,6) => -9524, (2,7) => -18824,
    (3,4) => 4112, (3,5) => 6408, (3,6) => -1768, (3,7) => -3412,
    (4,5) => 6032, (4,6) => -3904, (4,7) => -7400,
    (5,6) => -2248, (5,7) => -4552,
    (6,7) => 4
)
sigmas_base = Dict(1=>0, 2=>1, 3=>-1, 4=>x4_base, 5=>x5_base, 6=>x6_base, 7=>x7_base)

f4 = build_chy_poly(4, sigmas_base, s_base)
f5 = build_chy_poly(5, sigmas_base, s_base)
f6 = build_chy_poly(6, sigmas_base, s_base)
f7 = build_chy_poly(7, sigmas_base, s_base)

F_base = System([f4, f5, f6, f7])
result_base = solve(F_base)

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

physical_sols = filter(is_physical, solutions(result_base))
sorted_sols = sort(physical_sols, by = sol -> (real(sol[1]), imag(sol[1]), real(sol[2])))
println("Found $(length(sorted_sols)) physical solutions")

# Conjugation permutation
function find_conj_perm(sols)
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

τ = find_conj_perm(sorted_sols)
println("τ (conjugation) extracted")

function match_solutions(starts, ends; tol=1e-4)
    n = length(starts)
    m = length(ends)
    perm = zeros(Int, n)
    
    for i in 1:n
        best_j = 0
        best_dist = Inf
        
        for j in 1:m
            dist = norm(starts[i] .- ends[j])
            if dist < best_dist
                best_dist = dist
                best_j = j
            end
        end
        if best_dist < tol
            perm[i] = best_j
        else
            perm[i] = -1  # No match
        end
    end
    return perm
end

function cycle_decomposition(perm)
    n = length(perm)
    visited = falses(n)
    cycles = []
    
    for i in 1:n
        if visited[i] || perm[i] <= 0 continue end
        
        cycle = [i]
        visited[i] = true
        j = perm[i]
        
        while j != i && j > 0 && !visited[j]
            push!(cycle, j)
            visited[j] = true
            j = perm[j]
        end
        
        if length(cycle) > 1
            push!(cycles, cycle)
        end
    end
    
    return cycles
end

# Try many random loops
println("\n" * "="^70)
println("TRYING RANDOM TRIANGLE LOOPS")
println("="^70)

found_generator = false
best_σ = nothing
best_order = 1

for trial in 1:20
    Random.seed!(trial * 12345)
    
    # Random complex perturbations - larger variance
    scale = 1000.0 + trial * 500.0
    p_b = p0 .+ randn(ComplexF64, 3) .* scale
    p_c = p0 .+ randn(ComplexF64, 3) .* scale
    
    print("Trial $trial (scale=$scale): ")
    
    try
        # A→B→C→A
        res_AB = solve(G, sorted_sols; start_parameters=p0, target_parameters=p_b)
        sols_B = solutions(res_AB)
        
        if length(sols_B) < 24
            println("A→B: $(length(sols_B))/24")
            continue
        end
        
        res_BC = solve(G, sols_B; start_parameters=p_b, target_parameters=p_c)
        sols_C = solutions(res_BC)
        
        if length(sols_C) < 24
            println("B→C: $(length(sols_C))/24")
            continue
        end
        
        res_CA = solve(G, sols_C; start_parameters=p_c, target_parameters=p0)
        sols_final = solutions(res_CA)
        
        if length(sols_final) < 24
            println("C→A: $(length(sols_final))/24")
            continue
        end
        
        # Match permutation
        σ = match_solutions(sorted_sols, sols_final)
        
        if any(σ[i] <= 0 for i in 1:24)
            println("incomplete matching")
            continue
        end
        
        # Check if trivial
        is_identity = all(σ[i] == i for i in 1:24)
        is_tau = all(σ[i] == τ[i] for i in 1:24)
        
        if is_identity
            println("identity")
            continue
        elseif is_tau
            println("= τ")
            continue
        end
        
        # Compute order
        local order = 1
        local current = copy(σ)
        while !all(current[i] == i for i in 1:24) && order < 100
            local next = [current[σ[i]] for i in 1:24]
            current = next
            order += 1
        end
        
        cycles = cycle_decomposition(σ)
        cycle_lengths = sort([length(c) for c in cycles], rev=true)
        
        println("NEW! order=$order, cycles=$cycle_lengths")
        
        if order > best_order
            best_order = order
            best_σ = σ
            found_generator = true
        end
        
        # If we find order 24, we're done!
        if order == 24
            println("\n✓✓✓ Found ORDER 24 element! ✓✓✓")
            break
        end
        
    catch e
        println("error: $e")
    end
end

# Results
println("\n" * "="^70)
println("RESULTS")
println("="^70)

if found_generator && best_σ !== nothing
    println("Best generator found:")
    println("  Order: $best_order")
    
    cycles = cycle_decomposition(best_σ)
    println("  Cycle type: ", sort([length(c) for c in cycles], rev=true))
    
    # Check dihedral relation
    σ = best_σ
    σ_inv = zeros(Int, 24)
    for i in 1:24
        σ_inv[σ[i]] = i
    end
    
    τστ = [τ[σ[τ[i]]] for i in 1:24]
    is_dihedral = all(τστ[i] == σ_inv[i] for i in 1:24)
    
    println("  Dihedral relation τστ = σ⁻¹: ", is_dihedral)
    
    if is_dihedral
        println("\n✓✓✓ DIHEDRAL STRUCTURE CONFIRMED! ✓✓✓")
        println("G₇ has dihedral structure")
        
        if best_order == 24
            println("G₇ ≅ D₂₄ (EVIDENCE FOUND)")
        else
            println("G₇ contains D_{$(best_order)} as subgroup")
        end
    end
else
    println("No non-trivial generator found.")
    println("The loops may not encircle discriminant components.")
    println("\nAlternative: Use Frobenius cycle statistics at split primes")
    println("to verify D₂₄ via Chebotarev density.")
end

println("\n" * "="^70)
println("DONE")
println("="^70)
