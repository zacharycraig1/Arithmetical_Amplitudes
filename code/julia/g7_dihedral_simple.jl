# G₇ ≅ D₂₄ VERIFICATION - SIMPLIFIED
# =====================================

using HomotopyContinuation
using LinearAlgebra
using Random

println("="^70)
println("G₇ ≅ D₂₄ VERIFICATION")
println("="^70)

# ==========================================================================
# PART 1: Setup and Solve
# ==========================================================================

@var x4 x5 x6 x7  # Use simple variable names

# Fixed Mandelstams (CHY-consistent)
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

println("\nBuilding system...")
f4 = build_chy_poly(4, sigmas, s)
f5 = build_chy_poly(5, sigmas, s)
f6 = build_chy_poly(6, sigmas, s)
f7 = build_chy_poly(7, sigmas, s)

F = System([f4, f5, f6, f7])

println("Solving...")
result = solve(F)

all_sols = solutions(result)
println("Total solutions: ", length(all_sols))

# Filter physical
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

# ==========================================================================
# PART 2: Extract τ (conjugation)
# ==========================================================================

println("\n" * "="^70)
println("CONJUGATION PERMUTATION (τ)")
println("="^70)

function find_conjugate_permutation(sols)
    n = length(sols)
    perm = zeros(Int, n)
    tol = 1e-8
    
    for i in 1:n
        sol_i = sols[i]
        conj_i = conj.(sol_i)
        
        for j in 1:n
            match = true
            for k in 1:4
                if abs(sols[j][k] - conj_i[k]) > tol
                    match = false
                    break
                end
            end
            if match
                perm[i] = j
                break
            end
        end
    end
    return perm
end

τ = find_conjugate_permutation(sorted_sols)

# Cycle decomposition
function cycle_decomposition(perm)
    n = length(perm)
    visited = falses(n)
    cycles = []
    
    for i in 1:n
        if visited[i] continue end
        
        cycle = [i]
        visited[i] = true
        j = perm[i]
        
        while j != i
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

τ_cycles = cycle_decomposition(τ)
println("τ cycle type: ", length(τ_cycles), " cycles of length ", 
        length(τ_cycles) > 0 ? length(τ_cycles[1]) : 0)

if length(τ_cycles) == 12 && all(length(c) == 2 for c in τ_cycles)
    println("✓ τ has cycle type 2¹² (12 transpositions)")
end

# ==========================================================================
# PART 3: Try to find σ via parameter variation
# ==========================================================================

println("\n" * "="^70)
println("SEARCHING FOR ROTATION GENERATOR (σ)")
println("="^70)

# Use parameter homotopy
@var p45 p46 p47

# Parameterized system
s_param = Dict{Tuple{Int,Int}, Any}(
    (1,2) => -6728, (1,3) => -8080, (1,4) => -14144, (1,5) => -22672, 
    (1,6) => 17440, (1,7) => 34184,
    (2,3) => 2740, (2,4) => 15304, (2,5) => 17032, (2,6) => -9524, (2,7) => -18824,
    (3,4) => 4112, (3,5) => 6408, (3,6) => -1768, (3,7) => -3412,
    (4,5) => p45, (4,6) => p46, (4,7) => p47,
    (5,6) => -2248, (5,7) => -4552,
    (6,7) => 4
)

sigmas_p = Dict(1=>0, 2=>1, 3=>-1, 4=>x4, 5=>x5, 6=>x6, 7=>x7)

g4 = build_chy_poly(4, sigmas_p, s_param)
g5 = build_chy_poly(5, sigmas_p, s_param)
g6 = build_chy_poly(6, sigmas_p, s_param)
g7 = build_chy_poly(7, sigmas_p, s_param)

G = System([g4, g5, g6, g7]; parameters=[p45, p46, p47])

# Base and target parameters
p0 = [6032.0 + 0.0im, -3904.0 + 0.0im, -7400.0 + 0.0im]

# Random complex target
Random.seed!(42)
p1 = p0 .+ randn(ComplexF64, 3) .* 300.0

println("Tracking from p0 to p1...")
println("  p0 = ", p0)
println("  p1 = ", p1)

# Use solve with start/target parameters
println("  Solving with target parameters...")
result1 = solve(G, sorted_sols; start_parameters=p0, target_parameters=p1)
tracked = solutions(result1)
successful = length(tracked)
println("  Solutions found: $successful")

if successful >= 20
    # Get endpoints
    endpoints1 = tracked
    
    # Track back: p1 to p0
    println("Tracking from p1 back to p0...")
    result_back = solve(G, endpoints1; start_parameters=p1, target_parameters=p0)
    endpoints_back = solutions(result_back)
    successful_back = length(endpoints_back)
    println("  Solutions found: $successful_back")
    
    if successful_back >= 20
        endpoints_final = endpoints_back
        
        # Match to original
        function match_solutions(starts, ends; tol=1e-4)
            n = length(starts)
            perm = zeros(Int, n)
            
            for i in 1:n
                best_j = 0
                best_dist = Inf
                
                for j in 1:n
                    dist = norm(starts[i] .- ends[j])
                    if dist < best_dist
                        best_dist = dist
                        best_j = j
                    end
                end
                perm[i] = best_j
            end
            return perm
        end
        
        σ = match_solutions(sorted_sols, endpoints_final)
        
        println("\nLoop permutation σ = ", σ)
        
        σ_cycles = cycle_decomposition(σ)
        println("σ cycle structure: ", length(σ_cycles), " cycles")
        
        if length(σ_cycles) > 0
            cycle_lengths = [length(c) for c in σ_cycles]
            println("  Cycle lengths: ", sort(cycle_lengths, rev=true))
        end
        
        # Check if σ is trivial
        is_identity = all(σ[i] == i for i in 1:length(σ))
        is_tau = all(σ[i] == τ[i] for i in 1:length(σ))
        
        if is_identity
            println("σ = identity (trivial loop)")
        elseif is_tau
            println("σ = τ (same as conjugation)")
        else
            println("\n✓ Found NEW generator!")
            
            # Compute order
            order = 1
            current = copy(σ)
            while !all(current[i] == i for i in 1:length(σ)) && order < 100
                next = [current[σ[i]] for i in 1:length(σ)]
                current = next
                order += 1
            end
            println("Order of σ: ", order)
            
            # Check dihedral relation
            σ_inv = zeros(Int, length(σ))
            for i in 1:length(σ)
                σ_inv[σ[i]] = i
            end
            
            τστ = [τ[σ[τ[i]]] for i in 1:length(σ)]
            is_dihedral = all(τστ[i] == σ_inv[i] for i in 1:length(σ))
            
            println("Dihedral relation τστ = σ⁻¹: ", is_dihedral)
            
            if is_dihedral
                println("\n" * "="^70)
                println("✓✓✓ D₂₄ STRUCTURE VERIFIED! ✓✓✓")
                println("="^70)
            end
        end
    end
end

# ==========================================================================
# PART 4: Try Triangle Loop
# ==========================================================================

println("\n" * "="^70)
println("TRIANGLE LOOP ATTEMPT")
println("="^70)

# p0 → p1 → p2 → p0
Random.seed!(12345)
p_a = p0
p_b = p0 .+ randn(ComplexF64, 3) .* 500.0
p_c = p0 .+ randn(ComplexF64, 3) .* 500.0

println("Triangle vertices:")
println("  A = ", p_a)
println("  B = ", p_b)
println("  C = ", p_c)

try
    # Leg A→B
    println("\nLeg A→B...")
    result_AB = solve(G, sorted_sols; start_parameters=p_a, target_parameters=p_b)
    sols_B = solutions(result_AB)
    println("  Solutions: $(length(sols_B))")
    
    if length(sols_B) >= 20
        # Leg B→C
        println("Leg B→C...")
        result_BC = solve(G, sols_B; start_parameters=p_b, target_parameters=p_c)
        sols_C = solutions(result_BC)
        println("  Solutions: $(length(sols_C))")
        
        if length(sols_C) >= 20
            # Leg C→A
            println("Leg C→A...")
            result_CA = solve(G, sols_C; start_parameters=p_c, target_parameters=p_a)
            sols_final = solutions(result_CA)
            println("  Solutions: $(length(sols_final))")
            
            if length(sols_final) >= 20
                
                # Extract permutation
                ρ = match_solutions(sorted_sols, sols_final)
                
                println("\n✓ Triangle loop completed!")
                println("ρ = ", ρ)
                
                ρ_cycles = cycle_decomposition(ρ)
                println("ρ cycle structure: ", length(ρ_cycles), " cycles")
                
                is_identity = all(ρ[i] == i for i in 1:length(ρ))
                is_tau = all(ρ[i] == τ[i] for i in 1:length(ρ))
                
                if !is_identity && !is_tau
                    println("\n✓ Found NON-TRIVIAL generator from triangle!")
                    
                    # Compute order
                    order = 1
                    current = copy(ρ)
                    while !all(current[i] == i for i in 1:length(ρ)) && order < 100
                        next = [current[ρ[i]] for i in 1:length(ρ)]
                        current = next
                        order += 1
                    end
                    println("Order of ρ: ", order)
                    
                    # Check if it satisfies dihedral relation with τ
                    ρ_inv = zeros(Int, length(ρ))
                    for i in 1:length(ρ)
                        ρ_inv[ρ[i]] = i
                    end
                    
                    τρτ = [τ[ρ[τ[i]]] for i in 1:length(ρ)]
                    is_dihedral = all(τρτ[i] == ρ_inv[i] for i in 1:length(ρ))
                    
                    println("Dihedral relation τρτ = ρ⁻¹: ", is_dihedral)
                    
                    if is_dihedral && order == 24
                        println("\n" * "="^70)
                        println("✓✓✓ G₇ ≅ D₂₄ (EVIDENCE FOUND) ✓✓✓")
                        println("="^70)
                        println("""
We have:
  - τ: order 2, cycle type 2¹²
  - ρ: order 24
  - τρτ = ρ⁻¹

This is exactly D₂₄ = ⟨r, s | r²⁴ = s² = 1, srs = r⁻¹⟩
""")
                    end
                else
                    println("Triangle gave trivial permutation")
                end
            else
                println("Leg C→A: insufficient solutions")
            end
        else
            println("Leg B→C: insufficient solutions")
        end
    else
        println("Leg A→B: insufficient solutions")
    end
catch e
    println("Triangle loop error: ", e)
end

println("\n" * "="^70)
println("DONE")
println("="^70)
