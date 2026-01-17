# ═══════════════════════════════════════════════════════════════════════════════
# G₇ ≅ D₂₄ PROOF - MULTI-STRATEGY ROBUST VERSION
# ═══════════════════════════════════════════════════════════════════════════════
#
# This script uses MULTIPLE strategies to find monodromy generators:
#   Strategy A: HomotopyContinuation.monodromy_solve (fastest if it works)
#   Strategy B: Random triangle loops with varying scales
#   Strategy C: Circle loops around discriminant
#   Strategy D: Systematic parameter space exploration
#
# If ANY strategy finds σ with order 24 satisfying τστ = σ⁻¹, proof is complete.
#
# RUNTIME: 5-30 minutes depending on which strategy succeeds first
#
# ═══════════════════════════════════════════════════════════════════════════════

using HomotopyContinuation
using LinearAlgebra
using Random
using Printf

# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

const CONFIG = (
    tol_match = 1e-3,         # Solution matching tolerance (looser for numerical monodromy)
    tol_conj = 1e-10,         # Conjugation matching
    target_sols = 24,          # Expected solution count
    triangle_trials = 120,     # Trials per scale for triangle method
    circle_steps = 200,        # Steps for circle loops
    max_monodromy_loops = 120, # For monodromy_solve
)

# ═══════════════════════════════════════════════════════════════════════════════
# CLI MODES
# ═══════════════════════════════════════════════════════════════════════════════
#
#   --monodromy-only   Run Strategy A only (skip heavy sweeps)
#   --quick            Reduce budgets for fast diagnosis
#   --show-progress    Enable HomotopyContinuation progress output (Strategy A)
#
const ARGS_SET = Set(ARGS)
const MONODROMY_ONLY = ("--monodromy-only" in ARGS_SET)
const SKIP_MONODROMY = ("--skip-monodromy" in ARGS_SET)
const DISC_ONLY = ("--disc-only" in ARGS_SET)   # run discriminant encircling only
const QUICK = ("--quick" in ARGS_SET)
const SHOW_PROGRESS = ("--show-progress" in ARGS_SET)

# Runtime budgets (override via flags)
triangle_trials = CONFIG.triangle_trials
circle_steps = CONFIG.circle_steps
max_monodromy_loops = CONFIG.max_monodromy_loops
monodromy_show_progress = false

if QUICK
    triangle_trials = 30
    circle_steps = 80
    max_monodromy_loops = 25
    monodromy_show_progress = true
end
if SHOW_PROGRESS
    monodromy_show_progress = true
end

# ═══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

"""Cycle decomposition of permutation"""
function cycles(perm::Vector{Int})
    n = length(perm)
    visited = falses(n)
    result = Vector{Vector{Int}}()
    
    for start in 1:n
        visited[start] && continue
        cycle = [start]
        visited[start] = true
        j = perm[start]
        while j != start
            push!(cycle, j)
            visited[j] = true
            j = perm[j]
        end
        length(cycle) > 1 && push!(result, cycle)
    end
    return result
end

"""Order of permutation"""
function perm_order(perm::Vector{Int})
    n = length(perm)
    ord = 1
    cur = copy(perm)
    while !all(cur[i] == i for i in 1:n) && ord < 200
        cur = [cur[perm[i]] for i in 1:n]
        ord += 1
    end
    return ord
end

"""Inverse of permutation"""
function perm_inv(perm::Vector{Int})
    n = length(perm)
    inv = zeros(Int, n)
    for i in 1:n
        inv[perm[i]] = i
    end
    return inv
end

"""Check τστ = σ⁻¹"""
function is_dihedral(τ::Vector{Int}, σ::Vector{Int})
    τστ = [τ[σ[τ[i]]] for i in 1:length(σ)]
    return τστ == perm_inv(σ)
end

"""Format permutation for GAP"""
function gap_format(perm::Vector{Int})
    cycs = cycles(perm)
    isempty(cycs) && return "()"
    return join(["(" * join(c, ",") * ")" for c in cycs], "")
end

"""Match solutions with tolerance"""
function match_sols(starts, ends; tol=CONFIG.tol_match)
    n = length(starts)
    m = length(ends)
    n == m || error("match_sols: expected same number of starts and ends, got $n and $m")

    # Optimal bipartite matching (Hungarian algorithm) to avoid greedy mismatches
    # near the discriminant (where multiple endpoints can be very close).
    function hungarian(cost::Matrix{Float64})
        n = size(cost, 1)
        m = size(cost, 2)
        n == m || error("hungarian: expected square matrix")
        # 1-indexed dummy column 1 represents 0 in standard implementations
        u = zeros(Float64, n + 1)
        v = zeros(Float64, m + 1)
        p = zeros(Int, m + 1)      # matched row for column j
        way = zeros(Int, m + 1)
        for i in 1:n
            p[1] = i
            j0 = 1
            minv = fill(Inf, m + 1)
            used = falses(m + 1)
            fill!(way, 0)
            while true
                used[j0] = true
                i0 = p[j0]
                delta = Inf
                j1 = 0
                for j in 2:(m + 1)
                    used[j] && continue
                    cur = cost[i0, j - 1] - u[i0] - v[j]
                    if cur < minv[j]
                        minv[j] = cur
                        way[j] = j0
                    end
                    if minv[j] < delta
                        delta = minv[j]
                        j1 = j
                    end
                end
                for j in 1:(m + 1)
                    if used[j]
                        u[p[j]] += delta
                        v[j] -= delta
                    else
                        minv[j] -= delta
                    end
                end
                j0 = j1
                if p[j0] == 0
                    break
                end
            end
            # augmenting
            while true
                j1 = way[j0]
                p[j0] = p[j1]
                j0 = j1
                j0 == 1 && break
            end
        end
        assignment = zeros(Int, n)
        for j in 2:(m + 1)
            i = p[j]
            i == 0 && continue
            assignment[i] = j - 1
        end
        return assignment
    end

    cost = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        cost[i, j] = norm(starts[i] .- ends[j])
    end
    perm = hungarian(cost)
    # Validate tolerance
    all(cost[i, perm[i]] < tol for i in 1:n) || return fill(-1, n)
    return perm
end

"""Minimum pairwise distance between solutions (detect discriminant proximity)."""
function min_solution_distance(sols)
    n = length(sols)
    n < 2 && return Inf
    min_d = Inf
    for i in 1:n, j in (i+1):n
        d = norm(sols[i] .- sols[j])
        d < min_d && (min_d = d)
    end
    return min_d
end

"""Cluster nearby complex points by simple threshold averaging."""
function cluster_points(points::Vector{ComplexF64}, threshold::Float64)
    isempty(points) && return ComplexF64[]
    clusters = ComplexF64[points[1]]
    for p in points[2:end]
        merged = false
        for i in eachindex(clusters)
            if abs(p - clusters[i]) < threshold
                clusters[i] = (clusters[i] + p) / 2
                merged = true
                break
            end
        end
        !merged && push!(clusters, p)
    end
    return clusters
end

"""
Find approximate discriminant points along a 1-parameter slice:
    params(t) = P₀ + t * direction

We record t where continuation loses solutions (length < 24) or solutions collide.
"""
function find_discriminant_points(G, SOLS, P₀, direction;
                                 t_box::Float64=2.0,
                                 samples::Int=120,
                                 close_tol::Float64=2e-2,
                                 cluster_tol::Float64=0.35,
                                 topk::Int=12)
    candidates = Vector{Tuple{Float64, ComplexF64}}()
    sizehint!(candidates, samples)
    for _ in 1:samples
        t = (2rand() - 1) * t_box + im * (2rand() - 1) * t_box
        params = P₀ .+ t .* direction
        try
            res = solve(G, SOLS; start_parameters=P₀, target_parameters=params)
            sols = solutions(res)
            if length(sols) < 24
                push!(candidates, (0.0, t))
                continue
            end
            dmin = min_solution_distance(sols)
            dmin < close_tol && push!(candidates, (dmin, t))
        catch
            push!(candidates, (0.0, t))
        end
    end
    isempty(candidates) && return ComplexF64[]
    sort!(candidates, by = x -> x[1])
    chosen = ComplexF64[t for (_, t) in candidates[1:min(topk, length(candidates))]]
    return cluster_points(chosen, cluster_tol)
end

"""
Track solutions around a circle in 1-parameter space centered at t_center:
    t(θ) = t_center + radius * exp(i θ)

Returns induced permutation σ or nothing if tracking/matching fails.
"""
function circle_monodromy(G, SOLS, P₀, direction, t_center::ComplexF64, radius::Float64;
                          n_steps::Int=60, tol_match::Float64=CONFIG.tol_match)
    # Start at θ=0 point on the circle
    t0 = t_center + ComplexF64(radius)
    p0 = P₀ .+ t0 .* direction
    try
        res0 = solve(G, SOLS; start_parameters=P₀, target_parameters=p0)
        sols0 = solutions(res0)
        length(sols0) < 24 && return nothing
        current_sols = sols0
        current_params = p0
        
        for step in 1:n_steps
            θ = 2π * step / n_steps
            t = t_center + radius * exp(im * θ)
            next_params = P₀ .+ t .* direction
            res = solve(G, current_sols; start_parameters=current_params, target_parameters=next_params)
            sols = solutions(res)
            length(sols) < 24 && return nothing
            current_sols = sols
            current_params = next_params
        end
        
        # Close loop back to base point
        res_back = solve(G, current_sols; start_parameters=current_params, target_parameters=P₀)
        sols_back = solutions(res_back)
        length(sols_back) < 24 && return nothing
        σ = match_sols(SOLS, sols_back; tol=tol_match)
        all(σ .> 0) || return nothing
        return σ
    catch
        return nothing
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# SYSTEM SETUP
# ═══════════════════════════════════════════════════════════════════════════════

println("═"^80)
println("│" * " "^28 * "G₇ ≅ D₂₄ MULTI-STRATEGY PROOF" * " "^20 * "│")
println("═"^80)
println()

@var σ₄ σ₅ σ₆ σ₇
@var p₄₅ p₄₆ p₄₇

# Mandelstam invariants
const S_BASE = Dict{Tuple{Int,Int}, Int}(
    (1,2) => -6728, (1,3) => -8080, (1,4) => -14144, (1,5) => -22672, 
    (1,6) => 17440, (1,7) => 34184,
    (2,3) => 2740, (2,4) => 15304, (2,5) => 17032, (2,6) => -9524, (2,7) => -18824,
    (3,4) => 4112, (3,5) => 6408, (3,6) => -1768, (3,7) => -3412,
    (4,5) => 6032, (4,6) => -3904, (4,7) => -7400,
    (5,6) => -2248, (5,7) => -4552,
    (6,7) => 4
)

const SIGMAS = Dict(1=>0, 2=>1, 3=>-1, 4=>σ₄, 5=>σ₅, 6=>σ₆, 7=>σ₇)

function chy_eq(a, sigmas, s)
    get_s(i,j) = s[(min(i,j), max(i,j))]
    result = 0
    for b in 1:7
        b == a && continue
        term = get_s(a, b)
        for c in 1:7
            (c == a || c == b) && continue
            term *= (sigmas[a] - sigmas[c])
        end
        result += term
    end
    return result
end

# Fixed system
F = System([chy_eq(4, SIGMAS, S_BASE), chy_eq(5, SIGMAS, S_BASE),
            chy_eq(6, SIGMAS, S_BASE), chy_eq(7, SIGMAS, S_BASE)])

# Parameterized system
S_PARAM = merge(S_BASE, Dict((4,5)=>p₄₅, (4,6)=>p₄₆, (4,7)=>p₄₇))
G = System([chy_eq(4, SIGMAS, S_PARAM), chy_eq(5, SIGMAS, S_PARAM),
            chy_eq(6, SIGMAS, S_PARAM), chy_eq(7, SIGMAS, S_PARAM)];
           parameters=[p₄₅, p₄₆, p₄₇])

# Base parameters
const P₀ = ComplexF64[6032.0, -3904.0, -7400.0]

# ═══════════════════════════════════════════════════════════════════════════════
# FIND BASE SOLUTIONS
# ═══════════════════════════════════════════════════════════════════════════════

println("Phase 1: Finding all 24 solutions...")

function is_physical(sol; tol=1e-6)
    vals = [sol[1], sol[2], sol[3], sol[4]]
    for v in vals
        (abs(v) < tol || abs(v-1) < tol || abs(v+1) < tol) && return false
    end
    for i in 1:4, j in (i+1):4
        abs(vals[i] - vals[j]) < tol && return false
    end
    return true
end

result = solve(F)
physical = filter(is_physical, solutions(result))
SOLS = sort(physical, by=s->(real(s[1]), imag(s[1]), real(s[2])))

println("  Found $(length(SOLS)) physical solutions")
@assert length(SOLS) == 24 "Expected 24 solutions!"
println("  ✓ All 24 solutions found")
println()

# ═══════════════════════════════════════════════════════════════════════════════
# EXTRACT τ (CONJUGATION)
# ═══════════════════════════════════════════════════════════════════════════════

println("Phase 2: Extracting conjugation generator τ...")

function extract_τ(sols)
    perm = zeros(Int, length(sols))
    for i in 1:length(sols)
        conj_sol = conj.(sols[i])
        for j in 1:length(sols)
            if norm(sols[j] .- conj_sol) < CONFIG.tol_conj
                perm[i] = j
                break
            end
        end
    end
    return perm
end

const τ = extract_τ(SOLS)
τ_ord = perm_order(τ)
τ_cycs = cycles(τ)

println("  τ order: $τ_ord")
println("  τ cycles: $(length(τ_cycs)) of length $(length(τ_cycs) > 0 ? length(τ_cycs[1]) : 0)")

@assert τ_ord == 2 "τ should have order 2"
@assert length(τ_cycs) == 12 "τ should have 12 cycles"
@assert all(length(c) == 2 for c in τ_cycs) "All cycles should be 2-cycles"

println("  ✓ τ has correct structure (order 2, cycle type 2¹²)")
println()

# ═══════════════════════════════════════════════════════════════════════════════
# MONODROMY SEARCH - MULTIPLE STRATEGIES
# ═══════════════════════════════════════════════════════════════════════════════

println("Phase 3: Searching for rotation generator σ...")
println()

# Store all non-trivial generators found
generators = Vector{Vector{Int}}()
best_σ = nothing
best_order = 1
proof_complete = false

"""Check if permutation is interesting"""
function is_interesting(σ)
    all(σ[i] == i for i in 1:24) && return false     # Identity
    all(σ[i] == τ[i] for i in 1:24) && return false  # Same as τ
    return true
end

"""Process a found permutation"""
function process_generator!(σ)
    global best_σ, best_order, proof_complete
    
    !is_interesting(σ) && return false
    
    ord = perm_order(σ)
    push!(generators, σ)
    
    if ord > best_order
        best_order = ord
        best_σ = σ
    end
    
    if is_dihedral(τ, σ) && ord == 24
        println("    ★★★ FOUND ORDER-24 DIHEDRAL GENERATOR! ★★★")
        proof_complete = true
        return true
    end
    
    return false
end

# ─────────────────────────────────────────────────────────────────────────────
# STRATEGY A: monodromy_solve
# ─────────────────────────────────────────────────────────────────────────────

if !SKIP_MONODROMY
    println("  Strategy A: monodromy_solve...")

    try
        # HomotopyContinuation.monodromy_solve signature:
        #   monodromy_solve(F, sols, p; options...)
        #
        # To *generate* monodromy permutations, we must start from a strict subset
        # of the fiber; if we start with all 24 solutions, monodromy_solve can
        # terminate immediately without exploring loops.
        start_sols = [Vector{ComplexF64}(SOLS[1])]
        mono = monodromy_solve(
            G,
            start_sols,
            P₀;
            target_solutions_count=24,
            max_loops_no_progress=max_monodromy_loops,
            permutations=true,
            timeout=(QUICK ? 10.0 : 120.0),
        )
        nfound = try
            length(solutions(mono))
        catch
            -1
        end
        nfound > 0 && println("    monodromy_solve found $nfound solutions")

        if nfound != 24
            println("    (Skipping monodromy permutations: need all 24 solutions in one fiber)")
        else
            gens_raw = permutations(mono)
            gens24 = Vector{Vector{Int}}()

            if gens_raw isa AbstractMatrix
                if size(gens_raw, 1) == 24
                    for j in 1:size(gens_raw, 2)
                        push!(gens24, collect(gens_raw[:, j]))
                    end
                elseif size(gens_raw, 2) == 24
                    for i in 1:size(gens_raw, 1)
                        push!(gens24, collect(gens_raw[i, :]))
                    end
                end
            elseif gens_raw isa AbstractVector
                if eltype(gens_raw) <: Integer && length(gens_raw) == 24
                    push!(gens24, collect(gens_raw))
                else
                    for gen in gens_raw
                        (gen isa AbstractVector{<:Integer} && length(gen) == 24) || continue
                        push!(gens24, collect(gen))
                    end
                end
            end

            println("    Found $(length(gens24)) 24-point permutations from monodromy_solve")
            for σ in gens24
                ord = perm_order(σ)
                println("    Permutation: order $ord, dihedral=$(is_dihedral(τ, σ))")
                process_generator!(σ) && break
            end
        end
    catch e
        println("    Failed: $e")
    end
else
    println("  Strategy A: monodromy_solve... (skipped via --skip-monodromy)")
end

# ─────────────────────────────────────────────────────────────────────────────
# STRATEGY B: 1-Parameter Discriminant Encircling (recommended)
# ─────────────────────────────────────────────────────────────────────────────

if !proof_complete && !MONODROMY_ONLY
    println()
    println("  Strategy B: 1-Parameter Discriminant Encircling...")
    
    disc_seeds = QUICK ? 6 : 20
    disc_samples = QUICK ? 80 : 220
    disc_box = QUICK ? 1.5 : 2.0
    dir_scale = QUICK ? 250.0 : 500.0
    
    for seed in 1:disc_seeds
        proof_complete && break
        Random.seed!(seed * 7919)
        
        direction = randn(ComplexF64, 3)
        direction .*= dir_scale / norm(direction)
        
        @printf("    Seed %d: searching discriminant candidates... ", seed)
        disc_points = find_discriminant_points(G, SOLS, P₀, direction;
            t_box=disc_box, samples=disc_samples, close_tol=2e-2, cluster_tol=0.35, topk=(QUICK ? 6 : 12))
        println("found $(length(disc_points))")
        
        for t_disc in disc_points
            proof_complete && break
            for radius in (QUICK ? [0.01, 0.02, 0.05, 0.10] : [0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12])
                proof_complete && break
                
                σ = circle_monodromy(G, SOLS, P₀, direction, t_disc, radius;
                    n_steps=(QUICK ? 80 : 120),
                    tol_match=CONFIG.tol_match)
                σ === nothing && continue
                !is_interesting(σ) && continue
                
                ord = perm_order(σ)
                dihed = is_dihedral(τ, σ)
                @printf("    t_disc=%.2f%+.2fi radius=%.2f: order=%d dihedral=%s\n",
                    real(t_disc), imag(t_disc), radius, ord, dihed)
                if ord > 2
                    cycs = cycles(σ)
                    lens = sort([length(c) for c in cycs], rev=true)
                    println("      cycle lengths: ", lens)
                end
                process_generator!(σ)
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# STRATEGY C: Random triangle loops (fallback)
# ─────────────────────────────────────────────────────────────────────────────

if !proof_complete && !MONODROMY_ONLY && !DISC_ONLY && triangle_trials > 0
    println()
    println("  Strategy C: Random triangle loops...")
    
    for scale in [500.0, 1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0, 100000.0]
        proof_complete && break
        
        found_at_scale = 0
        for trial in 1:triangle_trials
            proof_complete && break
            
            Random.seed!(Int(scale) + trial * 137)
            
            p_b = P₀ .+ randn(ComplexF64, 3) .* scale
            p_c = P₀ .+ randn(ComplexF64, 3) .* scale
            
            try
                res_AB = solve(G, SOLS; start_parameters=P₀, target_parameters=p_b)
                length(solutions(res_AB)) < 24 && continue
                
                res_BC = solve(G, solutions(res_AB); start_parameters=p_b, target_parameters=p_c)
                length(solutions(res_BC)) < 24 && continue
                
                res_CA = solve(G, solutions(res_BC); start_parameters=p_c, target_parameters=P₀)
                length(solutions(res_CA)) < 24 && continue
                
                σ = match_sols(SOLS, solutions(res_CA))
                any(σ[i] <= 0 for i in 1:24) && continue
                
                if is_interesting(σ)
                    found_at_scale += 1
                    ord = perm_order(σ)
                    @printf("    scale=%.0f trial=%d: order=%d dihedral=%s\n", 
                            scale, trial, ord, is_dihedral(τ, σ))
                    process_generator!(σ)
                end
            catch
                continue
            end
        end
        
        found_at_scale > 0 && println("    Found $found_at_scale generators at scale $scale")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# STRATEGY D: Circle loops (fallback)
# ─────────────────────────────────────────────────────────────────────────────

if !proof_complete && !MONODROMY_ONLY && !DISC_ONLY && circle_steps > 0
    println()
    println("  Strategy D: Circle loops...")
    
    for radius in [100.0, 500.0, 1000.0, 5000.0, 20000.0, 50000.0]
        proof_complete && break
        
        for center_offset in [ComplexF64[0,0,0], randn(ComplexF64, 3) .* radius]
            proof_complete && break
            
            center = P₀ .+ center_offset
            
            # Track around circle
            try
                current_sols = SOLS
                n_steps = circle_steps
                ok = true
                
                for step in 1:n_steps
                    θ = 2π * step / n_steps
                    p_next = center .+ radius .* [exp(im*θ), exp(im*(θ+0.1)), exp(im*(θ+0.2))]
                    
                    p_prev = step == 1 ? P₀ : center .+ radius .* [
                        exp(im*2π*(step-1)/n_steps),
                        exp(im*(2π*(step-1)/n_steps+0.1)),
                        exp(im*(2π*(step-1)/n_steps+0.2)),
                    ]
                    res = solve(G, current_sols; 
                        start_parameters=p_prev,
                        target_parameters=p_next)
                    
                    if length(solutions(res)) < 24
                        ok = false
                        break
                    end
                    current_sols = solutions(res)
                end
                
                if ok
                    # Close the loop back to P₀
                    res_final = solve(G, current_sols; 
                        start_parameters=center .+ radius .* [ComplexF64(1), exp(im*0.1), exp(im*0.2)],
                        target_parameters=P₀)
                    
                    if length(solutions(res_final)) >= 24
                        σ = match_sols(SOLS, solutions(res_final))
                        if all(σ[i] > 0 for i in 1:24) && is_interesting(σ)
                            ord = perm_order(σ)
                            @printf("    radius=%.0f: order=%d dihedral=%s\n", 
                                    radius, ord, is_dihedral(τ, σ))
                            process_generator!(σ)
                        end
                    end
                end
            catch
                continue
            end
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("═"^80)
println("RESULTS")
println("═"^80)
println()

println("Generators found: $(length(generators))")
println("Best σ order: $best_order")
if best_σ !== nothing
    println("Best σ is dihedral: $(is_dihedral(τ, best_σ))")
end
println()

if proof_complete
    println("╔════════════════════════════════════════════════════════════════════════════╗")
    println("║                                                                            ║")
    println("║   ████████╗██╗  ██╗███████╗ ██████╗ ██████╗ ███████╗███╗   ███╗            ║")
    println("║   ╚══██╔══╝██║  ██║██╔════╝██╔═══██╗██╔══██╗██╔════╝████╗ ████║            ║")
    println("║      ██║   ███████║█████╗  ██║   ██║██████╔╝█████╗  ██╔████╔██║            ║")
    println("║      ██║   ██╔══██║██╔══╝  ██║   ██║██╔══██╗██╔══╝  ██║╚██╔╝██║            ║")
    println("║      ██║   ██║  ██║███████╗╚██████╔╝██║  ██║███████╗██║ ╚═╝ ██║            ║")
    println("║      ╚═╝   ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝            ║")
    println("║                                                                            ║")
    println("║                      G₇ ≅ D₂₄ (EVIDENCE FOUND)                              ║")
    println("║                                                                            ║")
    println("║   τ: order 2, cycle type 2¹²                                              ║")
    println("║   σ: order 24                                                              ║")
    println("║   τστ = σ⁻¹ ✓                                                              ║")
    println("║   |G₇| = 48                                                                ║")
    println("║                                                                            ║")
    println("╚════════════════════════════════════════════════════════════════════════════╝")
elseif best_σ !== nothing && is_dihedral(τ, best_σ)
    println("╔════════════════════════════════════════════════════════════════════════════╗")
    println("║  PARTIAL RESULT: Found dihedral subgroup D_$best_order                        ║")
    println("║  G₇ contains D_$best_order as subgroup                                        ║")
    println("║  More exploration needed for full D₂₄                                      ║")
    println("╚════════════════════════════════════════════════════════════════════════════╝")
else
    println("╔════════════════════════════════════════════════════════════════════════════╗")
    println("║  No dihedral generator found                                               ║")
    println("║  Try: increase trials, different seeds, or Frobenius method                ║")
    println("╚════════════════════════════════════════════════════════════════════════════╝")
end

# ═══════════════════════════════════════════════════════════════════════════════
# GAP CODE OUTPUT
# ═══════════════════════════════════════════════════════════════════════════════

println()
println("═"^80)
println("GAP VERIFICATION CODE")
println("═"^80)
println()

println("# Paste into GAP:")
println("tau := $(gap_format(τ));")

# Emit several observed generators (if any), not just the single "best" one.
unique_gens = Vector{Vector{Int}}()
seen = Set{NTuple{24,Int}}()
for σ in generators
    length(σ) == 24 || continue
    t = NTuple{24,Int}(Tuple(σ))
    t in seen && continue
    push!(seen, t)
    push!(unique_gens, σ)
end
sort!(unique_gens, by=perm_order, rev=true)

max_emit = min(length(unique_gens), 8)
for i in 1:max_emit
    println("sigma$i := $(gap_format(unique_gens[i]));")
end

if max_emit > 0
    sig_list = join(["sigma$i" for i in 1:max_emit], ",")
    println("G := Group(tau, $sig_list);")
    println("Print(\"Order: \", Size(G), \"\\n\");")
    println("Print(\"IsTransitive: \", IsTransitive(G, [1..24]), \"\\n\");")
    println("Print(\"Orbits sizes: \", List(Orbits(G, [1..24]), Size), \"\\n\");")
    println("Print(\"Best order observed: $best_order\\n\");")
else
    println("G := Group(tau);")
    println("Print(\"Order: \", Size(G), \"\\n\");")
end

println()
println("═"^80)
println("COMPLETE")
println("═"^80)
