# ═══════════════════════════════════════════════════════════════════════════════
# G₇ MONODROMY VIA FACTORIZATION CHANNEL LOOPS (n=7 CHY)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Motivation (from project handoffs):
# - Random loops in (s45,s46,s47) are unlikely to link the discriminant.
# - Instead, explicitly loop around *channel* discriminant components:
#     s123 := s12+s13+s23 = 0
#     s234 := s23+s24+s34 = 0
#   and compose loops to mix sectors (potentially yielding order-24 elements).
#
# This script is a diagnostic: it reports permutation orders / cycle types.
#
# Run:
#   julia --project=. prove_D24_channels.jl
#   julia --project=. prove_D24_channels.jl --quick
#
# ═══════════════════════════════════════════════════════════════════════════════

using HomotopyContinuation
using LinearAlgebra
using Random
using Printf

const ARGS_SET = Set(ARGS)
const QUICK = ("--quick" in ARGS_SET)

const N = 24
const TOL_CONJ = 1e-10
const TOL_MATCH = QUICK ? 3e-3 : 1e-3

println("═"^80)
println("G₇ MONODROMY - FACTORIZATION CHANNEL LOOPS")
println("═"^80)
println()

# ------------------------------------------------------------------------------
# Permutation utilities
# ------------------------------------------------------------------------------

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

function perm_order(perm::Vector{Int})
    cycs = cycles(perm)
    isempty(cycs) && return 1
    return lcm([length(c) for c in cycs]...)
end

function cycle_type(perm::Vector{Int})
    cycs = cycles(perm)
    n = length(perm)
    moved = sum(length.(cycs); init=0)
    fixed = n - moved
    lens = sort([length(c) for c in cycs], rev=true)
    fixed > 0 && append!(lens, fill(1, fixed))
    return Tuple(sort(lens, rev=true))
end

function perm_inv(perm::Vector{Int})
    inv = zeros(Int, length(perm))
    for i in eachindex(perm)
        inv[perm[i]] = i
    end
    return inv
end

function is_dihedral(τ::Vector{Int}, σ::Vector{Int})
    τστ = [τ[σ[τ[i]]] for i in eachindex(σ)]
    return τστ == perm_inv(σ)
end

function gap_format(perm::Vector{Int})
    cycs = cycles(perm)
    isempty(cycs) && return "()"
    return join(["(" * join(c, ",") * ")" for c in cycs], "")
end

# Hungarian algorithm for optimal matching (square cost matrix)
function hungarian(cost::Matrix{Float64})
    n = size(cost, 1)
    n == size(cost, 2) || error("hungarian expects square cost matrix")
    u = zeros(Float64, n + 1)
    v = zeros(Float64, n + 1)
    p = zeros(Int, n + 1)
    way = zeros(Int, n + 1)
    for i in 1:n
        p[1] = i
        j0 = 1
        minv = fill(Inf, n + 1)
        used = falses(n + 1)
        fill!(way, 0)
        while true
            used[j0] = true
            i0 = p[j0]
            delta = Inf
            j1 = 0
            for j in 2:(n + 1)
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
            for j in 1:(n + 1)
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
        while true
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            j0 == 1 && break
        end
    end
    assignment = zeros(Int, n)
    for j in 2:(n + 1)
        i = p[j]
        i == 0 && continue
        assignment[i] = j - 1
    end
    return assignment
end

function match_sols(starts, ends; tol::Float64=TOL_MATCH)
    n = length(starts)
    length(ends) == n || return fill(-1, n)
    cost = Matrix{Float64}(undef, n, n)
    for i in 1:n, j in 1:n
        cost[i, j] = norm(starts[i] .- ends[j])
    end
    perm = hungarian(cost)
    all(cost[i, perm[i]] < tol for i in 1:n) || return fill(-1, n)
    return perm
end

# ------------------------------------------------------------------------------
# CHY system setup
# ------------------------------------------------------------------------------

@var σ₄ σ₅ σ₆ σ₇
@var q_s12 q_s23

const S_BASE_FULL = Dict{Tuple{Int,Int}, Int}(
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

# Fixed system at the base point
F_fixed = System([
    chy_eq(4, SIGMAS, S_BASE_FULL),
    chy_eq(5, SIGMAS, S_BASE_FULL),
    chy_eq(6, SIGMAS, S_BASE_FULL),
    chy_eq(7, SIGMAS, S_BASE_FULL),
])

# Parameterized system: only s12 and s23 vary (channel loops)
S_LOOP = merge(
    Dict{Tuple{Int,Int}, Any}((k, v) for (k, v) in S_BASE_FULL),
    Dict{Tuple{Int,Int}, Any}((1,2) => q_s12, (2,3) => q_s23),
)

G = System([
    chy_eq(4, SIGMAS, S_LOOP),
    chy_eq(5, SIGMAS, S_LOOP),
    chy_eq(6, SIGMAS, S_LOOP),
    chy_eq(7, SIGMAS, S_LOOP),
]; parameters=[q_s12, q_s23])

const P₀ = ComplexF64[S_BASE_FULL[(1,2)], S_BASE_FULL[(2,3)]]
const S123_BASE = S_BASE_FULL[(1,2)] + S_BASE_FULL[(1,3)] + S_BASE_FULL[(2,3)]
const S234_BASE = S_BASE_FULL[(2,3)] + S_BASE_FULL[(2,4)] + S_BASE_FULL[(3,4)]

println("Base parameters:")
println("  s12 = $(P₀[1])")
println("  s23 = $(P₀[2])")
println("Derived channel invariants at base:")
println("  s123 = $S123_BASE")
println("  s234 = $S234_BASE")
println()

# ------------------------------------------------------------------------------
# Solve for 24 physical solutions at base
# ------------------------------------------------------------------------------

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

println("Phase 1: solving base system...")
res0 = solve(F_fixed)
physical = filter(is_physical, solutions(res0))
SOLS = sort(physical, by = s -> (real(s[1]), imag(s[1]), real(s[2])))
println("  Physical solutions: $(length(SOLS))")
@assert length(SOLS) == N "Expected 24 solutions"

# Conjugation τ
println("Phase 2: extracting τ (conjugation)...")
function extract_τ(sols)
    perm = zeros(Int, length(sols))
    for i in eachindex(sols)
        conj_sol = conj.(sols[i])
        for j in eachindex(sols)
            if norm(sols[j] .- conj_sol) < TOL_CONJ
                perm[i] = j
                break
            end
        end
        perm[i] == 0 && error("No conjugate match for solution $i")
    end
    return perm
end

τ = extract_τ(SOLS)
println("  τ order: $(perm_order(τ))")
println("  τ cycle type: $(cycle_type(τ))")
println()

# ------------------------------------------------------------------------------
# Loop parameterizations (start and end at base)
# ------------------------------------------------------------------------------

"""
Loop around s123 = 0, while starting/ending at the base point.

We vary s12(t) so that:
  s123(t) = s123_base*( (1-α) + α*exp(2π i t) )
For α=1 this is exactly s123(t)=s123_base*exp(2π i t), a circle around 0.
"""
function params_s123(t::Float64; α::Float64=1.0)
    shift = α * ComplexF64(S123_BASE) * (exp(2π * im * t) - 1)
    s12 = P₀[1] + shift
    return ComplexF64[s12, P₀[2]]
end

"""
Loop around s234 = 0, while starting/ending at base.

We vary s23(t) so that:
  s234(t) = s234_base*( (1-α) + α*exp(2π i t) )
"""
function params_s234(t::Float64; α::Float64=1.0)
    shift = α * ComplexF64(S234_BASE) * (exp(2π * im * t) - 1)
    s23 = P₀[2] + shift
    return ComplexF64[P₀[1], s23]
end

"""
Twisted (simultaneous) loop varying both channels with a rational frequency ratio.
For freq=2, exp(2π i * freq * t) returns to 1 at t=1.
"""
function params_twisted(t::Float64; α1::Float64=1.0, α2::Float64=1.0, freq::Int=2)
    shift1 = α1 * ComplexF64(S123_BASE) * (exp(2π * im * t) - 1)
    shift2 = α2 * ComplexF64(S234_BASE) * (exp(2π * im * freq * t) - 1)
    return ComplexF64[P₀[1] + shift1, P₀[2] + shift2]
end

function track_param_loop(loop_func; n_steps::Int=120)
    current_sols = SOLS
    current_params = P₀
    for step in 1:n_steps
        t = step / n_steps
        next_params = loop_func(t)
        res = solve(G, current_sols; start_parameters=current_params, target_parameters=next_params)
        sols = solutions(res)
        length(sols) < N && return nothing
        current_sols = sols
        current_params = next_params
    end
    # Since loop_func(1)=P₀, we're already back at base parameters.
    σ = match_sols(SOLS, current_sols; tol=TOL_MATCH)
    all(σ .> 0) || return nothing
    return σ
end

function compose_perms(σ2::Vector{Int}, σ1::Vector{Int})
    # (σ2 ∘ σ1)(i) = σ2(σ1(i))
    return [σ2[σ1[i]] for i in eachindex(σ1)]
end

function run_loop(name, loop_func; n_steps::Int=120)
    σ = track_param_loop(loop_func; n_steps=n_steps)
    if σ === nothing
        println("  $name: FAILED (lost solutions or matching)")
        return nothing
    end
    ord = perm_order(σ)
    ct = cycle_type(σ)
    dihed = is_dihedral(τ, σ)
    @printf("  %s: order=%d, dihedral=%s, cycle_type=%s\n", name, ord, dihed, ct)
    return σ
end

println("Phase 3: channel-loop experiments")
println()

n_steps = QUICK ? 80 : 160

println("Strategy A (simple closed loops through basepoint):")
σ123 = run_loop("Loop s123 (α=1)", t -> params_s123(t; α=1.0); n_steps=n_steps)
σ234 = run_loop("Loop s234 (α=1)", t -> params_s234(t; α=1.0); n_steps=n_steps)
σtw = run_loop("Twisted (freq=2)", t -> params_twisted(t; α1=1.0, α2=1.0, freq=2); n_steps=n_steps)

println()
println("Strategy B (basepoint → circle around channel → basepoint):")

"""
Encircle s123=0 by taking s123(t)=ε*exp(2πit) and transporting from basepoint.
This matches the loop construction in CURSOR_FIX_MONODROMY.md.
"""
function encircle_s123(; ε::Float64=500.0, n_steps::Int=120)
    # parameter path for s12 such that s123(t)=ε e^{2πit}
    params_at(t) = ComplexF64[P₀[1] + (ε * exp(2π * im * t) - ComplexF64(S123_BASE)), P₀[2]]
    p_start = params_at(0.0)
    current_sols = SOLS
    current_params = P₀
    # transport to circle start
    res = solve(G, current_sols; start_parameters=current_params, target_parameters=p_start)
    sols = solutions(res)
    length(sols) < N && return nothing
    current_sols = sols
    current_params = p_start
    # go around circle
    for step in 1:n_steps
        t = step / n_steps
        p_next = params_at(t)
        res = solve(G, current_sols; start_parameters=current_params, target_parameters=p_next)
        sols = solutions(res)
        length(sols) < N && return nothing
        current_sols = sols
        current_params = p_next
    end
    # back to basepoint
    res_back = solve(G, current_sols; start_parameters=current_params, target_parameters=P₀)
    sols_back = solutions(res_back)
    length(sols_back) < N && return nothing
    σ = match_sols(SOLS, sols_back; tol=TOL_MATCH)
    all(σ .> 0) || return nothing
    return σ
end

"""Encircle s234=0 by taking s234(t)=ε*exp(2πit) and transporting from basepoint."""
function encircle_s234(; ε::Float64=500.0, n_steps::Int=120)
    params_at(t) = ComplexF64[P₀[1], P₀[2] + (ε * exp(2π * im * t) - ComplexF64(S234_BASE))]
    p_start = params_at(0.0)
    current_sols = SOLS
    current_params = P₀
    res = solve(G, current_sols; start_parameters=current_params, target_parameters=p_start)
    sols = solutions(res)
    length(sols) < N && return nothing
    current_sols = sols
    current_params = p_start
    for step in 1:n_steps
        t = step / n_steps
        p_next = params_at(t)
        res = solve(G, current_sols; start_parameters=current_params, target_parameters=p_next)
        sols = solutions(res)
        length(sols) < N && return nothing
        current_sols = sols
        current_params = p_next
    end
    res_back = solve(G, current_sols; start_parameters=current_params, target_parameters=P₀)
    sols_back = solutions(res_back)
    length(sols_back) < N && return nothing
    σ = match_sols(SOLS, sols_back; tol=TOL_MATCH)
    all(σ .> 0) || return nothing
    return σ
end

for ε in (QUICK ? [200.0, 500.0, 1000.0] : [100.0, 200.0, 500.0, 1000.0, 2000.0])
    σ = encircle_s123(; ε=ε, n_steps=(QUICK ? 80 : 140))
    if σ === nothing
        @printf("  Encircle s123 ε=%.0f: FAILED\n", ε)
    else
        @printf("  Encircle s123 ε=%.0f: order=%d dihedral=%s cycle_type=%s\n",
            ε, perm_order(σ), is_dihedral(τ, σ), cycle_type(σ))
    end
end

for ε in (QUICK ? [200.0, 500.0, 1000.0] : [100.0, 200.0, 500.0, 1000.0, 2000.0])
    σ = encircle_s234(; ε=ε, n_steps=(QUICK ? 80 : 140))
    if σ === nothing
        @printf("  Encircle s234 ε=%.0f: FAILED\n", ε)
    else
        @printf("  Encircle s234 ε=%.0f: order=%d dihedral=%s cycle_type=%s\n",
            ε, perm_order(σ), is_dihedral(τ, σ), cycle_type(σ))
    end
end

println()
println("═"^80)
println("GAP SNIPPET (best-effort)")
println("═"^80)
println("tau := $(gap_format(τ));")
if σtw !== nothing
    println("sigma := $(gap_format(σtw));")
    println("G := Group(tau, sigma);")
    println("Size(G);")
    println("TransitiveIdentification(G);")
    println("tau * sigma * tau = sigma^(-1);")
end

println()
println("NOTE:")
println("  If only orders ≤ 6 are ever observed, that is evidence AGAINST D₂₄-as-monodromy.")
println("  If an order 8/12/24 element appears, that supports D₂₄-like dihedral structure.")
