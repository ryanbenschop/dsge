using LinearAlgebra
using Statistics
# using DSGE
# using Plots

println("Hello bestiepie")

# The following function takes in matrices A, B, C, D, and E, and returns the matrices Γ0, Γ1, C_sims, Ψ, and Π, where:
# Γ0 = [A B; I 0] has twice as many rows and columns as A.
# Γ1 = -[0 C; I 0] has twice as many rows and columns as A.
# C_sims = -[D; 0] has twice as many rows as A.
# Ψ = -[E; 0] has twice as many rows as A.
# Π = [0 I] has twice as many rows and columns as A.
function get_Sims_mats(A::Array{Float64}, B::Array{Float64}, C::Array{Float64}, D::Array{Float64}, E::Array{Float64})
    num_vars = size(A, 1)
    num_shocks = size(E, 2)

    Γ0 = zeros(Float64, num_vars * 2, num_vars * 2)
    Γ1 = zeros(Float64, num_vars * 2, num_vars * 2)
    C_sims = zeros(Float64, num_vars * 2)
    Ψ = zeros(Float64, num_vars * 2, num_shocks)
    Π = zeros(Float64, num_vars * 2, num_vars)

    Γ0[1:num_vars, 1:num_vars] = A
    Γ0[1:num_vars, num_vars + 1:num_vars * 2] = B
    Γ0[num_vars + 1:num_vars * 2, num_vars + 1:num_vars * 2] = Matrix{Float64}(I, num_vars, num_vars)

    Γ1[1:num_vars, num_vars + 1:num_vars * 2] = -C
    Γ1[num_vars + 1:num_vars * 2, 1:num_vars] = Matrix{Float64}(I, num_vars, num_vars)

    C_sims[1:num_vars] = -D

    Ψ[1:num_vars, 1:num_shocks] = -E

    Π[num_vars + 1:num_vars * 2, 1:num_vars] = Matrix{Float64}(I, num_vars, num_vars)

    return Γ0, Γ1, C_sims, Ψ, Π
end

#--------------------------------------------------------------------------------------------------------------
# Test
#--------------------------------------------------------------------------------------------------------------

β::Float64 = 0.96;
γ::Float64 = 2;
α::Float64 = 0.33;
δ::Float64 = 0.1;
φ::Float64 = 0.975;
ρ::Float64 = 0.9;
σ::Float64 = 0.013;
k_ss::Float64 = ((1 / α) * ((1 / β) - 1 + δ)) ^ (1 / (α - 1));
i_ss::Float64 = δ * k_ss;
c_ss::Float64 = k_ss^α - i_ss;

A = zeros(Float64, 3, 3)
B = zeros(Float64, 3, 3)
C = zeros(Float64, 3, 3)
D = zeros(Float64, 3, 1)
E = zeros(Float64, 3, 1)

# Equation 1: Euler equation
A[1, 1] = β * (α * k_ss^(α - 1) + 1 - δ) * γ
A[1, 3] = -β * α * k_ss^(α - 1)

B[1, 1] = -γ
B[1, 2] = -β * α * (α - 1) * k_ss^(α - 1)

# Equation 2: Capital LOM
B[2, 1] = c_ss
B[2, 2] = k_ss
B[2, 3] = -k_ss^α

C[2, 2] = -(α * k_ss^α  + (1 - δ) * k_ss)

# Equation 3: Technology LOM
B[3, 3] = 1

C[3, 3] = -ρ

E[3, 1] = -σ

Γ0, Γ1, C_sims, Ψ, Π = get_Sims_mats(A, B, C, D, E)

# sol = DSGE.gensys(Γ0, Γ1, C_sims, Ψ, Π)

println("Goodbye bestiepie")