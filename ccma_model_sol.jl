# using LinearAlgebra
# using ForwardDiff
# using Plots
using DSGE
# include("ccma_model_def.jl")

println("hello")

# ccma_ss_args = [ccma_log_ss..., ccma_log_ss..., ccma_log_ss...]

# J = ForwardDiff.jacobian(x -> ccma_model_00(x), ccma_ss_args)

# A = J[1:42, 1:42]
# B = J[1:42, 43:78]
# C = J[1:42, 79:117]

# function get_Sims_mats(A::Array{Float64}, B::Array{Float64}, C::Array{Float64}, D::Array{Float64}, E::Array{Float64})
#     num_vars = size(A, 1)
#     num_shocks = size(E, 2)

#     Γ0 = zeros(Float64, num_vars * 2, num_vars * 2)
#     Γ1 = zeros(Float64, num_vars * 2, num_vars * 2)
#     C_sims = zeros(Float64, num_vars * 2)
#     Ψ = zeros(Float64, num_vars * 2, num_shocks)
#     Π = zeros(Float64, num_vars * 2, num_vars)

#     Γ0[1:num_vars, 1:num_vars] = A
#     Γ0[1:num_vars, num_vars + 1:num_vars * 2] = B
#     Γ0[num_vars + 1:num_vars * 2, num_vars + 1:num_vars * 2] = Matrix{Float64}(I, num_vars, num_vars)

#     Γ1[1:num_vars, num_vars + 1:num_vars * 2] = -C
#     Γ1[num_vars + 1:num_vars * 2, 1:num_vars] = Matrix{Float64}(I, num_vars, num_vars)

#     C_sims[1:num_vars] = -D

#     Ψ[1:num_vars, 1:num_shocks] = E

#     Π[num_vars + 1:num_vars * 2, 1:num_vars] = Matrix{Float64}(I, num_vars, num_vars)

#     return Γ0, Γ1, C_sims, Ψ, Π
# end

# E = zeros(Float64, 39, 1)
# E[22, 1] = 1
# E[23, 1] = 1
# E[24, 1] = 1
# E[25, 1] = 1

# Γ0, Γ1, C_sims, Ψ, Π = get_Sims_mats(A, B, C, zeros(Float64, 39, 1), E)

# sol = DSGE.gensys(Γ0, Γ1, C_sims, Ψ, Π)

# P = sol[1]
# Q = sol[3]

# eps = 1

# s = [Q * eps]

# periods = 100

# for t in 2:periods
#     s_next = P * s[t - 1]
#     push!(s, s_next)
# end

# c = []
# k = []
# z = []
# a = []

# for t in 1:periods
#     push!(c, s[t][1])
#     push!(k, s[t][6])
#     push!(z, s[t][22])
#     push!(a, s[t][23])
# end

# c_plt = plot(c, label="c")
# k_plt = plot(k, label="k")
# z_plt = plot(z, label="z")
# a_plt = plot(a, label="a")

# # plt = plot(c_plt)
# # plt = plot(k_plt)
# # plt = plot(z_plt)
# plt = plot(a_plt)

# display(plt)

println("goodbye")