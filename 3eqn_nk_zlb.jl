using LinearAlgebra
using DSGE
using Plots

##################################################################################################################################
# This file computes the solution of a baseline new Keynesian DSGE model, as described in Galí (2015), Chapter 3. It incorporates 
# forward guidance via a commitment to the ZLB by the monetary authority, as implemented in Chen, Cúrdia, and Ferrero (2012).
##################################################################################################################################

##################################################################################################################################
# 1. Declare parameters
##################################################################################################################################

# Coefficient of relative risk aversion
σ = 1

# Discount factor
β = 0.99

# Calvo parameter
θ = 0.75

# Elasticity of substitution between goods/inputs
ε = 9

# Taylor rule inflation response coefficient
φ_π = 1.5

# Taylor rule output gap response coefficient
φ_x = 0.125

# Inverse Frisch elasticity of substitution
φ = 5

# Labour share of output
α = 0.25

# TFP shock persistence
ρ_a = 0.5

# Monetary policy shock persistence
ρ_v = 0.5

# Cost push shock persistence
ρ_u = 0.5

# Demand shock persistence
ρ_z = 0.5

# Auxiliary parameters
Θ = (1 - α) / (1 - α * (1 - ε))
κ = (((1 - α) * σ + α + φ) / (1 - α)) * (((1 - θ) * (1 - θ * β) * Θ) / (θ))
ψ_ya = ((1 + φ) / (σ * (1 - α) + α + φ))

# Steady state parameters
i_ss = 1 / β - 1

##################################################################################################################################
# 2. Set up matrices in canonical Sims form and partitioned form.
# Sims form: Γ_4^{s_t} z_t = \bar Γ_0^{s_t} + Γ_1^{s_t} z_{t-1} + Γ_2^{s_t} ε_t + Γ_3^{s_t} η_t
# Partitioned system:
# Γ_4^{s_t}(j) \mathbb E_t z_{t+1} = Γ_0^{s_t}(j) + Γ_1^{s_t}(j) z_t
# Γ_4^{s_t}(i) z_t = Γ_0^{s_t}(i) + Γ_1^{s_t}(i) z_{t-1} + Γ_2^{s_2}(i) ε_t
##################################################################################################################################

variables = ["x", "x(+1)", "π", "π(+1)", "i", "r", "r_n", "a", "v", "z", "u"]
shocks = ["ε_a", "ε_v", "ε_u", "ε_z"]
num_variables = length(variables)
num_shocks = length(shocks)
num_exp_errors = 2

# Initialise Sims matrices for normal times
Γ4_n = zeros(Float64, num_variables, num_variables)
Γ0_bar_n = zeros(Float64, num_variables)
Γ1_n = zeros(Float64, num_variables, num_variables)
Γ2_n = zeros(Float64, num_variables, num_shocks)
Γ3_n = zeros(Float64, num_variables, num_exp_errors)

# Equation 1: Dynamic IS
# x - x(+1) + 1/σ(i - π(+1) - r_n) = 0
Γ4_n[1, 1] = 1
Γ4_n[1, 2] = -1
Γ4_n[1, 5] = 1 / σ
Γ4_n[1, 4] = -1 / σ
Γ4_n[1, 7] = -1 / σ

# Equation 2: NKPC
# π - κx - βπ(+1) = 0
Γ4_n[2, 3] = 1
Γ4_n[2, 1] = -κ
Γ4_n[2, 4] = -β
Γ4_n[2, 11] = -1

# Equation 3: Policy rule in normal times
# i - φ_π π - φ_x x - v = 0
Γ4_n[3, 5] = 1
Γ4_n[3, 3] = -φ_π
Γ4_n[3, 1] = -φ_x
Γ4_n[3, 9] = -1

# Equation 4: Natural rate process
# r_n - σψ_ya(ρ_a - 1)a - (1 - ρ_z)z = 0
Γ4_n[4, 7] = 1
Γ4_n[4, 8] = -σ * ψ_ya * (ρ_a - 1)
Γ4_n[4, 10] = ρ_z - 1

# Equation 5: Real rate definition
# r - i + π(+1) = 0
Γ4_n[5, 6] = 1
Γ4_n[5, 5] = -1
Γ4_n[5, 4] = 1

# Equation 6: TFP shock process
# a = ρ_a a(-1) + ε_a
Γ4_n[6, 8] = 1
Γ1_n[6, 8] = ρ_a
Γ2_n[6, 1] = 1

# Equation 7: Monetary policy shock process
# v = ρ_v v(-1) + ε_v
Γ4_n[7, 9] = 1
Γ1_n[7, 9] = ρ_v
Γ2_n[7, 2] = 1

# Equation 8: Demand shock process
# z = ρ_z z(-1) + ε_z
Γ4_n[8, 10] = 1
Γ1_n[8, 10] = ρ_z
Γ2_n[8, 3] = 1

# Equation 9: Cost push shock process
# u = ρ_u u(-1) + ε_u
Γ4_n[9, 11] = 1
Γ1_n[9, 11] = ρ_u
Γ2_n[9, 4] = 1

# Equation 10: Output gap expectation error
Γ4_n[10, 1] = 1
Γ1_n[10, 2] = 1
Γ3_n[10, 1] = 1

# Equation 11: Inflation expectation error
Γ4_n[11, 3] = 1
Γ1_n[11, 4] = 1
Γ3_n[11, 2] = 1

# Set up matrices for ZLB state
Γ4_zlb = copy(Γ4_n)
Γ0_zlb = copy(Γ0_bar_n)
Γ1_zlb = copy(Γ1_n)
Γ2_zlb = copy(Γ2_n)
Γ3_zlb = copy(Γ3_n)

# Interest rate rule at ZLB:
# i + i_ss = 0
Γ4_zlb[3, 5] = 1
Γ4_zlb[3, 3] = 0
Γ4_zlb[3, 1] = 0
Γ4_zlb[3, 9] = 0
Γ0_zlb[3] = -i_ss

# Forward-looking block – this only contains the equations defining the expectational errors
Γ4_n_FL = Γ4_n[10:11, :]
Γ0_n_FL = zeros(Float64, 2)
Γ1_n_FL = Γ1_n[10:11, :]
Γ3_n_FL = Γ3_n[10:11, :]

# Backward-looking block – all the other equations
Γ4_n_BL = Γ4_n[1:9, :]
Γ0_n_BL = Γ0_bar_n[1:9]
Γ1_n_BL = Γ1_n[1:9, :]
Γ2_n_BL = Γ2_n[1:9, :]

Γ4_zlb_FL = copy(Γ4_n_FL)
Γ0_zlb_FL = copy(Γ0_n_FL)
Γ1_zlb_FL = copy(Γ1_n_FL)
Γ3_zlb_FL = copy(Γ3_n_FL)

Γ4_zlb_BL = Γ4_zlb[1:9, :]
Γ0_zlb_BL = Γ0_zlb[1:9]
Γ1_zlb_BL = Γ1_zlb[1:9, :]
Γ2_zlb_BL = Γ2_zlb[1:9, :]

#######################################################################################################################
# 3. Solve system for "normal" times. i.e., find the representation:
# z_t = Φ_0^n + Φ_1^n z_{t-1} + Φ_2^n ε_t
#######################################################################################################################

sol = DSGE.gensys(Γ4_n, Γ1_n, Γ0_bar_n, Γ2_n, Γ3_n)
Φ0_n = sol[2]
Φ1_n = sol[1]
Φ2_n = sol[3]

#######################################################################################################################
# Recursively solve backwards for zlb times. i.e., find the representation:
# \tilde Γ_4(t) z_t = \tilde Γ_0(t) + \tilde Γ_1(t) z_{t-1} + \tilde Γ_2(t) ε_t
#######################################################################################################################

# Number of periods the monetary authority commits to the ZLB (call this K to simplify notation)
commitment_length = 4

# Initialise vector of shocks – each element should be a vector of the shocks in each period
shocks = []

for k in 1:commitment_length
    push!(shocks, zeros(num_shocks))
end

# Set a technology shock in the initial period
shocks[1][1] = 0.25

# Matrices for last period of commitment
curlyΓ0_K = [
    Γ0_zlb_FL - Γ4_zlb_FL * (Φ0_n + Φ2_n * shocks[commitment_length]);
    Γ0_zlb_BL
]

curlyΓ1_K = [
    zeros(Float64, 2, 11);
    Γ1_zlb_BL
]

curlyΓ2_K = [
    zeros(Float64, 2, 4);
    Γ2_zlb_BL
]

curlyΓ4_K = [
    Γ4_zlb_FL * Φ1_n - Γ1_zlb_FL;
    Γ4_zlb_BL
]

Φ0_K = (pinv(curlyΓ4_K)) * curlyΓ0_K
Φ1_K = (pinv(curlyΓ4_K)) * curlyΓ1_K
Φ2_K = (pinv(curlyΓ4_K)) * curlyΓ2_K

Φ0 = [Φ0_K]
Φ1 = [Φ1_K]
Φ2 = [Φ2_K]

# Matrices for remaining periods
for k in 1:commitment_length-1
    curlyΓ0_k = [
        Γ0_zlb_FL - Γ4_zlb_FL * (Φ0[k] + Φ2[k] * shocks[commitment_length - k + 1]);
        Γ0_zlb_BL
    ]

    curlyΓ1_k = [
        zeros(Float64, 2, 11);
        Γ1_zlb_BL
    ]

    curlyΓ2_k = [
        zeros(Float64, 2, 4);
        Γ2_zlb_BL
    ]

    curlyΓ4_k = [
        Γ4_zlb_FL * Φ1[k] - Γ1_zlb_FL;
        Γ4_zlb_BL
    ]

    Φ0_k = (pinv(curlyΓ4_k)) * curlyΓ0_k
    Φ1_k = (pinv(curlyΓ4_k)) * curlyΓ1_k
    Φ2_k = (pinv(curlyΓ4_k)) * curlyΓ2_k
    
    push!(Φ0, Φ0_k)
    push!(Φ1, Φ1_k)
    push!(Φ2, Φ2_k)
end

#######################################################################################################################
# 4. Compute IRFs
#######################################################################################################################

s = [Φ0[commitment_length] + Φ2[commitment_length] * shocks[1]]
for k in 1:commitment_length-1
    s_next = Φ0[commitment_length - k] + Φ1[commitment_length - k] * s[k] + Φ2[commitment_length - k] * shocks[k+1]
    push!(s, s_next)
end

# Number of periods to compute
periods = 15
for t in 1:periods-commitment_length
    s_next = Φ0_n + Φ1_n * s[t + commitment_length - 1]
    push!(s, s_next)
end

#######################################################################################################################
# 5. Plot IRFs
#######################################################################################################################

x = []
π = []
i = []
r = []
rn = []
a = []
v = []
z = []
u = []

for k in 1:periods
    push!(x, s[k][1])
    push!(π, s[k][3] * 4) # Multiplied by 4 to annualise
    push!(i, s[k][5] * 4) # Multiplied by 4 to annualise
    push!(r, s[k][6] * 4) # Multiplied by 4 to annualise
    push!(rn, s[k][7] * 4) # Multiplied by 4 to annualise
    push!(a, s[k][8])
    push!(v, s[k][9])
    push!(z, s[k][10])
    push!(u, s[k][11])
end

function plot_irf(series, series_name, legend_position)
    plot(1:periods, series, label = "", title = series_name, grid = false, markershape = :circle, markersize = 2.5, markerstrokewidth = 0, legend = legend_position)
    plot!(1:periods, zeros(periods), linestyle = :dash, color = :grey, label = "")
end

x_irf = plot_irf(x, "Output gap", :topright)
π_irf = plot_irf(π, "Inflation", :bottomright)
i_irf = plot_irf(i, "Nominal rate", :bottomright)
r_irf = plot_irf(r, "Real rate", :bottomright)
rn_irf = plot_irf(rn, "Wicksellian rate", :bottomright)
a_irf = plot_irf(a, "Productivity shock", :bottomright)
v_irf = plot_irf(v, "Monetary policy shock", :bottomright)
z_irf = plot_irf(z, "Preference shock", :bottomright)
u_irf = plot_irf(u, "Cost push shock", :bottomright)

plt = plot(x_irf, π_irf, i_irf, r_irf, rn_irf, a_irf, layout = @layout([a b; c d; e f]))
display(plt)