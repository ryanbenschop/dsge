using LinearAlgebra
using DSGE
using Plots

##################################################################################################################################
# This file computes the solution of a baseline new Keynesian DSGE model, as described in Galí (2015), Chapter 3.
##################################################################################################################################

##################################################################################################################################
# 1. Declare parameters
##################################################################################################################################

# Coefficient of relative risk aversion
σ = 1

# Discount factor
β = 0.99

# Inverse Frisch elasticity of substitution
φ = 5

# Elasticity of substitution between goods/inputs 
ε = 9

# Calvo parameter
θ = 0.75

# Taylor rule inflation response coefficient
φ_π = 1.5

# Taylor rule output gap (/output deviation from steady state) response coefficient
φ_y = 0.125

# Labour share of output
α = 0.25

# Steady state markup
μ = ε / (ε - 1)

# Productivity shock persistence
ρ_a = 0.9

# Monetary policy shock persistence
ρ_v = 0.5

# Preference shock persistence
ρ_z = 0.5

# Steady state variables used in linearised equations
R_ss = 1 / β
mc_ss = 1 / μ
l_ss = ((1 - α) * mc_ss)^(1 / ((1 - σ) * α + φ + σ))
y_ss = l_ss^(1-α)
λ_ss = y_ss^(-σ)
xn_ss = (λ_ss * mc_ss * y_ss) / (1 - β * θ)
xd_ss = (λ_ss * y_ss) / (1 - β * θ)


##################################################################################################################################
# 2. Set up matrices Γ0. Γ1, C, Ψ, Π corresponding with the canonical Sims form:
# Γ0 s_t = Γ1 s_{t-1} + C + Ψ ε_t + Π η_t
##################################################################################################################################

s = ("λ", "λp1", "c", "w", "πstar", "xn", "xnp1", "xd", "xdp1", "π", "πp1", "R", "q", "r", "l", "y", "yn", "yg", "mc", "p", "a", "v", "z")

Γ0 = zeros(Float64, 23, 23)
Γ1 = zeros(Float64, 23, 23)
C = zeros(Float64, 23)
Ψ = zeros(Float64, 23, 3)
Π = zeros(Float64, 23, 4)

# Equation 1: Lagrange multiplier definition
Γ0_λ_1 = 1
Γ0_z_1 = -1
Γ0_c_1 = σ

Γ0[1, 1] = Γ0_λ_1
Γ0[1, 23] = Γ0_z_1
Γ0[1, 3] = Γ0_c_1

# Equation 2: Euler equation
Γ0_q_2 = 1
Γ0_λp1_2 = -1
Γ0_λ_2 = 1
Γ0_πp1_2 = 1

Γ0[2, 13] = Γ0_q_2
Γ0[2, 2] = Γ0_λp1_2
Γ0[2, 1] = Γ0_λ_2
Γ0[2, 11] = Γ0_πp1_2

# Equation 3: Labour supply
Γ0_w_3 = 1
Γ0_z_3 = -1
Γ0_l_3 = -φ
Γ0_λ_3 = 1

Γ0[3, 4] = Γ0_w_3
Γ0[3, 23] = Γ0_z_3
Γ0[3, 15] = Γ0_l_3
Γ0[3, 1] = Γ0_λ_3

# Equation 4: Price setting
Γ0_πstar_4 = (1 + ((ε * α) / (1 - α)))
Γ0_xn_4 = -1
Γ0_xd_4 = 1

Γ0[4, 5] = Γ0_πstar_4
Γ0[4, 6] = Γ0_xn_4
Γ0[4, 8] = Γ0_xd_4

# Equation 5: Price setting numerator recursion
Γ0_xn_5 = xn_ss
Γ0_λ_5 = -λ_ss * y_ss * mc_ss
Γ0_y_5 = -λ_ss * y_ss * mc_ss
Γ0_mc_5 = -λ_ss * y_ss * mc_ss
Γ0_πp1_5 = -β * θ * xn_ss * (ε + ((ε * α) / (1 - α)))
Γ0_xnp1_5 = -β * θ * xn_ss

Γ0[5, 6] = Γ0_xn_5
Γ0[5, 1] = Γ0_λ_5
Γ0[5, 16] = Γ0_y_5
Γ0[5, 19] = Γ0_mc_5
Γ0[5, 11] = Γ0_πp1_5
Γ0[5, 7] = Γ0_xnp1_5

# Equation 6: Price setting denominator recursion
Γ0_xd_6 = xd_ss
Γ0_λ_6 = -λ_ss * y_ss
Γ0_y_6 = -λ_ss * y_ss
Γ0_πp1_6 = -β * θ * xd_ss * (ε - 1)
Γ0_xdp1_6 = -β * θ * xd_ss

Γ0[6, 8] = Γ0_xd_6
Γ0[6, 1] = Γ0_λ_6
Γ0[6, 16] = Γ0_y_6
Γ0[6, 11] = Γ0_πp1_6
Γ0[6, 9] = Γ0_xdp1_6

# Equation 7: Price level LOM
Γ0_πstar_7 = 1
Γ0_π_7 = -θ/(1 - θ)

Γ0[7, 5] = Γ0_πstar_7
Γ0[7, 10] = Γ0_π_7

# Equation 8: Monetary policy rule
Γ0_R_8 = 1
Γ0_π_8 = -φ_π
#Γ0_yg_8 = -φ_y         # Standard rule: CB responds to output gap
Γ0_y_8 = -φ_y           # Galí rule: CB responds to output deviations from steady state
Γ0_v_8 = -1

Γ0[8, 12] = Γ0_R_8
Γ0[8, 10] = Γ0_π_8
#Γ0[8, 18] = Γ0_yg_8    # Standard rule
Γ0[8, 16] = Γ0_y_8      # Galí rule
Γ0[8, 22] = Γ0_v_8

# Equation 9: Nominal rate definition
Γ0_R_9 = 1
Γ0_q_9 = 1

Γ0[9, 12] = Γ0_R_9
Γ0[9, 13] = Γ0_q_9

# Equation 10: Real rate Fisher relation
Γ0_r_10 = 1
Γ0_R_10 = -1
Γ0_πp1_10 = 1

Γ0[10, 14] = Γ0_r_10
Γ0[10, 12] = Γ0_R_10
Γ0[10, 11] = Γ0_πp1_10

# Equation 11: Aggregate output
Γ0_y_11 = 1
Γ0_a_11 = -1
Γ0_l_11 = α - 1

Γ0[11, 16] = Γ0_y_11
Γ0[11, 21] = Γ0_a_11
Γ0[11, 15] = Γ0_l_11

# Equation 12: Market clearing
Γ0_c_12 = 1
Γ0_y_12 = -1

Γ0[12, 3] = Γ0_c_12
Γ0[12, 16] = Γ0_y_12

# Equation 13: Natural output
Γ0_yn_13 = 1
Γ0_a_13 = -(1 + φ) / (σ * (1 - α) + α + φ)

Γ0[13, 17] = Γ0_yn_13
Γ0[13, 21] = Γ0_a_13

# Equation 14: Output gap definition
Γ0_yg_14 = 1
Γ0_y_14 = -1
Γ0_yn_14 = 1

Γ0[14, 18] = Γ0_yg_14
Γ0[14, 16] = Γ0_y_14
Γ0[14, 17] = Γ0_yn_14

# Equation 15: Marginal cost
Γ0_mc_15 = 1
Γ0_w_15 = -1
Γ0_y_15 = -α / (1 - α)
Γ0_a_15 = 1 / (1 - α)

Γ0[15, 19] = Γ0_mc_15
Γ0[15, 4] = Γ0_w_15
Γ0[15, 16] = Γ0_y_15
Γ0[15, 21] = Γ0_a_15

# Equation 16: Inflation/price level relation
Γ0_π_16 = 1
Γ0_p_16 = -1

Γ0[16, 10] = Γ0_π_16
Γ0[16, 20] = Γ0_p_16

Γ1_p_16 = -1

Γ1[16, 20] = Γ1_p_16

# Equation 17: Technology shock process
Γ0_a_17 = 1

Γ0[17, 21] = Γ0_a_17

Γ1_a_17 = ρ_a

Γ1[17, 21] = Γ1_a_17

Ψ_eps_a_17 = 1

Ψ[17, 1] = Ψ_eps_a_17

# Equation 18: Monetary policy shock process
Γ0_v_18 = 1

Γ0[18, 22] = Γ0_v_18

Γ1_v_18 =ρ_v

Γ1[18, 22] = Γ1_v_18

Ψ_eps_v_18 = 1

Ψ[18, 2] = Ψ_eps_v_18

# Equation 19: Preference shock process
Γ0_z_19 = 1

Γ0[19, 23] = Γ0_z_19

Γ1_z_19 = ρ_z

Γ1[19, 23] = Γ1_z_19

Ψ_eps_z_19 = 1

Ψ[19, 3] = Ψ_eps_z_19

# Equation 20: λ expectational error
Γ0_λ_20 = 1

Γ0[20, 1] = Γ0_λ_20

Γ1_λp1_20 = 1

Γ1[20, 2] = Γ1_λp1_20

Π_λ_20 = 1
Π[20, 1] = Π_λ_20

# Equation 21: π expectational error
Γ0_π_21 = 1

Γ0[21, 10] = Γ0_π_21

Γ1_πp1_21 = 1

Γ1[21, 11] = Γ1_πp1_21

Π_π_21 = 1
Π[21, 2] = Π_π_21

# Equation 22: xn expectational error
Γ0_xn_22 = 1

Γ0[22, 6] = Γ0_xn_22

Γ1_xnp1_22 = 1

Γ1[22, 7] = Γ1_xnp1_22

Π_xn_22 = 1
Π[22, 3] = Π_xn_22

# Equation 23: xd expectational error
Γ0_xd_23 = 1

Γ0[23, 8] = Γ0_xd_23

Γ1_xdp1_23 = 1

Γ1[23, 9] = Γ1_xdp1_23

Π_xd_23 = 1
Π[23, 4] = Π_xd_23

##################################################################################################################################
# 3. Solve the system for the matrices T, R in the transition equation
# s_t = T s_{t-1} + R ε_t + C
##################################################################################################################################

sol = DSGE.gensys(Γ0, Γ1, C, Ψ, Π)
T = sol[1]
R = sol[3]

##################################################################################################################################
# 4. Set shocks
##################################################################################################################################
# One time technology shock
eps_a = 1

# One time monetary policy shock
eps_v = 0

# One time preference shock
eps_z = 0

shocks = [
    eps_a;
    eps_v;
    eps_z
]

##################################################################################################################################
# 5. Generate impulse responses to one time shocks
##################################################################################################################################

# Number of periods to compute
periods = 15

# Calculate initial responses to the shocks
# s_{-1} = 0 by construction, so s_0 = R ε_1
s = [R * shocks]

# Calculate future values of the endogenous variables
# ε_κ = 0 for k > 0, so s_t = T s_{t-1} for t > 0
for k in 2:15
    s_next = T * s[k - 1]
    push!(s, s_next)
end

##################################################################################################################################
# 6. Plot IRFs
##################################################################################################################################

y = []
yg = []
π = []
l = []
w = []
p = []
R = []
mc = []
r = []
a = []
v = []
z = []

for k in 1:periods
    push!(y, s[k][16])
    push!(yg, s[k][18])
    push!(l, s[k][15])
    push!(w, s[k][4])
    push!(p, s[k][20])
    push!(mc, s[k][19])
    push!(π, s[k][10] * 4)
    push!(R, s[k][12] * 4)
    push!(r, s[k][14] * 4)
    push!(a, s[k][21])
    push!(v, s[k][22])
    push!(z, s[k][23])
end

function plot_irf(series, series_name, legend_position)
    plot(1:periods, series, label = "", title = series_name, grid = false, markershape = :circle, markersize = 2.5, markerstrokewidth = 0, legend = legend_position)
end

yg_irf = plot_irf(yg, "Output gap", :topright)
π_irf = plot_irf(π, "Inflation", :bottomright)
y_irf = plot_irf(y, "Output", :bottomright)
l_irf = plot_irf(l, "Employment", :bottomright)
w_irf = plot_irf(w, "Real wage", :bottomright)
p_irf = plot_irf(p, "Price level", :bottomright)
mc_irf = plot_irf(mc, "Marginal cost", :bottomright)
a_irf = plot_irf(a, "Productivity", :topright)
#v_irf = plot_irf(v, "Monetary policy shock", :topright)
#z_irf = plot_irf(z, "Preference shock", :bottomright)

plt = plot(yg_irf, π_irf, y_irf, l_irf, w_irf, p_irf, mc_irf, a_irf, layout = @layout([a b; c d; e f; g h]), size = (500, 800))
display(plt)