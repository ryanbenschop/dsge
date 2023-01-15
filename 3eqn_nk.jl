using LinearAlgebra
using DSGE
using Plots

##################################################################################################################################
# 1. Declare parameters
##################################################################################################################################

# Coefficient of relative risk aversion
sigma = 1

# Discount factor
beta = 0.99

# Calvo parameter
theta = 0.75

# Elasticity of substitution between goods/inputs
varepsilon = 9

# Taylor rule inflation response coefficient
phi_pi = 1.5

# Taylor rule output gap response coefficient
phi_y = 0.125

# Inverse Frisch elasticity of substitution
varphi = 5

# Labour share of output
alpha = 0.25

# TFP shock persistence
rho_a = 0.5

# Monetary policy shock persistence
rho_v = 0.5

# Auxiliary parameters
bigtheta = (1 - alpha) / (1 - alpha * (1 - varepsilon))
kappa = (((1 - alpha) * sigma + alpha + varphi) / (1 - alpha)) * (((1 - theta) * (1 - theta * beta) * bigtheta) / (theta))
psi_yn = ((1 + varphi) / (sigma * (1 - alpha) + alpha + varphi))

##################################################################################################################################
# 2. Set up matrices Γ0. Γ1, C, Ψ, Π corresponding with the canonical Sims form:
# Γ0 s_t = Γ1 s_{t-1} + C + Ψ ε_t + Π η_t
##################################################################################################################################

Γ0 = [
#   y_g     y_g_p1      y_n     y_n_p1      π           π_p1        i           r       r_n         v       a
    1       -1          0       0           0           -1/sigma    1/sigma     0       -1/sigma    0       0;          # Equation 1: Dynamic IS
    -kappa  0           0       0           1           -beta       0           0       0           0       0;          # Equation 2: NKPC
    -phi_y  0           0       0           -phi_pi     0           1           0       0           -1      0;          # Equation 3: Monetary policy rule
    0       0           sigma   -sigma      0           0           0           0       1           0       0;          # Equation 4: Wicksellian real rate definition
    0       0           1       0           0           0           0           0       0           0       -psi_yn;    # Equation 5: Natural output production
    0       0           0       0           0           -1          1           -1      0           0       0;          # Equation 6: Fisher relation
    0       0           0       0           0           0           0           0       0           0       1;          # Equation 7: TFP shock process
    0       0           0       0           0           0           0           0       0           1       0;          # Equation 8: Monetary policy shock process
    1       0           0       0           0           0           0           0       0           0       0;          # Equation 9: Output gap expectational error
    0       0           0       0           1           0           0           0       0           0       0;          # Equation 10: Inflation expectational error
    0       0           1       0           0           0           0           0       0           0       0           # Equation 11: Natural output expectational error
]

Γ1 = [
#   y_g     y_g_p1      y_n     y_n_p1      π           π_p1        i           r       r_n         v       a
    0.0     0           0       0           0           0           0           0       0           0       0;      # Equation 1: Dynamic IS
    0       0           0       0           0           0           0           0       0           0       0;      # Equation 2: NKPC
    0       0           0       0           0           0           0           0       0           0       0;      # Equation 3: Monetary policy rule
    0       0           0       0           0           0           0           0       0           0       0;      # Equation 4: Wicksellian real rate definition
    0       0           0       0           0           0           0           0       0           0       0;      # Equation 5: Natural output production
    0       0           0       0           0           0           0           0       0           0       0;      # Equation 6: Fisher relation
    0       0           0       0           0           0           0           0       0           0       rho_a;  # Equation 7: TFP shock process
    0       0           0       0           0           0           0           0       0           rho_v   0;      # Equation 8: Monetary policy shock process
    0       1           0       0           0           0           0           0       0           0       0;      # Equation 9: Output gap expectational error
    0       0           0       0           0           1           0           0       0           0       0;      # Equation 10: Inflation expectational error
    0       0           0       1           0           0           0           0       0           0       0       # Equation 11: Natural output expectational error
]

C = [
    0.0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0
]

Ψ = [
#   eps_a   eps_v
    0.0     0;      # Equation 1
    0       0;      # Equation 2
    0       0;      # Equation 3
    0       0;      # Equation 4
    0       0;      # Equation 5
    0       0;      # Equation 6
    1       0;      # Equation 7
    0       1;      # Equation 8
    0       0;      # Equation 9
    0       0;      # Equation 10
    0       0      # Equation 11
]

Π = [
#   y_g     π       y_n
    0.0     0       0;  # Equation 1
    0       0       0;  # Equation 2
    0       0       0;  # Equation 3
    0       0       0;  # Equation 4
    0       0       0;  # Equation 5
    0       0       0;  # Equation 6
    0       0       0;  # Equation 7
    0       0       0;  # Equation 8
    1       0       0;  # Equation 9
    0       1       0;  # Equation 10
    0       0       1  # Equation 11
]

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
eps_a = 0

# One time monetary policy shock
eps_v = 0.25

shocks = [
    eps_a;
    eps_v
]

##################################################################################################################################
# 5. Generate impulse responses to one time shocks
##################################################################################################################################

# Number of periods to compute
periods = 15

# Calculate initial responses to the shocks
# s_0 = 0 by construction, so s1 = R ε_1
s = [R * shocks]

# Calculate future values of the endogenous variables
# ε_κ = 0 for k > 1, so s_t = T s_{t-1} for t > 1
for t in 2:periods
    s_next = T * s[t - 1]
    push!(s, s_next)
end

##################################################################################################################################
# Plot IRFs
##################################################################################################################################

x = []
π = []
i = []
v = []

for t in 1:periods
    push!(x, s[t][1])
    push!(i, s[t][7] * 4) # Multiplied by 4 to annualise
    push!(π, s[t][5] * 4) # Multiplied by 4 to annualise
    push!(v, s[t][10])
end

function plot_irf(series, series_name, legend_position)
    plot(1:periods, series, label = "", title = series_name, grid = false, markershape = :circle, markersize = 2.5, markerstrokewidth = 0, legend = legend_position)
end

x_irf = plot_irf(x, "Output Gap", :bottomright)
π_irf = plot_irf(π, "Inflation", :topright)
i_irf = plot_irf(i, "Nominal Rate", :topright)
v_irf = plot_irf(v, "Monetary Policy Shock", :topright)

irfs = plot(x_irf, π_irf, i_irf, v_irf, layout = @layout([a b; c d]))
display(irfs)
