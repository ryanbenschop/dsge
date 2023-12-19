using LinearAlgebra
using DSGE
using Plots

println("Hello bestiepie")

β::Float64 = 0.96;
γ::Float64 = 2;
α::Float64 = 0.33;
δ::Float64 = 0.1;
φ::Float64 = 0.975;
ρ::Float64 = 0.9;
σ::Float64 = 0.013;
σ::Float64 = 1;
k_ss::Float64 = ((1 / α) * ((1 / β) - 1 + δ)) ^ (1 / (α - 1));
i_ss::Float64 = δ * k_ss;
c_ss::Float64 = k_ss^α - i_ss;

s = ("c", "cp1", "i", "k", "a", "ap1", "λ")

Γ0 = zeros(Float64, 7, 7)
Γ1 = zeros(Float64, 7, 7)
C_sims = zeros(Float64, 7)
Ψ = zeros(Float64, 7, 1)
Π = zeros(Float64, 7, 2)

# Equation 1: Euler equation
Γ0[1, 1] = -γ
Γ0[1, 2] = γ * β * (1 - δ + α * k_ss^(α - 1))
Γ0[1, 4] = -α * (α - 1) * β * k_ss^(α - 1)
Γ0[1, 6] = -α * β * k_ss^(α - 1)

# Equation 2: Lagrange definition
Γ0[2, 7] = 1

# Equation 3: Resource constraint
Γ0[3, 1] = c_ss
Γ0[3, 3] = i_ss
Γ0[3, 5] = -k_ss^α

Γ1[3, 4] = α * k_ss^α

# Equation 4: Capital evolution
Γ0[4, 3] = -i_ss
Γ0[4, 4] = k_ss

Γ1[4, 4] = (1 - δ) * k_ss

# Equation 4: Technology evolution
Γ0[5, 5] = 1

Γ1[5, 5] = ρ

Ψ[5, 1] = σ

# Equation 5: Consumption expectation error definition
Γ0[6, 1] = 1

Γ1[6, 2] = 1

Π[6, 1] = 1

# Equation 6: Technology expectation error definition
Γ0[7, 5] = 1

Γ1[7, 6] = 1

Π[7, 2] = 1

sol = DSGE.gensys(Γ0, Γ1, C_sims, Ψ, Π)

P_sims = sol[1]
C_sims = sol[2]
Q_sims = sol[3]

################################################################################################
# Sims stuff done here. Now the OccBin stuff
################################################################################################

x = ("c", "i", "k", "a", "λ")

A = [
#   c                                       i           k           a                           λ
    γ * β * (1 - δ + α * k_ss^(α - 1))      0           0           -α * β * k_ss^(α - 1)       0;
    0                                       0           0           0                           0;
    0                                       0           0           0                           0;
    0                                       0           0           0                           0;
    0                                       0           0           0                           0;
]

B = [
#   c           i           k                                   a           λ
    -γ          0           -α * (α - 1) * β * k_ss^(α - 1)     0           c_ss^γ * (β * (1 - δ) - 1);
    0           0           0                                   0           1;
    c_ss        i_ss        0                                   -k_ss^α     0;
    0           -i_ss       k_ss                                0           0;
    0           0           0                                   1           0;
]

C = [
#   c           i           k               a           λ
    0           0           0               0           0;
    0           0           0               0           0;
    0           0           -α * k_ss^α     0           0;
    0           0           -(1 - δ) * k_ss 0           0;
    0           0           0               -ρ          0;
]

E = [
#   ε
    0;
    0;
    0;
    0;
    -σ;
]

A_star = copy(A)

B_star = copy(B)
B_star[2, 5] = 0
B_star[2, 2] = 1

C_star = copy(C)

# D_star = zeros(Float64, 5, 1)
# D_star[2, 1] = -(log(φ * i_ss) - log(i_ss))
D_star = [
    0;
    # 0;
    -(log(φ * i_ss) - log(i_ss));
    0;
    0;
    0;
]

println("Lower bound on Investment:")
display(D_star[2, 1])

E_star = copy(E)

P = [
    P_sims[1, 1]    P_sims[1, 3]    P_sims[1, 4]    P_sims[1, 5]    P_sims[1, 7];
    P_sims[3, 1]    P_sims[3, 3]    P_sims[3, 4]    P_sims[3, 5]    P_sims[3, 7];
    P_sims[4, 1]    P_sims[4, 3]    P_sims[4, 4]    P_sims[4, 5]    P_sims[4, 7];
    P_sims[5, 1]    P_sims[5, 3]    P_sims[5, 4]    P_sims[5, 5]    P_sims[5, 7];
    P_sims[6, 1]    P_sims[6, 3]    P_sims[6, 4]    P_sims[6, 5]    P_sims[6, 7];
]

R = [
    C_sims[1, 1];
    C_sims[2, 1];
    C_sims[3, 1];
    C_sims[4, 1];
    C_sims[5, 1];
]

T = 14  # Number of periods constraint is binding

P_mats = [P]
R_mats = [R]

for j in 1:T
    P_mat = -pinv(A_star * P_mats[j] + B_star) * C_star

    R_mat = -pinv(A_star * P_mats[j] + B_star) * (A_star * R_mats[j] + D_star)
    push!(P_mats, P_mat)
    push!(R_mats, R_mat)
end

Q1 = -pinv(A_star * P_mats[T] + B_star) * E_star
R1 = -pinv(A_star * P_mats[T] + B_star) * (A_star * R_mats[T] + D_star)

x = []

eps = -3

x1 = Q1 * eps + R1

println("Q1:")
display(Q1)

push!(x, x1)

for j in 1:T-1
    x_next = P_mats[T-j] * x[j] + R_mats[T-j]
    push!(x, x_next)
end

################################################################################################
# OccBin stuff done. Simulate the later periods.
################################################################################################

sim_periods = 50 - T + 1

for j in 1:sim_periods
    x_next = P * x[T+j-1] + R_mats[1]
    push!(x, x_next)
end

c = [x[t][1] for t in 1:T+sim_periods]
i = [x[t][2] for t in 1:T+sim_periods]
k = [x[t][3] for t in 1:T+sim_periods]
a = [x[t][4] for t in 1:T+sim_periods]

c_plt = plot(c, title="Consumption")
i_plt = plot(i, title="Investment")
k_plt = plot(k, title="Capital")
a_plt = plot(a, title="Technology")

plt = plot(c_plt, i_plt, k_plt, a_plt, layout=(2, 2), legend=false)
display(plt)

println("Goodbye bestiepie")