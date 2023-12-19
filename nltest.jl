using NLsolve

println("hello bestie")

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

# The model is a system of nonlinear equations in the variables:
# c, i, k, a, λ
# The equations are:
# 1. Euler equation
#C_t^(-γ) - Λ_t - β * E_t[C_{t+1}^(-γ) * (1 - δ + α * A_{t+1} * K_{t+1}^{α - 1}) - (1 - δ) Λ_t] = 0
# 2. Lagrange definition (contraint not binding)
# Λ_t = 0
# 3. Resource constraint
# C_t + I_t - A_t * K_t^α = 0
# 4. Capital evolution
# K_{t+1} - (1 - δ) K_t - I_t = 0
# 5. Technology evolution
# A_{t+1} - ρ A_t = 0

function f!(F, x)
    F[1] = x[1]^(-γ) - x[5] - β * (x[1]^(-γ) * (1 - δ + α * x[4] * x[3]^(α - 1)) - (1 - δ) * x[5])
    F[2] = x[5]
    F[3] = x[1] + x[2] - x[4] * x[3]^α
    F[4] = x[3] - (1 - δ) * x[3] - x[2]
    F[5] = log(x[4]) - ρ * log(x[4])
end

# x = nlsolve(f!, [ 0.1; 1.2], autodiff = :forward)
x_ss = nlsolve(f!, [0.1, 0.1, 0.1, 0.1, 0.1])
println(x_ss.zero)

println("goodbye bestie")