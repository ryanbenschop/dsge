using LinearAlgebra
using NLsolve

println("Hello bestie")

#--------------------------------------------------------------------------------------------------------------------
# Parameters
#--------------------------------------------------------------------------------------------------------------------
# Calibrated parameters
m::Float64 = 0.9;                   # Maximum LTV
η::Float64 = 1;                     # Labour disutility
β::Float64 = 0.995;                 # Discount factor for patient agents
π_ss::Float64 = 1.005;              # Steady state gross inflation
α::Float64 = 0.33;                  # Capital share of output
δ_k::Float64 = 0.025;               # Capital depreciation rate
j_bar::Float64 = 0.04;              # Housing weight in utility function
x_p_ss::Float64 = 1.2;              # Steady state price markup
x_w_ss::Float64 = 1.2;              # Steady state wage markup

# Estimated parameters (posterior mode)
βp::Float64 = 0.9922;               # Discount factor for impatient agents
ε_c::Float64 = 0.6842;              # Habit in consumption
ε_h::Float64 = 0.8799;              # Habit in housing
φ::Float64 = 4.1209;                # Investment adjustment cost
σ::Float64 = 0.5013;                # Wage share for impatient agents
r_π::Float64 = 1.5379;              # Inflation coefficient in Taylor rule
r_R::Float64 = 0.5509;              # Interest rate coefficient in Taylor rule
r_y::Float64 = 0.0944;              # Output coefficient in Taylor rule
θ_π::Float64 = 0.8913;              # Calvo parameter for prices
θ_w::Float64 = 0.9163;              # Calvo parameter for wages
γ::Float64 = 0.6945;                # Inertia in borrowing constraint
ρ_J::Float64 = 0.9835;              # Persistence in housing shock
ρ_K::Float64 = 0.7859;              # Persistence in investment shock
ρ_R::Float64 = 0.6232;              # Persistence in monetary policy shock
ρ_Z::Float64 = 0.7556;              # Persistence in preference shock
σ_J::Float64 = 0.0079;              # Standard deviation of housing shock
σ_K::Float64 = 0.0079;              # Standard deviation of investment shock
σ_P::Float64 = 0.0079;              # Standard deviation of price markup shock
σ_R::Float64 = 0.0079;              # Standard deviation of monetary policy shock
σ_W::Float64 = 0.0079;              # Standard deviation of wage markup shock
σ_Z::Float64 = 0.0079;              # Standard deviation of preference shock

# Auxiliary parameters
ε_π::Float64 = (1 - θ_π) * (1 - β * θ_π) / θ_π;
ε_w::Float64 = (1 - θ_w) * (1 - β * θ_w) / θ_w;
ε_wp::Float64 = (1 - θ_w) * (1 - βp * θ_w) / θ_w;
Γ_c::Float64 = (1 - ε_c) / (1 - β * ε_c);
Γ_c_p::Float64 = (1 - ε_c) / (1 - βp * ε_c);
Γ_h::Float64 = (1 - ε_h) / (1 - β * ε_h);
Γ_h_p::Float64 = (1 - ε_h) / (1 - βp * ε_c);

#--------------------------------------------------------------------------------------------------------------------
# Steady solution
#--------------------------------------------------------------------------------------------------------------------
var_names = [
    # Normal variables
    "c",        # 1  
    "c'",       # 2
    "h",        # 3
    "h'",       # 4
    "i",        # 5
    "k",        # 6
    "y",        # 7
    "b",        # 8
    "n",        # 9
    "n'",       # 10
    "w",        # 11
    "w'",       # 12
    "π",        # 13
    "q",        # 14
    "R",        # 15
    "λ",        # 16
    "x_p",      # 17
    "x_w",      # 18
    "x'_w",     # 19
    "r_k",      # 20
    "q_k",      # 21
    # Exogenous shock variables
    "a",        # 22
    "j",        # 23
    "z",        # 24
    "e_r",      # 25
    # Auxiliary variables
    "u_c",      # 26
    "u_c'",     # 27
    "u_h",      # 28
    "u_h'",     # 29
    "u_n",      # 30
    "u_n'",     # 31
    "div",      # 32
    "div'",     # 33
    "ω",        # 34
    "ω'",       # 35
    "πA"        # 36
]

function ccF!(F, x)
    F[1] = x[1] + x[5] - (x[15] * x[8]) / x[13] - (x[11] * x[9]) / x[18] - x[20] * x[6] + x[8] - x[32]
    F[2] = x[26] * x[21] - x[26]
    F[3] = x[26] * (x[21] / x[22]) - β * (x[26] * (x[20] + x[21] * ((1 - δ_k) / x[22])))
    F[4] = x[6] - x[22] * x[5] - (1 - δ_k) * x[6]
    F[5] = x[26] - β * (x[26] * x[15] / x[13])
    F[6] = (x[11] / x[18]) * x[26] - x[30]
    F[7] = x[14] * x[26] - x[28] - β * x[14] * x[27]
    F[8] = x[2] + (x[15] / x[13]) * x[8] - (x[12] / x[19]) * x[10] - x[8] - x[33]
    # F[9] //////
    F[9] = x[16]
    F[10] = (1 - x[16]) * x[27] - βp * (x[15] - γ * x[16]) * x[27] / x[13]
    F[11] = (x[12] / x[19]) * x[27] - x[31]
    F[12] = x[14] * x[27] - x[29] - βp * x[14] * x[27] - x[27] * x[16] * (1 - γ) * m * x[14]
    F[13] = x[7] - x[9]^((1 - σ) * (1 - α)) * x[10]^(σ * (1 - α)) * x[6]^α
    F[14] = (1 - α) * (1 - σ) * x[7] - x[17] * x[11] * x[9]
    F[15] = (1 - α) * σ * x[7] - x[17] * x[12] * x[10]
    F[16] = α * x[7] - x[17] * x[20] * x[6]
    F[17] = log(x[13] / π_ss) - β * (log(x[13] / π_ss)) + ε_π * (log(x[17] / x_p_ss))
    F[18] = log(x[34] / π_ss) - β * (log(x[34] / π_ss)) + ε_w * (log(x[18] / x_w_ss))
    F[19] = log(x[35] / π_ss) - βp * (log(x[35] / π_ss)) + ε_wp * (log(x[19] / x_w_ss))
    F[20] = (x[15]^r_R) * ((x[36] / π_ss)^((1 - r_R) * r_π)) * ((x[7] / x[7])^((1 - r_R) * r_y))
    F[21] = x[3] + x[4] - 1
    F[22] = log(x[22]) - ρ_K * log(x[22])
    F[23] = log(x[23]) - (1 - ρ_J) - ρ_J * log(x[23])
    F[24] = log(x[24]) - ρ_Z * log(x[24])
    F[25] = log(x[25]) - ρ_R * log(x[25])
    F[26] = x[26] - x[24] * Γ_c * (1 / (x[1] - ε_c * x[1]))
    F[27] = x[27] - x[24] * Γ_c_p * (1 / (x[2] - ε_c * x[2]))
    F[28] = x[28] - x[24] * x[23] * Γ_h * (1 / (x[3] - ε_c * x[3]))
    F[29] = x[29] - x[24] * x[23] * Γ_h_p * (1 / (x[4] - ε_c * x[4]))
    F[30] = x[30] + x[24] * x[9]^η
    F[31] = x[31] + x[24] * x[10]^η
    F[32] = x[32] - ((x[17] - 1) / x[17]) * x[7] - ((x[18] - 1) / x[17]) * x[18] * x[9]
    F[33] = x[33] - ((x[19] - 1) / x[19]) * x[12] * x[10]
    F[34] = x[34] - x[13]
    F[35] = x[35] - x[13]
    F[36] = x[36] - x[13]
end

ss_guess = ones(36)
ss = nlsolve(ccF!, ss_guess)

println("Goodbye bestie")