mutable struct CCMAParams
    # Calibrated parameters
    m::Float64;                     # Maximum LTV
    η::Float64;                     # Labour disutility
    β::Float64;                     # Discount factor for patient agents
    π_ss::Float64;                  # Steady state gross inflation
    α::Float64;                     # Capital share of output
    δ_k::Float64;                   # Capital depreciation rate
    j_bar::Float64;                 # Housing weight in utility function
    x_p_ss::Float64;                # Steady state price markup
    x_w_ss::Float64;                # Steady state wage markup

    # Estimated parameters (posterior mode)
    βp::Float64;                    # Discount factor for impatient agents
    ε_c::Float64;                   # Habit in consumption
    ε_h::Float64;                   # Habit in housing
    φ::Float64;                     # Investment adjustment cost
    σ::Float64;                     # Wage share for impatient agents
    r_π::Float64;                   # Inflation coefficient in Taylor rule
    r_R::Float64;                   # Interest rate coefficient in Taylor rule
    r_y::Float64;                   # Output coefficient in Taylor rule
    θ_π::Float64;                   # Calvo parameter for prices
    θ_w::Float64;                   # Calvo parameter for wages
    γ::Float64;                     # Inertia in borrowing constraint
    ρ_J::Float64;                   # Persistence in housing shock
    ρ_K::Float64;                   # Persistence in investment shock
    ρ_R::Float64;                   # Persistence in monetary policy shock
    ρ_Z::Float64;                   # Persistence in preference shock
    σ_J::Float64;                   # Standard deviation of housing shock
    σ_K::Float64;                   # Standard deviation of investment shock
    σ_P::Float64;                   # Standard deviation of price markup shock
    σ_R::Float64;                   # Standard deviation of monetary policy shock
    σ_W::Float64;                   # Standard deviation of wage markup shock
    σ_Z::Float64;                   # Standard deviation of preference shock

    # Auxiliary parameters
    ε_π::Float64;
    ε_w::Float64;
    ε_wp::Float64;
    Γ_c::Float64;
    Γ_c_p::Float64;
    Γ_h::Float64;
    Γ_h_p::Float64;
end

function initialise_ccma_params(
    m::Float64,
    η::Float64,
    β::Float64,
    π_ss::Float64,
    α::Float64,
    δ_k::Float64,
    j_bar::Float64,
    x_p_ss::Float64,
    x_w_ss::Float64,
    βp::Float64,
    ε_c::Float64,
    ε_h::Float64,
    φ::Float64,
    σ::Float64,
    r_π::Float64,
    r_R::Float64,
    r_y::Float64,
    θ_π::Float64,
    θ_w::Float64,
    γ::Float64,
    ρ_J::Float64,
    ρ_K::Float64,
    ρ_R::Float64,
    ρ_Z::Float64,
    σ_J::Float64,
    σ_K::Float64,
    σ_P::Float64,
    σ_R::Float64,
    σ_W::Float64,
    σ_Z::Float64
)
    ε_π::Float64 = (1 - θ_π) * (1 - β * θ_π) / θ_π;
    ε_w::Float64 = (1 - θ_w) * (1 - β * θ_w) / θ_w;
    ε_wp::Float64 = (1 - θ_w) * (1 - βp * θ_w) / θ_w;
    Γ_c::Float64 = (1 - ε_c) / (1 - β * ε_c);
    Γ_c_p::Float64 = (1 - ε_c) / (1 - βp * ε_c);
    Γ_h::Float64 = (1 - ε_h) / (1 - β * ε_h);
    Γ_h_p::Float64 = (1 - ε_h) / (1 - βp * ε_c);

    params = CCMAParams(
        m,
        η,
        β,
        π_ss,
        α,
        δ_k,
        j_bar,
        x_p_ss,
        x_w_ss,
        βp,
        ε_c,
        ε_h,
        φ,
        σ,
        r_π,
        r_R,
        r_y,
        θ_π,
        θ_w,
        γ,
        ρ_J,
        ρ_K,
        ρ_R,
        ρ_Z,
        σ_J,
        σ_K,
        σ_P,
        σ_R,
        σ_W,
        σ_Z,
        ε_π,
        ε_w,
        ε_wp,
        Γ_c,
        Γ_c_p,
        Γ_h,
        Γ_h_p
    )

    return params
end

ccma_params = initialise_ccma_params(
    0.9,                # m
    1.0,                # η
    0.995,              # β
    1.005,              # π_ss
    0.33,               # α
    0.025,              # δ_k
    0.04,               # j_bar
    1.2,                # x_p_ss
    1.2,                # x_w_ss
    0.9922,             # βp
    0.6842,             # ε_c
    0.8799,             # ε_h
    4.1209,             # φ
    0.5013,             # σ
    1.5379,             # r_π
    0.5509,             # r_R
    0.0944,             # r_y
    0.8913,             # θ_π
    0.9163,             # θ_w
    0.6945,             # γ
    0.9835,             # ρ_J
    0.7859,             # ρ_K
    0.6232,             # ρ_R
    0.7556,             # ρ_Z
    0.0079,             # σ_J
    0.0079,             # σ_P
    0.0079,             # σ_R
    0.0079,             # σ_W
    0.0079,             # σ_Z
    0.0079,             # σ_J
)