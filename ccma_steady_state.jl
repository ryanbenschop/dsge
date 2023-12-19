include("ccma_parameterisation.jl")

mutable struct CCMASteadyState
    c_ss::Float64
    c_p_ss::Float64
    h_ss::Float64
    h_p_ss::Float64
    i_ss::Float64
    k_ss::Float64
    y_ss::Float64
    b_ss::Float64
    n_ss::Float64
    n_p_ss::Float64
    w_ss::Float64
    w_p_ss::Float64
    π_ss::Float64
    q_ss::Float64
    R_ss::Float64
    λ_ss::Float64
    x_p_ss::Float64
    x_w_ss::Float64
    x_w_p_ss::Float64
    r_k_ss::Float64
    q_k_ss::Float64
    j_ss::Float64
    z_ss::Float64
    a_ss::Float64
    e_ss::Float64
    u_c_ss::Float64
    u_c_p_ss::Float64
    u_h_ss::Float64
    u_h_p_ss::Float64
    u_n_ss::Float64
    u_n_p_ss::Float64
    div_ss::Float64
    div_p_ss::Float64
    ω_ss::Float64
    ω_p_ss::Float64
    π_A_ss::Float64
    bnot_ss::Float64
    maxlev_ss::Float64
    rnot_ss::Float64
end

function compute_ccma_steady_state(params::CCMAParams)
    R_ss = params.π_ss / params.β
    r_k_ss = 1 / params.β - (1 - params.δ_k)
    λ_ss = (1 - params.βp / params.β) / (1 - params.βp * params.γ / params.π_ss)

    qhtoc = params.j_bar / (1 - params.β)
    qhtocp = params.j_bar / (1 - params.βp - λ_ss * (1 - params.γ) * params.m)
    ktoy = params.α / (params.x_p_ss * r_k_ss)
    btoqh1 = params.m * (1 - params.γ) / (1 - params.γ / params.π_ss)
    c1toy = (1 - params.α) * params.σ / (1 + (1 / params.β - 1) * btoqh1 * qhtocp) * (1 / params.x_p_ss)
    ctoy = 1 - c1toy - params.δ_k * ktoy

    n_ss = ((1 - params.σ) * (1 - params.α) / (params.x_p_ss * params.x_w_ss * ctoy))^(1 / (1 + params.η))
    n_p_ss = (params.σ * (1 - params.α) / (params.x_p_ss * params.x_w_ss * c1toy))^(1 / (1 + params.η))

    y_ss = ktoy^(params.α / (1 - params.α)) * n_ss^(1 - params.σ) * n_p_ss^params.σ
    c_ss = ctoy * y_ss
    c_p_ss = c1toy * y_ss
    k_ss = ktoy * y_ss
    i_ss = params.δ_k * k_ss

    w_ss = params.x_w_ss * c_ss * n_ss^params.η
    w_p_ss = params.x_w_ss * c_p_ss * n_p_ss^params.η

    q_ss = qhtoc * c_ss + qhtocp * c_p_ss

    h_ss = qhtoc * c_ss / q_ss
    h_p_ss = qhtocp * c_p_ss / q_ss

    b_ss = btoqh1 * q_ss * h_p_ss

    u_c_ss = 1 / c_ss
    u_c_p_ss = 1 / c_p_ss
    u_h_ss = params.j_bar / h_ss
    u_h_p_ss = params.j_bar / h_p_ss
    u_n_ss = n_ss^params.η
    u_n_p_ss = n_p_ss^params.η

    q_k_ss = 1

    ω_ss = params.π_ss
    ω_p_ss = params.π_ss

    div_ss = c_ss + i_ss - R_ss * b_ss / params.π_ss - w_ss * n_ss / params.x_w_ss - r_k_ss * k_ss + b_ss
    div_p_ss = c_p_ss + R_ss * b_ss / params.π_ss - w_p_ss * n_p_ss / params.x_w_ss - b_ss

    π_ss = params.π_ss
    x_p_ss = params.x_p_ss
    x_w_ss = params.x_w_ss
    x_w_p_ss = params.x_w_ss

    j_ss = params.j_bar
    z_ss = 1
    a_ss = 1
    e_ss = 1

    π_A_ss = π_ss

    bnot_ss = b_ss

    maxlev_ss = 0

    rnot_ss = R_ss

    ss = CCMASteadyState(
        c_ss,
        c_p_ss,
        h_ss,
        h_p_ss,
        i_ss,
        k_ss,
        y_ss,
        b_ss,
        n_ss,
        n_p_ss,
        w_ss,
        w_p_ss,
        π_ss,
        q_ss,
        R_ss,
        λ_ss,
        x_p_ss,
        x_w_ss,
        x_w_p_ss,
        r_k_ss,
        q_k_ss,
        j_ss,
        z_ss,
        a_ss,
        e_ss,
        u_c_ss,
        u_c_p_ss,
        u_h_ss,
        u_h_p_ss,
        u_n_ss,
        u_n_p_ss,
        div_ss,
        div_p_ss,
        ω_ss,
        ω_p_ss,
        π_A_ss,
        bnot_ss,
        maxlev_ss,
        rnot_ss
    )
end

ccma_ss = compute_ccma_steady_state(ccma_params)

function compute_log_ss(ss::CCMASteadyState)
    log_ss = [
        log(ss.c_ss),
        log(ss.c_p_ss),
        log(ss.h_ss),
        log(ss.h_p_ss),
        log(ss.i_ss),
        log(ss.k_ss),
        log(ss.y_ss),
        log(ss.b_ss),
        log(ss.n_ss),
        log(ss.n_p_ss),
        log(ss.w_ss),
        log(ss.w_p_ss),
        log(ss.π_ss),
        log(ss.q_ss),
        log(ss.R_ss),
        log(ss.λ_ss),
        log(ss.x_p_ss),
        log(ss.x_w_ss),
        log(ss.x_w_p_ss),
        log(ss.r_k_ss),
        log(ss.q_k_ss),
        log(ss.j_ss),
        log(ss.z_ss),
        log(ss.a_ss),
        log(ss.e_ss),
        log(ss.u_c_ss),
        log(ss.u_c_p_ss),
        log(ss.u_h_ss),
        log(ss.u_h_p_ss),
        log(ss.u_n_ss),
        log(ss.u_n_p_ss),
        log(ss.div_ss),
        log(ss.div_p_ss),
        log(ss.ω_ss),
        log(ss.ω_p_ss),
        log(ss.π_A_ss),
        log(ss.bnot_ss),
        ss.maxlev_ss,
        log(ss.rnot_ss)
    ]
end

ccma_log_ss = compute_log_ss(ccma_ss)