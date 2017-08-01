using PyPlot, Interpolations

type AiyagariEconomy

    # Define de parameters
    β::Float64
    γ::Float64
    α::Float64
    δ::Float64
    ρ::Float64
    r::Float64
    w::Float64
    σe::Float64
    σy::Float64

    # Define grids
    my::Float64
    ny::Int64
    ne::Int64
    na::Int64
    num_nodes::Int64
    ygrid::Vector{Float64}
    ϵgrid::Vector{Float64}
    agrid::Vector{Float64}
    a_ast_mat::Array{Float64}
    agrid_finer::Vector{Float64}
    snodes::Array{Float64, 2}

    # Probability matrix
    Πy::Array{Float64, 2}
    Πyinv::Vector{Float64}

    # Aggregate variables
    K::Float64
    H::Float64

    # Policy Functions
    cpol_mat::Array{Float64, 2}

    # Define options
    dist_tol::Float64
    maxit::Int
    lb::Float64
    ub::Float64
end


function aiyagari_fun(; β   = 0.96,
                        γ   = 2.0,
                        α   = 1/3,
                        δ   = 0.08,
                        ρ   = 0.9,
                        r   = 0.029,
                        w   = 0.0,
                        σe  = 0.1,
                        my  = 0.0,
                        ny  = 7,
                        ne  = ny,
                        na  = 100,
                        dist_tol = 1e-5,
                        maxit = 1000,
                        lb  = 0.037 )


    ub      = 1.0/β-1.0-0.0001
    my      = 0
    σy      = sqrt((σe^2)/(1-ρ^2))
    amax    = 60.0
    a_temp  = log(linspace(exp(0.0), na, na))
    agrid   = (cumsum(a_temp)/sum(a_temp) * amax)


    """ Rouwenhorst Method """
    p = (1.0+ρ)/2; q = p
    ψ = σy * sqrt(ny-1)

    ymax  = ψ; ymin = -ψ
    ygrid = linspace(ymin, ymax, ny)

    Πy2 = [[p 1.0-p];[1.0-q q]]
    Πyn1 = copy(Πy2)

    for jj=1:(ny-2)
        num_rows = size(Πyn1, 1)
        mat1     = zeros(num_rows+1, num_rows+1)
        mat2, mat3, mat4 = copy(mat1), copy(mat1), copy(mat1)

        mat1[1:end-1, 1:end-1]  = Πyn1
        mat2[1:end-1, 2:end]    = Πyn1
        mat3[2:end, 1:end-1]    = Πyn1
        mat4[2:end, 2:end]      = Πyn1

        Πyn1 = p*mat1 + (1-p)*mat2 + (1-q)*mat3 + q*mat4
        Πyn1[2:end-1, :] = Πyn1[2:end-1, :] / 2
    end

    Πy    = copy(Πyn1)

    vals, vecs  = eig(Πy')
    _, ind_val  = findmin(abs(vals-1.0))
    Πyinv       = vecs[:, ind_val]/sum(vecs[:, ind_val])

    sum(Πyinv.>=0.0) == ny || throw(error("Negative elements in invariant distribution"))

    ϵgrid   = exp(ygrid + my*ones(ny))
    H       = dot(ϵgrid, Πyinv)
    K       = 0.0

    cpol_mat = zeros(na, ny)

    a_ast_mat   = repmat(agrid, 1, ny)
    # agrid_finer = linspace(agrid[1], agrid[end], 2*na)

    a_finer_temp  = log(linspace(exp(0.0), 3*na, 3*na))
    agrid_finer   = (cumsum(a_temp)/sum(a_temp) * amax)


    snodes      = [repmat(agrid_finer, ny) kron(ygrid, agrid_finer)]
    num_nodes   = size(snodes, 1)

    return AiyagariEconomy( β, γ, α, δ, ρ, r, w, σe, σy, my, ny, ne, na, num_nodes, ygrid, ϵgrid, agrid, a_ast_mat, agrid_finer, snodes, Πy, Πyinv, K, H, cpol_mat, dist_tol, maxit, lb, ub )
end


### Utility Function ###
uprime(ay::AiyagariEconomy, c) = c.^(-ay.γ)
upinv(ay::AiyagariEconomy, vv) = vv.^(-1.0/ay.γ)


### Production Function and F.O.C. ###
prod_fun(ay::AiyagariEconomy, K) = K^ay.α*ay.H^(1.0-ay.α)                               # computes output
mpk_fun(ay::AiyagariEconomy, K)  = (ay.α*K^(ay.α-1.0)*ay.H^(1.0-ay.α))                  # computes r
mpl_fun(ay::AiyagariEconomy, K)  = (1.0-ay.α)*K^ay.α*ay.H^-ay.α                         # computes w
inv_mpk(ay::AiyagariEconomy, rr) = ay.H^(1.0-ay.α)*(ay.α/(rr+ay.δ))^(1.0/(1.0-ay.α))  # computes k


function endogrid!(ay::AiyagariEconomy)

    R = 1.0+ay.r
    b = ay.agrid[1]

    zgrid = ay.ϵgrid*ay.w
    
    cpol_out = Array(Float64, ay.na, ay.ne)

    B = (ay.β * R * ay.Πy * uprime(ay, ay.cpol_mat'))' # RHS of the Bellman Equation

    # Solve for the ctil policy
    ctil  = upinv(ay, B)

    # Compute the endoegenous grid. Save it for piecewise linear approximation of the invariant distribution
    a_ast = (repmat(ay.agrid, 1, ay.ne) + ctil - repmat(zgrid', ay.na, 1)) / R # value of assets today that induces aprime tomorrow
        
    # Update the consumption policy function 
    for (i_z, z_v) in enumerate(zgrid)
       for (i_a, a_v) in enumerate(ay.agrid)
            if a_v<= a_ast[1, i_z]; # case where you are constrained
              cpol_out[i_a, i_z] = R*a_v + b + z_v

            elseif  a_v >= a_ast[end, i_z] # out of the range (to the right), linearly extrapolate
              cpol_out[i_a,i_z] = ctil[end, i_z] + (a_v-a_ast[end, i_z])*(ctil[end, i_z] - ctil[end-1, i_z])/(a_ast[end, i_z]-a_ast[end-1, i_z])

            else # inside the range, linearly interpolate
              ind1 = searchsortedfirst(a_ast[:,i_z], a_v)
              ind  = ind1-1
              cpol_out[i_a, i_z] = ctil[ind, i_z] + (a_v-a_ast[ind, i_z])*(ctil[ind1, i_z] - ctil[ind, i_z])/(a_ast[ind1, i_z]-a_ast[ind, i_z])
            end    
          
       end
    end

    copy!(ay.cpol_mat, cpol_out)
    copy!(ay.a_ast_mat, a_ast)

    Void
end


function iterate_pol!(ay::AiyagariEconomy; do_plots::Bool = false)

    iter = 0
    dist = 10.0

    while iter < ay.maxit && dist > ay.dist_tol

        iter    += 1
        old_pol = copy(ay.cpol_mat)

        endogrid!(ay)    

        dist = maxabs(old_pol - ay.cpol_mat)

        if iter%100 == 0
            @printf("Finished EndoGdrid iteration %d with distance %.3g\n", iter, dist)
        end
    end

    if do_plots

        ap_egm = (1.0+ay.r) * repmat(ay.agrid, 1, ay.ny) + ay.w*repmat(ay.ϵgrid', ay.na, 1) - ay.cpol_mat

        fig, axis = subplots(1, 2)
        ax = axis[1]; ax[:plot](ay.agrid, ap_egm);      ax[:set_title]("Assets Policy at Equilibrium")
        ax = axis[2]; ax[:plot](ay.agrid, ay.cpol_mat); ax[:set_title]("Consumption Policy at Equilibrium")
    end

    Void
end

function compute_invariant_pwlinear!(ay::AiyagariEconomy, Λ_invariant)

    Λn_mat   = zeros(length(ay.agrid_finer), ay.ny)
    Λnm1_mat = copy(Λ_invariant)

    iter = 0
    dist = 10.0

    a_ast_itp = interpolate((ay.agrid, ay.ygrid), ay.a_ast_mat, Gridded(Linear())) # next period's assets, current y

    while iter<2000 && dist>1e-4

        iter += 1

        for i_y = 1:ay.ny # next period
            for (i_a, a_v) in enumerate(ay.agrid_finer) # next period

                vec_temp = zeros(ay.ny)

                for (i_y0, y0_v) in enumerate(ay.ygrid) # last period y
                    aval  = min(max(a_ast_itp[a_v, y0_v], ay.agrid[1]), ay.agrid[end]) # today's assets (endogenous grid)
                    ind_r = min(max(searchsortedfirst(ay.agrid_finer, aval), 2), length(ay.agrid_finer))

                    Λval = Λnm1_mat[ind_r-1, i_y0] + (Λnm1_mat[ind_r, i_y0]- Λnm1_mat[ind_r-1, i_y0]) / (ay.agrid_finer[ind_r] - ay.agrid_finer[ind_r-1]) * (aval - ay.agrid_finer[ind_r-1])
                    vec_temp[i_y0] = Λval
                end

                Λn_mat[i_a, i_y] = dot(ay.Πy[:, i_y], vec_temp)
            end
        end

        dist = maxabs(Λn_mat - Λnm1_mat)
        copy!(Λnm1_mat, Λn_mat)

        # if iter%200 == 0.0
        #     @printf("Iteration # %d on distribution with distance %.3g\n", iter, dist)
        # end
    end

    copy!(Λ_invariant, Λnm1_mat)

    Void
end


function compute_supply(ay::AiyagariEconomy, rr)

    r0_fixed = -0.1*ay.δ + 0.9*(1.0/ay.β - 1.0)

    """ Initialize consumption policy """
    KK0   = inv_mpk(ay, r0_fixed)
    ww0   = mpl_fun(ay, KK0)
    zgrid = ay.ϵgrid*ww0

    cpol_init = zeros(ay.na, ay.ny)

    for (i_z, z_v) in enumerate(zgrid)
        for (i_a, a_v) in enumerate(ay.agrid)
            c_max = (1.0 + r0_fixed)*a_v + ww0*z_v-0.05
            cpol_init[i_a, i_z] = c_max
        end
    end

    copy!(ay.cpol_mat, cpol_init)

    ay.r = rr
    ay.K = inv_mpk(ay, ay.r)
    ay.w = mpl_fun(ay, ay.K)

    copy!(ay.cpol_mat, cpol_init)

    """ Solve using EndoGrid """
    iterate_pol!(ay)


    """ Compute invariant distribution given r0 """
    Λ_invariant = zeros(length(ay.agrid_finer), ay.ny)

    for i_y=1:ay.ny
        for (i_a, a_v) in enumerate(ay.agrid_finer)
            Λ_invariant[i_a, i_y] = (a_v - ay.agrid[1]) / (ay.agrid[end] - ay.agrid[1]) * ay.Πyinv[i_y]
        end
    end

    compute_invariant_pwlinear!(ay, Λ_invariant)


    """ Compute supply of capital """

    A_supply = 0.0

    for i_y = 1:ay.ny

        sum_val_a = 0.0

        for i_a = 1:(length(ay.agrid_finer)-1)
            anp1 = ay.agrid_finer[i_a+1]
            an   = ay.agrid_finer[i_a]
            sum_val_a += 0.5*(Λ_invariant[i_a+1, i_y] - Λ_invariant[i_a, i_y]) * (anp1 + an)
        end

        A_supply += sum_val_a + Λ_invariant[1, i_y]*ay.agrid_finer[1]
    end

    return A_supply, Λ_invariant
end


ay = aiyagari_fun()

nnr   = 20
rvec  = linspace(0.03, 1/ay.β-1-0.001, nnr)
A_vec = zeros(nnr)
K_vec = zeros(nnr)

for i_r=1:nnr
    A_vec[i_r], _ = compute_supply(ay, rvec[i_r])
    K_vec[i_r] = inv_mpk(ay, rvec[i_r])
end

fig, ax = subplots()
ax[:plot](A_vec, [rvec ones(rvec)*(1/ay.β-1)]); ax[:set_title]("Supply")

fig, ax = subplots()
ax[:plot](K_vec, rvec); ax[:set_title]("Demand")

fig, ax = subplots()
ax[:plot](A_vec, rvec); ax[:plot](K_vec, rvec); ax[:set_xlim]([4; 10]); ax[:set_title]("Equilibrium")

AS, Λ_invar = compute_supply(ay, 0.036)

distrib_a = zeros(length(ay.agrid_finer))

for i_a=1:length(ay.agrid_finer)
    distrib_a[i_a] = sum(Λ_invar[i_a, :])
end


fig, ax = subplots()
ax[:plot](ay.agrid_finer, distrib_a); ax[:set_title]("Distribution of Savings")