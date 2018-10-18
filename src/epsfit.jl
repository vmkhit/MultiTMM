@inline function eps_drude(w, epsb, wp, gammap)
    epsd = epsb - wp^2/w/(w + 1.0im*gammap)
    return epsd
end

@inline function eps_lorentzian(w, f, wp, w0, gammap)
    return f*wp^2/(w0^2 - w^2 - im*w*gammap)
end

@inline function eps_debay1(w, eps0, epsinf, wT, gammap)
    return epsinf + (epsinf - eps0)/(1 - (w/wT)^2 - 1.0im*gammap*(w/wT))
end

@inline function eps_debay2(w, epsinf, wL, wT, gammaL, gammaT)
    return epsinf*(wL^2 - w^2 + 1.0im*w*gammaL)/(wT^2 - w^2 + 1.0im*w*gammaT)
end

mutable struct Drude
    epsb
    wp
    gammap
end

mutable struct Lorentzian
    f
    wp
    w0
    gammap
end

mutable struct Debay1
    eps0
    epsinf
    wT
    gammap
end

mutable struct Debay2
    epsinf
    wL
    wT
    gammaL
    gammaT
end

mutable struct File
    filename
end



function nonlinear_fit(x, fun, a0, eps=1e-8, maxiter=200)
    na = length(a0)
    np = size(x, 1)
    nv = size(x, 2)
    for i in 1:na
        if a0[i] == 0
            a0[i] = 0.01
        end
    end
    da = a0 / 10
    a1 = zeros(na)
    a = zeros(na)
    for k = 1:na
        a1[k] = a0[k] - da[k]
    end
    xp = zeros(nv)
    A = zeros(np, na)
    r0 = zeros(np)
    r1 = zeros(np)
    for p = 1:np
        for k = 1:nv
            xp[k] = x[p,k]
        end
        r0[p] = fun(xp, a0)
        r1[p] = fun(xp, a1)
    end
    maxerr = maximum(abs, [maximum(abs, r0), maximum(abs, r1)])
    iter = 1
    convergence = false
    for i = 1:maxiter
        iter = i
        for p = 1:np
            for k = 1:nv
                xp[k] = x[p,k]
            end
            for k = 1:na
                aa = a1[k]
                a1[k] = a1[k] + da[k]
                r = fun(xp, a1)
                a1[k] = aa
                A[p,k] = -(r1[p] - r) / da[k]
            end
        end
        da = A \ r1
        for k = 1:na
            a1[k] = a1[k] - da[k]
        end
        for p = 1:np
            r0[p] = r1[p]
            for k = 1:nv

                xp[k] = x[p,k]
            end
            r1[p] = fun(xp, a1)
        end
        if maximum(abs, r1) < eps * maxerr
            convergence = true
            break
        end
    end
    return a1, convergence, iter

end



tmm_ellipso(lambda::Real, kp::Vector{<:Real}, S::Stack)

function fit_ellips(wexp::Vector{Number}, TanPsi_exp::Vector{Number}, CosDelta_epx::Vector{Number}, theta_inc::Real, dlist, epslist, siglist; weight = 0.5, maxitter = 5000)
    local const ev2nm = 1239.84193
    lam0 = ev2nm/wexp
    k = 2.0*pi/lam0
    const  kp::Vector = k.*[sin(theta_inc), 0]

    mod = Model(solver = Ipopt.IpoptSolver(max_iter=maxitter))

    JuMP.register(mod, :eps_drude,       4, eps_drude,       autodiff=true)
    JuMP.register(mod, :eps_lorentzian,  5, eps_lorentzian,  autodiff=true)
    JuMP.register(mod, :eps_debay1,      5, eps_debay1,      autodiff=true)
    JuMP.register(mod, :eps_debay2,      6, eps_debay2,      autodiff=true)
    JuMP.register(mod, :tmm_ellipso,     3, tmm_ellipso,     autodiff=true)

    @variable(mod, 0.0 <=epsb[i=1:npeaks]<= 1, start=amp_guess[i])
    @variable(mod, 0.0 <=[i=1:npeaks]<= 1, start=pos_guess[i])
    @variable(mod, 0.0 <=gam[i=1:npeaks]<= 1, start=gam_guess[i])


    local const nlayers = length(dlist)
    Stk = Stack([Layer()],[0])
    for i = 1:nlayers
        mat = Material("eps", epslist[i])
        modify_stack_layers(Stk, 1, Layer(mat, dlist[i]), siglist[i])
    end


    @NLexpression(mod, shape[j=1:npoints, i=1:npeaks], Gaussian_model(xexp[j], amp[i], pos[i], gam[i]))


    @NLobjective(mod, Min, sum((sum(shape[j, i] for i = 1:npeaks) - yexp[j])^2/npoints for j = 1:npoints))
    solve(mod) # solving the model

    params = [getvalue(amp) getvalue(pos) getvalue(sig) getvalue(gam)]

    return params
end
