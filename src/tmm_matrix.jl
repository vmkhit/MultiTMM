using LinearAlgebra
function epston(epsilon)
    n = sqrt(abs(epsilon) + real(epsilon))/sqrt(2)
    k = sqrt(abs(epsilon) - real(epsilon))/sqrt(2)
    return n+im*k
end
mutable struct Material
    id::String
    eps::Number
    mu::Number
    nk::Number
    function Material(id::String = "nk", kwargs...)
        @assert id == "epsmu" || id == "nk"
        if isempty(kwargs)
            eps = mu = nk = 1.0
        else
            if id == "nk"
                nk = kwargs[1]
                eps = nk^2
                mu = 1.0
            else
                if length(kwargs) < 2
                    eps = kwargs[1]
                    mu = 1.0
                    nk = epston(eps)
                else
                    eps, mu = kwargs
                    nk = epston(eps*mu)
                end
            end
        end
        return new(id, eps, mu, nk)
    end
end
# Material can be specified as Material(eps, mu), or just Material(eps), this will by default take mu = 1.0

mutable struct Interface
    mat1::Material   # dielectric function of medium 1
    mat2::Material   # dielectric function of medium 2
    sig::Number      # conductivity at the interface
    function Interface(mat1 = Material(), mat2 = Material(), args...)
        if isempty(args)
            sig = 0.0
        else
            sig = args[1]
        end
        return new(mat1, mat2, sig)
     end
end

mutable struct Layer
    id::String  # coherent, incoherent
    mat::Material  # layer material
    d::Real        # layer thickness
    function Layer(id ="c", mat = Material(), kwargs...)
        if isempty(kwargs)
            d = 0.0
        else
            d = kwargs[1]
        end
        return new(id, mat, d)
    end
end

function extract_params(input::Union{Material, Interface, Layer})
    if typeof(input) == Material
        return input.eps, input.mu
    elseif typeof(input) == Interface
        return input.mat1, input.mat2, input.sig
    else
        return input.id, input.mat, input.d
    end
end

function update_interface_list(L::Vector{Layer}, sig::Vector{<:Number}, nl::Integer)
    Is =  Vector{Interface}(undef, nl -1)
    for i = 1:(nl -1)
        Is[i] = Interface(L[i].mat, L[i+1].mat, sig[i])
    end
    return Is
end

mutable struct Stack
    Layers::Vector{Layer} #List of layers
    Sigmas::Vector{Number}
    Interfs::Vector{Interface}  # List of Interfaces
    function Stack(Layers, Sigmas)
        nL = length(Layers)
        Interfs = update_interface_list(Layers, Sigmas, nL)
        return new(Layers, Sigmas, Interfs)
    end
end


function modify_stack_layers(S::Stack, r::UnitRange{<:Integer}, L = Layer[], sigs = Number[])
    #This isn't finished yet, need to modified to add the interface update as well
    splice!(S.Layers, r, L)
    splice!(S.Sigmas, r, sigs)
    S.Interfs = update_interface_list(S.Layers, S.Sigmas, length(S.Layers))
end

function betz(mat::Material, k::Number, q::Real)
    return sqrt(mat.eps*mat.mu*k^2 - q^2 + 0.0im)
end
function intface_rt(p::Integer, Iij::Interface, k0::Real, kp::Vector{<:Real})
    # Nottice: This is also coherent with Novotny-Hecht, since r, t
    # in both cases are just the ratios of the electric fields for given polarization
    # in contrast to other formulations with magnetic field ratio for p-polarized incidence
    @assert p == 1 || p == 2
    local q = norm(kp)
    mat1, mat2, σ = extract_params(Iij)
    kz1 = betz(mat1, k0, q)
    kz2 = betz(mat2, k0, q)
    gs = mat1.mu/mat2.mu
    if p == 1 # s-polarization
        g = gs
    else
        g = mat1.eps/mat2.eps
    end
    # these are only for ideal case i.e. for smooth interfaces with no conductive sheet
    # the roughfness and the conductive sheets or polarizable sheets are not implimeted yet
    r = (kz1 - g*kz2)/(kz1 + g*kz2)
    t = sqrt(g/gs)*(1 + r)
    trev = mat1.mu*kz2*t/(mat2.mu*kz1)
    rrev = -r
    return r, t, rrev, trev
end

function Matrix2x2_coh(r, t, rrev, trev, d)
    # I used here the general form, which will allow me to impliment the rough surface case as well
    p = cis(-d)
    return [p r*p; r/p (t*trev - r*rrev)/p]./t
end

function Matrix2x2_inc(r, t, rrev, trev, d)
    p = abs2(cis(-d))
    [p -abs2(r)*p; abs2(r)/p (abs2(t*trev) - abs2(r*rrev))/p]./abs2(t)
end

function rt_from_matrix(T)
    r = T[2,1]/T[1,1]
    t = 1.0/T[1,1]
    rrev = -T[1,2]/T[1,1]
    trev = (T[1,1]*T[2,2] - T[1,2]*T[2,1])/T[1,1]
    return r, t, rrev, trev
end

function tmm_matrix(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)
    local k0 = 2.0*pi/lambda;
    local q = norm(kp);
    Mg = Matrix{Complex{Float64}}(I, 2,2)
    Ml = zeros(Complex64, (2,2))
    L::Layer = Layer()
    Ifs::Interface = Interface()
    nL = length(S.Layers)
    for i = 1:(nL-1)
        L = S.Layers[i]
        Ifs = S.Interfs[i]
        r, t, rrev, trev = intface_rt(p, Ifs, k0, kp)
        di = betz(L.mat, k0, q)*L.d
        Ml = Matrix2x2_coh(r, t, rrev, trev, di)
        Mg  = Mg*Ml
    end
    r, t, rrev, trev = rt_from_matrix(Mg)
    return r, t
end

function tmm_RT(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)
    local k0 = 2.0*pi/lambda;
    local q = norm(kp);
    r, t = tmm_matrix(p, lambda, kp, S)
    mat1 = S.Layers[1].mat
    matN  = S.Layers[end].mat
    kz1 = betz(mat1, k0, q)
    kzN = betz(matN, k0, q)
    R = abs2(r)
    T = abs2(t)
    T = (mat1.mu/matN.mu)*real(kzN/kz1)*T
    return R, T
end

function tmm_ellipso(lambda::Real, kp::Vector{<:Real}, S::Stack)
    rs, ts =  tmm_matrix(1, lambda, kp, S)
    rp, tp =  tmm_matrix(2, lambda, kp, S)
    rr = rp/rs
    return abs(rr), cos(angle(rr))
end


function RT_matrix_inc(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)
    local k0 = 2.0*pi/lambda;
    local q = norm(kp);
    Mg = Matrix{Complex{Float64}}(I,2,2)
    Ml = zeros(Complex64, (2,2))
    L::Layer = Layer()
    Ifs::Interface = Interface()
    nL = length(S.Layers)
    for i = 1:(nL-2)
        L = S.Layers[i]
        Ifs = S.Interfs[i]
        r, t, rrev, trev = intface_rt(p, Ifs, k0, kp)
        di = betz(L.mat, k0, q)*L.d
        Ml = Matrix2x2_coh(r, t, rrev, trev, di)
        Mg  = Mg*Ml
    end
    # calculating substrate separatly
    L = S.Layers[nL-1]
    Ifs = S.Interfs[nL-1]
    r, t, rrev, trev = intface_rt(p, Ifs, k0, kp)
    di = betz(L.mat, k0, q)*L.d
    Ml = Matrix2x2_inc(r, t, rrev, trev, di)
    r, t, rrev, trev = rt_from_matrix(Mg)
    Mg = Matrix2x2_inc(r, t, rrev, trev, 0.0)*Ml
    R = Mg[2,1]/Mg[1,1]
    T = 1.0/Mg[1,1]
    return R, T
end
