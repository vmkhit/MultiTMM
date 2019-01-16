function kzfn(ee, k, q)
    return sqrt(ee*k^2 - q^2 + 0.0im)
end

function interface_rt(p::Integer, lam, e1, e2, q)
    @assert p == 1 || p == 2
    k0 = 2*pi/lam
    kz1 = kzfn(e1, k0, q)
    kz2 = kzfn(e1, k0, q)
    if p == 1 # s-polarized
        r = (kz1 - kz2)/(kz1 + kz2)
        t = 2*kz1/(kz1 + kz2)
    else
        r =  (e2*kz1 - e1*kz2)/(e2*kz1 + e1*kz2)
        t = 2*e2*kz1/(e2*kz1 + e1*kz2)
    end
    return r, t
 end

function rt_layer(p::Integer, lam, ee, d, q)
    #ee = [e1, e2, e3]
    k0 = 2*pi/lam
    kz2 = kzfn(ee[2], k0, q)
    r12, t12 = interface_rt(p, lam, ee[1], ee[2], q)
    r23, t23 = interface_rt(p, lam, ee[2], ee[3], q)
    phi = cis(-kz2*d)
    r = (r12 + r23*phi^2)/(1 - r12*r23*phi^2)
    t = t12*t23*phi/(1 - r12*r23*phi^2)
    return r, t
end

function layerellipso(lam, ee, d, q)
    rs = rt_layer(1, lam, ee, d, q)
    rp = rt_layer(2, lam, ee, d, q)
    rr = rp/rs
    return abs(rr), cos(angle(rr))
end
