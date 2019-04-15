using MultiTMM
using PyPlot
pygui(true)
matplotlib["rcParams"][:update](["font.size" => 18, "font.weight"=>"normal", "font.family" => "Arial", "text.usetex"=>true])

nm2eV = 1239.84193
nth = 50
th = range(25, 75, length = nth)

nw = 500
ww = range(1.5, 5, length = nw)

nAg_JC= nk_import("Ag(JC-eV)",  ww)
nAg_P= nk_import("Ag(Palik-eV)", ww)
nAg_B = nk_import("Ag(Barber-eV)", ww)
nAg_H = nk_import("Ag(Hageman-eV)", ww)
nAg_We = nk_import("Ag(Werner-eV)", ww)
nAg_S = nk_import("Ag(sopra-eV)", ww)
nSi = nk_import("Si111(sopra)", ww)
nSiO2 = nk_import("SiO2(sopra)", ww)


function drude_pole(w, epsb, wp, gam)
    return epsb - wp^2/(w*(w+im*gam))
end


begin
    blist = [0.5, 1, 2, 5, 10, 20, 50]
    TanPsi = Array{Float64}(undef, nw, nth)
    CosDel = Array{Float64}(undef, nw, nth)
    NML = 12
    tc = 2.0
    tML = 0.40853/sqrt(3)
    for j = 1:nth
        for i = 1:nw
            λ = nm2eV/ww[i]
            k = 2.0*pi/λ
            LAir = Layer()
            LSiO2 = Layer("c", Material("nk", nSiO2[i]), tc)
            LSi_bot = Layer("c", Material("nk", nSi[i]))
            #LAg = Layer("c", Material("nk", nAg_S(ww[j])))
            epsAg = nAg_S[i]^2  # - drude_pole(ww[i], 0, 9.17, 0.021) + drude_pole(ww[i], 0, 9.17, 0.021)
            LAg = Layer("c", Material("epsmu", epsAg), NML*tML)
            S = Stack([LAir; LSiO2; LAg; LSi_bot], zeros(3))
            TanPsi[i, j], CosDel[i, j] = tmm_ellipso(λ, [k*sin(pi*th[j]/180), 0], S)
        end
    end
end

plot(ww, TanPsi[:, 50])
