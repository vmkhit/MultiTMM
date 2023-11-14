include("../src/MultiTMM.jl")
include("../src/tmm_matrix.jl")
include("../src/eps.jl")
using DelimitedFiles
using MAT
using PyPlot
pygui(true)
matplotlib["rcParams"][:update](["font.size" => 18, "font.weight"=>"normal", "font.family" => "Arial", "text.usetex"=>true])


function drude_pole(w, epsb, wp, gam)
    return epsb - wp^2/(w*(w+im*gam))
end

#parameters and w-theta ranges
nm2eV = 1239.84193
nth = 4
th = range(60, 75, length = nth)

path_Saad = "C:/Users/vmkhi/OneDrive - ICFO/NPTExperiments/VIS/S122/Post customs issues/Ribbons/R5/Ref on Ag mirror/"
DSi_N = readdlm(join([path_Saad, "BARE_Si_N.txt"]), skipstart=1)
DSi_V = readdlm(join([path_Saad, "BARE_Si_V.txt"]), skipstart=1)

begin
    plot(DSi_V[:, 1], DSi_V[:, 2])
    plot(DSi_N[:, 1], DSi_N[:, 2])
    xlim([450, 1800])
    ylim([0, 100])
end

ww = range(1,  2, length = 1000)
#ww = range(1, 6, length = 500)
nw = length(ww)

nAg_JC= nk_import("Ag(JC-eV)",  ww)
nAg_P= nk_import("Ag(Palik-eV)", ww)
nAg_B = nk_import("Ag(Barber-eV)", ww)
nAg_H = nk_import("Ag(Hageman-eV)", ww)
nAg_We = nk_import("Ag(Werner-eV)", ww)
nAg_S = nk_import("Ag(sopra-eV)", ww)
nSi = nk_import("Si111(sopra)", ww)
nSiO2 = nk_import("SiO2(sopra)", ww)
tML = 0.40853/sqrt(3)
tc = 2

Rs = Array{Float64}(undef, nw, nth)
Ts = Array{Float64}(undef, nw, nth)
Rp = Array{Float64}(undef, nw, nth)
Tp = Array{Float64}(undef, nw, nth)

for j = 1:nth
    for i = 1:nw
        λ = nm2eV/ww[i]
        k = 2.0*pi/λ
        LAir = Layer()
        LSiO2 = Layer("c", Material("nk", nSiO2[i]), tc)
        LSi_bot = Layer("c", Material("nk", nSi[i]))
        #LAg = Layer("c", Material("nk", nAg_S(ww[j])))
        #epsAg = nAg_S[i]^2  - drude_pole(ww[i], 0, 9.17, 0.021) + drude_pole(ww[i], 0, 9.17, bl*0.021)
        #LAg = Layer("c", Material("epsmu", epsAg), NML*tML)
        S = Stack([LAir; LSiO2; LSi_bot], zeros(2))
        Rs[i, j], Ts[i, j] = tmm_RT(1, λ, [k * sin(pi * th[j] / 180), 0], S)
        Rp[i, j], Tp[i, j] = tmm_RT(2, λ, [k * sin(pi * th[j] / 180), 0], S)
    end
end



plot(nm2eV./ww, Rs) 