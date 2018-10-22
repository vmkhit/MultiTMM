using MultiTMM
using PyPlot

nAg_JC(x)= nk_import("Ag(JC-eV)", x)
nAg_P(x)= nk_import("Ag(Palik-eV)", x)
nAg_B(x) = nk_import("Ag(Barber-eV)", x)
nAg_H(x) = nk_import("Ag(Hageman-eV)", x)
nAg_Wi(x) = nk_import("Ag(film)", x)
nAg_We(x) = nk_import("Ag(Werner-eV)", x)
nAg_S(x) = nk_import("Ag(sopra-eV)", x)
nSi(x) = nk_import("Si111(sopra)", x)
nSiO2(x) = nk_import("SiO2(sopra)", x)
nInAs(x) = nk_import("InAs(film)", x)
nAl(x) = nk_import("Al(film)", x)

function drude_pole(w, epsb, wp, gam)
    return epsb - wp^2/(w*(w+im*gam))
end

#parameters and w-theta ranges
nm2eV = 1239.84193
deg = pi/180

nth = 4
th = linspace(60, 75, nth)


# import the experimental data
# import the experimental data
DATA1 = Array{Float64}(nth, 1103-79, 5)
DATA2 = Array{Float64}(nth, 1103-79, 5)
DATA3 = Array{Float64}(nth, 1103-79, 5)

path_ell = "C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Nanogune\\Ellipsometry\\"
for i = 1:nth
    theta = round(Int, th[i])
    file1 = readdlm(join([path_ell, "181015_Chip-Zaka-S58_Ag_10ML-15ML-20ML_an15_pos1_", string(theta),".pae"]))
    file2 = readdlm(join([path_ell, "181015_Chip-Zaka-S58_Ag_10ML-15ML-20ML_an15_pos2_", string(theta),".pae"]))
    file3 = readdlm(join([path_ell, "181015_Chip-Zaka-S58_Ag_10ML-15ML-20ML_an15_pos3_", string(theta),".pae"]))
    DATA1[i, :, :] = file1[79:end-1, 1:5]
    DATA2[i, :, :] = file2[79:end-1, 1:5]
    DATA3[i, :, :] = file3[79:end-1, 1:5]
end

ww = DATA1[1, :, 1]
ww = linspace(1, 6, 500)
nw = length(ww)



# Drude parameters
wp = 9.17
gam0 = 0.022

nb = 100
blist = linspace(0.05, 20, nb)



# start the calculation
TanPsi = Array{Float64}(nth, nw)
CosDel = Array{Float64}(nth, nw)

abulk = 0.40853
aa = abulk/sqrt(3)
d = 10
bb = 1.0
tc = 2
tSiO2 = 2

for j = 1:nw
    λ = nm2eV/ww[j]
    k = 2.0*pi/λ
    epsAg = nAg_JC(ww[j])^2 - drude_pole(ww[j], 0, wp, gam0) + drude_pole(ww[j], 0, wp, bb*gam0)
    LAir = Layer()
    LSiO2 = Layer("c", Material("nk", nSiO2(ww[j])), tSiO2)
    LSi_top = Layer("c", Material("nk", nSi(ww[j])), tc - tSiO2)
    LSi_bot = Layer("c", Material("nk", nSi(ww[j])))

    #epsAg = nAg(ww[j])^2  - drude_pole(ww[j], 0, wp, gam0) + drude_pole(ww[j], 0, wp, bb*gam0)
    LAg = Layer("c", Material("epsmu", epsAg, d*aa))
    S = Stack([LAir; LSiO2; LSi_top; LAg; LSi_bot], zeros(4))
    for i = 1:nth
        θ = th[i]*deg
        TanPsi[i, j], CosDel[i, j] = tmm_ellipso(λ, [k*sin(θ), 0], S)
        #Rs[i, j], Ts[i, j] = RT_calc(1, λ, [k*sin(θ), 0], S)
        #Rp[i, j], Tp[i, j] = RT_calc(2, λ, [k*sin(θ), 0], S)
    end
end




begin
    cmap1 = ColorMap("brg")
    cmap2 = ColorMap("jet")
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure("$d ML ellipsometric parameters", figsize=(12, 5))
    for i = 1:nth
        subplot(121)
        title("$d ML ellipsometric parameters")
        plot(DATA1[i, :, 1], DATA1[i, :, 2], linestyle="--", color = "black")
        plot(ww, TanPsi[i, :], linestyle="-", label = labels[i], color = cmap2(i/float(nth)))
        xlim([1.55, 5])
        ylim([0, 1.0])
        xlabel("Energy (eV)")
        ylabel(L"Tan(\Psi)")

        subplot(122)
        title("$d ML ellipsometric parameters")
        plot(DATA1[i, :, 1], DATA1[i, :, 3], linestyle="--", color = "black")
        plot(ww, CosDel[i, :], linestyle="-", label = labels[i], color = cmap2(i/float(nth)))
        xlim([1.55, 5])
        ylim([-1, 0.25])
        xlabel("Energy (eV)")
        ylabel(L"Cos(\Delta)")
    end
    legend()
    #savefig("C:\\Users\\vmkhi\\Desktop\\Ag_10ML\\8ML-16ML_Ellipsometric\\TanPsiCosDel_Sopra_$d ML_02.png")
end
