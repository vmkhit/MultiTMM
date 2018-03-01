using MultiTMM
using HDF5, JLD
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

nth = 8
th = linspace(40, 75, nth)


# import the experimental data
DATA = Array{Float64}(nth, 995, 5)
for i = 1:nth
    theta = round(Int, th[i])
    file = readdlm("C:\\Users\\vmkhi\\Desktop\\Ag_10ML\\171214_Chip-Zaka-1710-2\\20180215_chip-Zaka-1813-1-8ML_16ML_Plane_CDD_inc40to75_an45_MS_pos_8ML_Y119.5_"string(theta)".pae")
    DATA[i, :, :] = file[79:end-1, 1:5]
end
#=
begin
    for i = 1:nth
        plot(DATA[i, :, 1],  DATA[i, :, 2], linestyle = "-", label = "new", color = "black")
        xlim([1.55, 5])
        ylim([0, 1.0])
        xlabel("Energy (eV)")
        ylabel(L"Tan(\Psi)")
    end
    legend()
end
=#

ww = DATA[1, :, 1]
ww = linspace(1, 6, 500)
nw = length(ww)



# Drude parameters
wp = 9.01
gam0 = 0.022

nb = 100
blist = linspace(0.05, 20, nb)



# start the calculation
TanPsi = Array{Float64}(nth, nw)
CosDel = Array{Float64}(nth, nw)
#Rp = Array{Float64}(nth, nw)
#Tp = Array{Float64}(nth, nw)
#Rs = Array{Float64}(nth, nw)
#Ts = Array{Float64}(nth, nw)

abulk = 0.357
aa = abulk/sqrt(3)
d = 0.0
for j = 1:nw
    λ = nm2eV/ww[j]
    k = 2.0*pi/λ
    LAir = Layer()
    LSiO2 = Layer("c", Material("nk", nSiO2(ww[j])), 2.0)
    LSi = Layer("c", Material("nk", nSi(ww[j])))
    #epsAg = nAg(ww[j])^2  - drude_pole(ww[j], 0, wp, gam0) + drude_pole(ww[j], 0, wp, bb*gam0)
    LAg = Layer("c", Material("nk", nAg_S(ww[j])), d*aa)
    S = Stack([LAir; LSiO2; LAg; LSi], zeros(3))
    for i = 1:nth
        θ = th[i]*deg
        TanPsi[i, j], CosDel[i, j] = tmm_ellipso(λ, [k*sin(θ), 0], S)
        #Rs[i, j], Ts[i, j] = RT_calc(1, λ, [k*sin(θ), 0], S)
        #Rp[i, j], Tp[i, j] = RT_calc(2, λ, [k*sin(θ), 0], S)
    end
end


cmap1 = ColorMap("brg")
cmap2 = ColorMap("jet")
labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
begin
    fig = figure("$d ML ellipsometric parameters", figsize=(12, 5))
    for i = 1:nth
        subplot(121)
        title("$d ML ellipsometric parameters")
        plot(DATA[i, :, 1], DATA[i, :, 2], linestyle="--", color = "black")
        plot(ww, TanPsi[i, :], linestyle="-", label = labels[i], color = cmap2(i/float(nth)))
        xlim([1.55, 5])
        ylim([0, 1.0])
        xlabel("Energy (eV)")
        ylabel(L"Tan(\Psi)")

        subplot(122)
        title("$d ML ellipsometric parameters")
        plot(DATA[i, :, 1], DATA[i, :, 3], linestyle="--", color = "black")
        plot(ww, CosDel[i, :], linestyle="-", label = labels[i], color = cmap2(i/float(nth)))
        xlim([1.55, 5])
        ylim([-1, 0.25])
        xlabel("Energy (eV)")
        ylabel(L"Cos(\Delta)")
    end
    legend()
    #savefig("C:\\Users\\vmkhi\\Desktop\\Ag_10ML\\8ML-16ML_Ellipsometric\\TanPsiCosDel_Sopra_$d ML_02.png")
end
