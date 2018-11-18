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
    file1 = readdlm(join([path_ell, "181030_Chip-Zaka-GaSiAg_2ML_10ML_an15_Ag_2ML_", string(theta),".pae"]))
    file2 = readdlm(join([path_ell, "181030_Chip-Zaka-GaSiAg_2ML_10ML_an15_Ag_10ML_", string(theta),".pae"]))
    file3 = readdlm(join([path_ell, "181030_Chip-Zaka-GaSiAg_2ML_10ML_an15_Si_", string(theta),".pae"]))
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
tGa = 0.3136
d = 2.0
bb = 1.0
tc = 1.0
tSiO2 = 1.0

for j = 1:nw
    λ = nm2eV/ww[j]
    k = 2.0*pi/λ
    LAir = Layer()
    LSiO2 = Layer("c", Material("nk", nSiO2(ww[j])), tSiO2)
    LSi_top = Layer("c", Material("nk", nSi(ww[j])), tc - tSiO2)
    LSi_bot = Layer("c", Material("nk", nSi(ww[j])))
    #LAg = Layer("c", Material("nk", nAg_S(ww[j])))

    epsAg = nAg_S(ww[j])^2  + drude_pole(ww[j], 0, wp, gam0) - drude_pole(ww[j], 0, wp, bb*gam0)
    LAg = Layer("c", Material("epsmu", epsAg), d*aa)
    S = Stack([LAir; LSiO2; LSi_top; LAg; LSi_bot], zeros(4))
    for i = 1:nth
        θ = th[i]*deg
        TanPsi[i, j], CosDel[i, j] = tmm_ellipso(λ, [k*sin(θ), 0], S)
        #Rs[i, j], Ts[i, j] = RT_calc(1, λ, [k*sin(θ), 0], S)
        #Rp[i, j], Tp[i, j] = RT_calc(2, λ, [k*sin(θ), 0], S)
    end
end

begin
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("winter")
    cmap3 = ColorMap("summer")
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure("$d ML ellipsometric parameters", figsize=(12, 5))
    for i = 1:nth
        thi = round(Int, th[i])
        subplot(121)
        title("2ML and 10ML ellipsometric parameters")
        plot(DATA1[i, :, 1], DATA1[i, :, 2], linestyle="-", color = cmap1((i+2)/float(2*nth)), label = "Ag 2ML, θ =$thi\$^o\$")
        #plot(DATA2[i, :, 1], DATA2[i, :, 2], linestyle="-", color = cmap2((i+2)/float(2*nth)), label = "Ag 10ML, θ =$thi\$^o\$")
        #plot(DATA3[i, :, 1], DATA3[i, :, 2], linestyle="--", color = "black", label = "Si, θ =$thi\$^o\$")

        plot(ww, TanPsi[i, :], linestyle="--", color = "red", label = "theory 2ML, θ =$thi\$^o\$")
        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Energy (eV)")
        ylabel(L"Tan(\Psi)")

        subplot(122)
        title("2ML and 10ML ellipsometric parameters")
        plot(DATA1[i, :, 1], DATA1[i, :, 3], linestyle="-", color = cmap1((i+2)/float(2*nth)), label = "Ag 2ML, θ =$thi\$^o\$")
        #plot(DATA2[i, :, 1], DATA2[i, :, 3], linestyle="-", color = cmap2((i+2)/float(2*nth)), label = "Ag 10ML, θ =$thi\$^o\$")
        #plot(DATA3[i, :, 1], DATA3[i, :, 3], linestyle="--", color = "black", label = "Si, θ =$thi\$^o\$ ")

        plot(ww, CosDel[i, :], linestyle="--", color = "red", label = "theory 2ML, θ =$thi\$^o\$")
        xlim([1.55, 5.5])
        ylim([-1, 0.25])
        xlabel("Energy (eV)")
        ylabel(L"Cos(\Delta)")
    end

    legend(bbox_to_anchor=[1.01, 1.0], loc=2, borderaxespad=0, fontsize = 8)
    #ax[:set_position]([0.03,0.03,0.71,1.01])
    tight_layout(rect=[0,0,0.95,0.99])
    savefig("D:\\Projects\\Ribbons\\Chip-Zaka-Laura-S66-Ga.Si\\Ellipsometry Pictures\\S66_2ML_exp_theory.png")
    clf()
end
