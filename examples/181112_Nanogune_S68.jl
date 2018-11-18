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
const nm2eV = 1239.84193
const nth = 4
th = linspace(60, 75, nth)


# import the experimental data
# import the experimental data
DATA1 = Array{Float64}(1103-79, 5, nth)
DATA2 = Array{Float64}(1103-79, 5, nth)
DATA3 = Array{Float64}(1103-79, 5, nth)

path_ell = "C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Nanogune\\Ellipsometry\\"

for i = 1:nth
    theta = round(Int, th[i])
    file1 = readdlm(join([path_ell, "181109_Chip-Zaka-Laura-S68-10ML_20ML_an15_Ag10ML_", string(theta),".pae"]))
    file2 = readdlm(join([path_ell, "181109_Chip-Zaka-Laura-S68-10ML_20ML_an15_Ag20ML_", string(theta),".pae"]))
    file3 = readdlm(join([path_ell, "181109_Chip-Zaka-Laura-S68-10ML_20ML_an15_Si_", string(theta),".pae"]))
    DATA1[:, :, i] = file1[79:end-1, 1:5]
    DATA2[:, :, i] = file2[79:end-1, 1:5]
    DATA3[:, :, i] = file3[79:end-1, 1:5]
end

ww = DATA1[1, :, 1]
ww = linspace(1, 6, 500)
nw = length(ww)

begin
    TanPsi = Array{Float64}(nw, nth)
    CosDel = Array{Float64}(nw, nth)
    for j = 1:nth
        for i = 1:nw
            λ = nm2eV/ww[i]
            k = 2.0*pi/λ
            LAir = Layer()
            LSiO2 = Layer("c", Material("nk", nSiO2(ww[i])), 1.0)
            LSi_bot = Layer("c", Material("nk", nSi(ww[i])))
            #LAg = Layer("c", Material("nk", nAg_S(ww[j])))

            epsAg = nAg_S(ww[i])^2  - drude_pole(ww[i], 0, 9.17, 0.022) + drude_pole(ww[i], 0, 9.17, 50*0.022)
            LAg = Layer("c", Material("epsmu", epsAg), 20*0.40853/sqrt(3))
            S = Stack([LAir; LSiO2; LAg; LSi_bot], zeros(3))
            TanPsi[i, j], CosDel[i, j] = tmm_ellipso(λ, [k*sin(pi*th[j]/180), 0], S)
        end
    end
end



CosDel_sim = Array{Float64}(undef, 500, 4, 5, 2)
TanPsi_sim = Array{Float64}(undef, 500, 4, 5, 2)
bsim = [0.5, 1.0, 5.0, 10.0, 50.0]
ws = linspace(1.55, 5.3, nw)
for (i, bi) = enumerate(bsim)
    CosDel_sim[:, :, i, 1] = readdlm(join(["C:\\Users\\vmkhi\\Desktop\\Ag10ML_CosDel_b$bi",".txt"]))
    TanPsi_sim[:, :, i, 1] = readdlm(join(["C:\\Users\\vmkhi\\Desktop\\Ag10ML_tanPsi_b$bi",".txt"]))
    CosDel_sim[:, :, i, 2] = readdlm(join(["C:\\Users\\vmkhi\\Desktop\\Ag20ML_CosDel_b$bi",".txt"]))
    TanPsi_sim[:, :, i, 2] = readdlm(join(["C:\\Users\\vmkhi\\Desktop\\Ag20ML_tanPsi_b$bi",".txt"]))
end



begin
    clf()
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("winter")
    cmap3 = ColorMap("summer")
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure("20 ML ellipsometric parameters", figsize=(12, 5))

    for i = 1:nth
        thi = round(Int, th[i])
        subplot(121)
        #plot(DATA1[:, 1, i], DATA1[:, 2, i], linestyle="-", color = cmap1((i+2)/float(2*nth)), label = "Ag 10ML, θ =$thi\$^o\$")
        plot(DATA2[:, 1, i], DATA2[:, 2, i], linestyle="-", color = cmap2((i+2)/float(2*nth)), label = "Ag 20ML, θ =$thi\$^o\$")
        #plot(DATA3[:, 1, i], DATA3[:, 2, i], linestyle="--", color = "black", label = "Si, θ =$thi\$^o\$")

        plot(ws, TanPsi_sim[:, i, 1:5, 2],  linestyle="--", color = "black")
        #plot(ww, TanPsi[:, i], linestyle="--", color = "black", label = "theory 20 ML, θ =$thi\$^o\$")
        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Energy (eV)")
        ylabel(L"Tan(\Psi)")

        subplot(122)
        #plot(DATA1[:, 1, i], DATA1[:, 3, i], linestyle="-", color = cmap1((i+2)/float(2*nth)), label = "Ag 10 ML, θ =$thi\$^o\$")
        plot(DATA2[:, 1, i], DATA2[:, 3, i], linestyle="-", color = cmap2((i+2)/float(2*nth)), label = "Ag 20ML, θ =$thi\$^o\$")
        #plot(DATA3[:, 1, i], DATA3[:, 3, i], linestyle="--", color = "black", label = "Si, θ =$thi\$^o\$ ")

        #plot(ww, CosDel[:, i], linestyle ="--", color = "black", label = "theory 20 ML, θ =$thi\$^o\$")
        plot(ws, CosDel_sim[:, i, 1:5, 2],  linestyle="--", color = "black")
        xlim([1.55, 5.5])
        ylim([-1, 0.25])
        xlabel("Energy (eV)")
        ylabel(L"Cos(\Delta)")
    end
    #suptitle("20ML with \$ \\gamma \$ = 1.0 \$ \\gamma_{JC}\$", fontsize = 16)
    #legend(bbox_to_anchor=(1.05, 1), loc=2)
    #tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Chip-Zaka-Laura-S68\\Ellipsometry Pictures\\S68_20ML_exp_theory_SiO2_1nm.png")
end
