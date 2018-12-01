using MultiTMM
using DelimitedFiles
using PyPlot
pygui(true)

function drude_pole(w, epsb, wp, gam)
    return epsb - wp^2/(w*(w+im*gam))
end

#parameters and w-theta ranges
nm2eV = 1239.84193
nth = 4
th = range(60, 75, length = nth)


# import the experimental data
DATA1 = Array{Float64}(undef, 1103-79, 5, nth)
DATA2 = Array{Float64}(undef, 1103-79, 5, nth)
DATA3 = Array{Float64}(undef, 1103-79, 5, nth)
DATA4 = Array{Float64}(undef, 1103-79, 5, nth)


path_ell = "C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Nanogune\\Ellipsometry\\"

for i = 1:nth
    theta = round(Int, th[i])
    file1 = readdlm(join([path_ell, "181129_Chip-Zaka-Laura-S78-8-12ML_an15_Ag_pos1_", string(theta),".pae"]))
    file2 = readdlm(join([path_ell, "181129_Chip-Zaka-Laura-S78-8-12ML_an15_Ag_pos2_", string(theta),".pae"]))
    file3 = readdlm(join([path_ell, "181129_Chip-Zaka-Laura-S78-8-12ML_an15_Ag_pos3_", string(theta),".pae"]))
    file4 = readdlm(join([path_ell, "181129_Chip-Zaka-Laura-S78-8-12ML_an15_Si_pos1_", string(theta),".pae"]))

    DATA1[:, :, i] = file1[79:end-1, 1:5]
    DATA2[:, :, i] = file2[79:end-1, 1:5]
    DATA3[:, :, i] = file3[79:end-1, 1:5]
    DATA4[:, :, i] = file4[79:end-1, 1:5]

end

ww = DATA1[1, :, 1]
ww = range(1, 6, length = 500)
nw = length(ww)


nAg_JC= nk_import("Ag(JC-eV)",  ww)
nAg_P= nk_import("Ag(Palik-eV)", ww)
nAg_B = nk_import("Ag(Barber-eV)", ww)
nAg_H = nk_import("Ag(Hageman-eV)", ww)
nAg_We = nk_import("Ag(Werner-eV)", ww)
nAg_S = nk_import("Ag(sopra-eV)", ww)
nSi = nk_import("Si111(sopra)", ww)
nSiO2 = nk_import("SiO2(sopra)", ww)

begin
    blist = [0.5, 1, 2, 5, 10, 20, 50]
    TanPsi = Array{Float64}(undef, nw, nth, length(blist))
    CosDel = Array{Float64}(undef, nw, nth, length(blist))
    NML = 12
    tc = 1.0
    for (l, bl) = enumerate(blist)
        for j = 1:nth
            for i = 1:nw
                λ = nm2eV/ww[i]
                k = 2.0*pi/λ
                LAir = Layer()
                LSiO2 = Layer("c", Material("nk", nSiO2[i]), tc)
                LSi_bot = Layer("c", Material("nk", nSi[i]))
                #LAg = Layer("c", Material("nk", nAg_S(ww[j])))

                epsAg = nAg_S[i]^2  - drude_pole(ww[i], 0, 9.17, 0.021) + drude_pole(ww[i], 0, 9.17, bl*0.021)
                LAg = Layer("c", Material("epsmu", epsAg), NML*0.40853/sqrt(3))
                S = Stack([LAir; LSiO2; LAg; LSi_bot], zeros(3))
                TanPsi[i, j, l], CosDel[i, j, l] = tmm_ellipso(λ, [k*sin(pi*th[j]/180), 0], S)
            end
        end
    end
end


begin
    clf()
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("winter")
    cmap3 = ColorMap("summer")
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure("20 ML ellipsometric parameters", figsize=(14, 6))
    for i = 1:2:nth
        thi = round(Int, th[i])
        subplot(1, 2, 1)
        plot(DATA1[:, 1, i], DATA1[:, 2, i], linestyle="-", color = cmap1((i+2)/float(2*nth)), label = "Ag 12ML, θ =$thi\$^o\$")
        #plot(DATA2[:, 1, i], DATA2[:, 2, i], linestyle="-", color = cmap2((i+2)/float(2*nth)), label = "Ag 10ML, θ =$thi\$^o\$")
        #plot(DATA3[:, 1, i], DATA3[:, 2, i], linestyle="-", color = cmap3((i+2)/float(2*nth)), label = "Ag 8ML, θ =$thi\$^o\$")
        #plot(DATA4[:, 1, i], DATA4[:, 2, i], linestyle="--", color = "black", label = "Si, θ =$thi\$^o\$")

        for (l, bl) = enumerate(blist)
            plot(ww, TanPsi[:, i, l], linestyle="--", color = cmap2((l)/float(length(blist))), label = label = "γ/γ\$_{JC}\$ = $bl")
        end
        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Energy (eV)")
        ylabel(L"Tan(\Psi)")
        title("Ag 12ML + 1nm SiO2", fontsize = 16)

        subplot(1, 2, 2)
        #color = cmap2((i+2)/float(2*nth))
        plot(DATA1[:, 1, i], DATA1[:, 3, i], linestyle="-", color = cmap1((i+2)/float(2*nth)), label = "Ag 12ML, θ =$thi\$^o\$")
        #plot(DATA2[:, 1, i], DATA2[:, 3, i], linestyle="-", color = cmap2((i+2)/float(2*nth)), label = "Ag 10ML, θ =$thi\$^o\$")
        #plot(DATA3[:, 1, i], DATA3[:, 3, i], linestyle="-", color = cmap3((i+2)/float(2*nth)), label = "Ag 8ML, θ =$thi\$^o\$")
        #plot(DATA4[:, 1, i], DATA4[:, 3, i], linestyle="--", color = "black", label = "Si, θ =$thi\$^o\$ ")
        for (l, bl) = enumerate(blist)
            plot(ww, CosDel[:, i, l], linestyle ="--", color = cmap2((l)/float(length(blist))), label = "γ/γ\$_{JC}\$ = $bl")
        end
        xlim([1.55, 5.5])
        ylim([-1, 0.25])
        xlabel("Energy (eV)")
        ylabel(L"Cos(\Delta)")
        title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end
    #suptitle("Ag 12ML + 1nm SiO2", fontsize = 16)
    legend(bbox_to_anchor=(0.9999, 1), loc=2)
    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Chip-Zaka-Laura-S78\\Pictures\\S78_exp_fit_12ML.png")
end
