using MultiTMM
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

# import the experimental data
DS90_I_Ag1 = Array{Float64}(undef, 1103-79, 5, nth)
DS90_I_Ag2 = Array{Float64}(undef, 1103-79, 5, nth)
DS90_I_Si = Array{Float64}(undef, 1103-79, 5, nth)

DS90_II_Ag1 = Array{Float64}(undef, 995, 5, nth)
DS90_II_Si = Array{Float64}(undef, 995, 5, nth)

DS90_III_Ag1 = Array{Float64}(undef, 995, 5, nth)
DS90_III_Ag2 = Array{Float64}(undef, 995, 5, nth)



path_S90 = "C:\\Users\\vmkhitaryan\\Documents\\Projects\\Ag Ribbons\\Chip-Andrew-S90\\Ellipsometry\\"
for i = 1:nth
    theta = round(Int, th[i])
    fI1 = readdlm(join([path_S90, "191117_Chip-Andrew-S90-16MLplus_an15_Ag_pos1_", string(theta),".pae"]))
    fI2 = readdlm(join([path_S90, "191117_Chip-Andrew-S90-16MLplus_an15_Ag_pos2_", string(theta),".pae"]))
    fI3 = readdlm(join([path_S90, "191117_Chip-Andrew-S90-16MLplus_an15_Si_", string(theta),".pae"]))

    fII1 = readdlm(join([path_S90, "2019-11-23-Chip-Andrew-S90-inc60to75_an15_MS_Ag_", string(theta),".pae"]))
    fII2 = readdlm(join([path_S90, "2019-11-23-Chip-Andrew-S90-inc60to75_an15_MS_Si_", string(theta),".pae"]))


    fIII1 = readdlm(join([path_S90, "2019-11-25-Chip-Andrew-S90-inc60to75_an15_MS_Ag_pos1_", string(theta),".pae"]))
    fIII2 = readdlm(join([path_S90, "2019-11-25-Chip-Andrew-S90-inc60to75_an15_MS_Ag_pos2_", string(theta),".pae"]))


    DS90_I_Ag1[:, :, i] = fI1[79:end-1, 1:5]
    DS90_I_Ag2[:, :, i] = fI2[79:end-1, 1:5]
    DS90_I_Si[:, :, i] = fI3[79:end-1, 1:5]

    DS90_II_Ag1[:, :, i] = fII1[79:end-1, 1:5]
    DS90_II_Si[:, :, i] = fII2[79:end-1, 1:5]

    DS90_III_Ag1[:, :, i] = fIII1[79:end-1, 1:5]
    DS90_III_Ag2[:, :, i] = fIII2[79:end-1, 1:5]
end

w1 = DS90_I_Ag1[:, 1, 1]
w2 = DS90_II_Ag1[:, 1, 1]
nw1 = length(w1)
nw2 = length(w2)


nAg_JC= nk_import("Ag(JC-eV)",  w1)
nAg_P= nk_import("Ag(Palik-eV)", w1)
nAg_B = nk_import("Ag(Barber-eV)", w1)
nAg_H = nk_import("Ag(Hageman-eV)", w1)
nAg_We = nk_import("Ag(Werner-eV)", w1)
nAg_S = nk_import("Ag(sopra-eV)", w1)
nSi = nk_import("Si111(sopra)", w1)
nSiO2 = nk_import("SiO2(sopra)", w1)

begin
    blist = [5, 10, 20, 50]
    TanPsi = Array{Float64}(undef, nw1, nth, length(blist))
    CosDel = Array{Float64}(undef, nw1, nth, length(blist))
    NML = 15
    tc = 2.0
    tML = 0.40853/sqrt(3)
    for (l, bl) = enumerate(blist)
        for j = 1:nth
            for i = 1:nw1
                λ = nm2eV/w1[i]
                k = 2.0*pi/λ
                LAir = Layer()
                LSiO2 = Layer("c", Material("nk", nSiO2[i]), tc)
                LSi_bot = Layer("c", Material("nk", nSi[i]))
                #LAg = Layer("c", Material("nk", nAg_S(ww[j])))

                epsAg = nAg_S[i]^2  - drude_pole(w1[i], 0, 9.17, 0.021) + drude_pole(w1[i], 0, 9.17, bl*0.021)
                LAg = Layer("c", Material("epsmu", epsAg), NML*tML)
                S = Stack([LAir; LSiO2; LAg; LSi_bot], zeros(3))
                TanPsi[i, j, l], CosDel[i, j, l] = tmm_ellipso(λ, [k*sin(pi*th[j]/180), 0], S)
            end
        end
    end
end

# 1st measurements
begin
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("YlGnBu")
    cmap3 = ColorMap("summer")
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 20
    bs = 1

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure("18 ML ellipsometric parameters", figsize=(12, 6))
    clf()
    subplot(1, 2, 1)
    for i = 1:nth
        thi = round(Int, th[i])
        plot(DS90_I_Ag1[:, 1, i], DS90_I_Ag1[:, 2, i], linestyle="-", linewidth=2, color = "blue")
        plot(DS90_I_Ag2[:, 1, i], DS90_I_Ag2[:, 2, i], linestyle="-", linewidth=2, color = "red")
        plot(DS90_I_Si[:, 1, i], DS90_I_Si[:, 2, i], linestyle="-", linewidth=2, color = "green")
        #=
        for (l, bl) = enumerate(blist)
             plot(ww, TanPsi[:, i, l], linestyle="--", color = cmap2((l)/float(length(blist))))
        end
        =#
        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="-", linewidth=2, color = "blue", label = "17/11/19 - S90/Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "red", label = "17/11/19 - S90/Ag pos.2")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "17/11/19 - S90/Si")
    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DS90_I_Ag1[:, 1, i], DS90_I_Ag1[:, 3, i], linestyle="-", linewidth=2, color = "blue")
        plot(DS90_I_Ag2[:, 1, i], DS90_I_Ag2[:, 3, i], linestyle="-", linewidth=2, color = "red")
        plot(DS90_I_Si[:, 1, i], DS90_I_Si[:, 3, i], linestyle="-", linewidth=2, color = "green")
        #=
        for (l, bl) = enumerate(blist)
            plot(ww, CosDel[:, i, l], linestyle ="--", color = cmap2((l)/float(length(blist))))
        end
        =#
        xlim([1.55, 5.5])
        ylim([-1, 0.25])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="-", linewidth=2, color = "blue", label = "17/11/19 - S90/Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "red", label = "17/11/19 - S90/Ag pos.2")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "17/11/19 - S90/Si")

    legend(frameon=false, fontsize = 14)
    tight_layout()
    suptitle("Ag 15ML +  1 nm SiO2", fontsize = 16)
    savefig(join([path_S90,"\\Pictures\\S90_mes_2019-11-17_exp_only.png"]), transparent = true)
end
# 1st measurement vs fit
begin
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("YlGnBu")
    cmap3 = ColorMap("summer")
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 20
    bs = 1

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure("18 ML ellipsometric parameters", figsize=(12, 6))
    clf()
    subplot(1, 2, 1)
    for i = 1:nth
        thi = round(Int, th[i])
        plot(DS90_I_Ag1[:, 1, i], DS90_I_Ag1[:, 2, i], linestyle="-", linewidth=2, color = "green")
        for (l, bl) = enumerate(blist)
             plot(w1, TanPsi[:, i, l], linestyle="--", color = cm((l)/float(length(blist))))
        end
        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "17/11/19 - S90/Ag pos.1")
    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DS90_I_Ag1[:, 1, i], DS90_I_Ag1[:, 3, i], linestyle="-", linewidth=2, color = "green")
        for (l, bl) = enumerate(blist)
            plot(w1, CosDel[:, i, l], linestyle ="--", color = cm((l)/float(length(blist))))
        end
        xlim([1.55, 5.5])
        ylim([-1, 0.25])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "17/11/19 - S90/Ag pos.1")
    legend(frameon=false, fontsize = 14)
    tight_layout()
    suptitle("Ag 15ML/1nm SiO2 fited 15ML/2nm SiO2", fontsize = 16)
    savefig(join([path_S90,"\\Pictures\\S90_mes_2019-11-17_exp_fit_15ML.png"]), transparent = true)
end
# 1st 2nd and 3th measurements
begin
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("YlGnBu")
    cmap3 = ColorMap("summer")
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 20
    bs = 1

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure("18 ML ellipsometric parameters", figsize=(12, 6))
    clf()
    subplot(1, 2, 1)
    for i = 1:nth
        thi = round(Int, th[i])
        plot(DS90_I_Ag1[:, 1, i], DS90_I_Ag1[:, 2, i], linestyle="-", linewidth=2, color = "blue")
        plot(DS90_II_Ag1[:, 1, i], DS90_II_Ag1[:, 2, i], linestyle="-", linewidth=2, color = "red")
        plot(DS90_III_Ag1[:, 1, i], DS90_III_Ag1[:, 2, i], linestyle="-", linewidth=2, color = "green")
        #=
        for (l, bl) = enumerate(blist)
             plot(ww, TanPsi[:, i, l], linestyle="--", color = cmap2((l)/float(length(blist))))
        end
        =#
        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="-", linewidth=2, color = "blue", label = "17/11/19 - S90/Ag")
    plot([],[], linestyle="-", linewidth=2, color = "red", label = "23/11/19 - S90/Ag")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "25/11/19 - S90/Ag")
    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DS90_I_Ag1[:, 1, i], DS90_I_Ag1[:, 3, i], linestyle="-", linewidth=2, color = "blue")
        plot(DS90_II_Ag1[:, 1, i], DS90_II_Ag1[:, 3, i], linestyle="-", linewidth=2, color = "red")
        plot(DS90_III_Ag1[:, 1, i], DS90_III_Ag1[:, 3, i], linestyle="-", linewidth=2, color = "green")
        #=
        for (l, bl) = enumerate(blist)
            plot(ww, CosDel[:, i, l], linestyle ="--", color = cmap2((l)/float(length(blist))))
        end
        =#
        xlim([1.55, 5.5])
        ylim([-1, 0.25])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="-", linewidth=2, color = "blue", label = "17/11/19 - S90/Ag")
    plot([],[], linestyle="-", linewidth=2, color = "red", label = "23/11/19 - S90/Ag")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "25/11/19 - S90/Ag")

    legend(frameon=false, fontsize = 14)
    tight_layout()
    suptitle("Ag 15ML +  1 nm SiO2", fontsize = 16)
    savefig(join([path_S90,"\\Pictures\\S90_all_exps.png"]), transparent = true)
end
