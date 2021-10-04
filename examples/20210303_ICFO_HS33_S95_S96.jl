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
path0 = "C:/Users/vmkhitaryan/Documents/Projects/Ag Ribbons/Chip-Andrew-HS33/Ellipsometry/"

#HS30
DHS30_Si = Array{Float64}(undef, 1103-79, 5, nth)
DHS30_Ag_pos1 = Array{Float64}(undef, 1103-79, 5, nth)
DHS30_Ag_pos2 = Array{Float64}(undef, 1103-79, 5, nth)
# HS33
DHS33_Si = Array{Float64}(undef, 1103-79, 5, nth)
DHS33_Ag_pos1 = Array{Float64}(undef, 1103-79, 5, nth)
DHS33_Ag_pos2 = Array{Float64}(undef, 1103-79, 5, nth)
DHS33_Ag_pos3 = Array{Float64}(undef, 1103-79, 5, nth)

# S95
DS95_Si = Array{Float64}(undef, 1103-79, 5, nth)
DS95_Ag_pos1 = Array{Float64}(undef, 1103-79, 5, nth)
DS95_Ag_pos2 = Array{Float64}(undef, 1103-79, 5, nth)

# S96
DS96_Si = Array{Float64}(undef, 1103-79, 5, nth)
DS96_Ag20ML_pos1 = Array{Float64}(undef, 1103-79, 5, nth)
DS96_Ag20ML_pos2 = Array{Float64}(undef, 1103-79, 5, nth)
DS96_Ag12ML_pos1 = Array{Float64}(undef, 1103-79, 5, nth)
DS96_Ag12ML_pos2 = Array{Float64}(undef, 1103-79, 5, nth)

for i = 1:nth
    theta = round(Int, th[i])

    # HS30
    fHS30_Si = readdlm(join([path0, "2021-03-2-Chip-Andrew-HS30-an15_MS_Si_", string(theta),".pae"]))
    fHS30_Ag_pos1 = readdlm(join([path0, "2021-03-2-Chip-Andrew-HS30-an15_MS_Ag_pos1_", string(theta),".pae"]))
    fHS30_Ag_pos2 = readdlm(join([path0, "2021-03-2-Chip-Andrew-HS30-an15_MS_Ag_pos2_", string(theta),".pae"]))

    # HS33
    fHS33_Si = readdlm(join([path0, "2021-03-2-Chip-Andrew-HS33-an15_MS_Si_", string(theta),".pae"]))
    fHS33_Ag_pos1 = readdlm(join([path0, "2021-03-2-Chip-Andrew-HS33-an15_MS_Ag_pos1_", string(theta),".pae"]))
    fHS33_Ag_pos2 = readdlm(join([path0, "2021-03-2-Chip-Andrew-HS33-an25_MS_Ag_pos2_", string(theta),".pae"]))
    fHS33_Ag_pos3 = readdlm(join([path0, "2021-03-2-Chip-Andrew-HS33-an25_MS_Ag_pos3_", string(theta),".pae"]))

    # S95
    fS95_Si = readdlm(join([path0, "2021-03-2-Chip-Andrew-S95-16ML-an25_MS_Si_", string(theta),".pae"]))
    fS95_Ag_pos1 = readdlm(join([path0, "2021-03-2-Chip-Andrew-S95-16ML-an25_MS_Ag_pos1_", string(theta),".pae"]))
    fS95_Ag_pos2 = readdlm(join([path0, "2021-03-2-Chip-Andrew-S95-16ML-an25_MS_Ag_pos2_", string(theta),".pae"]))

    # S96
    fS96_Si = readdlm(join([path0, "2021-03-2-Chip-Andrew-S96-20ML-12ML-an25_MS_Si_", string(theta),".pae"]))
    fS96_Ag20ML_pos1 = readdlm(join([path0, "2021-03-2-Chip-Andrew-S96-20ML-12ML-an25_MS_Ag20ML_pos1_", string(theta),".pae"]))
    fS96_Ag20ML_pos2 = readdlm(join([path0, "2021-03-2-Chip-Andrew-S96-20ML-12ML-an25_MS_Ag20ML_pos2_", string(theta),".pae"]))
    fS96_Ag12ML_pos1 = readdlm(join([path0, "2021-03-2-Chip-Andrew-S96-20ML-12ML-an25_MS_Ag12ML_pos1_", string(theta),".pae"]))
    fS96_Ag12ML_pos2 = readdlm(join([path0, "2021-03-2-Chip-Andrew-S96-20ML-12ML-an25_MS_Ag12ML_pos2_", string(theta),".pae"]))


    DHS30_Si[:, :, i] = fHS30_Si[79:end-1, 1:5]
    DHS30_Ag_pos1[:, :, i] = fHS30_Ag_pos1[79:end-1, 1:5]
    DHS30_Ag_pos2[:, :, i] = fHS30_Ag_pos2[79:end-1, 1:5]

    DHS33_Si[:, :, i] = fHS30_Si[79:end-1, 1:5]
    DHS33_Ag_pos1[:, :, i] = fHS33_Ag_pos1[79:end-1, 1:5]
    DHS33_Ag_pos2[:, :, i] = fHS33_Ag_pos2[79:end-1, 1:5]
    DHS33_Ag_pos3[:, :, i] = fHS33_Ag_pos3[79:end-1, 1:5]

    DS95_Si[:, :, i] = fS95_Si[79:end-1, 1:5]
    DS95_Ag_pos1[:, :, i] = fS95_Ag_pos1[79:end-1, 1:5]
    DS95_Ag_pos2[:, :, i] = fS95_Ag_pos2[79:end-1, 1:5]

    DS96_Si[:, :, i] = fS95_Si[79:end-1, 1:5]
    DS96_Ag20ML_pos1[:, :, i] = fS96_Ag20ML_pos1[79:end-1, 1:5]
    DS96_Ag20ML_pos2[:, :, i] = fS96_Ag20ML_pos2[79:end-1, 1:5]
    DS96_Ag12ML_pos1[:, :, i] = fS96_Ag12ML_pos1[79:end-1, 1:5]
    DS96_Ag12ML_pos2[:, :, i] = fS96_Ag12ML_pos2[79:end-1, 1:5]
end

ww = DHS30_Si[:, 1, 1]
nw = length(ww)

#================================ Plot HS30 ===================================#
begin
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("YlGnBu")
    cmap3 = ColorMap("summer")
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 20
    bs = 2

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure("18 ML ellipsometric parameters", figsize=(12, 6))
    clf()
    subplot(1, 2, 1)
    for i = 1:nth
        thi = round(Int, th[i])
        plot(DHS30_Si[:, 1, i], DHS30_Si[:, 2, i], linestyle="--", linewidth=2, color = "blue")
        plot(DHS30_Ag_pos1[:, 1, i], DHS30_Ag_pos1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        plot(DHS30_Ag_pos2[:, 1, i], DHS30_Ag_pos2[:, 2, i], linestyle="-", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos. 1")
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag pos.2")

    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DHS30_Si[:, 1, i], DHS30_Si[:, 3, i], linestyle="--", linewidth=2, color = "blue")
        plot(DHS30_Ag_pos1[:, 1, i], DHS30_Ag_pos1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        plot(DHS30_Ag_pos2[:, 1, i], DHS30_Ag_pos2[:, 3, i], linestyle="-", linewidth=2, color = "orange")
        xlim([1.55, 5.5])
        ylim([-1, 1])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos. 1")
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag pos.2")
    legend(frameon=false, fontsize = 14)
    tight_layout()
    suptitle("Sample HS30 Ag 18ML + 2 nm SiO2", fontsize = 16)
    savefig(join([path0,"\\Pictures\\HS30_exp.png"]), transparent = true)
end

#================================ Plot HS33 ===================================#
begin
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("YlGnBu")
    cmap3 = ColorMap("summer")
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 20
    bs = 2

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure("18 ML ellipsometric parameters", figsize=(12, 6))
    clf()
    subplot(1, 2, 1)
    for i = 1:nth
        thi = round(Int, th[i])
        plot(DHS33_Si[:, 1, i], DHS33_Si[:, 2, i], linestyle="--", linewidth=2, color = "blue")
        plot(DHS33_Ag_pos1[:, 1, i], DHS33_Ag_pos1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        plot(DHS33_Ag_pos2[:, 1, i], DHS33_Ag_pos2[:, 2, i], linestyle="-", linewidth=2, color = "orange")
        plot(DHS33_Ag_pos3[:, 1, i], DHS33_Ag_pos3[:, 2, i], linestyle="-", linewidth=2, color = "green")

        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag pos.2")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.3")


    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DHS33_Si[:, 1, i], DHS33_Si[:, 3, i], linestyle="--", linewidth=2, color = "blue")
        plot(DHS33_Ag_pos1[:, 1, i], DHS33_Ag_pos1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        plot(DHS33_Ag_pos2[:, 1, i], DHS33_Ag_pos2[:, 3, i], linestyle="-", linewidth=2, color = "orange")
        plot(DHS33_Ag_pos3[:, 1, i], DHS33_Ag_pos3[:, 3, i], linestyle="-", linewidth=2, color = "green")

        xlim([1.55, 5.5])
        ylim([-1, 1])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag pos.2")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.3")

    legend(frameon=false, fontsize = 14)
    tight_layout()
    suptitle("Sample HS33 Ag 16ML + 1 nm SiO2", fontsize = 16)
    savefig(join([path0,"\\Pictures\\HS33_exp.png"]), transparent = true)
end

#================================ Plot S95 ===================================#
begin
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("YlGnBu")
    cmap3 = ColorMap("summer")
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 20
    bs = 2

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure("18 ML ellipsometric parameters", figsize=(12, 6))
    clf()
    subplot(1, 2, 1)
    for i = 1:nth
        thi = round(Int, th[i])
        plot(DS95_Si[:, 1, i], DS95_Si[:, 2, i], linestyle="--", linewidth=2, color = "blue")
        plot(DS95_Ag_pos1[:, 1, i], DS95_Ag_pos1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        plot(DS95_Ag_pos2[:, 1, i], DS95_Ag_pos2[:, 2, i], linestyle="-", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag pos.2")


    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DS95_Si[:, 1, i], DS95_Si[:, 3, i], linestyle="--", linewidth=2, color = "blue")
        plot(DS95_Ag_pos1[:, 1, i], DS95_Ag_pos1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        plot(DS95_Ag_pos2[:, 1, i], DS95_Ag_pos2[:, 3, i], linestyle="-", linewidth=2, color = "orange")


        xlim([1.55, 5.5])
        ylim([-1, 1])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag pos.2")


    legend(frameon=false, fontsize = 14)
    tight_layout()
    suptitle("Sample S95 Flat Ag 16ML + 1 nm SiO2", fontsize = 16)
    savefig(join([path0,"\\Pictures\\S95_exp.png"]), transparent = true)
end

#================================ Plot S96 ===================================#
begin
    cmap1 = ColorMap("autumn")
    cmap2 = ColorMap("YlGnBu")
    cmap3 = ColorMap("summer")
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 20
    bs = 2

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure("18 ML ellipsometric parameters", figsize=(12, 6))
    clf()
    subplot(1, 2, 1)
    for i = 1:nth
        thi = round(Int, th[i])
        plot(DS96_Si[:, 1, i], DS96_Si[:, 2, i], linestyle="--", linewidth=2, color = "blue")
        plot(DS96_Ag12ML_pos1[:, 1, i], DS96_Ag12ML_pos1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        plot(DS96_Ag12ML_pos2[:, 1, i], DS96_Ag12ML_pos2[:, 2, i], linestyle="-", linewidth=2, color = "orange")
        plot(DS96_Ag20ML_pos1[:, 1, i], DS96_Ag20ML_pos1[:, 2, i], linestyle="-", linewidth=2, color = "green")
        plot(DS96_Ag20ML_pos2[:, 1, i], DS96_Ag20ML_pos2[:, 2, i], linestyle="-", linewidth=2, color = "red")

        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag12ML pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag12ML pos.2")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag20ML pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "red", label = "Ag20ML pos.2")

    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DS96_Si[:, 1, i], DS96_Si[:, 3, i], linestyle="--", linewidth=2, color = "blue")
        plot(DS96_Ag12ML_pos1[:, 1, i], DS96_Ag12ML_pos1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        plot(DS96_Ag12ML_pos2[:, 1, i], DS96_Ag12ML_pos2[:, 3, i], linestyle="-", linewidth=2, color = "orange")
        plot(DS96_Ag20ML_pos1[:, 1, i], DS96_Ag20ML_pos1[:, 3, i], linestyle="-", linewidth=2, color = "green")
        plot(DS96_Ag20ML_pos2[:, 1, i], DS96_Ag20ML_pos2[:, 3, i], linestyle="-", linewidth=2, color = "red")

        xlim([1.55, 5.5])
        ylim([-1, 1])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag12ML pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag12ML pos.2")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag20ML pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "red", label = "Ag20ML pos.2")


    legend(frameon=false, fontsize = 14)
    tight_layout()
    suptitle("Sample S96 Flat Ag 12 ML and 20ML + 1 nm SiO2", fontsize = 16)
    savefig(join([path0,"\\Pictures\\S96_exp.png"]), transparent = true)
end
