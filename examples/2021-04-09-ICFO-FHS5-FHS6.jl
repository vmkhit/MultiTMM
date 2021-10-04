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
nth = 3
th = range(65, 75, length = nth)

# import the experimental data
path0 = "C:/Users/vmkhitaryan/Documents/Projects/Ag Ribbons/Chip-Andrew-FSH5-FHS6/Ellipsometry/"

#FHS5
DFHS5_Si = Array{Float64}(undef, 1074-79, 5, nth)
DFHS5_Ag_pos1 = Array{Float64}(undef, 1074-79, 5, nth)
DFHS5_Ag_pos2 = Array{Float64}(undef, 1074-79, 5, nth)

#FHS5
DFHS6_Si = Array{Float64}(undef, 1074-79, 5, nth)
DFHS6_Ag_pos1 = Array{Float64}(undef, 1074-79, 5, nth)
DFHS6_Ag_pos2 = Array{Float64}(undef, 1074-79, 5, nth)

fFHS5_Si = readdlm(join([path0, "2021-05-04-Chip-Andrew-FHS5_an25_MS_Si_", string(65),".pae"]))


for i = 1:nth
    theta = round(Int, th[i])
    # FHS4
    fFHS5_Si = readdlm(join([path0, "2021-05-04-Chip-Andrew-FHS5_an25_MS_Si_", string(theta),".pae"]))
    fFHS5_Ag_pos1 = readdlm(join([path0, "2021-05-04-Chip-Andrew-FHS5_an25_MS_Ag_pos1_", string(theta),".pae"]))
    fFHS5_Ag_pos2 = readdlm(join([path0, "2021-05-04-Chip-Andrew-FHS5_an25_MS_Ag_pos2_", string(theta),".pae"]))

    fFHS6_Si = readdlm(join([path0, "2021-05-04-Chip-Andrew-FHS6_an25_MS_Si_", string(theta),".pae"]))
    fFHS6_Ag_pos1 = readdlm(join([path0, "2021-05-04-Chip-Andrew-FHS6_an25_MS_Ag_pos1_", string(theta),".pae"]))
    fFHS6_Ag_pos2 = readdlm(join([path0, "2021-05-04-Chip-Andrew-FHS6_an25_MS_Ag_pos2_", string(theta),".pae"]))


    DFHS5_Si[:, 1:5, i] = fFHS5_Si[79:end-1, 1:5]
    DFHS5_Ag_pos1[:, 1:5, i] = fFHS5_Ag_pos1[79:end-1, 1:5]
    DFHS5_Ag_pos2[:, 1:5, i] = fFHS5_Ag_pos2[79:end-1, 1:5]

    DFHS6_Si[:, 1:5, i] = fFHS6_Si[79:end-1, 1:5]
    DFHS6_Ag_pos1[:, 1:5, i] = fFHS6_Ag_pos1[79:end-1, 1:5]
    DFHS6_Ag_pos2[:, 1:5, i] = fFHS6_Ag_pos2[79:end-1, 1:5]
end

ww = DFHS5_Si[:, 1, 1]
nw = length(ww)

#================================ Plot FHS5 ===================================#
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
        plot(DFHS5_Si[:, 1, i], DFHS5_Si[:, 2, i], linestyle="--", linewidth=2, color = "blue")
        plot(DFHS5_Ag_pos1[:, 1, i], DFHS5_Ag_pos1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        plot(DFHS5_Ag_pos2[:, 1, i], DFHS5_Ag_pos2[:, 2, i], linestyle="-", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black",  label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")


    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DFHS5_Si[:, 1, i], DFHS5_Si[:, 3, i], linestyle="--", linewidth=2, color = "blue")
        plot(DFHS5_Ag_pos1[:, 1, i], DFHS5_Ag_pos1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        plot(DFHS5_Ag_pos2[:, 1, i], DFHS5_Ag_pos2[:, 3, i], linestyle="-", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        ylim([-1, 1])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")

    legend(frameon=false, fontsize = 14)
    suptitle("Sample FHS5 Ag 5nm + 0.7 nm SiO\$_2\$", fontsize = 22, fontweight = "bold")
    tight_layout()
    savefig(join([path0,"\\Pictures\\FHS5_exp.png"]), transparent = true)
end
#================================ Plot FHS6 ===================================#
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
        plot(DFHS6_Si[:, 1, i], DFHS6_Si[:, 2, i], linestyle="--", linewidth=2, color = "blue")
        plot(DFHS6_Ag_pos1[:, 1, i], DFHS6_Ag_pos1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        plot(DFHS6_Ag_pos2[:, 1, i], DFHS6_Ag_pos2[:, 2, i], linestyle="-", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black",  label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")


    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DFHS6_Si[:, 1, i], DFHS6_Si[:, 3, i], linestyle="--", linewidth=2, color = "blue")
        plot(DFHS6_Ag_pos1[:, 1, i], DFHS6_Ag_pos1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        plot(DFHS6_Ag_pos2[:, 1, i], DFHS6_Ag_pos2[:, 3, i], linestyle="-", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        ylim([-1, 1])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")

    legend(frameon=false, fontsize = 14)
    suptitle("Sample FHS6 Ag 5nm + 0.7 nm SiO\$_2\$", fontsize = 22, fontweight = "bold")
    tight_layout()
    savefig(join([path0,"\\Pictures\\FHS6_exp.png"]), transparent = true)
end

#============================== Plot FHS5 - FHS6 ==============================#
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

        plot(DFHS5_Si[:, 1, i], DFHS5_Si[:, 2, i], linestyle="-", linewidth=2, color = "blue")
        plot(DFHS5_Ag_pos1[:, 1, i], DFHS5_Ag_pos1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        plot(DFHS5_Ag_pos2[:, 1, i], DFHS5_Ag_pos2[:, 2, i], linestyle="-", linewidth=2, color = "orange")

        plot(DFHS6_Si[:, 1, i], DFHS6_Si[:, 2, i], linestyle="--", linewidth=2, color = "blue")
        plot(DFHS6_Ag_pos1[:, 1, i], DFHS6_Ag_pos1[:, 2, i], linestyle="--", linewidth=2, color = "black")
        plot(DFHS6_Ag_pos2[:, 1, i], DFHS6_Ag_pos2[:, 2, i], linestyle="--", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="-", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black",  label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")


    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DFHS5_Si[:, 1, i], DFHS5_Si[:, 3, i], linestyle="-", linewidth=2, color = "blue")
        plot(DFHS5_Ag_pos1[:, 1, i], DFHS5_Ag_pos1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        plot(DFHS5_Ag_pos2[:, 1, i], DFHS5_Ag_pos2[:, 3, i], linestyle="-", linewidth=2, color = "orange")

        plot(DFHS6_Si[:, 1, i], DFHS6_Si[:, 3, i], linestyle="--", linewidth=2, color = "blue")
        plot(DFHS6_Ag_pos1[:, 1, i], DFHS6_Ag_pos1[:, 3, i], linestyle="--", linewidth=2, color = "black")
        plot(DFHS6_Ag_pos2[:, 1, i], DFHS6_Ag_pos2[:, 3, i], linestyle="--", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        ylim([-1, 1])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="-", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")

    legend(frameon=false, fontsize = 14)
    suptitle("Samples FHS5 and FHS6", fontsize = 22, fontweight = "bold")
    tight_layout()
    savefig(join([path0,"\\Pictures\\FHS5_FHS6_exp.png"]), transparent = true)
end
