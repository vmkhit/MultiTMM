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
path0 = "C:/Users/vmkhitaryan/Documents/Projects/Ag Ribbons/Chip-Andrew-FHS4/Ellipsometry/"

#HS30
DFHS4_Si = Array{Float64}(undef, 1103-79, 5, nth)
DFHS4_Ag_pos1 = Array{Float64}(undef, 1103-79, 5, nth)
DFHS4_Ag_pos2 = Array{Float64}(undef, 1103-79, 5, nth)
DFHS4_Ag_pos3 = Array{Float64}(undef, 1103-79, 5, nth)

for i = 1:nth
    theta = round(Int, th[i])
    # FHS4
    fFHS4_Si = readdlm(join([path0, "2021-04-08-Chip-FHS4_22ML_an25_Si_", string(theta),".pae"]))
    fFHS4_Ag_pos1 = readdlm(join([path0, "2021-04-08-Chip-FHS4_22ML_an25_Ag_pos1_", string(theta),".pae"]))
    #fFHS4_Ag_pos2 = readdlm(join([path0, "2021-04-08-Chip-FHS4_22ML_an25_Ag_pos2_", string(theta),".pae"]))
    fFHS4_Ag_pos3 = readdlm(join([path0, "2021-04-08-Chip-FHS4_22ML_an25_Ag_pos3_", string(theta),".pae"]))

    DFHS4_Si[:, :, i] = fFHS4_Si[79:end-1, 1:5]
    DFHS4_Ag_pos1[:, :, i] = fFHS4_Ag_pos1[79:end-1, 1:5]
    #DFHS4_Ag_pos2[:, :, i] = fFHS4_Ag_pos2[79:end-1, 1:5]
    DFHS4_Ag_pos3[:, :, i] = fFHS4_Ag_pos3[79:end-1, 1:5]
end

ww = DFHS4_Si[:, 1, 1]
nw = length(ww)

#================================ Plot FHS4 ===================================#
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
        plot(DFHS4_Si[:, 1, i], DFHS4_Si[:, 2, i], linestyle="--", linewidth=2, color = "blue")
        plot(DFHS4_Ag_pos1[:, 1, i], DFHS4_Ag_pos1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        #plot(DFHS4_Ag_pos2[:, 1, i], DFHS4_Ag_pos2[:, 2, i], linestyle="-", linewidth=2, color = "orange")
        plot(DFHS4_Ag_pos3[:, 1, i], DFHS4_Ag_pos3[:, 2, i], linestyle="-", linewidth=2, color = "green")

        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black",  label = "Ag pos.1")
    #plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag pos.2")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")


    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DFHS4_Si[:, 1, i], DFHS4_Si[:, 3, i], linestyle="--", linewidth=2, color = "blue")
        plot(DFHS4_Ag_pos1[:, 1, i], DFHS4_Ag_pos1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        #plot(DFHS4_Ag_pos2[:, 1, i], DFHS4_Ag_pos2[:, 3, i], linestyle="-", linewidth=2, color = "orange")
        plot(DFHS4_Ag_pos3[:, 1, i], DFHS4_Ag_pos3[:, 3, i], linestyle="-", linewidth=2, color = "green")

        xlim([1.55, 5.5])
        ylim([-1, 1])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "Si")
    plot([],[], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    #plot([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag pos.2")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")

    legend(frameon=false, fontsize = 14)
    suptitle("Sample FHS4 Ag 22ML + 0.7 nm SiO\$_2\$", fontsize = 22, fontweight = "bold")
    tight_layout()
    savefig(join([path0,"\\Pictures\\FHS4_exp.png"]), transparent = true)
end
