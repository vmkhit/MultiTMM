using MultiTMM
using DelimitedFiles
using PyPlot
using MAT
pygui(true)
matplotlib["rcParams"][:update](["font.size" => 18, "font.weight"=>"normal", "font.family" => "Arial", "text.usetex"=>false])


function drude_pole(w, epsb, wp, gam)
    return epsb - wp^2/(w*(w+im*gam))
end

file3 = matopen("C:\\Users\\vmkhi\\Desktop\\pymat_sopra1nm.mat")
energy = read(file3, "energy")
mNML = read(file3, "NML")
bb = read(file3, "gamb")
TanPsi_1nm = read(file3, "TanPsi")
CosDel_1nm = read(file3, "CosDel")
close(file3)

file4 = matopen("C:\\Users\\vmkhi\\Desktop\\pymat_sopra2nm.mat")
energy = read(file4, "energy")
mNML = read(file4, "NML")
bb = read(file4, "gamb")
TanPsi_2nm = read(file4, "TanPsi")
CosDel_2nm = read(file4, "CosDel")
close(file4)


file5 = matopen("C:\\Users\\vmkhi\\Desktop\\pymat_sopra3nm.mat")
energy = read(file5, "energy")
mNML = read(file5, "NML")
bb = read(file5, "gamb")
TanPsi_3nm = read(file5, "TanPsi")
CosDel_3nm = read(file5, "CosDel")
close(file5)


energy = energy[:]
#parameters and w-theta ranges
nm2eV = 1239.84193
nth = 4
th = range(60, 75, length = nth)


# import the experimental data
DATA1 = Array{Float64}(undef, 1103-79, 5, nth)
DATA2 = Array{Float64}(undef, 1103-79, 5, nth)
DATA3 = Array{Float64}(undef, 1103-79, 5, nth)
DATA4 = Array{Float64}(undef, 1103-79, 5, nth)
DATA5 = Array{Float64}(undef, 1103-79, 5, nth)
DATA6 = Array{Float64}(undef, 1103-79, 5, nth)


path_ell = "C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Nanogune\\Ellipsometry\\"



for i = 1:nth
    theta = round(Int, th[i])

    file1 = readdlm(join([path_ell, "190124_Chip-Inaki-Laura-HS5-12ML_an15_Ag_", string(theta),".pae"]))
    file2 = readdlm(join([path_ell, "190124_Chip-Inaki-Laura-HS5-12ML_an15_Ag_pos2_", string(theta),".pae"]))
    file3 = readdlm(join([path_ell, "190124_Chip-Inaki-Laura-HS5-12ML_an15_Si_", string(theta),".pae"]))

    file4 = readdlm(join([path_ell, "190120_Chip-Inaki-Laura-S80-10ML_an15_Ag_pos1_", string(theta),".pae"]))
    file5 = readdlm(join([path_ell, "190120_Chip-Inaki-Laura-S80-10ML_an15_Ag_pos2_", string(theta),".pae"]))
    file6 = readdlm(join([path_ell, "190120_Chip-Inaki-Laura-S80-10ML_an15_Si_", string(theta),".pae"]))

    DATA1[:, :, i] = file1[79:end-1, 1:5]
    DATA2[:, :, i] = file2[79:end-1, 1:5]
    DATA3[:, :, i] = file3[79:end-1, 1:5]
    DATA4[:, :, i] = file4[79:end-1, 1:5]
    DATA5[:, :, i] = file5[79:end-1, 1:5]
    DATA6[:, :, i] = file6[79:end-1, 1:5]
end

ww = DATA1[1, :, 1]
ww = range(1, 6, length = 500)
nw = length(ww)


begin
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 40
    bs = 1
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]

    fig = figure(figsize=(12, 6))
    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])

        for l = nb1:bs:nb2
            bl = bb[l]
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            ax1[:plot](energy, TanPsi_1nm[:, i, 5, l], linestyle="--", color = cm(cc), alpha=1.0) # Si
        end

        ax1[:plot](DATA1[:, 1, i], DATA1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        #ax1[:plot](DATA2[:, 1, i], DATA2[:, 2, i], linestyle="-", linewidth=2, color = "green")
        #ax1[:plot](DATA3[:, 1, i], DATA3[:, 2, i], linestyle="--", linewidth=2, color = "orange")

        #ax1[:plot](DATA4[:, 1, i], DATA4[:, 2, i], linestyle="-", linewidth=2, color = "black")
        #ax1[:plot](DATA5[:, 1, i], DATA5[:, 2, i], linestyle="-", linewidth=2, color = "green")
        #ax1[:plot](DATA6[:, 1, i], DATA6[:, 2, i], linestyle="--", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([0, 0.8])
        yticks([0.0, 0.2, 0.4, 0.6])
        title("Ag 15ML + 1nm SiO2: Prepatterned", fontsize = 16)
    end
    #ax1[:plot]([], [], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    #ax1[:plot]([], [], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")
    #ax1[:plot]([], [], linestyle="--", linewidth=2, color = "orange", label = "Si")

    ax1[:set_xlabel]("Photon energy (eV)", color = "black", alpha = 1.0)
    ax1[:set_ylabel]("tan(\$ \\Psi \$)", color = "black", alpha = 1.0)
    ax1[:tick_params](labelbottom=true)


    legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
                                         labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)

    ax2 = subplot(122)


    for i = 1:nth

        for l = nb1:bs:nb2
            bl = bb[l]
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            ax2[:plot](energy, CosDel_1nm[:, i, 5, l], linestyle="--", color = cm(cc), alpha=1.0)
        end

        ax2[:plot](DATA1[:, 1, i], DATA1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        #ax2[:plot](DATA2[:, 1, i], DATA2[:, 3, i], linestyle="-", linewidth=2, color = "green")
        #ax2[:plot](DATA3[:, 1, i], DATA3[:, 3, i], linestyle="--", linewidth=2, color = "orange")

        #ax2[:plot](DATA4[:, 1, i], DATA4[:, 3, i], linestyle="-", linewidth=2, color = "black")
        #ax2[:plot](DATA5[:, 1, i], DATA5[:, 3, i], linestyle="-", linewidth=2, color = "green")
        #ax2[:plot](DATA6[:, 1, i], DATA6[:, 3, i], linestyle="--", linewidth=2, color = "orange")

        ax2[:set_xlim]([2.0, 5.5])
        ax2[:set_xticks]([2.0, 3, 4, 5])
        ylim([-1, 0.25])
        yticks([-1, -0.8, -0.6, -0.4, -0.2, 0])
        title("Ag 15ML + 1nm SiO2: Prepatterned", fontsize = 16)
    end
    #ax2[:plot]([], [], linestyle="-", linewidth=2, color = "black", label = "Ag pos.1")
    #ax2[:plot]([], [], linestyle="-", linewidth=2, color = "green", label = "Ag pos.2")
    #ax2[:plot]([], [], linestyle="--", linewidth=2, color = "orange", label = "Si")

    ax2[:set_xlabel]("Photon energy (eV)", color = "black", alpha = 1.0)
    ax2[:set_ylabel]("cos(\$\\Delta\$)", color = "black", alpha = 1.0)

    ax2[:tick_params](labelbottom=true)
    legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
                                         labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)
    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Chip-Laura-Ignacio-HS5\\Pictures\\Ellipsometry_HS5_withFit_1nmSiO2.png", transparent = true, dpi = 600)
end
