using MultiTMM
#include("C:\\Users\\vmkhi\\Documents\\Github\\MultiTMM\\examples\\layer_fresnel.jl")
using DelimitedFiles
using MAT
using PyPlot
pygui(true)

matplotlib["rcParams"][:update](["font.size" => 18, "font.weight"=>"normal", "font.family" => "Arial", "text.usetex"=>false])

function drude_pole(w, epsb, wp, gam)
    return epsb - wp^2/(w*w + 1.0im*w*gam)
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


## read simulation from python

file1 = matopen("C:\\Users\\vmkhi\\Desktop\\pymat_sopra_Si.mat")
energy = read(file1, "energy")
mNML = read(file1, "NML")
bb = read(file1, "gamb")
m1TanPsi = read(file1, "TanPsi")
m1CosDel = read(file1, "CosDel")
close(file1)

file2 = matopen("C:\\Users\\vmkhi\\Desktop\\pymat_sopra.mat")
energy = read(file2, "energy")
mNML = read(file2, "NML")
bb = read(file2, "gamb")
m2TanPsi = read(file2, "TanPsi")
m2CosDel = read(file2, "CosDel")
close(file2)

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

nAg_JC= nk_import("Ag(JC-eV)",  ww)
nAg_P= nk_import("Ag(Palik-eV)", ww)
nAg_B = nk_import("Ag(Barber-eV)", ww)
nAg_H = nk_import("Ag(Hageman-eV)", ww)
nAg_We = nk_import("Ag(Werner-eV)", ww)
nAg_S = nk_import("Ag(sopra-eV)", ww)
nSi = nk_import("Si111(sopra)", ww)
nSiO2 = nk_import("SiO2(sopra)", ww)


############# Si #################################################
begin
    cm = ColorMap("cool")
    nb1 = 10
    nb2 = 10
    bs = 1
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure(figsize=(12, 4.685))

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])


        for l = nb1:bs:nb2
            bl = bb[l]
            ax1[:plot](energy, m2TanPsi[:, i, 1, l], linestyle="--", color = "black", alpha=1.0) # Si
        end
        ax1[:plot](DATA4[:, 1, i], DATA4[:, 2, i], linestyle="-", linewidth=2, color = "blue")  # Si

        xlim([1.55, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([0, 0.8])
        yticks([0.0, 0.2, 0.4, 0.6])
        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("tan(\$ \\Psi \$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end
    ax1[:plot]([],[], linestyle="-", linewidth=2, color = "blue", label = "bare Si") # Si
    ax1[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    ax1[:tick_params](labelbottom=false)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)

    ax2 = subplot(122)
    for i = 1:nth
        #color = cmap2((i+2)/float(2*nth))

        for l = nb1:bs:nb2
            bl = bb[l]
            ax2[:plot](energy, m2CosDel[:, i, 1, l], linestyle="--", color = "black", alpha=1.0)
        end
        ax2[:plot](DATA4[:, 1, i], DATA4[:, 3, i], linestyle="-", linewidth=2, color = "blue")  # Si

        xlim([2.0, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([-1, 0.25])
        yticks([-1, -0.8, -0.6, -0.4, -0.2, 0])

        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("cos(\$\\Delta\$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end

    ax2[:plot]([],[], linestyle="-", linewidth=2, color = "blue", label = "bare Si") #Si
    ax2[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    ax2[:tick_params](labelbottom=false)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)
    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_Si.eps", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_Si.pdf", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_Si.png", transparent = true, dpi = 600)
end

################ 8ML ############################################
begin
    cm = ColorMap("cool")
    nb1 = 10
    nb2 = 30
    bs = 1
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure(figsize=(12, 4.685))

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])

        for l = nb1:bs:nb2
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            bl = bb[l]
            ax1[:plot](energy, TanPsi_3nm[:, i, 2, l], linestyle="--", color = cm(cc), alpha=1.0) #8ML
        end
        ax1[:plot](DATA3[:, 1, i], DATA3[:, 2, i], linestyle="-", linewidth=2, color = "green") #8ML

        xlim([1.55, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([0, 0.8])
        yticks([0.0, 0.2, 0.4, 0.6])
        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("tan(\$ \\Psi \$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end

    ax1[:plot]([],[], linestyle="-", linewidth=2, color = "green", label = "8 ML Ag(111)")
    ax1[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    ax1[:tick_params](labelbottom=false)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)

    ax2 = subplot(122)
    for i = 1:nth
        #color = cmap2((i+2)/float(2*nth))
        for l = nb1:bs:nb2
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            bl = bb[l]
            ax2[:plot](energy, CosDel_3nm[:, i, 2, l], linestyle="--", color = cm(cc), alpha=1.0)
        end
        ax2[:plot](DATA3[:, 1, i], DATA3[:, 3, i], linestyle="-", linewidth=2, color = "green") #8ML

        xlim([2.0, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([-1, 0.25])
        yticks([-1, -0.8, -0.6, -0.4, -0.2, 0])
        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("cos(\$\\Delta\$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end


    ax2[:plot]([],[], linestyle="-", linewidth=2, color = "green", label = "8 ML Ag(111)")
    ax2[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    ax2[:tick_params](labelbottom=false)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)


    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_8ML.eps", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_8ML.pdf", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_8ML.png", transparent = true, dpi = 600)
end
############### 10ML ############################################
begin
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 10
    bs = 1
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure(figsize=(12, 4.685))

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])
        for l = nb1:bs:nb2
            bl = bb[l]
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            ax1[:plot](energy, TanPsi_2nm[:, i, 3, l], linestyle="--", color = cm(cc), alpha=1.0) #10ML
        end
        ax1[:plot](DATA2[:, 1, i], DATA2[:, 2, i], linestyle="-", linewidth=2, color = "red") #10ML

        xlim([1.55, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([0, 0.8])
        yticks([0.0, 0.2, 0.4, 0.6])
        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("tan(\$ \\Psi \$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end


    ax1[:plot]([],[], linestyle="-", linewidth=2, color = "red", label = "10 ML Ag(111)")
    ax1[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    ax1[:tick_params](labelbottom=false)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)

    ax2 = subplot(122)
    for i = 1:nth
        #color = cmap2((i+2)/float(2*nth))

        for l = nb1:bs:nb2
            bl = bb[l]
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            ax2[:plot](energy, CosDel_2nm[:, i, 3, l], linestyle="--", color = cm(cc), alpha=1.0)
        end
        ax2[:plot](DATA2[:, 1, i], DATA2[:, 3, i], linestyle="-", linewidth=2, color = "red") #10ML

        xlim([2.0, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([-1, 0.25])
        yticks([-1, -0.8, -0.6, -0.4, -0.2, 0])
        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("cos(\$\\Delta\$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end

    ax2[:plot]([],[], linestyle="-", linewidth=2, color = "red", label = "10 ML Ag(111)")
    ax2[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    ax2[:tick_params](labelbottom=false)

    #suptitle("Ag 12ML + 1nm SiO2", fontsize = 16)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)

    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_10ML.eps", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_10ML.pdf", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_10ML.png", transparent = true, dpi = 600)
end
############### 12 ML
begin
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 5
    bs = 1
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure(figsize=(12, 4.685))

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])
        for l = nb1:bs:nb2
            bl = bb[l]
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            ax1[:plot](energy, TanPsi_1nm[:, i, 4, l], linestyle="--", color = cm(cc), alpha=1.0) #12ML
        end
        ax1[:plot](DATA1[:, 1, i], DATA1[:, 2, i], linestyle="-",linewidth=2, color = "orange") #12ML

        xlim([1.55, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([0, 0.8])
        yticks([0.0, 0.2, 0.4, 0.6])
        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("tan(\$ \\Psi \$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end

    ax1[:plot]([],[], linestyle="-", linewidth=2, color = "orange", label = "10 ML Ag(111)")
    ax1[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    #suptitle("Ag 12ML + 1nm SiO2", fontsize = 16)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)

    ax2 = subplot(122)
    for i = 1:nth
        #color = cmap2((i+2)/float(2*nth))
        for l = nb1:bs:nb2
            bl = bb[l]
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            ax2[:plot](energy, CosDel_1nm[:, i, 4, l], linestyle="--", color = cm(cc), alpha=1.0)
        end
        ax2[:plot](DATA1[:, 1, i], DATA1[:, 3, i], linestyle="-", linewidth=2, color = "orange") #12ML

        xlim([2.0, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([-1, 0.25])
        yticks([-1, -0.8, -0.6, -0.4, -0.2, 0])
        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("cos(\$\\Delta\$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end

    ax2[:plot]([],[], linestyle="-", linewidth=2, color = "orange", label = "10 ML Ag(111)")
    ax2[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    #suptitle("Ag 12ML + 1nm SiO2", fontsize = 16)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)
    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_12ML.eps", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_12ML.pdf", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_12ML.png", transparent = true, dpi = 600)
end


################ All ########################################
#=
begin
    cm = ColorMap("cool")
    nb1 = 1
    nb2 = 10
    bs = 1
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure(figsize=(14, 6))

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])
        ax1[:plot](DATA1[:, 1, i], DATA1[:, 2, i], linestyle="-",linewidth=2, color = "orange") #12ML
        ax1[:plot](DATA2[:, 1, i], DATA2[:, 2, i], linestyle="-", linewidth=2, color = "red") #10ML
        #ax1[:plot](DATA3[:, 1, i], DATA3[:, 2, i], linestyle="-", linewidth=2, color = "green") #8ML
        #ax1[:plot](DATA4[:, 1, i], DATA4[:, 2, i], linestyle="-", linewidth=2, color = "black")  # Si

        for l = nb1:bs:nb2
            bl = bb[l]
            #ax1[:plot](energy, m2TanPsi[:, i, 1, l], linestyle="--", color = cm(0.7*(l-nb1)/(nb2 - nb1)), alpha=1.0) # Si
            #ax1[:plot](energy, TanPsi_3nm[:, i, 2, l], linestyle="--", color = cm(0.7*(l-nb1)/(nb2 - nb1)), alpha=1.0) #8ML
            ax1[:plot](energy, TanPsi_2nm[:, i, 3, l], linestyle="--", color = cm(0.7*(l-nb1)/(nb2 - nb1)), alpha=1.0) #10ML
            #ax1[:plot](energy, m2TanPsi[:, i, 4, l], linestyle="--", color = cm(0.7*(l-nb1)/(nb2 - nb1)), alpha=1.0) #12ML
        end

        xlim([1.55, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([0, 0.8])
        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("tan(\$ \\Psi \$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end

    #ax1[:plot]([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag 12ML")
    ax1[:plot]([],[], linestyle="-", linewidth=2, color = "red", label = "Ag 10ML")
    #ax1[:plot]([],[], linestyle="-", linewidth=2, color = "green", label = "Ag 8ML")
    #ax1[:plot]([],[], linestyle="-", linewidth=2, color = "black", label = "Si") # Si
    #ax1[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    #suptitle("Ag 12ML + 1nm SiO2", fontsize = 16)
    legend(frameon=false, fontsize = 16)

    ax2 = subplot(122)
    for i = 1:nth
        #color = cmap2((i+2)/float(2*nth))
        #ax2[:plot](DATA1[:, 1, i], DATA1[:, 3, i], linestyle="-", linewidth=2, color = "orange") #12ML
        ax2[:plot](DATA2[:, 1, i], DATA2[:, 3, i], linestyle="-", linewidth=2, color = "red") #10ML
        #ax2[:plot](DATA3[:, 1, i], DATA3[:, 3, i], linestyle="-", linewidth=2, color = "green") #8ML
        #ax2[:plot](DATA4[:, 1, i], DATA4[:, 3, i], linestyle="-", linewidth=2, color = "black")  # Si

        for l = nb1:bs:nb2
            bl = bb[l]
            #ax2[:plot](energy, m2CosDel[:, i, 1, l], linestyle="--", color = cm(0.7*(l-nb1)/(nb2 - nb1)), alpha=1.0)
            #ax2[:plot](energy, CosDel_3nm[:, i, 2, l], linestyle="--", color = cm(0.7*(l-nb1)/(nb2 - nb1)), alpha=1.0)
            ax2[:plot](energy, CosDel_2nm[:, i, 3, l], linestyle="--", color = cm(0.7*(l-nb1)/(nb2 - nb1)), alpha=1.0)
            #ax2[:plot](energy, m2CosDel[:, i, 4, l], linestyle="--", color = cm(0.7*(l-nb1)/(nb2 - nb1)), alpha=1.0)
        end

        xlim([2.0, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([-1, 0.25])
        xlabel("Photon energy (eV)", color = "white", alpha = 0)
        ylabel("cos(\$\\Delta\$)", color = "white", alpha = 0)
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end

    #ax2[:plot]([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag 12ML")
    ax2[:plot]([],[], linestyle="-", linewidth=2, color = "red", label = "Ag 10ML")
    #ax2[:plot]([],[], linestyle="-", linewidth=2, color = "green", label = "Ag 8ML")
    #ax2[:plot]([],[], linestyle="-", linewidth=2, color = "black", label = "Si") #Si
    #ax2[:plot]([],[], linestyle="--", linewidth=2, color = "black", label = "theory")

    #suptitle("Ag 12ML + 1nm SiO2", fontsize = 16)
    legend(frameon=false, fontsize = 16)
    tight_layout(pad = 0.5)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_10ML.eps", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_10ML.pdf", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_10ML.png", transparent = true, dpi = 600)
end
=#
