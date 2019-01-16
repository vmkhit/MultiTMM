#using MultiTMM
include("C:\\Users\\vmkhi\\Documents\\Github\\MultiTMM\\examples\\layer_fresnel.jl")
using DelimitedFiles
using MAT
using PyPlot
pygui(true)

matplotlib["rcParams"][:update](["font.size" => 18, "font.weight"=>"normal", "font.family" => "Arial", "text.usetex"=>false])

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


## read simulation from python

file1 = matopen("C:\\Users\\vmkhi\\Desktop\\pymat_JC.mat")
energy = read(file1, "energy")
mNML = read(file1, "NML")
bb = read(file1, "gamb")
m1TanPsi = read(file1, "TanPsi")
m1CosDel = read(file1, "CosDel")
close(file1)


file2 = matopen("C:\\Users\\vmkhi\\Desktop\\pymat_sopra.mat")
#energy = read(file1, "energy")
m2TanPsi = read(file2, "TanPsi")
m2CosDel = read(file2, "CosDel")
close(file2)

energy = energy[:]

nAg_JC= nk_import("Ag(JC-eV)",  ww)
nAg_P= nk_import("Ag(Palik-eV)", ww)
nAg_B = nk_import("Ag(Barber-eV)", ww)
nAg_H = nk_import("Ag(Hageman-eV)", ww)
nAg_We = nk_import("Ag(Werner-eV)", ww)
nAg_S = nk_import("Ag(sopra-eV)", ww)
nSi = nk_import("Si111(sopra)", ww)
nSiO2 = nk_import("SiO2(sopra)", ww)

begin
    TanPsi = Array{Float64}(undef, nw, nth, length(mNML), length(bb))
    CosDel = Array{Float64}(undef, nw, nth, length(mNML),  length(bb))
    tc = 1.0
    for (l, bl) = enumerate(bb)
        for (t, tNML) = enumerate(mNML)
            for j = 1:nth
                for i = 1:nw
                    λ = nm2eV/ww[i]
                    k = 2.0*pi/λ
                    LAir = Layer()
                    LSiO2 = Layer("c", Material("nk", nSiO2[i]), tc)
                    LSi_bot = Layer("c", Material("nk", nSi[i]))
                    #LAg = Layer("c", Material("nk", nAg_S(ww[j])))
                    epsAg = nAg_S[i]^2  - drude_pole(ww[i], 0, 9.17, 0.022) + drude_pole(ww[i], 0, 9.17, bl*0.022)
                    LAg = Layer("c", Material("epsmu", epsAg), tNML*0.40853/sqrt(3))
                    S = Stack([LAir; LSiO2; LAg; LSi_bot], zeros(3))
                    TanPsi[i, j, t, l], CosDel[i, j, t, l] = tmm_ellipso(λ, [k*sin(pi*th[j]/180), 0], S)
                end
            end
        end
    end
end

begin
    clf()
    cmap = ColorMap("autumn")
    cmap2 = ColorMap("winter")
    cmap3 = ColorMap("summer")
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure(figsize=(12, 6))
    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])
        ax1[:plot](DATA1[:, 1, i], DATA1[:, 2, i], linestyle="-",linewidth=2, color = "orange") #12ML
        #ax1[:plot](DATA2[:, 1, i], DATA2[:, 2, i], linestyle="-", linewidth=2, color = "red") #10ML
        #ax1[:plot](DATA3[:, 1, i], DATA3[:, 2, i], linestyle="-", linewidth=2, color = "green") #8ML
        #ax1[:plot](DATA4[:, 1, i], DATA4[:, 2, i], linestyle="-", linewidth=2, color = "blue")  # Si


        for l = 2
            bl = bb[l]
            #ax1[:plot](energy, m2TanPsi[:, i, 1, l], linestyle="--", color = "black") # Si
            #ax1[:plot](energy, m2TanPsi[:, i, 2, l], linestyle="--", color = "black") #10ML
            #ax1[:plot](energy, m2TanPsi[:, i, 3, l], linestyle="--", color = "black") #10ML
            ax1[:plot](energy, TanPsi[:, i, 4, l], linestyle="--", color = "black") #12ML

        end

        xlim([1.55, 5.5])
        ylim([0, 0.8])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end

    ax1[:plot]([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag 12ML")
    #ax1[:plot]([],[], linestyle="-", linewidth=2, color = "red", label = "Ag 10ML")
    #ax1[:plot]([],[], linestyle="-", linewidth=2, color = "green", label = "Ag 8ML")
    #ax1[:plot]([],[], linestyle="-", linewidth=2, color = "blue", label = "Si")
    #suptitle("Ag 12ML + 1nm SiO2", fontsize = 16)
    legend(frameon=false, fontsize = 14)

    ax2 = subplot(122)
    for i = 1:nth
        #color = cmap2((i+2)/float(2*nth))
        ax2[:plot](DATA1[:, 1, i], DATA1[:, 3, i], linestyle="-", linewidth=2, color = "orange") #12ML
        #ax2[:plot](DATA2[:, 1, i], DATA2[:, 3, i], linestyle="-", linewidth=2, color = "red") #10ML
        #ax2[:plot](DATA3[:, 1, i], DATA3[:, 3, i], linestyle="-", linewidth=2, color = "green") #8ML
        #ax2[:plot](DATA4[:, 1, i], DATA4[:, 3, i], linestyle="-", linewidth=2, color = "blue")  # Si


        for l = 2
            bl = bb[l]
            #ax2[:plot](energy, m2CosDel[:, i, 1, l], linestyle="--", color = cmap((l)/float(length(bb))), label = "8ML γ/γ\$_{JC}\$ = $bl")
            #ax2[:plot](energy, m2CosDel[:, i, 2, l], linestyle="--", color = cmap((l)/float(length(bb))), label = "10ML γ/γ\$_{JC}\$ = $bl")
            #ax2[:plot](energy, m2CosDel[:, i, 3, l], linestyle="--", color = cmap((l)/float(length(bb))), label = "12ML γ/γ\$_{JC}\$ = $bl")
            ax2[:plot](energy, CosDel[:, i, 4, l], linestyle="--", color = "black")
        end

        xlim([1.55, 5.5])
        ylim([-1, 0.25])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end
    ax2[:plot]([],[], linestyle="-", linewidth=2, color = "orange", label = "Ag 12ML")
    #ax2[:plot]([],[], linestyle="-", linewidth=2, color = "red", label = "Ag 10ML")
    #ax2[:plot]([],[], linestyle="-", linewidth=2, color = "green", label = "Ag 8ML")
    #ax2[:plot]([],[], linestyle="-", linewidth=2, color = "blue", label = "Si")
    #suptitle("Ag 12ML + 1nm SiO2", fontsize = 16)
    legend(frameon=false, fontsize = 14)

    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_12ML_wfit.png", transparent = true)
end
