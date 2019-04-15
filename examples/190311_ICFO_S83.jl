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

########## From Python #####################
file4 = matopen("C:\\Users\\vmkhi\\Desktop\\pymat_sopra2nm.mat")
energy = read(file4, "energy")
mNML = read(file4, "NML")
bb = read(file4, "gamb")
TanPsi_2nm = read(file4, "TanPsi")
CosDel_2nm = read(file4, "CosDel")
close(file4)
energy = energy[:]

# import the experimental data
DATA1 = Array{Float64}(undef, 1074-79, 5, nth)
DATA2 = Array{Float64}(undef, 1074-79, 5, nth)
DATA3 = Array{Float64}(undef, 1103-79, 5, nth)
DATA4 = Array{Float64}(undef, 1103-79, 5, nth)
#DATA5 = Array{Float64}(undef, 1103-79, 5, nth)


path_ell = "C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Nanogune\\Ellipsometry\\"

for i = 1:nth
    theta = round(Int, th[i])
    file1 = readdlm(join([path_ell, "20190301_chip-Zaka-S83_Plane_CDD_inc60o75_an15_MS_Ag_pos1_", string(theta),".pae"]))
    file2 = readdlm(join([path_ell, "20190301_chip-Zaka-S83_Plane_CDD_inc60o75_an15_MS_Ag_pos2_", string(theta),".pae"]))

    file3 = readdlm(join([path_ell, "190226_Chip-Inaki-Laura-S82-12ML_an15_Ag_pos1_", string(theta),".pae"]))
    file4 = readdlm(join([path_ell, "190226_Chip-Inaki-Laura-S82-12ML_an15_Ag_pos2_", string(theta),".pae"]))

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
    tc = 2.0
    tML = 0.40853/sqrt(3)
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
                LAg = Layer("c", Material("epsmu", epsAg), NML*tML)
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
    cm = ColorMap("cool")
    nb1 = 10
    nb2 = 30
    bs = 2

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure("20 ML ellipsometric parameters", figsize=(12, 6))
    for i = 1:nth
        thi = round(Int, th[i])
        ax1 = subplot(1, 2, 1)
        #=
        for l = nb1:bs:nb2
            bl = bb[l]
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            ax1[:plot](energy, TanPsi_2nm[:, i, 4, l], linestyle="--", color = cm(cc), alpha=1.0) # Si
        end
        =#

        plot(DATA1[:, 1, i], DATA1[:, 2, i], linestyle="-", linewidth=2, color = "orange")
        plot(DATA2[:, 1, i], DATA2[:, 2, i], linestyle="-", linewidth=2, color = "green")
        plot(DATA3[:, 1, i], DATA3[:, 2, i], linestyle="--", linewidth=2, color = "blue")
        plot(DATA4[:, 1, i], DATA4[:, 2, i], linestyle="--", linewidth=2, color = "red")


        #=
        for (l, bl) = enumerate(blist)
            plot(ww, TanPsi[:, i, l], linestyle="--", color = cmap2((l)/float(length(blist))))
        end
        =#
        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
        #title("Ag 15ML + 1.5m SiO2", fontsize = 16)
    end
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "S83, pos.1 /190231")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "S83, pos.2 / 190231")
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "S82, pos.1 /190226")
    plot([],[], linestyle="--", linewidth=2, color = "red", label = "S82, pos.2 /190226")

    legend(frameon=false, fontsize = 14)
    for i = 1:nth
        ax2 = subplot(1, 2, 2)
        #color = cmap2((i+2)/float(2*nth))
        #=
        for l = nb1:bs:nb2
            bl = bb[l]
            cc = 1.0 + (1 - 0.4)*(l - nb2)/(nb2 - nb1)
            ax2[:plot](energy, CosDel_2nm[:, i, 4, l], linestyle="--", color = cm(cc), alpha=1.0)
        end
        =#

        plot(DATA1[:, 1, i], DATA1[:, 3, i], linestyle="-", linewidth=2, color = "orange")
        plot(DATA2[:, 1, i], DATA2[:, 3, i], linestyle="-", linewidth=2, color = "green")
        plot(DATA3[:, 1, i], DATA3[:, 3, i], linestyle="--", linewidth=2, color = "blue")
        plot(DATA4[:, 1, i], DATA4[:, 3, i], linestyle="--", linewidth=2, color = "red")


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
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "S83, pos.1 /190231")
    plot([],[], linestyle="-", linewidth=2, color = "green", label = "S83, pos.2 /190231")
    plot([],[], linestyle="--", linewidth=2, color = "blue", label = "S82, pos.1 /190226")
    plot([],[], linestyle="--", linewidth=2, color = "red", label = "S82, pos.2 /190226")



    suptitle("Ag 12ML + 2nm SiO2", fontsize = 16)
    legend(frameon=false, fontsize = 14)
    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Chip-Laura-S83\\Ellipsometry Pictures\\S83_S82)Compare_Ellipsometry_12ML_S83_gam1_10.png", transparent = true)
end
