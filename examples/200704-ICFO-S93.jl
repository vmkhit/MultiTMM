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
DATA1 = Array{Float64}(undef, 1103-79, 5, nth, 4)
DATA2 = Array{Float64}(undef, 1103-79, 5, nth, 9)

path_ell = "C:\\Users\\vmkhitaryan\\Documents\\Projects\\Ag Ribbons\\Chip-Andrew-S93\\Ellipsometry\\"
for i = 1:nth
    theta = round(Int, th[i])

    file1 = readdlm(join([path_ell, "2020-06-04-Chip-Andrew-S93_an25_MS_Ag15ML_", string(theta),".pae"]))
    file2 = readdlm(join([path_ell, "2020-06-04-Chip-Andrew-S93_an25_MS_Ag15ML_again_", string(theta),".pae"]))
    file3 = readdlm(join([path_ell, "2020-06-04-Chip-Andrew-S93_an25_MS_Ag20ML_", string(theta),".pae"]))
    file4 = readdlm(join([path_ell, "2020-06-04-Chip-Andrew-S93_an25_MS_Si_", string(theta),".pae"]))

    DATA1[:, :, i, 1] = file1[79:end-1, 1:5]
    DATA1[:, :, i, 2] = file2[79:end-1, 1:5]
    DATA1[:, :, i, 3] = file3[79:end-1, 1:5]
    DATA1[:, :, i, 4] = file4[79:end-1, 1:5]
    for j = 1:9
        jj = j-1
        filex = readdlm(join([path_ell, "2020-06-29-Chip-S93-inc60to_an15_SiStart_$(jj)mm_", string(theta),".pae"]))
        DATA2[:, :, i, j] = filex[79:end-1, 1:5]
    end
end



w = DATA1[:, 1, 1, 1]
nw = length(w)

nAg_JC= nk_import("Ag(JC-eV)",  w)
nAg_P= nk_import("Ag(Palik-eV)", w)
nAg_B = nk_import("Ag(Barber-eV)", w)
nAg_H = nk_import("Ag(Hageman-eV)", w)
nAg_We = nk_import("Ag(Werner-eV)", w)
nAg_S = nk_import("Ag(sopra-eV)", w)
nSi = nk_import("Si111(sopra)", w)
nSiO2 = nk_import("SiO2(sopra)", w)

begin
    blist = [5, 10, 20, 50]
    TanPsi = Array{Float64}(undef, nw, nth, length(blist))
    CosDel = Array{Float64}(undef, nw, nth, length(blist))
    NML = 15
    tc = 2.0
    tML = 0.40853/sqrt(3)
    for (l, bl) = enumerate(blist)
        for j = 1:nth
            for i = 1:nw
                λ = nm2eV/w[i]
                k = 2.0*pi/λ
                LAir = Layer()
                LSiO2 = Layer("c", Material("nk", nSiO2[i]), tc)
                LSi_bot = Layer("c", Material("nk", nSi[i]))
                #LAg = Layer("c", Material("nk", nAg_S(ww[j])))

                epsAg = nAg_S[i]^2  - drude_pole(w[i], 0, 9.17, 0.021) + drude_pole(w[i], 0, 9.17, bl*0.021)
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
    cmap2 = ColorMap("gist_heat")
    nb1 = 1
    nb2 = 20
    bs = 1

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure("S92 ellipsometric parameters", figsize=(12, 6))
    clf()

    subplot(1, 2, 1)
    for i = 1:nth
        thi = round(Int, th[i])
        for i = 2:2
            color = cmap2((i+2)/float(2*nth))
            for j = 1:4
                plot(DATA1[:, 1, i, j], DATA1[:, 2, i, j], linestyle="-", linewidth=2, color =  cmap1(j/4))
            end
            for j = 1:9
                plot(DATA2[:, 1, i, j], DATA2[:, 2, i, j], linestyle="--", linewidth=2, color = cmap2(j/8))
            end
            #=
            for (l, bl) = enumerate(blist)
                plot(ww, CosDel[:, i, l], linestyle ="--", color = cmap2((l)/float(length(blist))))
            end
            =#
        end
    end
    xlim([1.55, 5.5])
    ylim([0, 1.0])
    xlabel("Photon energy (eV)")
    ylabel("tan(\$ \\Psi \$)")

    #legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 2:2
        color = cmap2((i+2)/float(2*nth))
        for j = 1:4
            plot(DATA1[:, 1, i, j], DATA1[:, 3, i,j], linestyle="-", linewidth=2, color =  cmap1(j/4))
        end
        for j = 1:9
            plot(DATA2[:, 1, i, j], DATA2[:, 3, i, j], linestyle="--", linewidth=2, color =  cmap2(j/8))
        end
        #=
        for (l, bl) = enumerate(blist)
            plot(ww, CosDel[:, i, l], linestyle ="--", color = cmap2((l)/float(length(blist))))
        end
        =#
    end
    xlim([1.55, 5.5])
    ylim([-1, 0.25])
    xlabel("Photon energy (eV)")
    ylabel("cos(\$\\Delta\$)")

    #legend(frameon=false, fontsize = 14)
    tight_layout()
    suptitle("S93", fontsize = 16)
    savefig(join([path_ell,"/Pictures/S93.png"]), transparent = false)
end
