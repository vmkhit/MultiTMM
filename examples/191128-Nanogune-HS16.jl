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
DHS16 = Array{Float64}(undef, 995, 5, nth)


path_HS16 = "C:\\Users\\vmkhitaryan\\Documents\\Projects\\Ag Ribbons\\Chip-Andrew-HS16\\Ellipsometery\\"
for i = 1:nth
    theta = round(Int, th[i])
    fHS16 = readdlm(join([path_HS16, "2019-11-23-Chip-Andrew-HS16-inc60to75_an15_MS_Ag_", string(theta),".pae"]))
    DHS16[:, :, i] = fHS16[79:end-1, 1:5]
end

ww = DHS16[:, 1, 1]
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
    blist = [5, 10, 20, 50]
    TanPsi = Array{Float64}(undef, nw, nth, length(blist))
    CosDel = Array{Float64}(undef, nw, nth, length(blist))
    NML = 18
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
        plot(DHS16[:, 1, i], DHS16[:, 2, i], linestyle="-", linewidth=2, color = "orange")
        for (l, bl) = enumerate(blist)
             plot(ww, TanPsi[:, i, l], linestyle="--", color = cmap2((l)/float(length(blist))))
        end
        xlim([1.55, 5.5])
        ylim([0, 1.0])
        xlabel("Photon energy (eV)")
        ylabel("tan(\$ \\Psi \$)")
    end
    plot([],[], linestyle="-", linewidth=2, color = "orange", label = "HS16")
    legend(frameon=false, fontsize = 14)
    tight_layout()

    subplot(1, 2, 2)
    for i = 1:nth
        color = cmap2((i+2)/float(2*nth))
        plot(DHS16[:, 1, i], DHS16[:, 3, i], linestyle="-", linewidth=2, color = "blue")
        for (l, bl) = enumerate(blist)
            plot(ww, CosDel[:, i, l], linestyle ="--", color = cmap2((l)/float(length(blist))))
        end
        xlim([1.55, 5.5])
        ylim([-1, 0.25])
        xlabel("Photon energy (eV)")
        ylabel("cos(\$\\Delta\$)")
        #title("Ag 12ML + 2nm SiO2", fontsize = 16)
    end
    plot([],[], linestyle="-", linewidth=2, color = "blue", label = "HS16")
    legend(frameon=false, fontsize = 14)
    tight_layout()
    suptitle("Ag 18ML + ? nm SiO2", fontsize = 16)
    savefig(join([path_HS16,"\\Pictures\\HS16_exp_fit.png"]), transparent = true)
end
