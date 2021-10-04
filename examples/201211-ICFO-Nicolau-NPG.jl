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
DATA1_VIS = Array{Float64}(undef, 1103-79, 5, nth)
DATA2_VIS = Array{Float64}(undef, 1103-79, 5, nth)
DATA3_VIS = Array{Float64}(undef, 1103-79, 5, nth)
DATA4_VIS = Array{Float64}(undef, 1103-79, 5, nth)

DATA1_NIR = Array{Float64}(undef, 335-79, 5, nth)
DATA2_NIR = Array{Float64}(undef, 335-79, 5, nth)
DATA3_NIR = Array{Float64}(undef, 335-79, 5, nth)
DATA4_NIR = Array{Float64}(undef, 335-79, 5, nth)

path_ell = "C:\\Users\\vmkhitaryan\\Documents\\Projects\\GNR Au Gratings\\Ellipsometry\\Nicolau\\"

for i = 1:nth
    theta = round(Int, th[i])

    file1_VIS = readdlm(join([path_ell, "2020-12-09-Anlyser45_45_75deg_step5deg_SOI_", string(theta),".pae"]))
    file2_VIS = readdlm(join([path_ell, "2020-12-09-Anlyser45_45_75deg_step5deg_AuMica_", string(theta),".pae"]))
    file3_VIS = readdlm(join([path_ell, "2020-12-09-Anlyser45_45_75deg_step5deg_SOI_NPG_VIS_", string(theta),".pae"]))
    file4_VIS = readdlm(join([path_ell, "2020-12-09-Anlyser45_45_75deg_step5deg_AuMicaNPG_VIS_", string(theta),".pae"]))

    file1_NIR = readdlm(join([path_ell, "2020-12-09-Anlyser45_45_75deg_step5deg_SOI_NIR_", string(theta),".pae"]))
    #file2_NIR = readdlm(join([path_ell, "2020-12-09-Anlyser45_45_75deg_step5deg_AuMica_", string(theta),".pae"])) # this needs to be added
    file3_NIR = readdlm(join([path_ell, "2020-12-09-Anlyser45_45_75deg_step5deg_SOI_NPG_NIR_", string(theta),".pae"]))
    file4_NIR = readdlm(join([path_ell, "2020-12-09-Anlyser45_45_75deg_step5deg_AuMicaNPG_NIR_", string(theta),".pae"]))

    DATA1_VIS[:, :, i] = file1_VIS[79:end-1, 1:5]
    DATA2_VIS[:, :, i] = file2_VIS[79:end-1, 1:5]
    DATA3_VIS[:, :, i] = file3_VIS[79:end-1, 1:5]
    DATA4_VIS[:, :, i] = file4_VIS[79:end-1, 1:5]

    DATA1_NIR[:, :, i] = file1_NIR[79:end-1, 1:5]
    #DATA2_NIR[:, :, i] = file2_NIR[79:end-1, 1:5]
    DATA3_NIR[:, :, i] = file3_NIR[79:end-1, 1:5]
    DATA4_NIR[:, :, i] = file4_NIR[79:end-1, 1:5]
end

w_VIS = DATA1_VIS[:, 1, 1, 1]
w_NIR = DATA1_NIR[:, 1, 1, 1]

nw1 = length(w_VIS)
nw2 = length(w_NIR)

#=============================== Plot UV-VIS ==================================#
begin
    cmap1 = ColorMap("tab10")
    cmap2 = ColorMap("gist_heat")
    nb1 = 1
    nb2 = 20
    bs = 1

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure("S92 ellipsometric parameters", figsize=(12, 6))
    clf()

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])
        ax1[:plot](w_VIS, DATA1_VIS[:, 2, i], linestyle = "--", linewidth=2, color =  "black")
        ax1[:plot](w_VIS, DATA3_VIS[:, 2, i], linestyle = "-", linewidth=2, color =  cmap1(i), label = "\$ \\theta = $(thi) ^{\\rm o}\$")
    end
    xlim([1.0, 5.5])
    ylim([0, 9.0])
    xlabel("Photon energy (eV)")
    ylabel("tan(\$ \\Psi \$)")
    legend(frameon=false, fontsize = 14, loc = 1)

    ax2 = subplot(122)
    for i = 1:nth
        thi = round(Int, th[i])
        ax2[:plot](w_VIS, DATA1_VIS[:, 3, i], linestyle = "--", linewidth=2, color =  "black")
        ax2[:plot](w_VIS, DATA3_VIS[:, 3, i], linestyle = "-", linewidth=2, color =  cmap1(i))
    end

    plot([], [], linestyle = "--", linewidth=2, color =  "black", label = "SOI")
    plot([], [], linestyle = "-", linewidth=2, color = "black", label = "NPG on SOI")

    xlim([1.0, 5.5])
    #ylim([-1, 0.25])
    xlabel("Photon energy (eV)")
    ylabel("cos(\$\\Delta\$)")
    legend(frameon=false, fontsize = 14)

    suptitle("NPG on SOI", fontsize = 16)

    tight_layout()

    savefig(join([path_ell,"/Pictures/SOI_UV-VIS.png"]), transparent = false)
end
#------------------ Au Mica ------------------------------#
begin
    cmap1 = ColorMap("tab10")
    cmap2 = ColorMap("gist_heat")
    nb1 = 1
    nb2 = 20
    bs = 1

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure("S92 ellipsometric parameters", figsize=(12, 6))
    clf()

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])
        ax1[:plot](w_VIS, DATA2_VIS[:, 2, i], linestyle = "--", linewidth=2, color =  "black")
        ax1[:plot](w_VIS, DATA4_VIS[:, 2, i], linestyle = "-", linewidth=2, color =  cmap1(i), label = "\$ \\theta = $(thi) ^{\\rm o}\$")
    end
    xlim([1.0, 5.5])
    ylim([0.4, 1.0])
    xlabel("Photon energy (eV)")
    ylabel("tan(\$ \\Psi \$)")
    legend(frameon=false, fontsize = 14, loc = 1)

    ax2 = subplot(122)
    for i = 1:nth
        thi = round(Int, th[i])
        ax2[:plot](w_VIS, DATA2_VIS[:, 3, i], linestyle = "--", linewidth=2, color =  "black")
        ax2[:plot](w_VIS, DATA4_VIS[:, 3, i], linestyle = "-", linewidth=2, color =  cmap1(i))
    end

    plot([], [], linestyle = "--", linewidth=2, color =  "black", label = "Au on Mica")
    plot([], [], linestyle = "-", linewidth=2, color =  "black", label = "NPG on AuMica")

    xlim([1.0, 5.5])
    #ylim([-1, 0.25])
    xlabel("Photon energy (eV)")
    ylabel("cos(\$\\Delta\$)")
    legend(frameon=false, fontsize = 14)
    suptitle("NPG on Au/Mica substrate", fontsize = 16)
    tight_layout()
    savefig(join([path_ell,"/Pictures/MICA_samples_UV-VIS.png"]), transparent = false)
end

#=============================== Plot UV-NIR ==================================#
begin
    cmap1 = ColorMap("tab10")
    cmap2 = ColorMap("gist_heat")
    nb1 = 1
    nb2 = 20
    bs = 1

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure("S92 ellipsometric parameters", figsize=(12, 6))
    clf()

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])
        ax1[:plot](w_NIR, DATA1_NIR[:, 2, i], linestyle = "--", linewidth=2, color =  "black")
        ax1[:plot](w_NIR, DATA3_NIR[:, 2, i], linestyle = "-", linewidth=2, color =  cmap1(i), label = "\$ \\theta = $(thi) ^{\\rm o}\$")
    end
    xlim([0.7, 1.3])
    #ylim([0, 9.0])
    xlabel("Photon energy (eV)")
    ylabel("tan(\$ \\Psi \$)")
    legend(frameon=false, fontsize = 14, loc = 1)

    ax2 = subplot(122)
    for i = 1:nth
        thi = round(Int, th[i])
        ax2[:plot](w_NIR, DATA1_NIR[:, 3, i], linestyle = "--", linewidth=2, color =  "black")
        ax2[:plot](w_NIR, DATA3_NIR[:, 3, i], linestyle = "-", linewidth=2, color =  cmap1(i))
    end

    plot([], [], linestyle = "--", linewidth=2, color =  "black", label = "SOI")
    plot([], [], linestyle = "-", linewidth=2, color = "black", label = "NPG on SOI")

    xlim([0.7, 1.3])
    #ylim([-1, 0.25])
    xlabel("Photon energy (eV)")
    ylabel("cos(\$\\Delta\$)")
    legend(frameon=false, fontsize = 14)

    suptitle("NPG on SOI", fontsize = 16)

    tight_layout()

    savefig(join([path_ell,"/Pictures/SOI_NIR.png"]), transparent = false)
end
#------------------ Au Mica ------------------------------#
begin
    cmap1 = ColorMap("tab10")
    cmap2 = ColorMap("gist_heat")
    nb1 = 1
    nb2 = 20
    bs = 1

    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure("S92 ellipsometric parameters", figsize=(12, 6))
    clf()

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])
        #ax1[:plot](w_NIR, DATA2_NIR[:, 2, i], linestyle = "--", linewidth=2, color =  "black")
        ax1[:plot](w_NIR, DATA4_NIR[:, 2, i], linestyle = "-", linewidth=2, color =  cmap1(i), label = "\$ \\theta = $(thi) ^{\\rm o}\$")
    end
    xlim([0.7, 1.3])
    ylim([0.5, 1.3])
    xlabel("Photon energy (eV)")
    ylabel("tan(\$ \\Psi \$)")
    legend(frameon=false, fontsize = 14, loc = 1)

    ax2 = subplot(122)
    for i = 1:nth
        thi = round(Int, th[i])
        #ax2[:plot](w_NIR, DATA2_NIR[:, 3, i], linestyle = "--", linewidth=2, color =  "black")
        ax2[:plot](w_NIR, DATA4_NIR[:, 3, i], linestyle = "-", linewidth=2, color =  cmap1(i))
    end

    plot([], [], linestyle = "--", linewidth=2, color =  "black", label = "Au on Mica")
    plot([], [], linestyle = "-", linewidth=2, color =  "black", label = "NPG on AuMica")

    xlim([0.7, 1.3])
    #ylim([-1, 0.25])
    xlabel("Photon energy (eV)")
    ylabel("cos(\$\\Delta\$)")
    legend(frameon=false, fontsize = 14)
    suptitle("NPG on Au/Mica substrate", fontsize = 16)
    tight_layout()
    savefig(join([path_ell,"/Pictures/MICA_samples_NIR.png"]), transparent = false)
end
