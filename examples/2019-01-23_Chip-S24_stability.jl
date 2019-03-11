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

path_ell = "C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Ellipsometer\\171214_Chip-Zaka-1710-2\\171214_Chip-Zaka-1710-2\\"

DATA1 = Array{Float64}(undef, 1103-79, 5, nth)
DATA2 = Array{Float64}(undef, 1103-79, 5, nth)
DATA3 = Array{Float64}(undef, 1103-79, 5, nth)

for i = 1:nth
    theta = round(Int, th[i])
    file1 = readdlm(join([path_ell, "20180408_chip-Zaka-S24_10ML_Plane_CDD_inc60to75_an15_MS_Ag10ML_pos1_", string(theta),".pae"]))
    file2 = readdlm(join([path_ell, "20180419_chip-Zaka-S24_10ML_Plane_CDD_inc25to75_an15_MS_UV-VIS_Ag10ML_pos1_", string(theta),".pae"]))
    file3 = readdlm(join([path_ell, "20180423_chip-Zaka-S24_10ML_Plane_CDD_inc60to75_an15_MS_Ag_area_pos1_", string(theta),".pae"]))
    DATA1[:, :, i] = file1[79:end-1, 1:5]
    DATA2[:, :, i] = file2[79:end-1, 1:5]
    DATA3[:, :, i] = file3[79:end-1, 1:5]
end

ww = DATA1[:, 1, 1]
nw = length(ww)

begin
    cm = ColorMap("cool")
    nb1 = 10
    nb2 = 10
    bs = 1
    labels = [L"40^o", L"45^o", L"50^o", L"55^o", L"60^o", L"65^o", L"70^o", L"75^o"]
    fig = figure(figsize=(12, 6))

    ax1 = subplot(121)
    for i = 1:nth
        thi = round(Int, th[i])
        ax1[:plot](DATA1[:, 1, i], DATA1[:, 2, i], linestyle="-", linewidth=2, color = "black")
        ax1[:plot](DATA2[:, 1, i], DATA2[:, 2, i], linestyle="-", linewidth=2, color = "green")
        ax1[:plot](DATA3[:, 1, i], DATA3[:, 2, i], linestyle="-", linewidth=2, color = "orange")
        xlim([1.55, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([0, 0.8])
        yticks([0.0, 0.2, 0.4, 0.6, 0.8])
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end
    ax1[:set_xlabel]("Photon energy (eV)", color = "white", alpha = 0)
    ax1[:set_ylabel]("tan(\$ \\Psi \$)", color = "white", alpha = 0)
    #ax1[:tick_params](labelbottom=true)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)
    ax1[:minorticks_on]()
    ax1[:tick_params](axis="x", which="minor", bottom ="on")
    ax1[:tick_params](axis="y", which="minor", left = "off")

    ax2 = subplot(122)
    for i = 1:nth
        ax2[:plot](DATA1[:, 1, i], DATA1[:, 3, i], linestyle="-", linewidth=2, color = "black")
        ax2[:plot](DATA2[:, 1, i], DATA2[:, 3, i], linestyle="-", linewidth=2, color = "green")
        ax2[:plot](DATA3[:, 1, i], DATA3[:, 3, i], linestyle="-", linewidth=2, color = "orange")

        xlim([1.55, 5.5])
        xticks([2.0, 3, 4, 5])
        ylim([-1, 0.25])
        yticks([-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2])
        #title("Ag 12ML + 1nm SiO2", fontsize = 16)
    end
    ax2[:set_xlabel]("Photon energy (eV)", color = "white", alpha = 0)
    ax2[:set_ylabel]("cos(\$\\Delta\$)", color = "white", alpha = 0)
    ax2[:minorticks_on]()
    ax2[:tick_params](axis="x", which="minor", bottom ="on")
    ax2[:tick_params](axis="y", which="minor", left = "off")
    #ax2[:tick_params](labelbottom=true)
    #legend(bbox_to_anchor=(-0.01, 1.01), borderaxespad=0., loc = 2, frameon=false, fontsize = 24,
    #                                     labelspacing = 0.3, handlelength =1.25, handletextpad=0.3)
    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_S24_stability.eps", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_S24_stability.pdf", transparent = true, dpi = 600)
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Ribbons\\Cleaned experimental data\\Ellipsometry_S24_stability.png", transparent = true, dpi = 600)
end
