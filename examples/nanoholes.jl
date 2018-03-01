using MultiTMM
#using HDF5, JLD
using PyPlot

nAg_JC(x)= nk_import("Ag(JC-eV)", x)
nAg_P(x)= nk_import("Ag(Palik-eV)", x)
nAg_B(x) = nk_import("Ag(Barber-eV)", x)
nAg_H(x) = nk_import("Ag(Hageman-eV)", x)
nAg_Wi(x) = nk_import("Ag(film)", x)
nAg_We(x) = nk_import("Ag(Werner-eV)", x)
nAg_S(x) = nk_import("Ag(sopra-eV)", x)
nSi(x) = nk_import("Si111(sopra)", x)
nSiO2(x) = nk_import("SiO2(sopra)", x)
nPMMA(x) = nk_import("PMMA", x)



mutable struct Lattice2D
    a::Real
    phi::Real
    A::Real
    function Lattice2D(a, phi)
        A = a*a*sin(phi)
        return new(a, phi, A)
    end
end


function make_hole(nlayers::Integer, mat_hole::Material, mat_layer::Material, thickness::Real, radius::Real, lattice::Lattice2D)
    h = thickness/nlayers
    ff = pi*radius^2/lattice.A
    eps_eff = ff*mat_hole.eps +(1.0 - ff)*mat_layer.eps
    mat = Material("epsmu", eps_eff)
    Ls = Vector{Layer}(nlayers)
    for i = 1:nlayers
        Ls[i] = Layer("c", mat, h)
    end
    return Ls
end

function make_cone(nlayers::Integer, mat_cone::Material, mat_host::Material, height::Real, rbot::Real, rtop::Real, lattice::Lattice2D)
    h = height/nlayers
    Ls = Vector{Layer}(nlayers)
    for i = 0:(nlayers-1)
        ri = rtop + i*(rbot - rtop)/nlayers
        ff = pi*ri^2/lattice.A
        eps_eff = ff*mat_cone.eps +(1.0 - ff)*mat_host.eps
        mat = Material("epsmu", eps_eff)
        Ls[i+1] = Layer("c", mat, h)
    end
    return Ls
end

function make_cone_vol(mat_cone::Material, mat_host::Material, height::Real, rbot::Real, rtop::Real, lattice::Lattice2D)
    Vcone = pi*(rtop^2 + rbot*rtop + rbot^2 )*height/3.0
    Vcell = height*lattice.A
    ff = Vcone/Vcell
    eps_eff = ff*mat_cone.eps +(1.0 - ff)*mat_host.eps
    mat = Material("epsmu", eps_eff)
    return Layer("c", mat, height)
end

mat1 = Material()
mat2 = Material("nk", 1.5)

make_hole(5, mat1, mat2, 100, 150, Lattice2D(500, pi/2))

make_cone(10, mat2, mat1, 150, 200, 20, Lattice2D(500, pi/2) )
make_cone_vol( mat2, mat1, 150, 200, 20, Lattice2D(500, pi/2))


nwl = 500
wl = linspace(200, 900, nwl)
R = Array{Float64}(nwl)
T = Array{Float64}(nwl)
Ri = Array{Float64}(nwl)
Ti = Array{Float64}(nwl)
mlat = Lattice2D(500, pi/2)
for i = 1:nwl
    matSiO2 =  Material("nk", nSiO2(wl[i]))
    matPMMA = Material("nk", 6)
    L0  = Layer()
    Lhole =  make_hole(1, Material(), matPMMA, 40, 10, mlat)
    L1 = Layer("c", matPMMA, 100)
    L2 = Layer("i", matSiO2, 1e6)
    Ls = [L0; Lhole; L2; L0]
    S = Stack(Ls, zeros(length(Ls)-1))
    R[i], T[i] = tmm_RT(1, wl[i], zeros(2), S)
    Ri[i], Ti[i] = RT_matrix_inc(1, wl[i], zeros(2), S)
end

begin
    plot(wl, T, label = "T")
    plot(wl, R, label = "R")
    plot(wl, Ti, label = "Tinc")
    plot(wl, Ri, label = "Rinc")
    xlim([300,  900])
    ylim([0, 1])
    xlabel("Wavelength (nm)")
    ylabel("Transmittance")
    legend()
end
