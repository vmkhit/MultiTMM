using DataFrames
using DelimitedFiles
using Interpolations
using Dierckx
path_eps = "../eps_data"
function nk_import(fname, x)
    for f in filter(x -> ismatch(r"\.csv|\.txt", x), readdir(path_eps))
        for f in filter(x -> startswith(x, fname), readdir(path_eps))
            n_data = readdlm(joinpath(path_eps, f))
            n_re = Spline1D(n_data[:,1], n_data[:, 2])
            n_im = Spline1D(n_data[:, 1], n_data[:, 3])
            return complex(n_re(x), n_im(x))
        end
    end
end


#nk_import("PMMA", 800)


function eps_Au(w)
    hartree = 27.2116;              #  2 * Rydberg in eV
    tunit = 0.66 / hartree         #  time unit in fs
    rs = 3                   #  electron gas parameter
    eps0 = 10                #  background dielectric constant
    gammad = tunit/10         #  Drude relaxation rate
    #  density in atomic units
    density = 3/(4*pi *rs^3)
    #  plasmon energy
    wp = sqrt(4*pi*density)

    gammad = gammad*hartree
    wp = wp*hartree
    eps = eps0 - wp^2 ./(w.*(w + 1im*gammad))
    return wp, gammad, eps
end
