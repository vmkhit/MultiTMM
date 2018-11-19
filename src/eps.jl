using DataFrames
using DelimitedFiles
using Interpolations
using Dierckx
path_eps = "C:\\Users\\vmkhi\\Documents\\Github\\MultiTMM\\eps_data"
function nk_import(filename, x; path = path_eps)
     file = filter(f-> startswith(f, filename) && endswith(f, ".txt"), readdir(eps_path))
     if !isempty(file)
         data = readdlm(joinpath(eps_path, file[1]))
         if size(data, 2) < 2
             throw("The data file must have at least two columns")
         elseif size(data, 2) == 2
             n_re = Spline1D(data[:,1], data[:, 2])
             return complex.(n_re(x), 0.0)
         else
             n_re = Spline1D(data[:,1], data[:, 2])
             n_im = Spline1D(data[:, 1], data[:, 3])
             return complex.(n_re(x), n_im(x))
         end
     else
         throw("The specified file is not in the directory")
     end
end

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
