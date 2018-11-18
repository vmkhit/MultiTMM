using MultiTMM
using BenchmarkTools
lam = linspace(500, 2000, 1000)

M1 = Material("nk", 10.0)
M2 = Material("nk", 1.4)
I12 = Interface(M1, M2)

r = Array{ComplexF32}(undef, length(lam))
t = similar(r)

for (i, lami) =enumerate(lam)
    r[i], t[i] = intface_rt(1, I12, 2*pi/lami, zeros(2))
end

@benchmark intface_rt(1, I12, 2*pi/500, zeros(2))
