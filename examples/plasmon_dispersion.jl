using MultiTMM
using DelimitedFiles
using MAT
using PyPlot
pygui(true)
matplotlib["rcParams"][:update](["font.size" => 18, "font.weight"=>"normal", "font.family" => "Arial", "text.usetex"=>true])


function drude_pole(w, epsb, wp, gam)
    return epsb - wp^2/(w*(w+im*gam))
end




function rp_sigma(q, w, wp, gam, d, е1, е2)
    еа = 0.5*(е1 + е2)
    sigw = (1.0im*d/(4*pi*еа))*wp^2/(w*(w+im*gam))
    r = 1.0/(1 - 0.5im/(pi*q*sigw))
    return r
end





begin
    nw = 500
    nq = 500
    NML = 20
    tML = 0.40853/sqrt(3)
    d = tML*NML
    wp = 9.17
    gp = 0.021
    qq = range(0.0000000001, pi/d, length = nq);
    ww = range(0.0001, 5, length = nw);

    rp = Array{ComplexF64}(undef, nw, nq)

    for j = 1:nq
        for i = 1:nw
            rp[j, i] = rp_sigma(qq[j], ww[i], wp, gp, NML*tML, 2.0, 12.0)
        end
    end
end

begin
    extent = [qq[1]*d, qq[end]*d, ww[1]/9.17, ww[end]/9.17]
    imshow(log.(abs.(rp')), vmin = 0, vmax = 5, origin = "lower", aspect = "auto", cmap = "Blues", extent = extent)
    xlabel("\$k_{\\parallel} d\$")
    ylabel("\$\\omega/\\omega_{p}\$")
    tight_layout()
    savefig("C:\\Users\\vmkhi\\Documents\\Projects\\Conferences\\ICTON2019\\Random Pictures\\Dispersion.png")
end
