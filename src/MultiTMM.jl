module MultiTMM

include("tmm_matrix.jl")
include("eps.jl")

export nk_import
export epston, Material, Interface, Layer, Stack
export intface_rt, tmm_matrix, tmm_RT, tmm_ellipso, RT_matrix_inc

end
