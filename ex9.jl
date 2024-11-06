#!/usr/local/bin/julia

using BenchmarkTools

include("genreach.jl")

@variables x[1:4]
QE9=([2.0+2.0*x[1]+x[2]+3.0*x[3]+x[4], -1.0-x[1]-x[2]+x[3]+5.0*x[4]],["exists", 1, "forall", 2, "exists", 3, "exists", 4],4,2)
print_QE(QE9)

println("Results for ex9 - order 1")
@btime (global res_QE9=QEapprox_o1(QE9[1], QE9[2], QE9[3], QE9[4]))
println(res_QE9)
print_vol(res_QE9)
println

println("Results for ex9 - order 0")
@btime (global res_QE9_o0=QEapprox_o0(QE9[1], QE9[2], [["exists", 1, "forall", 2, "forall", 4, "exists", 3],
                                        ["forall", 1, "forall", 2, "forall", 3, "exists", 4]], QE9[3], QE9[4]))
println(res_QE9_o0)
print_vol_o0(res_QE9_o0, 2)
print_rel_vol_o1_o0(res_QE9, res_QE9_o0, 2)
println