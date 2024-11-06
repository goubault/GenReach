#!/usr/local/bin/julia

using BenchmarkTools

include("genreach.jl")

@variables x[1:4]
QE6=([0.1*x[1]+(1+0.01*x[2])*(x[4]+1)/4+0.000000131*x[3]*((x[4]+1)/4)^2],["exists", 1, "forall", 2, "exists", 3, "exists", 4],4,1)
print_QE(QE6)

println("Results for ex6 - order 1")
@btime (global res_QE6=QEapprox_o1(QE6[1], QE6[2], QE6[3], QE6[4]))
println(res_QE6)
print_vol(res_QE6)

println("Results for ex6 -order 0")
@btime (global res_QE6_o0=QEapprox_o0(QE6[1], QE6[2], [["exists", 1, "forall", 2, "exists", 3, "exists", 4]], QE6[3], QE6[4]))
println(res_QE6_o0)
print_vol_o0(res_QE6_o0,1)
print_rel_vol_o1_o0(res_QE6, res_QE6_o0, 1)
println