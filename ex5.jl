#!/usr/local/bin/julia

using BenchmarkTools

include("genreach.jl")

@variables x[1:3]
QE5=([x[1]^2/4.0+(x[2]+1.0)*(x[3]+2.0)+(x[3]+3.0)^2],["exists", 1, "forall", 2, "exists", 3],3,1)
print_QE(QE5)

println("Results for ex5 - order 1")
println("for ex5 -order 1")
@btime (global res_QE5=QEapprox_o1(QE5[1], QE5[2], QE5[3], QE5[4]))
println(res_QE5)
print_vol(res_QE5)

println("Results for ex5 -order 0")
@btime (global res_QE5_o0=QEapprox_o0(QE5[1], QE5[2], [["exists", 1, "forall", 2, "exists", 3]], QE5[3], QE5[4]))
println(res_QE5_o0)
print_vol_o0(res_QE5_o0,1)
print_rel_vol_o1_o0(res_QE5, res_QE5_o0, 1)
println