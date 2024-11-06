#!/usr/local/bin/julia

using BenchmarkTools

include("genreach.jl")

@variables x[1:8]
QE6xy=([0.1*x[1]+(1+0.01*x[2])*(x[8]+1)/4+0.000000131*x[3]*((x[8]+1)/4)^2, 
                                 0.1*x[4]+0.25*(0.01*x[6]+0.0025*x[7]*(1+x[8]))*(1+x[8])+0.005*x[5]*(0.25*(1+x[8]))^2],
                                 ["exists", 7, "exists", 1, "exists", 4, "exists", 6, "forall", 2, "exists", 3, "exists", 5, "exists", 8],8,2)
print_QE(QE6xy)

println("Results for ex6xy - order 1")
@btime (global res_QE6xy=QEapprox_o1(QE6xy[1], QE6xy[2], QE6xy[3], QE6xy[4]))
println(res_QE6xy)
print_vol(res_QE6xy)

println("Results for ex6xy -order 0")
@btime (global res_QE6xy_o0=QEapprox_o0(QE6xy[1], QE6xy[2], 
        [["forall", 7, "forall", 4, "forall", 6, "exists", 1, "forall", 2, "forall", 5, "forall", 8, "exists", 3],
         ["forall", 1, "exists", 7, "exists", 4, "exists", 6, "forall", 2, "forall", 3, "exists", 5, "exists", 8]], QE6xy[3], QE6xy[4]))
println(res_QE6xy_o0)
print_vol_o0(res_QE6xy_o0,1)
print_rel_vol_o1_o0(res_QE6xy, res_QE6xy_o0, 1)
println