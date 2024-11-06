#!/usr/local/bin/julia

using BenchmarkTools
using Random
using Distributions

include("genreach.jl")

@variables x[1:2000]

# Examples "linear", 1D
steps = 10 # 1000
p = 2*steps
n = 1

V = rand(Uniform(-1.0,1.0),p+1)

g = [sum(V[i]*x[i] for i in 1:p)+V[p+1]]

# encoding of the list of quantifiers
# forall x1, exists x2, forall x3, exists x4...

quantifiers = []
for i in 1:steps
  push!(quantifiers, "forall")
  push!(quantifiers, 2*i-1)
  push!(quantifiers, "exists")
  push!(quantifiers, 2*i)
end

q = [quantifiers]

QElin = (g, quantifiers, p, n)
print_QE(QElin)

print("Results for linear ")
print(steps)
println(" - order 0")
@btime (global QElin_o0=QEapprox_o0(g,quantifiers,[quantifiers],p,n))
println(QElin_o0)
print_vol_o0(QElin_o0,n)
println

print("Results for linear ")
print(steps)
println(" - order 1")
@btime (global QElin_o1=QEapprox_o0(g,quantifiers,[quantifiers],p,n))
println(QElin_o1)
print_vol(QElin_o1)
print_rel_vol_o1_o0(QElin_o1,QElin_o0,n)
println