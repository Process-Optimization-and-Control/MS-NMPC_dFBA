# FBA_solution_space
#Formulating the discrete dynamics
using JuMP
using Cbc
using LinearAlgebra
#using DifferentialEquations#, ForwardDiff, DiffEqSensitivity
#using OrdinaryDiffEq
using Sundials
using Plots
using DelimitedFiles
using Ipopt
#using TickTock
ENV["GKSwstype"] = "100"

S = readdlm("S.csv", ',');
eth=48  #ethanol reaction
obj=25  #objective function index
glu=52  #glucose uptake
o2=60  #oxygen uptake
ac=44
ATP=16

#E.coli iJR904
#obj=269  #objective function index
#glu=398  #glucose uptake
#o2=500  #oxygen uptake
#ac=289
#ATP=212

nv=size(S,2)
nm=size(S,1)
#c = readdlm("c.csv", ',');
#b = readdlm("b.csv", ',');
vlb = readdlm("lb.csv", ',');
vub = readdlm("ub.csv", ',');
vlb_max=vlb
vub_max=vub


include("FBA.jl")
#include("pFBA.jl")
include("FVA.jl")
include("FVA_min.jl")
include("KKT.jl")

vg_max= 10.0 #ok
vo_max= 15.0 #ok
vlb[o2]= -vo_max
vlb[ac] = 0.0
#vlb[ATP] = 10.0

# KKT
w_max = Vector{Float64}(undef,nv)  
w_min = Vector{Float64}(undef,nv)  
w_FBA = Vector{Float64}(undef,nv)  
w_max[ac] = -1e-4
w_min[ac] = 2.55e-4
d=0.0*Vector{Float64}(undef,nv) 
d[obj]=-1

#FVA
sub_FVA = 0.00  # sub-optimality

vo_val=[15 14 13 12 11 10 9 8 7 6 5 3 2 1 0] 
vg_val = [10 9.5 9 8.5 8 7.5 7 6.5 6 5.5 5.2 5 4.5 4 3.5]


v_FBA = Matrix{Float64}(undef,15,15)
v_max = Matrix{Float64}(undef,15,15)
v_max_OF = Matrix{Float64}(undef,15,15)
vk_max = Matrix{Float64}(undef,15,15)
vk_max_OF = Matrix{Float64}(undef,15,15)
flag_k_max = Matrix{String}(undef,15,15)
v_min = Matrix{Float64}(undef,15,15)
vk_min = Matrix{Float64}(undef,15,15)
vk_min_OF = Matrix{Float64}(undef,15,15)
vk_FBA = Matrix{Float64}(undef,15,15)
flag_k_min = Matrix{String}(undef,15,15)
v_min_OF = Matrix{Float64}(undef,15,15)

for i in 1:15
for j in 1:15
vg = vg_val[i]    
vo = vo_val[j]  

#println(vg,vo)

mi, vac = FBA(vg,vo)
v_FBA[i,j]= mi

mi,vac = FVA(mi,vg,vo)
v_max[i,j]= vac
#v_max_OF[i,j]= mi

mi,vac = FVA_min(mi,vg,vo)
v_min[i,j]= vac
#v_min_OF[i,j]= mi

#mi,vac = KKT(vg,vo,w_FBA)
#vk_FBA[i,j]= mi

#mi,vac, flag_max= KKT(vg,vo,w_max)
#vk_max[i,j]= vac
#vk_max_OF[i,j] = mi
#flag_k_max[i,j] =  string(flag_max)

#mi,vac, flag_min= KKT(vg,vo,w_min)
#vk_min[i,j]= vac
#vk_min_OF[i,j] = mi
#flag_k_min[i,j] =  string(flag_min)


#mi,vac = FBA(vg,vo)
#mi,vac = FVA_min(mi,vg,vo)
#v_min[i,j]= vac
end
end

writedlm("vg_val.csv",vg_val)
writedlm("vo_val.csv",vo_val)  
#writedlm("v_FBA.csv",v_FBA)  
writedlm("v_max_0.csv",v_max) 
#writedlm("v_max_1_OF.csv",v_max_OF) 
#writedlm("vk_max.csv",vk_max)  
#writedlm("vk_max_OF.csv",vk_max_OF)  
writedlm("v_min_0.csv",v_min)  
#writedlm("vk_min.csv",vk_min)  
#writedlm("v_min_1_OF.csv",v_min_OF)  
#writedlm("vk_min_OF.csv",vk_min_OF)  
#writedlm("vk_FBA.csv",vk_FBA)  