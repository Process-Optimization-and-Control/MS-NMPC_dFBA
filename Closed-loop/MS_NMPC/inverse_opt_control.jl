# inverse optimal control
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
using BlackBoxOptim
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

v_max = readdlm("v_max_1.csv") 
v_max_OF = readdlm("v_max_1_OF.csv") 

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
w_max[ac] = 0.0
w_min[ac] = 0.0
d=0.0*Vector{Float64}(undef,nv) 
d[obj]=-1

#FVA
sub_FVA = 0.01  # sub-optimality

vo_val=[15 14 13 12 11 10 9 8 7 6 5 3 2 1 0] 
vg_val = [10 9.5 9 8.5 8 7.5 7 6.5 6 5.5 5.2 5 4.5 4 3.5]


v_FBA = Matrix{Float64}(undef,15,15)
v_max = Matrix{Float64}(undef,15,15)
vk_max = Matrix{Float64}(undef,15,15)
vk_max_OF = Matrix{Float64}(undef,15,15)
flag_k_max = Matrix{String}(undef,15,15)
v_min = Matrix{Float64}(undef,15,15)
vk_min = Matrix{Float64}(undef,15,15)
vk_FBA = Matrix{Float64}(undef,15,15)
flag_k_min = Matrix{String}(undef,15,1615)


function min_inv(w_ac)
    #w_max[ac] = -w_ac

    setindex!(w_max, -w_ac[1], ac)
    
    for i in 1:15
    for j in 1:15
    vg = vg_val[i]    
    vo = vo_val[j]  


    mi,vac, flag_max= KKT(vg,vo,w_max)
    vk_max[i,j]= vac
    vk_max_OF[i,j] = mi

    end
    end

    residual1 = norm(vk_max - v_max)
    residual2 = norm(vk_max_OF - v_max_OF)
    residual = residual1 + residual2
    return residual
end

res = bboptimize(min_inv; SearchRange = (1e-8, 1e-2), NumDimensions = 1,MaxFuncEvals = 50)