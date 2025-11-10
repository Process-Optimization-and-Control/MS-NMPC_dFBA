using JuMP
using Ipopt, Cbc
using LinearAlgebra
#using Plots
#reading GSM model
using DelimitedFiles


include("FBA.jl")
include("FVA.jl")
include("FVA_min.jl")
#include("pFBA.jl")
#include("KKT.jl")
#include("KKT_r.jl")
#include("KKT_MPEC.jl")
#include("KKT_MPEC2.jl")


#stoichiometric matrix
S = readdlm("S.csv", ',');
#Sa,Vd,Dd = svd(S'*S)
#K = Dd[:,68:95]
#K = nullspace(S)
#ns = size(K,2)

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

nc = 3
nv=size(S,2)
nm=size(S,1)
#fluxes constraints
vlb2 = readdlm("lb.csv", ',');
vlb=vlb2[:,1]
vub2 = readdlm("ub.csv", ',');
vub=vub2[:,1]

#...Model fixed Parameters
#Kg = 0.015         #
#Ko=0.01
#vg_max= 10.5
vlb[o2]= -10.0
#vlb[ac] =-2.5
#vlb[ATP] = 0.0

#vlb[glu] = 0.0
vg=10.0
#vub[glu] = 0.0
#d=0.0*Vector{Float64}(undef,nv) 
#d[obj]=-1
#d[obj]=0
#phi=0.1
#w=1e-12

#relx = 1e-5
#relx2 = 1e-8

#initial conditions
#v0=0.0*Vector{Float64}(undef,nv) 
#alStar=0.0*Vector{Float64}(undef,nv) 
#lStar=0.0*Vector{Float64}(undef,nm) 

mi,vac = FBA(vg)
#mi,vac = pFBA(vg,mi)
#vStar2, lStar, alStar = KKT(v0)
#vStar3, alStar, auStar = KKT_r(vall)
#vStar2 = K*vStar3
#vStar, lStar, alStar = KKT_MPEC(v0,relx)
#vStar, lStar, alStar = KKT_MPEC2(v0,relx2,vStar, lStar, alStar)

println([mi vac])
#println(vStar2[obj])
#println(vStar2[ac])
#println(vStar[obj])
#println(vStar[ac])

#FVA
mi,vac_max = FVA(mi)
mi,vac_min = FVA_min(mi)

println([vac vac_max vac_min])
#[vStar[1:10] vStar[11:20] vStar[21:30] vStar[31:40] vStar[41:50] vStar[51:60] vStar[61:70] vStar[71:80]]
#[vStar[81:88] vStar[88:end] ]

#[vlb[1:10] vlb[11:20] vlb[21:30] vlb[31:40] vlb[41:50] vlb[51:60] vlb[61:70] vlb[71:80]]
#[vlb[81:88] vlb[88:end] ]

#[alStar[1:10] alStar[11:20] alStar[21:30] alStar[31:40] alStar[41:50] alStar[51:60] alStar[61:70] alStar[71:80]]
#[alStar[81:88] alStar[88:end] ]
 
# [auStar[1:10] auStar[11:20] auStar[21:30] auStar[31:40] auStar[41:50] auStar[51:60] auStar[61:70] auStar[71:80]]
#[auStar[81:88] auStar[88:end] ]

#[lStar[1:10] lStar[11:20] lStar[21:30] lStar[31:40] lStar[41:50] lStar[51:60] lStar[61:70]]
#lStar[71:72]