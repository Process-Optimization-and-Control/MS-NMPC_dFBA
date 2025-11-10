#=
      simulating the closed loop performance of a NMPC
=#

#−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
#Formulating the discrete dynamics
using JuMP
using Ipopt
using Cbc
using LinearAlgebra
#using DifferentialEquations#, ForwardDiff, DiffEqSensitivity
#using OrdinaryDiffEq
using Sundials
using Plots
using DelimitedFiles
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
nc = 5
nmv = 1 # number of manipulated variables
ns = 2   #number of scenarios

include("FBA.jl")
include("FVA.jl")
include("FVA_min.jl")
include("pFBA_KKT_flux.jl")

# FVA
ac_id = readdlm("ac_id_600.csv", ',')
FVA_ind_ini = "min"       #no, max or min
sub_FVA = 0.01
FVA_ac = ac_id[10,:]


#...Model fixed Parameters
Kg = 1.0        # g/L ok
Ki=100.0  # g/L   ok
Kcb= 100.0  # g/L         overflow switch on high glucose concentrations
Ko = 0.003 # mmol/L
Sf = 100.0   # g/L ok

vg_max= 10.0 #ok
vo_max= 15.0 #ok
vlb[o2]= -vo_max
vlb[ac] = 0.0
vlb[ATP] = 10.0

vlb_F = deepcopy(vlb)
vub_F = deepcopy(vub)


#initial conditions   
b0=0.1       # g/L                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
g0=10.0        # g/L  
ac0=0.0         # g/L 
O0=0.21         # mmol
v0=1.0
F0 = 0.0         # L/h
klo = 7.5  # 1/h 
O_sat = 0.21 # mmol 

#Process contraints
u_n = 0.0         # L/h
u_max_rate = 0.05 
u_max = 0.5  #0.5
u_min = 0.0
V_max = 3.0 # 3.0
Ac_max = 5.0 # 4.0

#flux rate change constraint inside the elements
v_dot = 0.1 

d=0.0*Vector{Float64}(undef,nv) 
d[obj]=-1
up=0.0*Vector{Float64}(undef,nv) 
up[glu] = 1
up2=0.0*Vector{Float64}(undef,nv) 
up2[o2] = 1
v0_ini=zeros(nv)


#scaling
cs=Vector{Float64}(undef,nc) 
cs[1]=1.0
cs[2]=1.0
cs[3]=1.0
cs[4]=1.0
cs[5]=1.0

vs=Vector{Float64}(undef,nv) 
for i in 1:nv

    if vlb[i] == 0.0 && vub[i] == 0.0

         vs[i] = 1.0
    else
        #vs[i] = max(abs(vlb[i]),abs(vub[i]))

        vs[i] = 1.0
      
    end

end


#integration parameters
nc = 5 # number of ODEs
nfe    =  600     # number of control intervals
nfe_r =  10          # finite elements into the controller
th_r = 0.1666           #Prediction horizon 
ncp    = 3           # number of collocation points
th     = 12.0        # time horizon
h = th/nfe           # length of one finite element on the time horizon
h_r = th_r/nfe_r 
hm=h_r*ones(nfe)'
var_h=0.1


w = Matrix{Float64}(undef,nv,ns)  
w[ac,1] = -1e-1   # max -1e-4
w[ac,2] = 1e0  # min 2.55e-4 
phi = 1.0
phi_all = 1e-6
phi1=phi_all
phi2=phi_all
phi3=phi_all
phi4 = 1e-4  # slacks


ts1     = Vector{Float64}(undef,nfe+1) # time series for plotting
xk1 = Matrix{Float64}(undef,nc,nfe+1)
vk1 = Matrix{Float64}(undef,4,nfe+1)
flags = Vector{String}(undef,nfe) # time series for plotting
time_opt = Vector{Float64}(undef,nfe) # time series for plotting
var_opt = Vector{Float64}(undef,nfe) # time series for plotting
u = ones(2,nfe+1) #F and klo

step_index=Int(nfe/2)
F0=0.0
kl0=klo
u[1,1:nfe-step_index]=F0*u[1,1:nfe-step_index]
u[2,1:nfe-step_index]=klo*u[2,1:nfe-step_index]

Fs=0.0
klos=klo
u[1,nfe-step_index+1:end]=Fs*u[1,nfe-step_index+1:end]
u[2,nfe-step_index+1:end]=klos*u[2,nfe-step_index+1:end]

       
          function f(xdot,x,p,t)
            

            #show("x1")
            #println(x[1])
            vg = vg_max*  (x[1])/( Kg + (x[1]) + ((x[1])^2/Ki) )
            vo = vo_max*  (x[4])/( Ko + (x[4]) )
  
          
            if x[1] < 0.1
                vg=0.0
                vo=0.0
                mi=0.0
                vac=0.0

            else
               mi,vac, vg, vo = FBA(vg,vo)
               
               if p[3] == "max"
                  mi,vac, vg, vo = FVA(mi,-vg,-vo)
               elseif p[3] == "min"
                  mi,vac, vg, vo = FVA_min(mi,-vg,-vo)
               else 
               end
               

            end
           
      

            xdot[1] = 0.18*vg*x[3] + (Sf - x[1])*(p[1]/x[5])
            xdot[2] = 0.06*vac*x[3] - x[2]*(p[1]/x[5])
            xdot[3] = mi*x[3] - x[3]*(p[1]/x[5])
            xdot[4] = vo*x[4]*x[5] + p[2]*(O_sat - x[4])
            xdot[5] = p[1]

          end
        

        x0 = [g0,ac0,b0,O0,v0]

        
      
        xk1[:,1]=x0
        ts1[1] = 0.0
        # Initial conditions
        vg = vg_max*  (xk1[1,1])/( Kg + (xk1[1,1]) + ((xk1[1,1])^2/Ki) )
        vo = vo_max*  (xk1[4,1])/( Ko + (xk1[4,1])  )
        mi,vac, vg, vo = FBA(vg,vo)
        if FVA_ind_ini == "max"
         mi,vac, vg, vo = FVA(mi,-vg,-vo)
        elseif FVA_ind_ini == "min"
         mi,vac, vg, vo = FVA_min(mi,-vg,-vo)
        else 
        end
        vk1[:,1] = [mi,vac, vg, vo]

       # initializing the NMPC solution
       u_ini, x_ini, v_ini, solf, soltime, solvar = pFBA_KKT_flux(x0,x0,0.1,v0_ini,u_n,nfe_r,th_r)

