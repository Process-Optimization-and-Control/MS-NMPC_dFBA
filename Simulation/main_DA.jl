#=
      simulation of dFBA by the direct approach
=#

#−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
#Formulating the discrete dynamics
using JuMP
using Cbc
using LinearAlgebra
#using DifferentialEquations#, ForwardDiff, DiffEqSensitivity
#using OrdinaryDiffEq
using Sundials
using Plots
using DelimitedFiles
using TickTock


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

include("FBA.jl")
#include("pFBA.jl")
include("FVA.jl")
include("FVA_min.jl")


#...Model fixed Parameters
Kg = 0.015         #
Ko=0.01
vg_max= 10.5
vo_max= 19.0
vlb[o2]= -19.0
vlb[ac] = 0.0
#vlb[ac] = -2.5
vlb[ATP] = 10.0


#initial conditions
b0=0.003        #                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
g0=10.4        #
ac0=0.3         #
O0=0.21         #
#v0=[0 1 0 0 9.46 12.92 0 1]
v0=0.0*Vector{Float64}(undef,nv) 


#integration parameters
nc = 3
nfe    = 3      # number of control intervals
th     = 6.7         # time horizon
#th     = 9.0    # time horizon
h = th/nfe           # length of one finite element on the time horizon


ts     = Vector{Float64}(undef,nfe) # time series for plotting
xkf = Matrix{Float64}(undef,nc,nfe+1)

for i in 1:nfe
    ts[i] = h*i
end
ts = pushfirst!(ts,0)

       
          function f(xdot,x,p,t)
            
            vg = vg_max*  (x[1])/(Kg+(x[1]))
          
            if x[1] < 0.0001
                vg=0
            else
            end

            mi,vac = FBA(vg)
            # println(mi)
            # println(vac)
           # mi,vac = pFBA(vg,mi)
           #FVA
            # mi,vac = FVA(mi,vg)
            mi,vac = FVA_min(mi,vg)

            if x[2] < 0.0001
                vac=0
                mi=0
            else
            end

            xdot[1] = -vg*x[3]
            xdot[2] = vac*x[3]
            xdot[3] = mi*x[3]

          end
        

        #...Time
        #t = (0.0,h)
        #t = (0.0,1.0)

        x0 = [g0,ac0,b0]
      
        xkf[:,1]=x0
tick()
   # for i in 1:nfe
    #    t = (0.0,h)
    #    println(i)

        
        #vo = vo_max*(xkf[2,i])/(Ko+(xkf[2,i]))
        
        
      


       # if xkf[1,i] < 0.01 && xkf[2,i] < 0.01
      #    xkf[:,i+1] = xkf[:,i]
      #  else
           u=[]
           t=(0,th)
           prob = ODEProblem(f,xkf[:,1],t,u)
           #sol = DifferentialEquations.solve(prob,Tsit5(),saveat=0.1)
           sol = Sundials.solve(prob,CVODE_Adams(),saveat=0.1)
           ts1=sol.t
           xk1=sol.u
           #xkf[:,i+1] = sol.u[end]
       # end
  
  #end
tock()
writedlm("xk_DA.csv",xk1)
writedlm("ts_DA.csv",ts1)  

      # return xk
#end

xk = readdlm("xk_DA.csv")
ts = readdlm("ts_DA.csv")


# p11 = plot(ts,xk[:,3],linewidth=2,xaxis="time [h]",yaxis="biomass [g]",legend=false)
# p11 = plot!(ts_max,xk_max[:,3],linewidth=2,xaxis="time [h]",yaxis="biomass [g]",legend=false)
# #C_C
# p12 = plot(ts,xk[:,2],linewidth=2,xaxis="time [h]",yaxis="Ac [mM]",legend=false)
# p12 = plot!(ts_max,xk_max[:,2],linewidth=2,xaxis="time [h]",yaxis="Ac [mM]",legend=false)

# # #C_D
# p14 = plot(ts,xk[:,1],linewidth=2,xaxis="time [h]",yaxis="Glucose [mM]",legend=false)
# p14 = plot!(ts_max,xk_max[:,1],linewidth=2,xaxis="time [h]",yaxis="Glucose [mM]",legend=false)

# g1 = plot(p11,p12,p14,layout=(3,1),reuse = true)
# #display(g1)

# savefig(g1,"test.pdf")

