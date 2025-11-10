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
#using JLD
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

u_optimal = readdlm("uk.csv")
nfe = readdlm("nfe.csv")
nfe = Int(reshape(nfe, 1)[1])
th = readdlm("th.csv")
th = reshape(th, 1)[1]

ac_id = readdlm("ac_id.csv", ',')

include("FBA.jl")
include("FVA.jl")
include("FVA_min.jl")


# FVA
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


#initial conditions   
b0=0.1       # g/L                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
g0=10.0        # g/L  
ac0=0.0         # g/L 
O0=0.21         # mmol
v0=1.0
F0 = 0.0         # L/h
klo = 7.5  # 1/h 
#klo = 0.8*7.5  # 1/h
O_sat = 0.21 # mmol 

#integration parameters
nc = 5 # number of ODEs
#nfe    = 30      # number of control intervals
#th     = 12.0        # time horizon
h = th/nfe           # length of one finite element on the time horizon


ts1     = Vector{Float64}(undef,nfe+1) # time series for plotting
xk1 = Matrix{Float64}(undef,nc,nfe+1)
vk1 = Matrix{Float64}(undef,4,nfe+1)

u = ones(2,nfe+1) #F and klo
u[1,:]= u_optimal
step_index=Int(nfe/2)
F0=0.0
kl0=klo
#u[1,1:nfe-step_index]=F0*u[1,1:nfe-step_index]
u[2,1:nfe-step_index]=klo*u[2,1:nfe-step_index]

Fs=0.0
klos=klo
#u[1,nfe-step_index+1:end]=Fs*u[1,nfe-step_index+1:end]
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
#tick()

    for i in 1:nfe
    #    t = (0.0,h)
    #    println(i)
           
       if FVA_ac[i] == 1
         FVA_ind = "max"
       else 
         FVA_ind = "min"
       end
       
        
           t=(0,h)
           prob = ODEProblem(f,xk1[:,i],t,[u[:,i];FVA_ind])
           #sol = DifferentialEquations.solve(prob,Tsit5(),saveat=0.1)
           sol = Sundials.solve(prob,CVODE_Adams(),saveat=0.1)
           #ts1=sol.t
           #xk1=sol.u
           ts1[i+1]=ts1[i]+h
           xk1[:,i+1] = sol.u[end]

          # Generating data for plotting
           if xk1[1,i+1] < 0.1
            vg=0.0
            vo=0.0
            mi=0.0
            vac=0.0
           else
            vg = vg_max*  (xk1[1,i+1])/( Kg + (xk1[1,i+1]) + ((xk1[1,i+1])^2/Ki) )
            vo = vo_max*  (xk1[4,i+1])/( Ko + (xk1[4,i+1])   )
            mi,vac, vg, vo = FBA(vg,vo)
             if FVA_ind == "max"
                mi,vac, vg, vo = FVA(mi,-vg,-vo)
             elseif FVA_ind == "min"
                mi,vac, vg, vo = FVA_min(mi,-vg,-vo)
             else 
             end
           end

           vk1[:,i+1] = [mi,vac, vg, vo]
    end
  
  #end
#tock()
writedlm("xk_DA.csv",xk1)
writedlm("ts_DA.csv",ts1)  
writedlm("vk_DA.csv",vk1)  

      # return xk
#end

ts = readdlm("ts_DA.csv")
xk = readdlm("xk_DA.csv")
xk = xk'
vk = readdlm("vk_DA.csv")
vk = vk'

 p11 = plot(ts,xk[:,3],linewidth=2,xaxis="time [h]",yaxis="biomass [g/L]",legend=false)
 #p11 = plot!(ts_min,xk_min[:,3],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="biomass [g/L]",label="Min")
 #p11 = plot!(ts_max,xk_max[:,3],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="biomass [g/L]",label="Max")

 p12 = plot(ts,xk[:,2],linewidth=2,xaxis="time [h]",yaxis="Ac [g/L]",legend=false)
 #p12 = plot!(ts_min,xk_min[:,2],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Ac [g/L]",label="Min")
 #p12 = plot!(ts_max,xk_max[:,2],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Ac [g/L]",label="Max")

 p14 = plot(ts,xk[:,1],linewidth=2,xaxis="time [h]",yaxis="Glucose [g/L]",legend=false)
 #p14 = plot!(ts_min,xk_min[:,1],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Glucose [g/L]",label="Min")
 #p14 = plot!(ts_max,xk_max[:,1],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max")

 p13 = plot(ts,xk[:,4],linewidth=2,xaxis="time [h]",yaxis="Oxygen [mmol/L]",legend=false)
 #p13 = plot!(ts_min,xk_min[:,4],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Min")
 #p13 = plot!(ts_max,xk_max[:,4],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max")


 g1 = plot(p11,p12,p14,p13,layout=(2,2),reuse = true)


 p21 = plot(ts,u[1,:],linewidth=2,xaxis="time [h]",yaxis="F [L/h]", linetype=:steppre,legend=false)
 p22 = plot(ts,u[2,:],linewidth=2,xaxis="time [h]",yaxis="Klo [1/h]",legend=false)
 p23 = plot(ts,xk[:,5],linewidth=2,xaxis="time [h]",yaxis="Volume [L]",legend=false)
 g2 = plot(p21,p22,p23,layout=(3,1),reuse = true)

 p31 = plot(ts,vk[:,1],linewidth=2,xaxis="time [h]",yaxis="Growth rate [1/h]",legend=false)
 #p31 = plot!(ts_min,vk_min[:,1],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Growth rate [1/h]",label="Min")
 #p31 = plot!(ts_max,vk_max[:,1],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Growth rate [1/h]",label="Max")

 p32 = plot(ts,vk[:,2],linewidth=2,xaxis="time [h]",yaxis="Acetate [mmol/gDW L]",legend=false)
 #p32 = plot!(ts_min,vk_min[:,2],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Acetate [mmol/gDW L]",label="Min")
 #p32 = plot!(ts_max,vk_max[:,2],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Acetate [mmol/gDW L]",label="Max")

 p33 = plot(ts,-vk[:,3],linewidth=2,xaxis="time [h]",yaxis="Glucose [mmol/gDW L]",legend=false)
 #p33 = plot!(ts_min,-vk_min[:,3],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Glucose [mmol/gDW L]",label="Min")
 #p33 = plot!(ts_max,-vk_max[:,3],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Glucose [mmol/gDW L]",label="Max")

 p34 = plot(ts,-vk[:,4],linewidth=2,xaxis="time [h]",yaxis="Oxygen [mmol/gDW L]",legend=false)
 #p34 = plot!(ts_min,-vk_min[:,4],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Oxygen [mmol/gDW L]",label="Min")
 #p34 = plot!(ts_max,-vk_max[:,4],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Oxygen [mmol/gDW L]",label="Max")

 g3= plot(p31,p32,p33,p34,layout=(2,2),reuse = true)

l = @layout [
    a    
    b  c 
]

lj = @layout [
    a{0.80w} _ b{0.01w}
    c{0.80w} _ d{0.01w}
]

lk = @layout [
    a{0.80w} _ b{0.01w}
]

p19=plot()

savefig(g1,"test.pdf")
savefig(g2,"test2.pdf")
savefig(g3,"test3.pdf")

