#=
      simulation of multi-stage dynamic optimization of dFBA by orthogonal collocation
=#

#−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
#Formulating the discrete dynamics
#using TickTock
using JuMP
using Ipopt
using LinearAlgebra
using Plots
ENV["GKSwstype"] = "100"
#reading GSM model
using DelimitedFiles

include("pFBA_KKT_flux.jl")

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
vlb = readdlm("lb.csv", ',');
vub = readdlm("ub.csv", ',');
nc = 5
nmv = 1 # number of manipulated variables
ns = 2   #number of scenarios


##...Model fixed Parameters
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
vl0=1.0
klo = 7.5  # 1/h 
O_sat = 0.21 # mmol 

#Process contraints
u_ini = 0.0         # L/h
u_max_rate = 0.05 
u_max = 0.5
u_min = 0.0
V_max = 3.0
Ac_max = 4.0

#flux rate change constraint inside the elements
v_dot = 0.1 

d=0.0*Vector{Float64}(undef,nv) 
d[obj]=-1
up=0.0*Vector{Float64}(undef,nv) 
up[glu] = 1
up2=0.0*Vector{Float64}(undef,nv) 
up2[o2] = 1
v0=zeros(nv)


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
nfe    = 24  # number of control intervals
ncp    = 3           # number of collocation points
th     = 12.0 # time horizon 
h = th/nfe           # length of one finite element on the time horizon
hm=h*ones(nfe)'
var_h=0.1

#w = Vector{Float64}(undef,nv)  
w = Matrix{Float64}(undef,nv,ns)  
w[ac,1] = -1e-1   # max -1e-4
w[ac,2] = 1e0  # min 2.55e-4
phi = 1.0
phi_all = 1e-4
phi1=phi_all
phi2=phi_all
phi3=phi_all

ts     = Matrix{Float64}(undef,nfe+1,ns) # time series for plotting
tsn     = Matrix{Float64}(undef,nfe*ncp+1,ns)
xk = Array{Float64}(undef,nc,nfe*ncp+1,ns)
uk = Array{Float64}(undef,nmv,nfe+1,ns)
vk = Array{Float64}(undef,nv,nfe+1,ns)
ts = setindex!(ts,0,1,1)
ts = setindex!(ts,0,1,2)

    x0 = [g0,ac0,b0,O0,vl0]

    x_ini = x0
    v_ini = v0

    #x_ini = sol
    #v_ini = solv

         solu, sol, solv, solcd, soll, solal, solag, solh = pFBA_KKT_flux(x0,x_ini,v_ini)
        
      for s in 1:ns   
        for i in 2:nfe+1
            #ts[i] = h*i
            ts[i,s] = ts[i-1,s] + solh[i-1,s]
            #ts[i] = ts[i-1] + hm[i-1]
        end
      end
        

#building vectors
#xk[:,1]= x0.*cs
for s in 1:ns
setindex!(xk, x0, :, 1,s)
setindex!(uk, u_ini, 1, 1,s)
tsn[1,s] = 0.0
end


radau  = [0.15505 0.64495 1.00000]
for s in 1:ns
   global kk = 2
for i in 1:nfe


   for j in 1:ncp
    xk[:,kk,s]=sol[:,i,j,s]
    
    #vk[:,kk]=K*solv[:,i,j]
    

     #tsn[kk] = ts[i] + radau[j]*h
     #tsn[kk] = ts[i] + radau[j]*hm[i]
     tsn[kk,s] = ts[i,s] + radau[j]*solh[i,s]
    global kk=kk+1
   end
   
   vk[:,i+1,s]=solv[:,i,s]
   uk[:,i+1,s]=solu[:,i,s]
end
end

for s in 1:ns
  vk[:,1,s]=vk[:,2,s]
end

writedlm("xk.csv",xk)
writedlm("tsn.csv",tsn)
#writedlm("../DA_Plant/uk.csv",uk)
#writedlm("../DA_Plant/nfe.csv",nfe)
#writedlm("../DA_Plant/th.csv",th)

xk1 = xk[:,:,1]'
xk2 = xk[:,:,2]'
vk1 = vk[:,:,1]'
vk2 = vk[:,:,2]'


p11 = plot(tsn[:,1],xk1[:,3],linewidth=2,xaxis="time [h]",yaxis="biomass [g/L]",label="Max")
p11 = plot!(tsn[:,1],xk2[:,3],linewidth=2,xaxis="time [h]",yaxis="biomass [g/L]",label="Min")

p12 = plot(tsn[:,1],xk1[:,2],linewidth=2,xaxis="time [h]",yaxis="Ac [g/L]",label="Max")
p12 = plot!(tsn[:,1],xk2[:,2],linewidth=2,xaxis="time [h]",yaxis="Ac [g/L]",label="Min")

p14 = plot(tsn[:,1],xk1[:,1],linewidth=2,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max")
p14 = plot!(tsn[:,1],xk2[:,1],linewidth=2,xaxis="time [h]",yaxis="Glucose [g/L]",label="Min")

p13 = plot(tsn[:,1],xk1[:,4],linewidth=2,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max")
p13 = plot!(tsn[:,1],xk2[:,4],linewidth=2,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Min")

g1 = plot(p11,p12,p14,p13,layout=(2,2),reuse = true)


p21 = plot(ts[:,1],uk[1,:,1],linewidth=2,xaxis="time [h]",yaxis="F [L/h]", linetype=:steppre,label="Max")
p21 = plot!(ts[:,1],uk[1,:,2],linewidth=2,xaxis="time [h]",yaxis="F [L/h]", linetype=:steppre,label="Min")
p23 = plot(tsn[:,1],xk1[:,5],linewidth=2,xaxis="time [h]",yaxis="Volume [L]",label="Max")
p23 = plot!(tsn[:,1],xk2[:,5],linewidth=2,xaxis="time [h]",yaxis="Volume [L]",label="Min")
g2 = plot(p21,p23,layout=(2,1),reuse = true)

p31 = plot(ts[:,1],vk1[:,obj],linewidth=2,xaxis="time [h]",yaxis="Growth rate [1/h]",label="Max")
p31 = plot!(ts[:,1],vk2[:,obj],linewidth=2,xaxis="time [h]",yaxis="Growth rate [1/h]",label="Min")

p32 = plot(ts[:,1],vk1[:,ac],linewidth=2,xaxis="time [h]",yaxis="Acetate [mmol/gDW L]",label="Max")
p32 = plot!(ts[:,1],vk2[:,ac],linewidth=2,xaxis="time [h]",yaxis="Acetate [mmol/gDW L]",label="Min")

p33 = plot(ts[:,1],vk1[:,glu],linewidth=2,xaxis="time [h]",yaxis="Glucose [mmol/gDW L]",label="Max")
p33 = plot!(ts[:,1],vk2[:,glu],linewidth=2,xaxis="time [h]",yaxis="Glucose [mmol/gDW L]",label="Min")

p34 = plot(ts[:,1],vk1[:,o2],linewidth=2,xaxis="time [h]",yaxis="Oxygen [mmol/gDW L]",label="Max")
p34 = plot!(ts[:,1],vk2[:,o2],linewidth=2,xaxis="time [h]",yaxis="Oxygen [mmol/gDW L]",label="Min")


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
