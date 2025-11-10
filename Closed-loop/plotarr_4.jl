using Plots
using DelimitedFiles

nfe = 70
ts_shr = readdlm("ts_shr.csv")
u_shr = readdlm("u_shr.csv")
x_shr = readdlm("x_shr.csv")
#x_shr = x_shr'
vk_shr = readdlm("vk_shr.csv")

ts_rec = readdlm("ts_rec.csv")
u_rec = readdlm("u_rec.csv")
x_rec = readdlm("xk_rec.csv")
#x_rec = x_rec'
vk_rec = readdlm("vk_rec.csv")

p11 = plot(ts_shr,x_shr[3,:],linewidth=2,xaxis="time [h]",yaxis="biomass [g/L]",label="FBA")
p11 = plot!(ts_rec,x_rec[3,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="biomass [g/L]",label="Max",legend=false,grid = false)

p12 = plot(ts_shr,x_shr[2,:],linewidth=2,xaxis="time [h]",yaxis="Acetate [g/L]",label="FBA")
p12 = plot!(ts_rec,x_rec[2,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)

p14 = plot(ts_shr,x_shr[1,:],linewidth=2,xaxis="time [h]",yaxis="Glucose [g/L]",label="FBA")
p14 = plot!(ts_rec,x_rec[1,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)

p13 = plot(ts_shr,x_shr[4,:],linewidth=2,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="FBA")
p13 = plot!(ts_rec,x_rec[4,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,ylims = (0,0.22),grid = false)


 g1 = plot(p11,p12,p14,p13,layout=(2,2),reuse = true)


# p21 = plot(ts,u[1,:],linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false)
# p22 = plot(ts,u[2,:],linewidth=2,xaxis="time [h]",yaxis="Klo [1/h]",legend=false)
# p23 = plot(ts,xk[:,5],linewidth=2,xaxis="time [h]",yaxis="Volume [L]"),label="Min"
# g2 = plot(p21,p22,p23,layout=(3,1),reuse = true)

 p31 = plot(ts_shr,vk_shr[1,:],linewidth=2,xaxis="time [h]",yaxis="Growth rate [1/h]",label="FBA")
 p31 = plot!(ts_rec,vk_rec[1,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="Growth rate [1/h]",label="Max",grid = false,legend=false)

 p32 = plot(ts_shr,vk_shr[2,:],linewidth=2,xaxis="time [h]",yaxis="Acetate [mmol/gDW L]",label="FBA")
 p32 = plot!(ts_rec,vk_rec[2,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="Ac. [mmol/gDW L]",label="Max",grid = false,legend=false)

 p33 = plot(ts_shr,vk_shr[3,:],linewidth=2,xaxis="time [h]",yaxis="Glucose [mmol/gDW L]",label="FBA")
 p33 = plot!(ts_rec,vk_rec[3,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="Glucose [mmol/gDW L]",label="Max",grid = false,legend=false)

 p34 = plot(ts_shr,vk_shr[4,:],linewidth=2,xaxis="time [h]",yaxis="Oxygen [mmol/gDW L]",label="FBA")
 p34 = plot!(ts_rec,vk_rec[4,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="Oxygen [mmol/gDW L]",label="Max",grid = false,legend=false)

 g3= plot(p31,p32,p33,p34,layout=(2,2),reuse = true)


p21 = plot(ts_shr,u_shr[1,:],linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false, c=:blue, linetype=:steppre,grid = false) 
p21 = plot!(ts_shr,0.5*ones(1,nfe+1)',linestyle = :dash, c=:red,linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false, linetype=:steppre,grid = false)     
p21 = plot!(ts_rec,u_rec[1,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="F [L/h]",legend=false, linetype=:steppre,grid = false) 

p22 = plot(ts_shr,x_shr[5,:],linewidth=2,xaxis="time [h]",yaxis="Volume [L]",label="FBA")
p22 = plot!(ts_shr,3.0*ones(1,nfe+1)',linestyle = :dash, c=:red,linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false, linetype=:steppre,grid = false) 
p22 = plot!(ts_rec,x_rec[5,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="Volume [L]",label="Max",legend=false,grid = false)

g2 = plot(p22,p21,layout=(2,1),reuse = true)

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

savefig(g1,"Comp_states.pdf")
savefig(g2,"Comp_control.pdf")
savefig(g3,"Comp_fluxes.pdf")