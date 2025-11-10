
using DelimitedFiles
using Plots
using JLD
ENV["GKSwstype"] = "100"

ts = readdlm("ts_open_loop.csv")
u = readdlm("u_open_loop.csv")
#save("x_all.jld", "data", x_all)
x_all= load("x_all.jld")["data"]

nfe    = 40  # number of control intervals

p11 = plot(ts,x_all[3,:,1],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)

p11 = plot!(ts,x_all[3,:,2],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)
p11 = plot!(ts,x_all[3,:,3],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)
p11 = plot!(ts,x_all[3,:,4],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)
p11 = plot!(ts,x_all[3,:,5],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)
p11 = plot!(ts,x_all[3,:,6],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)
p11 = plot!(ts,x_all[3,:,7],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)
p11 = plot!(ts,x_all[3,:,8],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)
p11 = plot!(ts,x_all[3,:,9],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)
p11 = plot!(ts,x_all[3,:,10],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass [g/L]",label="Max",legend=false,grid = false)



p12 = plot(ts,x_all[2,:,1],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,2],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,3],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,4],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,5],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,6],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,7],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,8],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,9],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,10],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,5.0*ones(1,nfe+1)',linestyle = :dash, c=:red,linewidth=2,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)


p13 = plot(ts,x_all[1,:,1],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)

p13 = plot!(ts,x_all[1,:,2],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)
p13 = plot!(ts,x_all[1,:,3],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)
p13 = plot!(ts,x_all[1,:,4],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)
p13 = plot!(ts,x_all[1,:,5],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)
p13 = plot!(ts,x_all[1,:,6],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)
p13 = plot!(ts,x_all[1,:,7],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)
p13 = plot!(ts,x_all[1,:,8],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)
p13 = plot!(ts,x_all[1,:,9],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)
p13 = plot!(ts,x_all[1,:,10],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Glucose [g/L]",label="Max",legend=false,grid = false)


p14 = plot(ts,x_all[4,:,1],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)

p14 = plot!(ts,x_all[4,:,2],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)
p14 = plot!(ts,x_all[4,:,3],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)
p14 = plot!(ts,x_all[4,:,4],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)
p14 = plot!(ts,x_all[4,:,5],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)
p14 = plot!(ts,x_all[4,:,6],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)
p14 = plot!(ts,x_all[4,:,7],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)
p14 = plot!(ts,x_all[4,:,8],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)
p14 = plot!(ts,x_all[4,:,9],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)
p14 = plot!(ts,x_all[4,:,10],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Oxygen [mmol/L]",label="Max",legend=false,grid = false)


p21 = plot(ts,x_all[5,:,1],linewidth=2, c=:blue,xaxis="time [h]",yaxis="Volume [L]",label="Max",legend=false,grid = false)
p21 = plot!(ts,3.0*ones(1,nfe+1)',linestyle = :dash, c=:red,linewidth=2,xaxis="time [h]",yaxis="Volume [L]",label="Max",legend=false,grid = false)

p23 = plot(ts,u[1,:],linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false, c=:blue, linetype=:steppre,grid = false) 
p23 = plot!(ts,0.5*ones(1,nfe+1)',linestyle = :dash, c=:red,linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false, linetype=:steppre,grid = false)     



g2 = plot(p21,p23,layout=(2,1),reuse = true)
 g1 = plot(p11,p12,p14,p13,layout=(2,2),reuse = true)
 #g3= plot(p31,p32,p33,p34,layout=(2,2),reuse = true)

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

savefig(g1,"open_loop_scn_states.pdf")
savefig(g2,"open_loop_scn_control.pdf")
#savefig(g3,"test3.pdf")

