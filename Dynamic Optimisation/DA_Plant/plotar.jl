
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
p11 = plot!(ts,x_all[3,:,10],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)

p12 = plot(ts,x_all[2,:,1],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,5.0*ones(1,nfe+1)',linestyle = :dash, c=:red,linewidth=2,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)


p12 = plot!(ts,x_all[2,:,2],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,3],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,4],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,5],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,6],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,7],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,8],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,9],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts,x_all[2,:,10],linewidth=0.5, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)

p13 = plot(ts,u[1,:],linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false, c=:blue, linetype=:steppre,grid = false) 
p13 = plot!(ts,0.5*ones(1,nfe+1)',linestyle = :dash, c=:red,linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false, linetype=:steppre,grid = false)     


 g1 = plot(p11,p12,p13,layout=(3,1),reuse = true)


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

savefig(g1,"open_loop.pdf")

