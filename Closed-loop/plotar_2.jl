using Plots
using DelimitedFiles

nfe = 70
ts_shr = readdlm("ts_shr.csv")
u_shr = readdlm("u_shr.csv")
x_shr = readdlm("x_shr.csv")
#x_shr = x_shr'
var_shr = readdlm("var_shr.csv")
time_shr = readdlm("time_shr.csv")

ts_rec = readdlm("ts_rec.csv")
u_rec = readdlm("u_rec.csv")
x_rec = readdlm("xk_rec.csv")
#x_rec = x_rec'
var_rec = readdlm("var_rec.csv")
time_rec = readdlm("time_rec.csv")


p11 = plot(ts_shr,x_shr[3,:],linewidth=2, c=:blue,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)

p11 = plot!(ts_rec,x_rec[3,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="biomass g/L",label="Max",legend=false,grid = false)


p12 = plot(ts_shr,x_shr[2,:],linewidth=2.0, c=:blue,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts_shr,5.0*ones(1,nfe+1)',linestyle = :dash, c=:red,linewidth=2,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)
p12 = plot!(ts_rec,x_rec[2,:],linewidth=2.0,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="Acetate [g/L]",label="Max",legend=false,grid = false)


p13 = plot(ts_shr,u_shr[1,:],linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false, c=:blue, linetype=:steppre,grid = false) 
p13 = plot!(ts_shr,0.5*ones(1,nfe+1)',linestyle = :dash, c=:red,linewidth=2,xaxis="time [h]",yaxis="F [L/h]",legend=false, linetype=:steppre,grid = false)     
p13 = plot!(ts_rec,u_rec[1,:],linewidth=2,linestyle = :dot, c=:green,xaxis="time [h]",yaxis="F [L/h]",legend=false, linetype=:steppre,grid = false) 

g1 = plot(p11,p12,p13,layout=(3,1),reuse = true)


p14 = plot(ts_shr[2:end],time_shr[1:end], seriestype=:scatter,xaxis="time [h]",yaxis="Comp. time [s]",legend=false, markersize = 3,markercolor = :white, markerstrokecolor=:blue,grid = false) 

p14 = plot!(ts_shr[2:end],time_rec[1:end], seriestype=:scatter,xaxis="time [h]",yaxis="Comp. time [s]",legend=false, markersize = 3,markercolor = :white, markerstrokecolor=:green,grid = false) 

p15 = plot(ts_shr[2:end],var_shr[1:end], seriestype=:scatter,xaxis="time [h]",yaxis="Opt. variables",legend=false, markersize = 3,markercolor = :white, markerstrokecolor=:blue,grid = false) 
p15 = plot!(ts_shr[2:end],var_rec[1:end], seriestype=:scatter,xaxis="time [h]",yaxis="Opt. variables",legend=false, markersize = 3,markercolor = :white, markerstrokecolor=:green,grid = false) 

g2 = plot(p14,p15,layout=(2,1),reuse = true)



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

savefig(g1,"horizon_1.pdf")
savefig(g2,"horizon_2.pdf")

