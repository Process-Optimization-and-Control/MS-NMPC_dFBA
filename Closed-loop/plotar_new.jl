#Ra32353235
#  daviddeo@sobol.biozone.utoronto.ca:9833
#  scp daviddeo@sobol.biozone.utoronto.ca:/nfs/homes/daviddeo/ts.csv ./
#  cd("Documents/UT/dFBA/code/ecoli/")

include("interpolation.jl")


using Plots
ENV["GKSwstype"] = "100"
using DelimitedFiles
#xk = readdlm("xk.csv")
#tsn = readdlm("tsn.csv")

xk_DA = readdlm("xk_DA.csv")
xk_DA=xk_DA'
xk_DA_=xk_DA[:,1:69]
ts_DA = readdlm("ts_DA.csv")

ts_DA_ = ts_DA[1:69]

xk_lab = readdlm("dFBAlab_core.csv",',')
xk_lab=xk_lab'
xk_lab_=xk_lab[:,1:69]
ts_lab = readdlm("dFBAlab_core_ts.csv")
ts_lab_ = ts_lab[1:69]

#vk = readdlm("vk.csv")
#lStar = readdlm("lStar.csv")
#auStar = readdlm("auStar.csv")
#alStar = readdlm("alStar.csv")
#agStar = readdlm("agStar.csv")




#p13 = plot(tsn_,xk_[2,:],linewidth=2,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF")
#p13 = plot!(ts_lab_,xk_lab_[2,:],linewidth=2,linestyle = :dot,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF",legend=:false)
#p13 = plot!(ts_DA_,xk_DA_[2,:],linewidth=2,color=:green,linestyle = :dash,xaxis="time [h]",yaxis="Acetate [mM]",label="DA",legend=:false)


# #C_D
#p13 = plot!(tsn_,xk_[1,:],linewidth=2,xaxis="time [h]",color=:red,yaxis="Glucose [mM]",label="NLP-PF",legend=:false)
#p13 = plot!(ts_lab_,xk_lab_[1,:],linewidth=2,linestyle = :dot,color=:red,xaxis="time [h]",yaxis="Concentration [mmol/L]",label="DA",legend=:false)
#p13 = plot!(ts_DA_,xk_DA_[1,:],linewidth=2,linestyle = :dash,color=:red,xaxis="time [h]",yaxis="Conc. [mmol/L]",label="DA",legend=:false,grid=:off)



#p13 = plot!(twinx(),tsn_,xk_[3,:],linewidth=2,color=:blue,legend=:false,yticks = 0:0.1:0.8,ylims=(0,0.8))
#p13 = plot!(twinx(),ts_lab_,xk_lab_[3,:],linewidth=2,color=:blue,linestyle = :dot,legend=:false,yticks = 0:0.2:0.8,ylims=(0,0.8))
#p13 = plot!(twinx(),ts_DA_,xk_DA_[3,:],linewidth=2,color=:blue,linestyle = :dash,xaxis="time [h]",yaxis="Biomass [g/L]",label="DA",legend=:false,yticks = 0:0.1:0.8,ylims=(0,0.8))





p14 = plot(tsn_i,xk_i[2,:],linewidth=2,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF")

#p14 = plot!(ts,zeros(size(ts)),seriestype = :scatter,linewidth=2,color=:black,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF")
#p14 = plot!(tsn_ll,xk_ll[2,:],seriestype = :scatter,linewidth=2,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF")

p14 = plot!(ts_lab,xk_lab[2,:],linewidth=2,linestyle = :dot,color=:green,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF",legend=:false)
p14 = plot!(ts_DA,xk_DA[2,:],linewidth=2,color=:green,linestyle = :dash,xaxis="time [h]",yaxis="Acetate [mM]",label="DA",legend=:false)


# #C_D
p14 = plot!(tsn_i,xk_i[1,:],linewidth=2,xaxis="time [h]",color=:red,yaxis="Glucose [mM]",label="NLP-PF",legend=:false)

#p14 = plot!(tsn_ll,xk_ll[1,:],seriestype = :scatter,linewidth=2,xaxis="time [h]",color=:red,yaxis="Glucose [mM]",label="NLP-PF",legend=:false)

p14 = plot!(ts_lab,xk_lab[1,:],linewidth=2,linestyle = :dot,color=:red,xaxis="time [h]",yaxis="Concentration [mmol/L]",label="DA",legend=:false)
p14 = plot!(ts_DA,xk_DA[1,:],linewidth=2,linestyle = :dash,color=:red,xaxis="time [h]",yaxis="Conc. [mmol/L]",label="DA",legend=:false,grid=:off,ylims=(0,10.5))



#p14 = plot!(twinx(),tsn_ll,xk_ll[3,:],seriestype = :scatter,linewidth=2,color=:blue,legend=:false,yticks = 0:0.2:0.9,ylims=(0,0.9))

p14 = plot!(twinx(),tsn_i,xk_i[3,:],linewidth=2,color=:blue,legend=:false,yticks = 0:0.2:0.9,ylims=(0,0.9))
p14 = plot!(twinx(),ts_lab,xk_lab[3,:],linewidth=2,color=:blue,linestyle = :dot,legend=:false,yticks = 0:0.2:0.9,ylims=(0,0.9))
p14 = plot!(twinx(),ts_DA,xk_DA[3,:],linewidth=2,color=:blue,linestyle = :dash,xaxis="time [h]",yaxis="Biomass [g/L]",label="DA",legend=:false,yticks = 0:0.2:0.9,ylims=(0,0.9))

#p12 = plot!(tsn_8,xk_8[3,:],linewidth=2,xaxis="time [h]",yaxis="biomass [g]",label="DOA-8.0h",legend=:topleft)
#p12 = plot!(tsn_9,xk_9[3,:],linewidth=2,xaxis="time [h]",yaxis="biomass [g]",label="DOA-9.0h",legend=:topleft)


#p15 = plot(tsn,vk[ac,:],linewidth=2,seriestype = :scatter,xaxis="time [h]",yaxis="fluxes",label="Ac",legend=:bottomleft)
#p16 = plot(tsn,vk[obj,:],linewidth=2,seriestype = :scatter,xaxis="time [h]",yaxis="fluxes",label="mi",legend=:bottomleft)
#p17 = plot(tsn,vk[glu,:],linewidth=2,seriestype = :scatter,xaxis="time [h]",yaxis="fluxes",label="glu",legend=:bottomleft)
#p18 = plot(tsn,vk[o2,:],linewidth=2,seriestype = :scatter,xaxis="time [h]",yaxis="fluxes",label="o2",legend=:bottomleft)

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

#g2 = plot(p13,p19,p14,p19,layout=lj,reuse = true)
g2 = plot(p14,p19,layout=lk,reuse = true)
#g2 = plot(p14,layout=(2,2),reuse = true)
#g3 = plot(p15,p16,p17,p18,layout=(2,2),reuse = true)
#display(g2)

savefig(g2,"test4.pdf")

#savefig(g3,"test2.pdf")