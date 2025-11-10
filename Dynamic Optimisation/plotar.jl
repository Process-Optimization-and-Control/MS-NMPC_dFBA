#Ra32353235
#  daviddeo@sobol.biozone.utoronto.ca:9833
#  scp daviddeo@sobol.biozone.utoronto.ca:/nfs/homes/daviddeo/ts.csv ./
#  cd("Documents/UT/dFBA/code/ecoli/")

using Plots
ENV["GKSwstype"] = "100"
using DelimitedFiles
#xk = readdlm("xk.csv")
#tsn = readdlm("tsn.csv")

xk_DA = readdlm("xk_DA.csv")
xk_DA=xk_DA'
#xk_DA = readdlm("dFBAlab_core.csv",',')
#xk_DA=xk_DA'
xk_DA=xk_DA[:,1:69]
ts_DA = readdlm("ts_DA.csv")
#ts_DA = readdlm("dFBAlab_core_ts.csv")
ts_DA = ts_DA[1:69]
#vk = readdlm("vk.csv")
#lStar = readdlm("lStar.csv")
#auStar = readdlm("auStar.csv")
#alStar = readdlm("alStar.csv")
#agStar = readdlm("agStar.csv")


p11 = plot(tsn,0.21*ones(size(tsn)),linewidth=2,xaxis="time [h]",yaxis="O2 [mM]",ylims=[0.0,0.3],legend=false)
#C_C
p12 = plot(tsn,xk[3,:],linewidth=2,xaxis="time [h]",yaxis="biomass [g]",label="NLP-PF",legend=:topleft)
p12 = plot!(ts_DA,xk_DA[3,:],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Biomass [g]",label="DA",legend=:topleft)
#p12 = plot!(tsn_8,xk_8[3,:],linewidth=2,xaxis="time [h]",yaxis="biomass [g]",label="DOA-8.0h",legend=:topleft)
#p12 = plot!(tsn_9,xk_9[3,:],linewidth=2,xaxis="time [h]",yaxis="biomass [g]",label="DOA-9.0h",legend=:topleft)

p13 = plot(tsn,xk[2,:],linewidth=2,xaxis="time [h]",yaxis="Ac [mM]",label="NLP-PF",legend=:topleft)
p13 = plot!(ts_DA,xk_DA[2,:],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Acetate [mM]",label="DA",legend=:topleft)
#p13 = plot!(tsn_8,xk_8[2,:],linewidth=2,xaxis="time [h]",yaxis="Ac [mM]",label="DOA-8.0h",legend=:topleft)
#p13 = plot!(tsn_9,xk_9[2,:],linewidth=2,xaxis="time [h]",yaxis="Ac [mM]",label="DOA-9.0h",legend=:topleft)

# #C_D
p14 = plot(tsn,xk[1,:],linewidth=2,xaxis="time [h]",yaxis="Glucose [mM]",label="NLP-PF",legend=:bottomleft)
p14 = plot!(twinx(),ts_DA,xk_DA[1,:],linewidth=2,linestyle = :dash,xaxis="time [h]",yaxis="Glucose [mM]",label="DA",legend=:bottomleft)
#p14 = plot!(tsn_8,xk_8[1,:],linewidth=2,xaxis="time [h]",yaxis="Glucose [mM]",label="DOA-8.0h",legend=:bottomleft)
#p14 = plot!(tsn_9,xk_9[1,:],linewidth=2,xaxis="time [h]",yaxis="Glucose [mM]",label="DOA-9.0h",legend=:bottomleft)

#p15 = plot(tsn,vk[ac,:],linewidth=2,seriestype = :scatter,xaxis="time [h]",yaxis="fluxes",label="Ac",legend=:bottomleft)
#p16 = plot(tsn,vk[obj,:],linewidth=2,seriestype = :scatter,xaxis="time [h]",yaxis="fluxes",label="mi",legend=:bottomleft)
#p17 = plot(tsn,vk[glu,:],linewidth=2,seriestype = :scatter,xaxis="time [h]",yaxis="fluxes",label="glu",legend=:bottomleft)
#p18 = plot(tsn,vk[o2,:],linewidth=2,seriestype = :scatter,xaxis="time [h]",yaxis="fluxes",label="o2",legend=:bottomleft)

l = @layout [
    a    
    b  c 
]

g2 = plot(p12,p13,p14,layout=l,reuse = true)
#g2 = plot(p11,p12,p13,p14,layout=(2,2),reuse = true)
#g3 = plot(p15,p16,p17,p18,layout=(2,2),reuse = true)
#display(g2)

savefig(g2,"test2.pdf")

#savefig(g3,"test2.pdf")