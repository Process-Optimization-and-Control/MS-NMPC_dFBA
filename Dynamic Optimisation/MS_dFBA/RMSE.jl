#RMSE
res_lab = Matrix{Float64}(undef,size(xk,1),size(ts,1))
id_lab = Vector{Float64}(undef,size(ts,1))
RMSE_lab = Vector{Float64}(undef,size(xk,1))

res_DA = Matrix{Float64}(undef,size(xk,1),size(ts,1))
id_DA = Vector{Float64}(undef,size(ts,1))
RMSE_DA = Vector{Float64}(undef,size(xk,1))

norma=[xk_i[1,1],maximum(xk_i[2,:]),xk_i[3,end]]

id_lab[1]=1
id_DA[1]=1

for i in 2:size(ts,1)
    
    for j in 2:size(ts_lab,1)
      
       
       if abs(ts_lab[j] - ts[i]) > abs(ts_lab[j-1] - ts[i]) 
         break
       end
      
       id_lab[i] = j

    end

    for j in 2:size(ts_DA,1)
      
       
      if abs(ts_DA[j] - ts[i]) > abs(ts_DA[j-1] - ts[i]) 
        break
    end
     
      id_DA[i] = j

   end

end

id_lab=Int.(id_lab)

res_lab[:,1] = xk[:,1] - xk_lab[:,1] 

id_DA=Int.(id_DA)

res_DA[:,1] = xk[:,1] - xk_DA[:,1] 

for i in 1:size(ts,1)-1
    #println(i)
    res_lab[:,i+1] = xk[:,(3*i)+1] - xk_lab[:,id_lab[i+1]] 

    res_DA[:,i+1] = xk[:,(3*i)+1] - xk_DA[:,id_DA[i+1]] 

    #println((3*i)+1)
end





for i in 1:size(xk,1)

RMSE_lab[i] = sqrt( sum(res_lab[i,:].^2)/size(ts,1) )

RMSE_DA[i] = sqrt( sum(res_DA[i,:].^2)/size(ts,1) )

end

RMSE_lab[:] = 100*(RMSE_lab[:]./norma[:])
RMSE_DA[:] = 100*(RMSE_DA[:]./norma[:])

#println(``RMSE_lab:``)
println(RMSE_lab)
#println(``RMSE_DA:``)
println(RMSE_DA)