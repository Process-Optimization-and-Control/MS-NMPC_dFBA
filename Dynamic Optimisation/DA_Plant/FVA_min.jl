#FBA
function FVA_min(mu_min,vg_min,vo_min)

     #index
  

     #input
     vlb[glu]=-vg_min
      #vub[glu]=vlb[glu]
     vub[glu]=0.0
      vlb[obj]=(1.0-sub_FVA)*mu_min
      vub[obj]=1000
      vlb[o2]=-vo_min
     #lb[o2]=-vo
     #ub[o2]=lb[o2]
      #ub[o2]=0.0


     # declare the JuMP model
     #m = Model(with_optimizer(Cbc.Optimizer,print_level = 0))
    # m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27") )
     #m = Model(with_optimizer(Ipopt.Optimizer,logLevel = 0))
     # m = Model(Cbc.Optimizer)
      #set_optimizer_attributes(m, "print_level" => 1)
     # m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 1, max_iter = 200000) )
     m_min = Model(Cbc.Optimizer)
      set_attribute(m_min, "logLevel", 0)

     # variables representing x1 and x2
     @variables(m_min, begin
     v_i[1:nv]
     end)

     @constraints(m_min, begin
     M1[i=1:nm],  sum(S[i,k]*v_i[k] for k in 1:nv) == 0
     M2[i=1:nv], v_i[i] <= vub[i]
     M3[i=1:nv], v_i[i] >= vlb[i]

     #M2[i=act_ub], v[i] <= vub[i]
     #M3[i=act_lb], v[i] >= vlb[i]
     end)

    # and we declare the objective function
    @objective(m_min, Min, v_i[ac])
    # now, we solve the model
    solveNLP = JuMP.optimize!
    status = solveNLP(m_min)

    mi=JuMP.value.(v_i[obj])
    vac=JuMP.value.(v_i[ac])
    vo=JuMP.value.(v_i[o2])
    vg=JuMP.value.(v_i[glu])

return  mi,vac, vg, vo
end
#end
