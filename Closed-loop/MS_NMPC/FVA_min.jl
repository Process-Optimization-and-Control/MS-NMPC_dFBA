#FBA
function FVA_min(mu_min,vg_min,vo_min)

     #index
  

     #input
     vlbi[iglu]=-vg_min
      #vub[glu]=vlb[glu]
     vubi[iglu]=0.0
      vlbi[iobj]=(1.0-sub_FVA)*mu_min
      vubi[iobj]=1000
      vlbi[io2]=-vo_min
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
     #m_min = Model(Cbc.Optimizer)
      #set_attribute(m_min, "logLevel", 0)
     # set_attribute(m_min, "print_level", 0)

     m_min = Model(Ipopt.Optimizer)
     set_attribute(m_min,"print_level", 0)
     set_attribute(m_min,"linear_solver", "ma27")

     # variables representing x1 and x2
     @variables(m_min, begin
     v_i[1:nvi]
     end)

     @constraints(m_min, begin
     M1[i=1:nmi],  sum(Si[i,k]*v_i[k] for k in 1:nvi) == 0
     M2[i=1:nvi], v_i[i] <= vubi[i]
     M3[i=1:nvi], v_i[i] >= vlbi[i]

     #M2[i=act_ub], v[i] <= vub[i]
     #M3[i=act_lb], v[i] >= vlb[i]
     end)

    # and we declare the objective function
    @objective(m_min, Min, v_i[iac])
    # now, we solve the model
    solveNLP = JuMP.optimize!
    status = solveNLP(m_min)

    mi=JuMP.value.(v_i[iobj])
    vac=JuMP.value.(v_i[iac])
    vo=JuMP.value.(v_i[io2])
    vg=JuMP.value.(v_i[iglu])

return  mi,vac, vg, vo
end
#end
