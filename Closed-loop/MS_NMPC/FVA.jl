#FBA
function FVA(mu_max,vg_max,vo_max)

     #index
  

     #input
     vlbi[iglu]=-vg_max
      #vub[glu]=vlb[glu]
      vubi[iglu]=0.0
      vlbi[iobj]=(1.0-sub_FVA)*mu_max
      vubi[iobj]=1000
     vlbi[io2]=-vo_max
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
     #m_max= Model(Cbc.Optimizer)
     #set_attribute(m_max, "logLevel", 0)
     #set_attribute(m_max, "print_level", 0)

     m_max = Model(Ipopt.Optimizer)
     set_attribute(m_max,"print_level", 0)
     set_attribute(m_max,"linear_solver", "ma27")

     # variables representing x1 and x2
     @variables(m_max, begin
     v_m[1:nvi]
     end)

     @constraints(m_max, begin
     M1[i=1:nmi],  sum(Si[i,k]*v_m[k] for k in 1:nvi) == 0
     M2[i=1:nvi], v_m[i] <= vubi[i]
     M3[i=1:nvi], v_m[i] >= vlbi[i]

     #M2[i=act_ub], v[i] <= vub[i]
     #M3[i=act_lb], v[i] >= vlb[i]
     end)

    # and we declare the objective function
    @objective(m_max, Max, v_m[iac])
    # now, we solve the model
    solveNLP = JuMP.optimize!
    status = solveNLP(m_max)

    mi=JuMP.value.(v_m[iobj])
    vac=JuMP.value.(v_m[iac])
    vo=JuMP.value.(v_m[io2])
    vg=JuMP.value.(v_m[iglu])
    vall=JuMP.value.(v_m[:])

return  mi,vac, vg, vo
end
#end
