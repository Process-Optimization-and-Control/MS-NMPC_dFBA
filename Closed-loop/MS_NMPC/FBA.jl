#FBA
function FBA(vg_f,vo_f)

     #index

     #input
     vlbi[iglu]=-vg_f
     # ub[glu]=lb[glu]
      vubi[iglu]=0.0
      vlbi[iobj]=0
      vubi[iobj]=1000
      vlbi[io2]=-vo_f
     #ub[o2]=lb[o2]
      #ub[o2]=0.0


     # declare the JuMP model
    # m = Model(with_optimizer(Cbc.Optimizer,print_level = 0))
    # m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27") )
     #m = Model(with_optimizer(Ipopt.Optimizer,logLevel = 0))
      #m = Model(Cbc.Optimizer)
     # set_attribute(m, "logLevel", 0)
     # set_attribute(m, "print_level", 0)
      m = Model(Ipopt.Optimizer)
    set_attribute(m,"print_level", 0)
    set_attribute(m,"linear_solver", "ma27")
     # m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 1, max_iter = 200000) )


     # variables representing x1 and x2
     @variables(m, begin
     v[1:nvi]
     end)

     @constraints(m, begin
     M1[i=1:nmi],  sum(Si[i,k]*v[k] for k in 1:nvi) == 0
     M2[i=1:nvi], v[i] <= vubi[i]
     M3[i=1:nvi], v[i] >= vlbi[i]

     #M2[i=act_ub], v[i] <= vub[i]
     #M3[i=act_lb], v[i] >= vlb[i]
     end)

    # and we declare the objective function
    @objective(m, Max, v[iobj])
    # now, we solve the model
    solveNLP = JuMP.optimize!
    status = solveNLP(m)

    mi_out=JuMP.value.(v[iobj])
    vac_out=JuMP.value.(v[iac])
    vo_out=JuMP.value.(v[io2])
    vg_out=JuMP.value.(v[iglu])
    vall=JuMP.value.(v[:])

return mi_out,vac_out, vg_out, vo_out
end
#end
