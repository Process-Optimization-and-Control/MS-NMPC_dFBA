#FBA
function FBA(vg_f,vo_f)

     #index

     #input
     vlb[glu]=-vg_f
     # ub[glu]=lb[glu]
      vub[glu]=0.0
      vlb[obj]=0
      vub[obj]=1000
      vlb[o2]=-vo_f
     #ub[o2]=lb[o2]
      #ub[o2]=0.0


     # declare the JuMP model
    # m = Model(with_optimizer(Cbc.Optimizer,print_level = 0))
    # m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27") )
     #m = Model(with_optimizer(Ipopt.Optimizer,logLevel = 0))
      m = Model(Cbc.Optimizer)
      set_attribute(m, "logLevel", 0)
     # m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 1, max_iter = 200000) )


     # variables representing x1 and x2
     @variables(m, begin
     v[1:nv]
     end)

     @constraints(m, begin
     M1[i=1:nm],  sum(S[i,k]*v[k] for k in 1:nv) == 0
     M2[i=1:nv], v[i] <= vub[i]
     M3[i=1:nv], v[i] >= vlb[i]

     #M2[i=act_ub], v[i] <= vub[i]
     #M3[i=act_lb], v[i] >= vlb[i]
     end)

    # and we declare the objective function
    @objective(m, Max, v[obj])
    # now, we solve the model
    solveNLP = JuMP.optimize!
    status = solveNLP(m)

    mi_out=JuMP.value.(v[obj])
    vac_out=JuMP.value.(v[ac])
    vo_out=JuMP.value.(v[o2])
    vg_out=JuMP.value.(v[glu])
    vall=JuMP.value.(v[:])

return mi_out,vac_out, vg_out, vo_out
end
#end
