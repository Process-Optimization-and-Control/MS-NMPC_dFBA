#FBA
function FBA(vg)

     #index
  

     #input
     vlb[glu]=-vg
     # ub[glu]=lb[glu]
      vub[glu]=0.0
      vlb[obj]=0
     #lb[o2]=-vo
     #ub[o2]=lb[o2]
      #ub[o2]=0.0


     # declare the JuMP model
     m = Model(with_optimizer(Cbc.Optimizer,print_level = 0))
    # m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27") )
     #m = Model(with_optimizer(Ipopt.Optimizer,logLevel = 0))
     # m = Model(Cbc.Optimizer)
      #set_optimizer_attributes(m, "print_level" => 1)
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

    mi=JuMP.value.(v[obj])
    vac=JuMP.value.(v[ac])
    vall=JuMP.value.(v[:])

return mi,vac
end
#end
