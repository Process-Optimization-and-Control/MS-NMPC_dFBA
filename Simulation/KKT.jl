function KKT(vo)

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # Collocation parameters and radau time series
    #N = M^-1 =
    colmat = [0.19681547722366  -0.06553542585020 0.02377097434822;
              0.39442431473909  0.29207341166523 -0.04154875212600;
              0.37640306270047  0.51248582618842 0.11111111111111]
    radau  = [0.15505 0.64495 1.00000]

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # JuMP model
    m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma57") )
                                           #
                                           #mu_init = 1e-3,
                                           #replace_bounds = "yes"

    # Set up variables for JuMP
    @variables(m, begin                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
        #dv[1:nv, 1:nfe, 1:ncp]
        v[1:nv]
        lambda[1:nm]  #multipliers
       alpha_U[1:nv]
        alpha_L[1:nv]
        #alpha_upt[1:nv]
   
    end)

    # Set up initial guesses for solver

            for k in 1:nv
            #set_start_value(v[k], vo[k])      #v
             
              #set_start_value(alpha_U[k], auo[k])  
           # set_start_value(alpha_L[k], alo[k])  
            end
            
            for k in 1:nm
            #set_start_value(lambda[k], lo[k]) 
             end

    # Set up objective function
    #@NLobjective(m, Max, c[3,end,end])
    #@NLobjective(m, Max, sum(
     #                       sum(c[3,mc,i]/(c0[3]*exp(misc*ts[mc+1])) for mc in 1:nfe) for i in 1:ncp))
    #@NLobjective(m, Max, sum(
     #   sum(c[3,mc,i] for mc in 1:nfe) for i in 1:ncp))

    #@NLobjective(m, Min, sum(((v[mc]-vlb[mc])*alpha_L[mc])^2 + ((v[mc]-vub[mc])*alpha_U[mc])^2 for mc in 1:nv))

    @NLobjective(m, Min, sum( -((v[mc]-vlb[mc])*alpha_L[mc]) - ((v[mc]-vub[mc])*alpha_U[mc]) for mc in 1:nv))

    #@NLobjective(m, Min, sum(((v[mc]-vlb[mc])*alpha_L[mc])^2 for mc in 1:nv))

   #@objective(m, Min, 1)
    
    #@NLobjective(m, Min,  -v[obj] + phi*sum(((v[mc]-vlb[mc])*alpha_L[mc])^2 for mc in 1:nv))
    #@NLobjective(m, Max, v[obj])

    #Set up the constraints
    @constraints(m, begin
        #------------- system constraints --------------#
        Sc[mc=1:nm],  sum(S[mc,k]*v[k] for k in 1:nv) == 0
        v_UB[mc=1:nv], v[mc] - vub[mc] <= 0
        v_LB[mc=1:nv], -v[mc] + vlb[mc] <= 0
        #M4[i=1:nv], FO[i] >= v[i] 
        #M5[i=1:nv], FO[i] >= -v[i]

        #Lagr[mc=1:nv], +d[mc] + alpha_L[mc] + alpha_U[mc] + sum(S[k,mc]*lambda[k] for k in 1:nm) == 0
        #Lagr[mc=1:nv], +d[mc] + alpha_L[mc] + sum(S[k,mc]*lambda[k] for k in 1:nm) == 0

        alpha1_LB[mc=1:nv], alpha_L[mc] <= 0.0
        alpha1_UB[mc=1:nv], alpha_U[mc] >= 0.0
        
    
    end)

    @NLconstraints(m, begin

    Lagr[mc=1:nv], +d[mc] + w*v[mc]/(sqrt(v[mc]^2)+ 1e-12) + alpha_L[mc] + alpha_U[mc] + sum(S[k,mc]*lambda[k] for k in 1:nm) == 0
    
    end)

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # Solve the model
    solveNLP = JuMP.optimize!
    status = solveNLP(m)

    # Get values for plotting
    vStar = JuMP.value.(v[:])
    lStar = JuMP.value.(lambda[:])
    #auStar = JuMP.value.(alpha_U[:])
    alStar = JuMP.value.(alpha_L[:])
    #agStar = JuMP.value.(alpha_upt[:,:,3])

return   vStar, lStar, alStar

end