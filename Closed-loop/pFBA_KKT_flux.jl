function pFBA_KKT_flux(c0,c_ini,u_ini,v0,u_0,nfe_n,th_n)
    

    

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # Collocation parameters and radau time series
    #N = M^-1 =
    colmat = [0.19681547722366  -0.06553542585020 0.02377097434822;
              0.39442431473909  0.29207341166523 -0.04154875212600;
              0.37640306270047  0.51248582618842 0.11111111111111]
    radau  = [0.15505 0.64495 1.00000]

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # JuMP model
    #m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27") )

    m = Model(Ipopt.Optimizer)
    set_attribute(m,"print_level", 0)
    set_attribute(m,"linear_solver", "ma27")
    set_attribute(m,"warm_start_init_point", "yes")
   # set_attribute(m,"mu_init", 1e-8)
   # set_attribute(m,"warm_start_bound_push", 1e-8)
   # set_attribute(m,"warm_start_slack_bound_push", 1e-8)
   # set_attribute(m,"warm_start_mult_bound_push", 1e-8)
    

    #m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27",tol=1e-4,dual_inf_tol=1e20,acceptable_dual_inf_tol=1e20,acceptable_iter=5,acceptable_tol=1e-2,mu_init=1e-8,
    #warm_start_bound_push=1e-8,warm_start_slack_bound_push=1e-8,warm_start_mult_bound_push=1e-8) )

                                           #
                                           #mu_init = 1e-3,
                                           #replace_bounds = "yes"

    # Set up variables for JuMP
    @variables(m, begin
        c[1:nc, 1:nfe_n, 1:ncp]     #States: 
        cdot[1:nc, 1:nfe_n, 1:ncp]  #Differential DifferentialEquations
        u_o[1:nmv,1:nfe_n]
        v[1:nv, 1:nfe_n] #fluxes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
        lambda[1:nm, 1:nfe_n]  #multipliers
        alpha_U[1:nv, 1:nfe_n]
        alpha_L[1:nv, 1:nfe_n]
        alpha_upt[1:nfe_n,1:2]

        FO_U[1:nv, 1:nfe_n]
        FO_L[1:nv, 1:nfe_n]
        FO_upt[1:nfe_n,1:2]

        s_UB[1:nfe_n]
        s_V[1:nfe_n, 1:ncp]
        s_Ac[i=1:nfe_n, 1:ncp]

        hv[1:nfe_n]
   
    end)

    
    #set_start_value.(u_o, u_ini) 
    set_start_value.(c, c_ini)     
    set_start_value.(v, v0)   
    set_start_value.(hv, hm)

     #scaling
   for i in 1:nc
    c0[i]=c0[i]/cs[i]

  end

    # Set up objective function

    @expressions(m, begin
    
    vg_o[i=1:nfe_n], vg_max*( (cs[1]*c[1,i,3])/(Kg+(cs[1]*c[1,i,3]) + ((cs[1]*c[1,i,3])^2/Ki) ))
    vo_o[i=1:nfe_n], vo_max* ( (cs[4]*c[4,i,3])/(Ko + (cs[4]*c[4,i,3])))
    
    end)
  
  
 # @objective(m, Min, -phi*c[3,nfe_n,3] + sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + sum(phi2*FO_upt[i,j] + phi4*s_V[i,j] + phi4*s_Ac[i,j] for j in 1:2) + phi4*s_UB[i] for i in 1:nfe_n))

 @objective(m, Min, -phi*c[3,nfe_n,3] + sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + sum(phi2*FO_upt[i,j] + phi4*s_Ac[i,j] for j in 1:2)  for i in 1:nfe_n))

  #@objective(m, Min, -phi*c[3,nfe_n,3] + sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + sum(phi2*FO_upt[i,j] for j in 1:2) for i in 1:nfe_n))

  #@objective(m, Min, sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + sum(phi2*FO_upt[i,j] for j in 1:2) for i in 1:nfe_n))
  
  #@objective(m, Min, 1)
  
  
    #Set up the constraints
    @constraints(m, begin
        # set up differential equations
        m1[i=1:nfe_n, j=1:ncp], cdot[1,i,j] == (0.18*vs[glu]*v[glu,i]*c[3,i,j]*cs[3] + (Sf - c[1,i,j]*cs[1])*(u_o[1,i]/c[5,i,j]*cs[5]))/cs[1]
        m2[i=1:nfe_n, j=1:ncp], cdot[2,i,j] == ( 0.06*vs[ac]*v[ac,i]*c[3,i,j]*cs[3] - c[2,i,j]*cs[2]*(u_o[1,i]/(c[5,i,j]*cs[5])))/cs[2]
        m3[i=1:nfe_n, j=1:ncp], cdot[3,i,j] == vs[obj]*v[obj,i]*c[3,i,j] - c[3,i,j]*(u_o[1,i]/c[5,i,j])
        m4[i=1:nfe_n, j=1:ncp], cdot[4,i,j] ==  (vs[o2]*v[o2,i]*c[4,i,j]*cs[4]*c[5,i,j]*cs[5] + klo*(O_sat - c[4,i,j]*cs[4]))/cs[4]
        m5[i=1:nfe_n, j=1:ncp], cdot[5,i,j] ==  u_o[1,i]

        # set up collocation equations - 2nd-to-nth point
        coll_c_n[l=1:nc, i=2:nfe_n, j=1:ncp], c[l,i,j] == c[l,i-1,ncp] + hv[i]*sum(colmat[j,k]*cdot[l,i,k] for k in 1:ncp)
        
        # set up collocation equations - 1st point
        coll_c_0[l=1:nc, j=1:ncp], c[l,1,j] == c0[l] + hv[1]*sum(colmat[j,k]*cdot[l,1,k] for k in 1:ncp)

        #------------- system constraints --------------#
        Sc[mc=1:nm,i=1:nfe_n],  sum(S[mc,k]*v[k,i]*vs[k] for k in 1:nv) == 0
        v_UB[mc=1:nv,i=1:nfe_n], v[mc,i]*vs[mc] - vub_F[mc]  <= 0
        v_LB[mc=1:nv,i=1:nfe_n], -v[mc,i]*vs[mc]  + vlb_F[mc] <= 0

        c_LB[mc=1:nc,i=1:nfe_n, j=1:ncp], -c[mc,i,j]  <= 0.0
        
        #Process constriants
        #V_UB[mc=5:5,i=1:nfe_n, j=1:ncp], c[mc,i,j]  <= V_max + s_V[i,j]
        Ac_UB[mc=2:2,i=1:nfe_n, j=1:ncp], c[mc,i,j]  <= Ac_max + s_Ac[i,j]
        V_UB[mc=5:5,i=1:nfe_n, j=1:ncp], c[mc,i,j]  <= V_max 
        #Ac_UB[mc=2:2,i=1:nfe_n, j=1:ncp], c[mc,i,j]  <= Ac_max 
        u_LB[i=1:nfe_n], -u_o[1,i]  <= u_min
       # u_UB[i=1:nfe_n],  u_o[1,i]  <= u_max + s_UB[i]
        u_UB[i=1:nfe_n],  u_o[1,i]  <= u_max 
        u_rate01,  u_o[1,1] - u_0  <= u_max_rate
        u_rate02,  u_o[1,1] - u_0  >= -u_max_rate
        u_rate_UB[i=2:nfe_n],  u_o[1,i] - u_o[1,i-1]  <= u_max_rate
        u_rate_LB[i=2:nfe_n],  u_o[1,i] - u_o[1,i-1]  >= -u_max_rate

        #s_V_UB[i=1:nfe_n, j=1:ncp],  -s_V[i,j] <= 0.0 
       s_Ac_UB[i=1:nfe_n, j=1:ncp],  -s_Ac[i,j] <= 0.0
        #s_UB_UB[i=1:nfe_n],  -s_UB[i] <= 0.0

        MFE1, sum(hv[i] for i in 1:nfe_n) == th_n
        MFE3[i=1:nfe_n], hv[i]  >= 0.0
        MFE4[i=1:nfe_n], hv[i]  >= (1-var_h)*hm[1]
        MFE5[i=1:nfe_n], hv[i]  <= (1+var_h)*hm[1]

        #pFBA
        # Lagr[mc=1:nv,i=1:nfe_n], +d[mc] + w[mc]*v[mc,i]*vs[mc]  + alpha_L[mc,i] + alpha_U[mc,i] + up[mc]*alpha_upt[i] + sum(S[k,mc]*lambda[k,i] for k in 1:nm) == 0
        # Lagr[mc=1:nv,i=1:nfe_n], +d[mc] + w[mc] + alpha_L[mc,i] + alpha_U[mc,i] + up[mc]*alpha_upt[i] + sum(S[k,mc]*lambda[k,i] for k in 1:nm) == 0
        #FBA
        Lagr[mc=1:nv,i=1:nfe_n], +d[mc] + alpha_L[mc,i] + alpha_U[mc,i] + up[mc]*alpha_upt[i,1] + up2[mc]*alpha_upt[i,2] + sum(S[k,mc]*lambda[k,i] for k in 1:nm) == 0

        alpha1_LB[mc=1:nv,i=1:nfe_n], alpha_L[mc,i] <= 0
        alpha4_LB[i=1:nfe_n,j=1:2], alpha_upt[i,j] <= 0
        alpha1_UB[mc=1:nv,i=1:nfe_n], alpha_U[mc,i] >= 0    
        
    end)

    @constraints(m, begin
    
      v_LB_g[i=1:nfe_n], -v[glu,i]*vs[glu]  - vg_o[i] <= 0  
      v_LB_o[i=1:nfe_n], -v[o2,i]*vs[o2]  - vo_o[i] <= 0  
      
      FO1[mc=1:nv,i=1:nfe_n], FO_L[mc,i] == (v[mc,i]*vs[mc] -vlb_F[mc])*alpha_L[mc,i]
      FO2[mc=1:nv,i=1:nfe_n], FO_U[mc,i] == (v[mc,i]*vs[mc] -vub_F[mc])*alpha_U[mc,i]
      FO3_upt[i=1:nfe_n], FO_upt[i,1] == (-v[glu,i]*vs[glu] -vg_o[i])*alpha_upt[i,1]
      FO4_upt[i=1:nfe_n], FO_upt[i,2] == (-v[o2,i]*vs[o2] -vo_o[i])*alpha_upt[i,2]

    end)

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # Solve the model
    #tick()
    solveNLP = JuMP.optimize!
    status = solveNLP(m)
    #tock()

    # Get values for plotting
    uStar = JuMP.value.(u_o[:,:])
    cStar = JuMP.value.(c[:,:,:]).*cs # time series for plotting
    vStar = JuMP.value.(v[:,:]).*vs
    cdStar = JuMP.value.(cdot[:,:,:])
    lStar = JuMP.value.(lambda[:,:])
    alStar = JuMP.value.(alpha_L[:,:])
    agStar = JuMP.value.(alpha_upt[:,:])
    hStar = JuMP.value.(hv[:])
    flag_star = termination_status(m)
    time_star = solve_time(m)
    num_star = num_variables(m)

return  uStar, cStar, vStar, flag_star, time_star, num_star

end