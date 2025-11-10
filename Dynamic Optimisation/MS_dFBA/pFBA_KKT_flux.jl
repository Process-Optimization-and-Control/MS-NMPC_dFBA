function pFBA_KKT_flux(c0,c_ini,v0)

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
    set_attribute(m,"print_level", 5)
    set_attribute(m,"linear_solver", "ma27")
    set_attribute(m,"warm_start_init_point", "yes")

    #m = Model(with_optimizer(Ipopt.Optimizer, warm_start_init_point = "yes", print_level = 5, linear_solver = "ma27",tol=1e-4,dual_inf_tol=1e20,acceptable_dual_inf_tol=1e20,acceptable_iter=5,acceptable_tol=1e-2,mu_init=1e-8,
    #warm_start_bound_push=1e-8,warm_start_slack_bound_push=1e-8,warm_start_mult_bound_push=1e-8) )

                                           #
                                           #mu_init = 1e-3,
                                           #replace_bounds = "yes"

    # Set up variables for JuMP
    @variables(m, begin
        c[1:nc, 1:nfe, 1:ncp, 1:ns]     #States: 
        cdot[1:nc, 1:nfe, 1:ncp, 1:ns]  #Differential DifferentialEquations
        u[1:nmv,1:nfe, 1:ns]
        v[1:nv, 1:nfe, 1:ns] #fluxes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
        lambda[1:nm, 1:nfe, 1:ns]  #multipliers
        alpha_U[1:nv, 1:nfe, 1:ns]
        alpha_L[1:nv, 1:nfe, 1:ns]
        alpha_upt[1:nfe,1:2, 1:ns]

        FO_U[1:nv, 1:nfe, 1:ns]
        FO_L[1:nv, 1:nfe, 1:ns]
        FO_upt[1:nfe,1:2, 1:ns]

        hv[1:nfe, 1:ns]
   
    end)

    set_start_value.(c[:,:,:,1], c_ini)     
    set_start_value.(c[:,:,:,2], c_ini) 
    set_start_value.(v[:,:,1], v0)   
    set_start_value.(v[:,:,2], v0)  
    set_start_value.(hv[:,1], hm)
    set_start_value.(hv[:,2], hm)

     #scaling
   for i in 1:nc
    c0[i]=c0[i]/cs[i]

  end

    # Set up objective function

    @expressions(m, begin
    
    vg[i=1:nfe,s=1:ns], vg_max*( (cs[1]*c[1,i,3,s])/(Kg+(cs[1]*c[1,i,3,s]) + ((cs[1]*c[1,i,3,s])^2/Ki) ))
    vo[i=1:nfe,s=1:ns], vo_max* ( (cs[4]*c[4,i,3,s])/(Ko + (cs[4]*c[4,i,3,s])))
    
    end)
  
  @objective(m, Min, sum(-phi*c[3,nfe,3,s] + sum( sum(-phi1*FO_L[mc,i,s] + -phi3*FO_U[mc,i,s] for mc in 1:nv) + sum(phi2*FO_upt[i,j,s] for j in 1:2) for i in 1:nfe) for s in 1:ns))
  #@objective(m, Min, sum( sum( sum(-phi1*FO_L[mc,i,s] + -phi3*FO_U[mc,i,s] for mc in 1:nv) + sum(phi2*FO_upt[i,j,s] for j in 1:2) for i in 1:nfe) for s in 1:ns))

  #@objective(m, Min, sum( sum(-phi1*FO_L[mc,i] + -phi3*FO_U[mc,i] for mc in 1:nv) + sum(phi2*FO_upt[i,j] for j in 1:2) for i in 1:nfe))
  
  #@objective(m, Min, 1)
  
  
    #Set up the constraints
    @constraints(m, begin
        # set up differential equations
        m1[i=1:nfe, j=1:ncp,s=1:ns], cdot[1,i,j,s] == (0.18*vs[glu]*v[glu,i,s]*c[3,i,j,s]*cs[3] + (Sf - c[1,i,j,s]*cs[1])*(u[1,i,s]/c[5,i,j,s]*cs[5]))/cs[1]
        m2[i=1:nfe, j=1:ncp,s=1:ns], cdot[2,i,j,s] == ( 0.06*vs[ac]*v[ac,i,s]*c[3,i,j,s]*cs[3] - c[2,i,j,s]*cs[2]*(u[1,i,s]/(c[5,i,j,s]*cs[5])))/cs[2]
        m3[i=1:nfe, j=1:ncp,s=1:ns], cdot[3,i,j,s] == vs[obj]*v[obj,i,s]*c[3,i,j,s] - c[3,i,j,s]*(u[1,i,s]/c[5,i,j,s])
        m4[i=1:nfe, j=1:ncp,s=1:ns], cdot[4,i,j,s] ==  (vs[o2]*v[o2,i,s]*c[4,i,j,s]*cs[4]*c[5,i,j,s]*cs[5] + klo*(O_sat - c[4,i,j,s]*cs[4]))/cs[4]
        m5[i=1:nfe, j=1:ncp,s=1:ns], cdot[5,i,j,s] ==  u[1,i,s]

        # set up collocation equations - 2nd-to-nth point
        coll_c_n[l=1:nc, i=2:nfe, j=1:ncp,s=1:ns], c[l,i,j,s] == c[l,i-1,ncp,s] + hv[i,s]*sum(colmat[j,k]*cdot[l,i,k,s] for k in 1:ncp)
        
        # set up collocation equations - 1st point
        #coll_c_0[l=1:nc, i=1, j=1:ncp], c[l,i,j] == c0[l] + hv[i]*sum(colmat[j,k]*cdot[l,i,k] for k in 1:ncp)
        coll_c_0[l=1:nc, j=1:ncp,s=1:ns], c[l,1,j,s] == c0[l] + hv[1,s]*sum(colmat[j,k]*cdot[l,1,k,s] for k in 1:ncp)

        #------------- system constraints --------------#
        Sc[mc=1:nm,i=1:nfe,s=1:ns],  sum(S[mc,k]*v[k,i,s]*vs[k] for k in 1:nv) == 0
        v_UB[mc=1:nv,i=1:nfe,s=1:ns], v[mc,i,s]*vs[mc] - vub[mc]  <= 0
        v_LB[mc=1:nv,i=1:nfe,s=1:ns], -v[mc,i,s]*vs[mc]  + vlb[mc] <= 0

        c_LB[mc=1:nc,i=1:nfe, j=1:ncp,s=1:ns], -c[mc,i,j,s]  <= 0.0
        
        #Process constriants
        V_UB[mc=5:5,i=1:nfe, j=1:ncp,s=1:ns], c[mc,i,j,s]  <= V_max
        Ac_UB[mc=2:2,i=1:nfe, j=1:ncp,s=1:ns], c[mc,i,j,s]  <= Ac_max
        u_LB[i=1:nfe,s=1:ns], -u[1,i,s]  <= u_min
        u_UB[i=1:nfe,s=1:ns],  u[1,i,s]  <= u_max
        u_rate01[s=1:ns],  u[1,1,s] - u_ini  <= u_max_rate
        u_rate02[s=1:ns],  u[1,1,s] - u_ini  >= -u_max_rate
        u_rate_UB[i=2:nfe,s=1:ns],  u[1,i,s] - u[1,i-1,s]  <= u_max_rate
        u_rate_LB[i=2:nfe,s=1:ns],  u[1,i,s] - u[1,i-1,s]  >= -u_max_rate

        #non-anticipativy constraints
        u_NAC1,  u[1,1,1] - u[1,1,2]  <= 1e-4    
        u_NAC2,  u[1,1,1] - u[1,1,2]  >= -1e-4                                                                                                                                 

        MFE1[s=1:ns], sum(hv[i,s] for i in 1:nfe) == th
        MFE3[i=1:nfe,s=1:ns], hv[i,s]  >= 0.0
        MFE4[i=1:nfe,s=1:ns], hv[i,s]  >= (1-var_h)*hm[1]
        MFE5[i=1:nfe,s=1:ns], hv[i,s]  <= (1+var_h)*hm[1]

        #pFBA
        Lagr[mc=1:nv,i=1:nfe,s=1:ns], +d[mc] + w[mc,s]*v[mc,i,s]*vs[mc]  + alpha_L[mc,i,s] + alpha_U[mc,i,s] + up[mc]*alpha_upt[i,1,s] + up2[mc]*alpha_upt[i,2,s] + sum(S[k,mc]*lambda[k,i,s] for k in 1:nm) == 0
        # Lagr[mc=1:nv,i=1:nfe], +d[mc] + w[mc] + alpha_L[mc,i] + alpha_U[mc,i] + up[mc]*alpha_upt[i] + sum(S[k,mc]*lambda[k,i] for k in 1:nm) == 0
        #FBA
        #Lagr[mc=1:nv,i=1:nfe], +d[mc] + alpha_L[mc,i] + alpha_U[mc,i] + up[mc]*alpha_upt[i,1] + up2[mc]*alpha_upt[i,2] + sum(S[k,mc]*lambda[k,i] for k in 1:nm) == 0

        alpha1_LB[mc=1:nv,i=1:nfe,s=1:ns], alpha_L[mc,i,s] <= 0
        alpha4_LB[i=1:nfe,j=1:2,s=1:ns], alpha_upt[i,j,s] <= 0
        alpha1_UB[mc=1:nv,i=1:nfe,s=1:ns], alpha_U[mc,i,s] >= 0    
        
    end)

    @constraints(m, begin
    
      v_LB_g[i=1:nfe,s=1:ns], -v[glu,i,s]*vs[glu]  - vg[i,s] <= 0  
      v_LB_o[i=1:nfe,s=1:ns], -v[o2,i,s]*vs[o2]  - vo[i,s] <= 0  
      
      FO1[mc=1:nv,i=1:nfe,s=1:ns], FO_L[mc,i,s] == (v[mc,i,s]*vs[mc] -vlb[mc])*alpha_L[mc,i,s]
      FO2[mc=1:nv,i=1:nfe,s=1:ns], FO_U[mc,i,s] == (v[mc,i,s]*vs[mc] -vub[mc])*alpha_U[mc,i,s]
      FO3_upt[i=1:nfe,s=1:ns], FO_upt[i,1,s] == (-v[glu,i,s]*vs[glu] -vg[i,s])*alpha_upt[i,1,s]
      FO4_upt[i=1:nfe,s=1:ns], FO_upt[i,2,s] == (-v[o2,i,s]*vs[o2] -vo[i,s])*alpha_upt[i,2,s]

    end)

    #−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−#
    # Solve the model
    #tick()
    solveNLP = JuMP.optimize!
    status = solveNLP(m)
    #tock()

    # Get values for plotting
    uStar = JuMP.value.(u[:,:,:])
    cStar = JuMP.value.(c[:,:,:,:]).*cs # time series for plotting
    vStar = JuMP.value.(v[:,:,:]).*vs
    cdStar = JuMP.value.(cdot[:,:,:,:])
    lStar = JuMP.value.(lambda[:,:,:])
    alStar = JuMP.value.(alpha_L[:,:,:])
    agStar = JuMP.value.(alpha_upt[:,:,:])
    hStar = JuMP.value.(hv[:,:])

return  uStar, cStar, vStar, cdStar, lStar, alStar, agStar, hStar

end