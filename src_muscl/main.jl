using ProgressMeter
using Dates

function main()
    
    # Start of time measurement
    start_t = now()

    out_dir  = "result"         # output dir
    PARAMDAT = "PARAMDAT.json"  # directory
    fwrite   = "write"          # write file 
    
    nval = 4         # number of conserved variables
    Rd   = 287.0     # gas constant of air, J/(kg K)
    R    = 8.314     # gas constant, J/(K mol)
    icell = 2        # imaginary cell
    
    # read grids and parameter
    xmax, ymax, nodes, vecAx, vecAy = read_allgrid()
    out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
    time_integ, nt, dt, every_outnum, in_nt, dtau, cfl, ad_scheme,
    init_rho, init_u, init_v, init_p, init_T, specific_heat_ratio, Rd, bdcon = input_para(PARAMDAT)
    
    # number of cells
    cellxmax = xmax - 1
    cellymax = ymax - 1
    
    # allocation
    Qbase, Qbase_ave, volume, cellcenter, wally, yplus, dx, dy, Qcon, Qcon_hat, mu, mut, mut_bd, lambda, 
    E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, RHS, QbaseU, QbaseD, QbaseL, QbaseR,
    QconU, QconD, QconL, QconR = common_allocation(cellxmax, cellymax, nval)
    
    # set initial condition
    Qbase, restartnum = set_initQbase(Qbase, cellxmax, cellymax, restart_file, init_rho, init_u, init_v, init_p, init_T,
                                      specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, Rd, nval, icell)
    
    # set initial condition for imaginary cell
    Qbase    = set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, nval, icell)

    # set volume, dx and dy
    volume = set_volume(nodes, cellxmax, cellymax, volume)
    cellcenter = set_cellcenter(cellcenter, nodes, cellxmax, cellymax)
    dx, dy = set_dx_lts(dx, dy, nodes, cellxmax, cellymax, icell)
    reset_write(fwrite)

    # wally
    wally, swith_wall = set_wally(nodes, bdcon, wally, cellcenter, cellxmax, cellymax, icell)

    # AUSM+up
    Minf = 0.0
    if ad_scheme == 2
        Minf = set_Minf(bdcon, specific_heat_ratio, Rd, nval)
    end

    # write number of threads
    print("threads num : ")
    println(Threads.nthreads())

    # check boundary condition
    check_bd(bdcon)

    output_result_yplus(0, Qbase, wally, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval, icell)
   
        
    # main loop
    loop_ite = 0
    if time_integ == "1"
        # exlicit scheme
        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)

            loop_ite += 1
            
            # step number
            evalnum = t + restartnum
            
            # set conserved variables in the general coordinate system
            Qbase    = set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, nval, icell)
            Qcon     = base_to_conservative(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
            Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
            
            # muscl
            QbaseU, QbaseD, QbaseL, QbaseR, 
            QconU, QconD, QconL, QconR = muscl(Qbase, QbaseU, QbaseD, QbaseL, QbaseR, 
                                                QconU, QconD, QconL, QconR, cellxmax, cellymax, vecAx, vecAy,
                                                nval, icell, specific_heat_ratio, volume, nodes)
            
            # set viscosity and thermal Conductivity
            mu     = set_mu(mu, Qbase, cellxmax, cellymax, specific_heat_ratio, Rd)
            lambda = set_lambda(lambda, Qbase, cellxmax, cellymax, mu, specific_heat_ratio, Rd)
            
            # yplus
            yplus = cal_yplus(yplus, Qbase, wally, swith_wall, mu, cellxmax, cellymax, vecAx, vecAy, volume, icell)
            #yplus = ones(cellxmax, cellymax)*100                  # yplus
            
            # advection_term
            E_adv_hat, F_adv_hat = AUSM(E_adv_hat, F_adv_hat, QbaseU, QbaseD, QbaseL, QbaseR, 
                                    QconU, QconD, QconL, QconR, cellxmax, cellymax, vecAx, vecAy, 
                                    specific_heat_ratio, volume, nval, Minf, ad_scheme, icell)
                        
            # viscos_term
            E_vis_hat, F_vis_hat, mut = central_diff(E_vis_hat, F_vis_hat, QbaseU, QbaseD, QbaseL, QbaseR, 
                                                QconU, QconD, QconL, QconR, cellxmax, cellymax, mu, mut, mut_bd,lambda,
                                                vecAx, vecAy, specific_heat_ratio, volume, Rd, nval, yplus, swith_wall, icell)
            
            #println(" fff ")
            #println(E_adv_hat[150,150,:])
            #println(E_adv_hat[20,20,:])
            #println(F_adv_hat[23,3,:])
            #println(E_vis_hat[150,150,:])
            #println(E_vis_hat[20,20,:])
            #println(F_vis_hat[23,3,:])
            

            #throw(UndefVarError(:x))
            
            #println(yplus[:,1])
            #println(" yplus ")
            #println(yplus[1000,2])
            #println(yplus[1000,3])
            #println(mu[1000,2] / Qbase[1000,2,1])
            #println(yplus[:,4])
            #println(yplus[:,5])
            #println(wally[1000,2])
            #println(wally[1000,3])
            #throw(UndefVarError(:x))
            

            # RHS
            RHS = setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume, icell)
            
            # time integral
            Qcon_hat = time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, nval, icell)
            
            # calculate primitive variables
            Qcon  = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
            Qbase = conservative_to_base(Qbase, Qcon, cellxmax, cellymax, specific_heat_ratio)
            Qbase_ave = cal_Qave(Qbase, Qbase_ave, cellxmax, cellymax, nval)
            
            # output
            if round(evalnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(evalnum)))
                #output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval, icell)
                output_result_yplus(evalnum, Qbase, yplus, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval, icell)
                output_ave(Qbase_ave, cellxmax, cellymax, out_file_front, out_ext, out_dir, Rd, nval, loop_ite, icell)
            end
            
            # Find out if the results were divergent
            check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite, icell)
        end
    elseif time_integ == "2"
        # implicit

        # allocation for implicit
        Qbasen, Qconn, Qconn_hat, Qbasem, dtau, lambda_facex, lambda_facey,
        A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig,
        jalphaP, jbetaP, delta_Q, delta_Q_temp, D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp, res,
        norm2, I = allocation_implicit(cellxmax, cellymax, nval)

        norm_rho = 0.0

        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)
            
            # step number
            evalnum = t + restartnum
            
            # write physical time 
            output_physicaltime(fwrite, t, dt)
            
            # copy
            for l in 1:nval
                for j in 1:cellymax
                    for i in 1:cellxmax
                        Qbasen[i,j,l] = Qbase[i,j,l]
                        Qbasem[i,j,l] = Qbase[i,j,l]
                    end
                end
            end

            # set conserved variables in the general coordinate system for inner iteration
            Qconn     = base_to_conservative(Qbasen, Qconn, cellxmax, cellymax, specific_heat_ratio)
            Qconn_hat = setup_Qcon_hat(Qconn, Qconn_hat, cellxmax, cellymax, volume, nval)     

            # start inner iteration
            for tau in 1:in_nt
                # set conserved variables in the general coordinate system
                Qbasem   = set_boundary(Qbasem, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, specific_heat_ratio, nval, icell)
                Qcon     = base_to_conservative(Qbasem, Qcon, cellxmax, cellymax, specific_heat_ratio)
                Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)

                # muscl
                QbaseU, QbaseD, QbaseL, QbaseR, 
                QconU, QconD, QconL, QconR = muscl(Qbasem, QbaseU, QbaseD, QbaseL, QbaseR, 
                                                    QconU, QconD, QconL, QconR, cellxmax, cellymax, vecAx, vecAy,
                                                    nval, icell, specific_heat_ratio, volume, nodes)

                # set viscosity and thermal Conductivity
                mu     = set_mu(mu, Qbasem, cellxmax, cellymax, specific_heat_ratio, Rd)
                lambda = set_lambda(lambda, Qbasem, cellxmax, cellymax, mu, specific_heat_ratio, Rd)
                
                yplus = cal_yplus(yplus, Qbasem, wally, swith_wall, mu, cellxmax, cellymax, vecAx, vecAy, volume, icell)
            
                #throw(UndefVarError(:x))
                
                # advection_term
                E_adv_hat, F_adv_hat = AUSM(E_adv_hat, F_adv_hat, QbaseU, QbaseD, QbaseL, QbaseR, 
                                            QconU, QconD, QconL, QconR, cellxmax, cellymax, vecAx, vecAy, 
                                            specific_heat_ratio, volume, nval, Minf, ad_scheme, icell)

                # viscos_term
                E_vis_hat, F_vis_hat, mut = central_diff(E_vis_hat, F_vis_hat, QbaseU, QbaseD, QbaseL, QbaseR, 
                                        QconU, QconD, QconL, QconR, cellxmax, cellymax, mu, mut, mut_bd,lambda,
                                        vecAx, vecAy, specific_heat_ratio, volume, Rd, nval, yplus, swith_wall, wally, icell)
                
                # RHS
                RHS = setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume, icell)
                
                #println(E_adv_hat[3,3,:])
                #println(F_adv_hat[3,3,:])
                #println(E_vis_hat[3,3,:])
                #println(F_vis_hat[3,3,:])
                #println(RHS[50,3,:])
                #throw(UndefVarError(:x))
                
                #throw(UndefVarError(:x))


                #println(E_adv_hat[150,150,:])
                #println(E_adv_hat[20,20,:])
                #println(F_adv_hat[23,3,:])
                #println(E_vis_hat[150,150,:])
                #println(E_vis_hat[20,20,:])
                #println(F_vis_hat[23,3,:])
                
                
                # 粘性の更新
                for i in 1:cellxmax
                    for j in 1:cellymax
                        mu[i,j] = mu[i,j] + mut[i,j]
                    end
                end

                # set inner time step by local time stepping
                dtau   = set_lts(dtau, lambda_facex, lambda_facey, Qbasem, cellxmax, cellymax, mu, dx, dy,
                                vecAx, vecAy, volume, specific_heat_ratio, cfl, icell)
                
                # lusgs_advection_term
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, 
                A_beta_shig, B_beta_shig = one_wave(A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, A_beta_shig, B_beta_shig, I,
                                                    Qbasem, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
                # lusgs_viscos_term
                jalphaP, jbetaP = central_diff_jacobian(jalphaP, jbetaP, Qbasem, Qcon, cellxmax, cellymax, mu, lambda,
                                                        vecAx, vecAy, specific_heat_ratio, volume, nval)

                #println(RHS[2,:,1])
                #println(A_adv_hat_m[5,5,:,:])
                #println(jalphaP[5,5])

                #throw(UndefVarError(:x))
                
                # LUSGS
                ite = 0
                while true
                    # copy delta_Q
                    for l in 1:nval
                        for j in 1:cellymax
                            for i in 1:cellxmax
                                delta_Q_temp[i,j,l] = delta_Q[i,j,l]
                            end
                        end
                    end
                                        
                    # Reversing the left-hand side by lusgs
                    delta_Q = lusgs(D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp, I, dt, dtau, Qcon_hat, Qconn_hat, delta_Q,
                                    A_adv_hat_p,  A_adv_hat_m,  B_adv_hat_p,  B_adv_hat_m,  A_beta_shig,  B_beta_shig,
                                     jalphaP,  jbetaP, RHS, cellxmax, cellymax, volume, nval, icell)
                    #println(delta_Q[5,5,:])
                    #println(delta_Q[3,3,:])
                    #println(dtau[3,3])
                    #throw(UndefVarError(:x))
                    
                    # cal Residuals by norm-2
                    res   = set_res(res, delta_Q, delta_Q_temp, cellxmax, cellymax, nval, icell)
                    norm2 = check_converge(res, RHS, cellxmax, cellymax, init_small, nval, icell)
                    
                    # Find out if the results converged
                    if norm2[1] < norm_ok && norm2[2] < norm_ok && norm2[3] < norm_ok && norm2[4] < norm_ok
                        break
                    end
                    if ite % 100 == 0
                        println(" now cal norm2 ")
                        println(norm2)
                        if isequal(norm2[1], NaN) == true
                            println("  ")
                            println(" norm2 = NaN ")
                            println(" stop cal ")
                            throw(UndefVarError(:x))
                        end
                    end

                    ite += 1
                end

                # output inner time
                output_innertime(fwrite, tau, norm2, nval)
                
                # Updating
                for l in 1:nval
                    for j in 2:cellymax-1
                        for i in 2:cellxmax-1
                            Qcon_hat[i,j,l] = Qcon_hat[i,j,l] + delta_Q[i,j,l]
                        end
                    end
                end
                
                # calculate primitive variables
                Qcon = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, volume, nval)
                Qbasem = conservative_to_base(Qbasem, Qcon, cellxmax, cellymax, specific_heat_ratio)

                # Find out if the results were divergent
                check_divrege(Qbasem, cellxmax, cellymax, Rd, fwrite, icell)

                # hantei inner ite 
                if tau == 1
                    norm_rho = norm2[1]
                elseif abs(norm2[1]/norm_rho) < 1.0e-3
                    break
                end
            end
            # End of the inner iteration

            # Updating in physical time
            for l in 1:nval
                for j in 1:cellymax
                    for i in 1:cellxmax
                        Qbase[i,j,l] = Qbasem[i,j,l]
                    end
                end
            end

            Qbase_ave = cal_Qave(Qbase_ave, Qbase, cellxmax, cellymax, nval, evalnum)
            
            # output
            if round(evalnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(evalnum)))
                # output_result(evalnum, Qbase, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval, icell)
                output_result_yplus(evalnum, Qbase, yplus, cellxmax, cellymax, specific_heat_ratio, out_file_front, out_ext, out_dir, Rd, nval, icell)

                output_ave(Qbase_ave, cellxmax, cellymax, out_file_front, out_ext, out_dir, Rd, nval, loop_ite, icell)
            end

            # Find out if the results were divergent
            check_divrege(Qbase, cellxmax, cellymax, Rd, fwrite, icell)
        end
    end
    
    # end of time measurement
    end_t = now()

    # output of calculation time
    output_fin(fwrite, start_t, end_t, nt, dt, in_nt, cellxmax, cellymax)
end


# -- main --
main()
#throw(UndefVarError(:x))