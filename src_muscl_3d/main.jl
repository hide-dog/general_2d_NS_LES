using ProgressMeter
using Dates

function main()
    
    # Start of time measurement
    start_t = now()

    out_dir  = "result"         # output dir
    PARAMDAT = "PARAMDAT.json"  # directory
    fwrite   = "write"          # write file 
    
    nval = 5         # number of conserved variables
    Rd   = 287.0     # gas constant of air, J/(kg K)
    R    = 8.314     # gas constant, J/(K mol)
    icell = 2        # imaginary cell
    
    # read grids and parameter
    xmax, ymax, zmax, nodes, vecAx, vecAy, vecAz = read_allgrid()
    out_file_front, out_ext, restartnum, restart_file, init_small, norm_ok,
    time_integ, nt, dt, every_outnum, in_nt, dtau, cfl, ad_scheme,
    init_rho, init_u, init_v, init_w, init_p, init_T, specific_heat_ratio, Rd, bdcon = input_para(PARAMDAT)
    
    # number of cells
    cellxmax = xmax - 1
    cellymax = ymax - 1
    cellzmax = zmax - 1
    
    # allocation
    Qbase, Qbase_ave, volume, cellcenter, wally, yplus, Qcon, Qcon_hat, mu, mut, mut_bd, lambda, 
    E_adv_hat, F_adv_hat, G_adv_hat, E_vis_hat, F_vis_hat, G_vis_hat, RHS, QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB,
    QconU, QconD, QconL, QconR, QconF, QconB = common_allocation(cellxmax, cellymax, cellzmax, nval)
    
    # set initial condition
    Qbase, restartnum = set_initQbase(Qbase, cellxmax, cellymax, cellzmax, restart_file, 
                                        init_rho, init_u, init_v, init_w, init_p, init_T,
                                        specific_heat_ratio, out_file_front, out_ext, out_dir, restartnum, Rd, nval, icell)
   
    # set initial condition for imaginary cell
    Qbase    = set_boundary(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, specific_heat_ratio, nval, icell)

    # set volume
    volume = set_volume(nodes, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, volume)
    cellcenter = set_cellcenter(cellcenter, nodes, cellxmax, cellymax, cellzmax)
    reset_write(fwrite)

    # wally
    wally, swith_wall = set_wally(nodes, bdcon, wally, cellcenter, cellxmax, cellymax, cellzmax, icell)

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
        
    # main loop
    loop_ite = 0
    if time_integ == "1"
        # exlicit scheme
        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)

            loop_ite += 1
            
            # step number
            stepnum = t + restartnum
            
            # set conserved variables in the general coordinate system
            Qbase    = set_boundary(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, 
                                    bdcon, Rd, specific_heat_ratio, nval, icell)
            Qcon     = base_to_conservative(Qbase, Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio)
            Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, cellzmax, volume, nval)
            
            # muscl
            QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB, 
            QconU, QconD, QconL, QconR, QconF, QconB = muscl(Qbase, QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB, 
                                                            QconU, QconD, QconL, QconR, QconF, QconB, 
                                                            cellxmax, cellymax, cellzmax, nval, icell, specific_heat_ratio, volume, nodes)
            
            # set viscosity and thermal Conductivity
            mu     = set_mu(mu, Qbase, cellxmax, cellymax, cellzmax, specific_heat_ratio, Rd)
            lambda = set_lambda(lambda, Qbase, cellxmax, cellymax, cellzmax, mu, specific_heat_ratio, Rd)
            
            # yplus
            yplus = cal_yplus(yplus, Qbase, wally, swith_wall, mu, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, volume, icell)
            
            # advection_term
            E_adv_hat, F_adv_hat, G_adv_hat = AUSM(E_adv_hat, F_adv_hat, G_adv_hat, QbaseF, QbaseB, QbaseU, QbaseD, QbaseL, QbaseR, 
                                                QconF, QconB, QconU, QconD, QconL, QconR, cellxmax, cellymax, cellzmax, 
                                                vecAx, vecAy, vecAz, specific_heat_ratio, volume, nval, Minf, ad_scheme, icell)
                        
            # viscos_term
            E_vis_hat, F_vis_hat, G_vis_hat, mut = central_diff(E_vis_hat, F_vis_hat, G_vis_hat, QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB, 
                                                            QconU, QconD, QconL, QconR, QconF, QconB, cellxmax, cellymax, cellzmax, mu, mut, mut_bd, lambda,
                                                            vecAx, vecAy, vecAz, specific_heat_ratio, volume, Rd, nval, yplus, swith_wall, icell)

            # RHS
            RHS = setup_RHS(RHS, cellxmax, cellymax, cellzmax, E_adv_hat, F_adv_hat, G_adv_hat, 
                            E_vis_hat, F_vis_hat, G_vis_hat, nval, volume, icell)
            
            # time integral
            Qcon_hat = time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, cellzmax, nval, icell)
            
            # calculate primitive variables
            Qcon  = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, cellzmax, volume, nval)
            Qbase = conservative_to_base(Qbase, Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio)
            Qbase_ave = cal_Qave(Qbase, Qbase_ave, cellxmax, cellymax, cellzmax, nval)
            
            # output
            if round(stepnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(stepnum)))
                output_result(stepnum, Qbase, cellxmax, cellymax, cellzmax, specific_heat_ratio, 
                                out_file_front, out_ext, out_dir, Rd, nval, icell)
                output_ave(Qbase_ave, cellxmax, cellymax, cellzmax, out_file_front, out_ext, out_dir, Rd, nval, loop_ite, icell)
            end
            
            # Find out if the results were divergent
            check_divrege(Qbase, cellxmax, cellymax, cellzmax, Rd, fwrite, icell)
        end
    elseif time_integ == "2"
        # implicit

        # allocation for implicit
        Qbasen, Qconn, Qconn_hat, Qbasem, dtau, lambda_facex, lambda_facey, lambda_facez, dx, dy, dz,
        A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, C_adv_hat_p, C_adv_hat_m, A_beta_shig, B_beta_shig, C_beta_shig,
        jalphaP, jbetaP, jgammaP, delta_Q, delta_Q_temp, D, Lx, Ly, Lz, Ux, Uy, Uz, LdQ, UdQ, RHS_temp, res,
        norm2, I = allocation_implicit(cellxmax, cellymax, cellzmax, nval)

        dx, dy, dz = set_dx_lts(dx, dy, dz, nodes, cellxmax, cellymax, cellzmax, icell)

        prog = Progress(nt,1)
        @time for t in 1:nt
            next!(prog)
            
            # step number
            loop_ite += 1
            stepnum = t + restartnum
            
            # write physical time 
            output_physicaltime(fwrite, t, dt)
            
            # copy
            for l in 1:nval
                for k in 1:cellzmax
                    for j in 1:cellymax
                        for i in 1:cellxmax
                            Qbasen[i,j,k,l] = Qbase[i,j,k,l]
                            Qbasem[i,j,k,l] = Qbase[i,j,k,l]
                        end
                    end
                end
            end

            # set conserved variables in the general coordinate system for inner iteration
            Qconn     = base_to_conservative(Qbasen, Qconn, cellxmax, cellymax, cellzmax, specific_heat_ratio)
            Qconn_hat = setup_Qcon_hat(Qconn, Qconn_hat, cellxmax, cellymax, cellzmax, volume, nval)
            

            # start inner iteration
            for tau in 1:in_nt
                # set conserved variables in the general coordinate system
                Qbasem    = set_boundary(Qbasem, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, 
                bdcon, Rd, specific_heat_ratio, nval, icell)
                Qcon     = base_to_conservative(Qbasem, Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio)
                Qcon_hat = setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, cellzmax, volume, nval)

                # muscl
                QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB, 
                QconU, QconD, QconL, QconR, QconF, QconB = muscl(Qbasem, QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB, 
                                                        QconU, QconD, QconL, QconR, QconF, QconB, 
                                                        cellxmax, cellymax, cellzmax, nval, icell, specific_heat_ratio, volume, nodes)

                # set viscosity and thermal Conductivity
                mu     = set_mu(mu, Qbasem, cellxmax, cellymax, cellzmax, specific_heat_ratio, Rd)
                lambda = set_lambda(lambda, Qbasem, cellxmax, cellymax, cellzmax, mu, specific_heat_ratio, Rd)
                
                # yplus
                yplus = cal_yplus(yplus, Qbasem, wally, swith_wall, mu, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, volume, icell)
                
                # advection_term
                E_adv_hat, F_adv_hat, G_adv_hat = AUSM(E_adv_hat, F_adv_hat, G_adv_hat, QbaseF, QbaseB, QbaseU, QbaseD, QbaseL, QbaseR, 
                                                    QconF, QconB, QconU, QconD, QconL, QconR, cellxmax, cellymax, cellzmax, 
                                                    vecAx, vecAy, vecAz, specific_heat_ratio, volume, nval, Minf, ad_scheme, icell)
                            
                # viscos_term
                E_vis_hat, F_vis_hat, G_vis_hat, mut = central_diff(E_vis_hat, F_vis_hat, G_vis_hat, QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB, 
                                                                QconU, QconD, QconL, QconR, QconF, QconB, cellxmax, cellymax, cellzmax, mu, mut, mut_bd, lambda,
                                                                vecAx, vecAy, vecAz, specific_heat_ratio, volume, Rd, nval, yplus, swith_wall, icell)

                # RHS
                RHS = setup_RHS(RHS, cellxmax, cellymax, cellzmax, E_adv_hat, F_adv_hat, G_adv_hat, 
                                E_vis_hat, F_vis_hat, G_vis_hat, nval, volume, icell)
            
                
                # 粘性の更新
                for k in 1:cellzmax
                    for j in 1:cellymax
                        for i in 1:cellxmax
                            mu[i,j,k] = mu[i,j,k] + mut[i,j,k]
                        end
                    end
                end

                # set inner time step by local time stepping
                dtau   = set_lts(dtau, lambda_facex, lambda_facey, lambda_facez, Qbasem, cellxmax, cellymax, cellzmax, 
                                    mu, dx, dy, dz, vecAx, vecAy, vecAz, volume, specific_heat_ratio, cfl, icell)
                
                # lusgs_advection_term
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, C_adv_hat_p, C_adv_hat_m, 
                A_beta_shig, B_beta_shig, C_beta_shig = one_wave(A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, C_adv_hat_m, C_adv_hat_p,
                                                                A_beta_shig, B_beta_shig, C_beta_shig, I, Qbasem, Qcon, cellxmax, cellymax, cellzmax,
                                                                vecAx, vecAy, vecAz, specific_heat_ratio, volume, nval)
                # lusgs_viscos_term
                jalphaP, jbetaP, jgammaP = central_diff_jacobian(jalphaP, jbetaP, jgammaP, Qbasem, Qcon, cellxmax, cellymax, cellzmax, mu, lambda,
                                                        vecAx, vecAy, vecAz, specific_heat_ratio, volume, nval)

                #println(RHS[2,:,1])
                #println(A_adv_hat_m[2,:,1])
                #println(jalphaP[2,:])
                
                # LUSGS
                ite = 0
                while true
                    # copy delta_Q
                    for l in 1:nval
                        for k in 1:cellzmax
                            for j in 1:cellymax
                                for i in 1:cellxmax
                                    delta_Q_temp[i,j,k,l] = delta_Q[i,j,k,l]
                                end
                            end
                        end
                    end
                                        
                    # Reversing the left-hand side by lusgs
                    delta_Q = lusgs(D, Lx, Ly, Lz, Ux, Uy, Uz, LdQ, UdQ, RHS_temp, I, dt, dtau, Qcon_hat, Qconn_hat, delta_Q,
                    A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, C_adv_hat_p, C_adv_hat_m, A_beta_shig, B_beta_shig, C_beta_shig,
                    jalphaP, jbetaP, jgammaP, RHS, cellxmax, cellymax, cellzmax, volume, nval, icell)
                    
                    # cal Residuals by norm-2
                    res   = set_res(res, delta_Q, delta_Q_temp, cellxmax, cellymax, cellzmax, nval, icell)
                    norm2 = check_converge(res, RHS, cellxmax, cellymax, cellzmax, init_small, nval, icell)
                    
                    # Find out if the results converged
                    if norm2[1] < norm_ok && norm2[2] < norm_ok && norm2[3] < norm_ok && norm2[4] < norm_ok && norm2[5] < norm_ok
                        break
                    end
                    if ite % 100 ==0
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
                    for k in 1+icell:cellzmax-icell
                        for j in 1+icell:cellymax-icell
                            for i in 1+icell:cellxmax-icell
                                Qcon_hat[i,j,k,l] = Qcon_hat[i,j,k,l] + delta_Q[i,j,k,l]
                            end
                        end
                    end
                end
                # calculate primitive variables
                Qcon  = Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, cellzmax, volume, nval)
                Qbase = conservative_to_base(Qbasem, Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio)
                Qbase_ave = cal_Qave(Qbasem, Qbase_ave, cellxmax, cellymax, cellzmax, nval)
                
                # Find out if the results were divergent
                check_divrege(Qbase, cellxmax, cellymax, cellzmax, Rd, fwrite, icell)
            end
            # End of the inner iteration

            # Updating in physical time
            for l in 1:nval
                for k in 1:cellzmax
                    for j in 1:cellymax
                        for i in 1:cellxmax
                            Qbase[i,j,k,l] = Qbasem[i,j,k,l]
                        end
                    end
                end
            end
            
            # output
            if round(stepnum) % every_outnum == 0
                println("\n")
                println("nt_______________________________"*string(round(stepnum)))
                output_result(stepnum, Qbase, cellxmax, cellymax, cellzmax, specific_heat_ratio, 
                                out_file_front, out_ext, out_dir, Rd, nval, icell)
                output_ave(Qbase_ave, cellxmax, cellymax, cellzmax, out_file_front, out_ext, out_dir, Rd, nval, loop_ite, icell)
            end

            # Find out if the results were divergent
            check_divrege(Qbase, cellxmax, cellymax, cellzmax, Rd, fwrite, icell)
        end
    end
    
    # end of time measurement
    end_t = now()

    # output of calculation time
    output_fin(fwrite, start_t, end_t, nt, dt, in_nt, cellxmax, cellymax, cellzmax)
end


# -- main --
main()
#throw(UndefVarError(:x))