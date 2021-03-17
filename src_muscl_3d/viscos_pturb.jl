# ------------------------------------
# calcucation of the viscosity term by central difference
# ------------------------------------
function central_diff(E_vis_hat, F_vis_hat, G_vis_hat, QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB, 
                    QconU, QconD, QconL, QconR, QconF, QconB, cellxmax, cellymax, mu, mut, mut_bd, lambda,
                    vecAx, vecAy, vecAz, specific_heat_ratio, volume, Rd, nval, yplus, swith_wall, icell)
    
    for k in 1+icell:cellzmax -icell
        for j in 1+icell:cellymax -icell
            for i in 1+icell:cellxmax+1 -icell
                rho_av = 0.5*(QbaseR[i,j,k,1] + QbaseL[i,j,k,1]) 
                u_av   = 0.5*(QbaseR[i,j,k,2] + QbaseL[i,j,k,2])
                v_av   = 0.5*(QbaseR[i,j,k,3] + QbaseL[i,j,k,3])
                w_av   = 0.5*(QbaseR[i,j,k,4] + QbaseL[i,j,k,4])
                
                mu_av     = 0.5*(mu[i-1,j,k] + mu[i,j,k]) 
                lambda_av = (0.5*(1/lambda[i-1,j,k] + 1/lambda[i,j,k]))^(-1)

                volume_av = 0.5*(volume[i,j,k] + volume[i-1,j,k])
        
                dudxi = QbaseR[i,j,k,2] - QbaseL[i,j,k,2]
                dvdxi = QbaseR[i,j,k,3] - QbaseL[i,j,k,3]
                dwdxi = QbaseR[i,j,k,4] - QbaseL[i,j,k,4]
                dTdxi = QbaseR[i,j,k,5]/(QbaseR[i,j,k,1]*Rd) - QbaseL[i,j,k,5]/(QbaseL[i,j,k,1]*Rd)
        
                dudeta = 0.25 * (QbaseU[i,  j+1,k,2] - QbaseD[i,  j+1,k,2] + QbaseU[i,  j,k,2] - QbaseD[i,  j,k,2]
                               + QbaseU[i-1,j+1,k,2] - QbaseD[i-1,j+1,k,2] + QbaseU[i-1,j,k,2] - QbaseD[i-1,j,k,2])
                dvdeta = 0.25 * (QbaseU[i,  j+1,k,3] - QbaseD[i,  j+1,k,3] + QbaseU[i,  j,k,3] - QbaseD[i,  j,k,3]
                               + QbaseU[i-1,j+1,k,3] - QbaseD[i-1,j+1,k,3] + QbaseU[i-1,j,k,3] - QbaseD[i-1,j,k,3])
                dwdeta = 0.25 * (QbaseU[i,  j+1,k,4] - QbaseD[i,  j+1,k,4] + QbaseU[i,  j,k,4] - QbaseD[i,  j,k,4]
                               + QbaseU[i-1,j+1,k,4] - QbaseD[i-1,j+1,k,4] + QbaseU[i-1,j,k,4] - QbaseD[i-1,j,k,4])
                
                dT1 = QbaseU[  i,j+1,k,5]/(QbaseU[  i,j+1,k,1]*Rd) - QbaseD[  i,j+1,k,5]/(QbaseD[  i,j+1,k,1]*Rd)
                dT2 = QbaseU[  i,  j,k,5]/(QbaseU[  i,  j,k,1]*Rd) - QbaseD[  i,  j,k,5]/(QbaseD[  i,  j,k,1]*Rd)
                dT3 = QbaseU[i-1,j+1,k,5]/(QbaseU[i-1,j+1,k,1]*Rd) - QbaseD[i-1,j+1,k,5]/(QbaseD[i-1,j+1,k,1]*Rd)
                dT4 = QbaseU[i-1,  j,k,5]/(QbaseU[i-1,  j,k,1]*Rd) - QbaseD[i-1,  j,k,5]/(QbaseD[i-1,  j,k,1]*Rd)
                dTdeta = 0.25 * (dT1 + dT2 + dT3 + dT4)

                dudzeta = 0.25 * (QbaseF[i,  j,k+1,2] - QbaseB[i,  j,k+1,2] + QbaseF[i,  j,k,2] - QbaseB[i,  j,k,2]
                                + QbaseF[i-1,j,k+1,2] - QbaseB[i-1,j,k+1,2] + QbaseF[i-1,j,k,2] - QbaseB[i-1,j,k,2])
                dvdzeta = 0.25 * (QbaseF[i,  j,k+1,3] - QbaseB[i,  j,k+1,3] + QbaseF[i,  j,k,3] - QbaseB[i,  j,k,3]
                                + QbaseF[i-1,j,k+1,3] - QbaseB[i-1,j,k+1,3] + QbaseF[i-1,j,k,3] - QbaseB[i-1,j,k,3])
                dwdzeta = 0.25 * (QbaseF[i,  j,k+1,4] - QbaseB[i,  j,k+1,4] + QbaseF[i,  j,k,4] - QbaseB[i,  j,k,4]
                                + QbaseF[i-1,j,k+1,4] - QbaseB[i-1,j,k+1,4] + QbaseF[i-1,j,k,4] - QbaseB[i-1,j,k,4])
                
                dT1 = QbaseF[  i,  j,k+1,5]/(QbaseF[  i,  j,k+1,1]*Rd) - QbaseB[  i,  j,k+1,5]/(QbaseB[  i,  j,k+1,1]*Rd)
                dT2 = QbaseF[  i,  j,  k,5]/(QbaseF[  i,  j,  k,1]*Rd) - QbaseB[  i,  j,  k,5]/(QbaseB[  i,  j,  k,1]*Rd)
                dT3 = QbaseF[i-1,  j,k+1,5]/(QbaseF[i-1,  j,k+1,1]*Rd) - QbaseB[i-1,  j,k+1,5]/(QbaseB[i-1,  j,k+1,1]*Rd)
                dT4 = QbaseF[i-1,  j,  k,5]/(QbaseF[i-1,  j,  k,1]*Rd) - QbaseB[i-1,  j,  k,5]/(QbaseB[i-1,  j,  k,1]*Rd)
                dTdzeta = 0.25 * (dT1 + dT2 + dT3 + dT4)

                vecAy_xav    = 0.25*( vecAy[i,j,k,1] + vecAy[i,j+1,k,1] + vecAy[i-1,j,k,1] + vecAy[i-1,j+1,k,1] )
                vecAy_yav    = 0.25*( vecAy[i,j,k,2] + vecAy[i,j+1,k,2] + vecAy[i-1,j,k,2] + vecAy[i-1,j+1,k,2] )
                vecAy_zav    = 0.25*( vecAy[i,j,k,3] + vecAy[i,j+1,k,3] + vecAy[i-1,j,k,3] + vecAy[i-1,j+1,k,3] )

                vecAz_xav    = 0.25*( vecAz[i,j,k,1] + vecAz[i,j,k+1,1] + vecAz[i-1,j,k,1] + vecAz[i-1,j,k+1,1] )
                vecAz_yav    = 0.25*( vecAz[i,j,k,2] + vecAz[i,j,k+1,2] + vecAz[i-1,j,k,2] + vecAz[i-1,j,k+1,2] )
                vecAz_zav    = 0.25*( vecAz[i,j,k,3] + vecAz[i,j,k+1,3] + vecAz[i-1,j,k,3] + vecAz[i-1,j,k+1,3] )
                

                dudx = vecAx[i,j,k,1] * dudxi + vecAy_xav * dudeta + vecAz_xav * dudzeta
                dudy = vecAx[i,j,k,2] * dudxi + vecAy_yav * dudeta + vecAz_yav * dudzeta
                dudz = vecAx[i,j,k,3] * dudxi + vecAy_zav * dudeta + vecAz_zav * dudzeta

                dvdx = vecAx[i,j,k,1] * dvdxi + vecAy_xav * dvdeta + vecAz_xav * dvdzeta
                dvdy = vecAx[i,j,k,2] * dvdxi + vecAy_yav * dvdeta + vecAz_yav * dvdzeta
                dvdz = vecAx[i,j,k,3] * dvdxi + vecAy_zav * dvdeta + vecAz_zav * dvdzeta

                dwdx = vecAx[i,j,k,1] * dwdxi + vecAy_xav * dwdeta + vecAz_xav * dwdzeta
                dwdy = vecAx[i,j,k,2] * dwdxi + vecAy_yav * dwdeta + vecAz_yav * dwdzeta
                dwdz = vecAx[i,j,k,3] * dwdxi + vecAy_zav * dwdeta + vecAz_zav * dwdzeta
                                
                sigma_xx = 2/3*mu_av*( 2*dudx -   dvdy -   dwdz )
                sigma_yy = 2/3*mu_av*( - dudx + 2*dvdy -   dwdz )
                sigma_zz = 2/3*mu_av*( - dudx -   dvdy + 2*dwdz )
                        
                sigma_xy = mu_av*( dudy + dvdx )
                sigma_xz = mu_av*( dvdz + dwdy )
                sigma_yz = mu_av*( dwdx + dudz )
                
                dTdx = vecAx[i,j,k,1]*dTdxi + vecAy_xav*dTdeta + vecAz_xav*dTdzeta
                dTdy = vecAx[i,j,k,2]*dTdxi + vecAy_yav*dTdeta + vecAz_yav*dTdzeta
                dTdz = vecAx[i,j,k,3]*dTdxi + vecAy_zav*dTdeta + vecAz_zav*dTdzeta
        
                betax = u_av*sigma_xx + v_av*sigma_xy + w_av*sigma_xz + lambda_av * dTdx
                betay = u_av*sigma_xy + v_av*sigma_yy + w_av*sigma_yz + lambda_av * dTdy
                betaz = u_av*sigma_xz + v_av*sigma_zy + w_av*sigma_zz + lambda_av * dTdz
                
                yplus_av = 0.5 *(yplus[i-1,j,k] + yplus[i,j,k])

                if swith_wall[1] == 1 && i == 1+icell
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, 
                    e_sgs_x, e_sgs_y, e_sgs_z, mut_bd[i,j,1] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                elseif swith_wall[2] == 1 && i == cellxmax - icell
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, 
                    e_sgs_x, e_sgs_y, e_sgs_z, mut_bd[i,j,1] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                else
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, 
                    e_sgs_x, e_sgs_y, e_sgs_z, mut_bd[i,j,1] = 
                    Smagorinsky_model(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dTdx, dTdy, dTdz,
                                    rho_av, u_av, v_av, w_av, volume_av, yplus_av, i, j)
                    
                    #tau_xx, tau_xy, tau_yy, e_sgs_x, e_sgs_y, mut_bd[i,j,1] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

                #=    
                    if i == 51 && j ==3
                        println(" tur tau ")
                        println(tau_xx)
                        println(tau_xy)
                        println(tau_yy)
                        println(e_sgs_x)
                        println(e_sgs_y)
                        println(" kiyo ")
                        println(vecAx[i,j,1]*sigma_xx + vecAx[i,j,2]*sigma_xy)
                        println(vecAx[i,j,1]*tau_xx + vecAx[i,j,2]*tau_xy)
                        println(vecAx[i,j,1]*tau_xx)
                        println(vecAx[i,j,2]*tau_xy)
                    end
                    =#
                end
                #=
                if i == 20 && j == 20
                    println("(20,20,du)")
                    println(dudxi)
                    println(dudeta)
                    println(vecAx[i,j,1])
                    println(vecAy_xav)
                    println(volume_av)
                    println((vecAx[i,j,1] * dudxi + vecAy_xav * dudeta)/volume_av)
                    println(tau_xx)
                    println(tau_xy)
                    println(tau_yy)
                end
                if i == 150 && j == 150
                    println("(150,150,du)")
                    println(dudxi)
                    println(dudeta)
                    println(vecAx[i,j,1])
                    println(vecAy_xav)
                    println(volume_av)
                    println((vecAx[i,j,1] * dudxi + vecAy_xav * dudeta)/volume_av)
                    println(tau_xx)
                    println(tau_xy)
                    println(tau_yy)
                end
                =#
                
                E_vis_hat[i,j,k,1] = 0.0
                E_vis_hat[i,j,k,2] = ((vecAx[i,j,k,1]*sigma_xx + vecAx[i,j,k,2]*sigma_xy + vecAx[i,j,k,3]*sigma_xz) -
                                      (vecAx[i,j,k,1]*tau_xx   + vecAx[i,j,k,2]*tau_xy   + vecAx[i,j,k,3]*tau_xz)) / volume_av
                E_vis_hat[i,j,k,3] = ((vecAx[i,j,k,1]*sigma_xy + vecAx[i,j,k,2]*sigma_yy + vecAx[i,j,k,3]*sigma_yz) - 
                                      (vecAx[i,j,k,1]*tau_xy   + vecAx[i,j,k,2]*tau_yy   + vecAx[i,j,k,3]*tau_yz)) / volume_av
                E_vis_hat[i,j,k,4] = ((vecAx[i,j,k,1]*sigma_xz + vecAx[i,j,k,2]*sigma_yz + vecAx[i,j,k,3]*sigma_zz) -
                                      (vecAx[i,j,k,1]*tau_xy   + vecAx[i,j,k,2]*tau_yy   + vecAx[i,j,k,3]*tau_zz)) / volume_av
                E_vis_hat[i,j,k,5] = ((vecAx[i,j,k,1]*betax    + vecAx[i,j,k,2]*betay    + vecAx[i,j,k,3]*betaz) -
                                      (vecAx[i,j,k,1]*e_sgs_x  + vecAx[i,j,k,2]*e_sgs_y  + vecAx[i,j,k,3]*e_sgs_z))/ volume_av
            end
        end
    end

    for k in 1+icell:cellzmax -icell
        for j in 1+icell:cellymax+1 -icell
            for i in 1+icell:cellxmax -icell
                rho_av = 0.5*(QbaseU[i,j,k,1] + QbaseD[i,j,k,1]) 
                u_av   = 0.5*(QbaseU[i,j,k,2] + QbaseD[i,j,k,2])
                v_av   = 0.5*(QbaseU[i,j,k,3] + QbaseD[i,j,k,3])
                w_av   = 0.5*(QbaseU[i,j,k,4] + QbaseD[i,j,k,4])

                mu_av     = 0.5*(mu[i,j-1,k] + mu[i,j,k])
                lambda_av = (0.5*(1/lambda[i,j-1,k] + 1/lambda[i,j,k]))^(-1)

                volume_av = 0.5*(volume[i,j,k] + volume[i,j-1,k])
        
                dudxi = 0.25 * (QbaseR[i+1,  j,k,2] - QbaseL[i+1,  j,k,2] + QbaseR[i,  j,k,2] - QbaseL[i,  j,k,2]
                              + QbaseR[i+1,j-1,k,2] - QbaseL[i+1,j-1,k,2] + QbaseR[i,j-1,k,2] - QbaseL[i,j-1,k,2])
                dvdxi = 0.25 * (QbaseR[i+1,  j,k,3] - QbaseL[i+1,  j,k,3] + QbaseR[i,  j,k,3] - QbaseL[i,  j,k,3]
                              + QbaseR[i+1,j-1,k,3] - QbaseL[i+1,j-1,k,3] + QbaseR[i,j-1,k,3] - QbaseL[i,j-1,k,3])
                dwdxi = 0.25 * (QbaseR[i+1,  j,k,4] - QbaseL[i+1,  j,k,4] + QbaseR[i,  j,k,4] - QbaseL[i,  j,k,4]
                              + QbaseR[i+1,j-1,k,4] - QbaseL[i+1,j-1,k,4] + QbaseR[i,j-1,k,4] - QbaseL[i,j-1,k,4])
                
                dT1 = QbaseR[i+1,  j,k,5]/(QbaseR[i+1,  j,k,1]*Rd) - QbaseL[i+1,  j,k,5]/(QbaseL[i+1,  j,k,1]*Rd)
                dT2 = QbaseR[  i,  j,k,5]/(QbaseR[  i,  j,k,1]*Rd) - QbaseL[  i,  j,k,5]/(QbaseL[  i,  j,k,1]*Rd)
                dT3 = QbaseR[i+1,j-1,k,5]/(QbaseR[i+1,j-1,k,1]*Rd) - QbaseL[i+1,j-1,k,5]/(QbaseL[i+1,j-1,k,1]*Rd)
                dT4 = QbaseR[  i,j-1,k,5]/(QbaseR[  i,j-1,k,1]*Rd) - QbaseL[  i,j-1,k,5]/(QbaseL[  i,j-1,k,1]*Rd)
                dTdxi = 0.25 * (dT1 + dT2 + dT3 + dT4)

                dudeta = QbaseU[i,j,k,2] - QbaseD[i,j,k,2]
                dvdeta = QbaseU[i,j,k,3] - QbaseD[i,j,k,3]
                dwdeta = QbaseU[i,j,k,4] - QbaseD[i,j,k,4]
                dTdeta = (QbaseU[i,j,k,5]/(QbaseU[i,j,k,1]*Rd) - QbaseD[i,j,k,5]/(QbaseD[i,j,k,1]*Rd))

                
                dudzeta = 0.25 * (QbaseF[i,  j,k+1,2] - QbaseB[i,  j,k+1,2] + QbaseF[i,  j,k,2] - QbaseB[i,  j,k,2]
                                + QbaseF[i,j-1,k+1,2] - QbaseB[i,j-1,k+1,2] + QbaseF[i,j-1,k,2] - QbaseB[i,j-1,k,2])
                dvdzeta = 0.25 * (QbaseF[i,  j,k+1,3] - QbaseB[i,  j,k+1,3] + QbaseF[i,  j,k,3] - QbaseB[i,  j,k,3]
                                + QbaseF[i,j-1,k+1,3] - QbaseB[i,j-1,k+1,3] + QbaseF[i,j-1,k,3] - QbaseB[i,j-1,k,3])
                dwdzeta = 0.25 * (QbaseF[i,  j,k+1,4] - QbaseB[i,  j,k+1,4] + QbaseF[i,  j,k,4] - QbaseB[i,  j,k,4]
                                + QbaseF[i,j-1,k+1,4] - QbaseB[i,j-1,k+1,4] + QbaseF[i,j-1,k,4] - QbaseB[i,j-1,k,4])
                
                dT1 = QbaseF[  i,  j,k+1,5]/(QbaseF[  i,  j,k+1,1]*Rd) - QbaseB[  i,  j,k+1,5]/(QbaseB[  i,  j,k+1,1]*Rd)
                dT2 = QbaseF[  i,  j,  k,5]/(QbaseF[  i,  j,  k,1]*Rd) - QbaseB[  i,  j,  k,5]/(QbaseB[  i,  j,  k,1]*Rd)
                dT3 = QbaseF[  i,j-1,k+1,5]/(QbaseF[  i,j-1,k+1,1]*Rd) - QbaseB[  i,j-1,k+1,5]/(QbaseB[  i,j-1,k+1,1]*Rd)
                dT4 = QbaseF[  i,j-1,  k,5]/(QbaseF[  i,j-1,  k,1]*Rd) - QbaseB[  i,j-1,  k,5]/(QbaseB[  i,j-1,  k,1]*Rd)
                dTdzeta = 0.25 * (dT1 + dT2 + dT3 + dT4)
            
                vecAx_xav    = 0.25*( vecAx[i,j,k,1] + vecAx[i+1,j,k,1] + vecAx[i,j-1,k,1] + vecAx[i+1,j-1,k,1] )
                vecAx_yav    = 0.25*( vecAx[i,j,k,2] + vecAx[i+1,j,k,2] + vecAx[i,j-1,k,2] + vecAx[i+1,j-1,k,2] )
                vecAx_zav    = 0.25*( vecAx[i,j,k,3] + vecAx[i+1,j,k,3] + vecAx[i,j-1,k,3] + vecAx[i+1,j-1,k,3] )

                vecAz_xav    = 0.25*( vecAz[i,j,k,1] + vecAz[i,j,k+1,1] + vecAz[i,j-1,k,1] + vecAz[i,j-1,k+1,1] )
                vecAz_yav    = 0.25*( vecAz[i,j,k,2] + vecAz[i,j,k+1,2] + vecAz[i,j-1,k,2] + vecAz[i,j-1,k+1,2] )
                vecAz_zav    = 0.25*( vecAz[i,j,k,3] + vecAz[i,j,k+1,3] + vecAz[i,j-1,k,3] + vecAz[i,j-1,k+1,3] )



                dudx = vecAx_xav * dudxi + vecAy[i,j,1] * dudeta
                dvdy = vecAx_yav * dvdxi + vecAy[i,j,2] * dvdeta
                dudy = vecAx_yav * dudxi + vecAy[i,j,2] * dudeta
                dvdx = vecAx_xav * dvdxi + vecAy[i,j,1] * dvdeta
            
                dudx = vecAx_xav * dudxi + vecAy[i,j,1] * dudeta + vecAz_xav * dudzeta
                dudy = vecAx_yav * dudxi + vecAy[i,j,2] * dudeta + vecAz_yav * dudzeta
                dudz = vecAx_zav * dudxi + vecAy[i,j,3] * dudeta + vecAz_zav * dudzeta

                dvdx = vecAx_xav * dvdxi + vecAy[i,j,1] * dvdeta + vecAz_xav * dvdzeta
                dvdy = vecAx_yav * dvdxi + vecAy[i,j,2] * dvdeta + vecAz_yav * dvdzeta
                dvdz = vecAx_zav * dvdxi + vecAy[i,j,3] * dvdeta + vecAz_zav * dvdzeta

                dwdx = vecAx_xav * dwdxi + vecAy[i,j,1] * dwdeta + vecAz_xav * dwdzeta
                dwdy = vecAx_yav * dwdxi + vecAy[i,j,2] * dwdeta + vecAz_yav * dwdzeta
                dwdz = vecAx_zav * dwdxi + vecAy[i,j,3] * dwdeta + vecAz_zav * dwdzeta
                                
                sigma_xx = 2/3*mu_av*( 2*dudx -   dvdy -   dwdz )
                sigma_yy = 2/3*mu_av*( - dudx + 2*dvdy -   dwdz )
                sigma_zz = 2/3*mu_av*( - dudx -   dvdy + 2*dwdz )
                        
                sigma_xy = mu_av*( dudy + dvdx )
                sigma_xz = mu_av*( dvdz + dwdy )
                sigma_yz = mu_av*( dwdx + dudz )
                
                dTdx = vecAx_xav*dTdxi + vecAy[i,j,1]*dTdeta + vecAz_xav*dTdzeta
                dTdy = vecAx_yav*dTdxi + vecAy[i,j,2]*dTdeta + vecAz_yav*dTdzeta
                dTdz = vecAx_zav*dTdxi + vecAy[i,j,3]*dTdeta + vecAz_zav*dTdzeta
        
                betax = u_av*sigma_xx + v_av*sigma_xy + w_av*sigma_xz + lambda_av * dTdx
                betay = u_av*sigma_xy + v_av*sigma_yy + w_av*sigma_yz + lambda_av * dTdy
                betaz = u_av*sigma_xz + v_av*sigma_zy + w_av*sigma_zz + lambda_av * dTdz
                        
                yplus_av = 0.5 *(yplus[i,j-1,k] + yplus[i,j,k])

                if swith_wall[3] == 1 && j == 1+icell
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, 
                    e_sgs_x, e_sgs_y, e_sgs_z, mut_bd[i,j,2] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                elseif swith_wall[3] == 1 && j == cellymax - icell
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, 
                    e_sgs_x, e_sgs_y, e_sgs_z, mut_bd[i,j,2] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                else
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, 
                    e_sgs_x, e_sgs_y, e_sgs_z, mut_bd[i,j,2] = 
                    Smagorinsky_model(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dTdx, dTdy, dTdz,
                                    rho_av, u_av, v_av, w_av, volume_av, yplus_av, i, j)
                end

                F_vis_hat[i,j,k,1] = 0.0
                F_vis_hat[i,j,k,2] = ((vecAy[i,j,k,1]*sigma_xx + vecAy[i,j,k,2]*sigma_xy + vecAy[i,j,k,3]*sigma_xz) -
                                      (vecAy[i,j,k,1]*tau_xx   + vecAy[i,j,k,2]*tau_xy   + vecAy[i,j,k,3]*tau_xz)) / volume_av
                F_vis_hat[i,j,k,3] = ((vecAy[i,j,k,1]*sigma_xy + vecAy[i,j,k,2]*sigma_yy + vecAy[i,j,k,3]*sigma_yz) - 
                                      (vecAy[i,j,k,1]*tau_xy   + vecAy[i,j,k,2]*tau_yy   + vecAy[i,j,k,3]*tau_yz)) / volume_av
                F_vis_hat[i,j,k,4] = ((vecAy[i,j,k,1]*sigma_xz + vecAy[i,j,k,2]*sigma_yz + vecAy[i,j,k,3]*sigma_zz) -
                                      (vecAy[i,j,k,1]*tau_xy   + vecAy[i,j,k,2]*tau_yy   + vecAy[i,j,k,3]*tau_zz)) / volume_av
                F_vis_hat[i,j,k,5] = ((vecAy[i,j,k,1]*betax    + vecAy[i,j,k,2]*betay    + vecAy[i,j,k,3]*betaz) -
                                      (vecAy[i,j,k,1]*e_sgs_x  + vecAy[i,j,k,2]*e_sgs_y  + vecAy[i,j,k,3]*e_sgs_z))/ volume_av
            end
        end
    end
    
    for k in 1+icell:cellzmax+1 -icell
        for j in 1+icell:cellymax -icell
            for i in 1+icell:cellxmax -icell
                rho_av = 0.5*(QbaseF[i,j,k,1] + QbaseB[i,j,k,1]) 
                u_av   = 0.5*(QbaseF[i,j,k,2] + QbaseB[i,j,k,2])
                v_av   = 0.5*(QbaseF[i,j,k,3] + QbaseB[i,j,k,3])
                w_av   = 0.5*(QbaseF[i,j,k,4] + QbaseB[i,j,k,4])
                
                mu_av     = 0.5*(mu[i,j,k-1] + mu[i,j,k]) 
                lambda_av = (0.5*(1/lambda[i,j,k-1] + 1/lambda[i,j,k]))^(-1)

                volume_av = 0.5*(volume[i,j,k] + volume[i,j,k-1])
        
                dudxi = 0.25 * (QbaseR[i+1,  j,  k,2] - QbaseL[i+1,  j,  k,2] + QbaseR[i,  j,  k,2] - QbaseL[i,  j,  k,2]
                              + QbaseR[i+1,  j,k-1,2] - QbaseL[i+1,  j,k-1,2] + QbaseR[i,  j,k-1,2] - QbaseL[i,  j,k-1,2])
                dvdxi = 0.25 * (QbaseR[i+1,  j,  k,3] - QbaseL[i+1,  j,  k,3] + QbaseR[i,  j,  k,3] - QbaseL[i,  j,  k,3]
                              + QbaseR[i+1,  j,k-1,3] - QbaseL[i+1,  j,k-1,3] + QbaseR[i,  j,k-1,3] - QbaseL[i,  j,k-1,3])
                dwdxi = 0.25 * (QbaseR[i+1,  j,  k,4] - QbaseL[i+1,  j,  k,4] + QbaseR[i,  j,  k,4] - QbaseL[i,  j,  k,4]
                              + QbaseR[i+1,  j,k-1,4] - QbaseL[i+1,  j,k-1,4] + QbaseR[i,  j,k-1,4] - QbaseL[i,  j,k-1,4])
                
                dT1 = QbaseR[i+1,  j,  k,5]/(QbaseR[i+1,  j,  k,1]*Rd) - QbaseL[i+1,  j,  k,5]/(QbaseL[i+1,  j,  k,1]*Rd)
                dT2 = QbaseR[  i,  j,  k,5]/(QbaseR[  i,  j,  k,1]*Rd) - QbaseL[  i,  j,  k,5]/(QbaseL[  i,  j,  k,1]*Rd)
                dT3 = QbaseR[i+1,  j,k-1,5]/(QbaseR[i+1,  j,k-1,1]*Rd) - QbaseL[i+1,  j,k-1,5]/(QbaseL[i+1,  j,k-1,1]*Rd)
                dT4 = QbaseR[  i,  j,k-1,5]/(QbaseR[  i,  j,k-1,1]*Rd) - QbaseL[  i,  j,k-1,5]/(QbaseL[  i,  j,k-1,1]*Rd)
                dTdxi = 0.25 * (dT1 + dT2 + dT3 + dT4)
        
                dudeta = 0.25 * (QbaseU[  i,j+1,  k,2] - QbaseD[  i,j+1,  k,2] + QbaseU[i,  j,  k,2] - QbaseD[i,  j,  k,2]
                               + QbaseU[  i,j+1,k-1,2] - QbaseD[  i,j+1,k-1,2] + QbaseU[i,  j,k-1,2] - QbaseD[i,  j,k-1,2])
                dvdeta = 0.25 * (QbaseU[  i,j+1,  k,3] - QbaseD[  i,j+1,  k,3] + QbaseU[i,  j,  k,3] - QbaseD[i,  j,  k,3]
                               + QbaseU[  i,j+1,k-1,3] - QbaseD[  i,j+1,k-1,3] + QbaseU[i,  j,k-1,3] - QbaseD[i,  j,k-1,3])
                dwdeta = 0.25 * (QbaseU[  i,j+1,  k,4] - QbaseD[  i,j+1,  k,4] + QbaseU[i,  j,  k,4] - QbaseD[i,  j,  k,4]
                               + QbaseU[  i,j+1,k-1,4] - QbaseD[  i,j+1,k-1,4] + QbaseU[i,  j,k-1,4] - QbaseD[i,  j,k-1,4])
                
                dT1 = QbaseU[  i,j+1,  k,5]/(QbaseU[  i,j+1,  k,1]*Rd) - QbaseD[  i,j+1,  k,5]/(QbaseD[  i,j+1,  k,1]*Rd)
                dT2 = QbaseU[  i,  j,  k,5]/(QbaseU[  i,  j,  k,1]*Rd) - QbaseD[  i,  j,  k,5]/(QbaseD[  i,  j,  k,1]*Rd)
                dT3 = QbaseU[  i,j+1,k-1,5]/(QbaseU[  i,j+1,k-1,1]*Rd) - QbaseD[  i,j+1,k-1,5]/(QbaseD[  i,j+1,k-1,1]*Rd)
                dT4 = QbaseU[  i,  j,k-1,5]/(QbaseU[  i,  j,k-1,1]*Rd) - QbaseD[  i,  j,k-1,5]/(QbaseD[  i,  j,k-1,1]*Rd)
                dTdeta = 0.25 * (dT1 + dT2 + dT3 + dT4)

                dudzeta = QbaseF[i,j,k,2] - QbaseB[i,j,k,2]
                dvdzeta = QbaseF[i,j,k,3] - QbaseB[i,j,k,3]
                dwdzeta = QbaseF[i,j,k,4] - QbaseB[i,j,k,4]
                dTdzeta = (QbaseF[i,j,k,5]/(QbaseF[i,j,k,1]*Rd) - QbaseB[i,j,k,5]/(QbaseB[i,j,k,1]*Rd))
                
                vecAx_xav    = 0.25*( vecAx[i,j,k,1] + vecAx[i+1,j,k,1] + vecAx[i,j,k-1,1] + vecAx[i+1,j,k-1,1] )
                vecAx_yav    = 0.25*( vecAx[i,j,k,2] + vecAx[i+1,j,k,2] + vecAx[i,j,k-1,2] + vecAx[i+1,j,k-1,2] )
                vecAx_zav    = 0.25*( vecAx[i,j,k,3] + vecAx[i+1,j,k,3] + vecAx[i,j,k-1,3] + vecAx[i+1,j,k-1,3] )

                vecAy_xav    = 0.25*( vecAy[i,j,k,1] + vecAy[i,j+1,k,1] + vecAy[i,j,k-1,1] + vecAy[i,j+1,k-1,1] )
                vecAy_yav    = 0.25*( vecAy[i,j,k,2] + vecAy[i,j+1,k,2] + vecAy[i,j,k-1,2] + vecAy[i,j+1,k-1,2] )
                vecAy_zav    = 0.25*( vecAy[i,j,k,3] + vecAy[i,j+1,k,3] + vecAy[i,j,k-1,3] + vecAy[i,j+1,k-1,3] )


                dudx = vecAx_xav * dudxi + vecAy_xav * dudeta + vecAz[i,j,k,1] * dudzeta
                dudy = vecAx_yav * dudxi + vecAy_yav * dudeta + vecAz[i,j,k,2] * dudzeta
                dudz = vecAx_zav * dudxi + vecAy_zav * dudeta + vecAz[i,j,k,3] * dudzeta

                dvdx = vecAx_xav * dvdxi + vecAy_xav * dvdeta + vecAz[i,j,k,1] * dvdzeta
                dvdy = vecAx_yav * dvdxi + vecAy_yav * dvdeta + vecAz[i,j,k,2] * dvdzeta
                dvdz = vecAx_zav * dvdxi + vecAy_zav * dvdeta + vecAz[i,j,k,3] * dvdzeta

                dwdx = vecAx_xav * dwdxi + vecAy_xav * dwdeta + vecAz[i,j,k,1] * dwdzeta
                dwdy = vecAx_yav * dwdxi + vecAy_yav * dwdeta + vecAz[i,j,k,2] * dwdzeta
                dwdz = vecAx_zav * dwdxi + vecAy_zav * dwdeta + vecAz[i,j,k,3] * dwdzeta
                                
                sigma_xx = 2/3*mu_av*( 2*dudx -   dvdy -   dwdz )
                sigma_yy = 2/3*mu_av*( - dudx + 2*dvdy -   dwdz )
                sigma_zz = 2/3*mu_av*( - dudx -   dvdy + 2*dwdz )
                        
                sigma_xy = mu_av*( dudy + dvdx )
                sigma_xz = mu_av*( dvdz + dwdy )
                sigma_yz = mu_av*( dwdx + dudz )
                
                dTdx = vecAx_xav*dTdxi + vecAy_xav*dTdeta + vecAz[i,j,k,1]*dTdzeta
                dTdy = vecAx_yav*dTdxi + vecAy_yav*dTdeta + vecAz[i,j,k,2]*dTdzeta
                dTdz = vecAx_zav*dTdxi + vecAy_zav*dTdeta + vecAz[i,j,k,3]*dTdzeta
        
                betax = u_av*sigma_xx + v_av*sigma_xy + w_av*sigma_xz + lambda_av * dTdx
                betay = u_av*sigma_xy + v_av*sigma_yy + w_av*sigma_yz + lambda_av * dTdy
                betaz = u_av*sigma_xz + v_av*sigma_zy + w_av*sigma_zz + lambda_av * dTdz
                
                yplus_av = 0.5 *(yplus[i,j,k] + yplus[i,j,k-1])

                if swith_wall[5] == 1 && i == 1+icell
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, 
                    e_sgs_x, e_sgs_y, e_sgs_z, mut_bd[i,j,3] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                elseif swith_wall[6] == 1 && i == cellzmax - icell
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, 
                    e_sgs_x, e_sgs_y, e_sgs_z, mut_bd[i,j,3] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                else
                    tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, 
                    e_sgs_x, e_sgs_y, e_sgs_z, mut_bd[i,j,3] = 
                    Smagorinsky_model(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dTdx, dTdy, dTdz,
                                    rho_av, u_av, v_av, w_av, volume_av, yplus_av, i, j)
                end
                
                G_vis_hat[i,j,k,1] = 0.0
                G_vis_hat[i,j,k,2] = ((vecAz[i,j,k,1]*sigma_xx + vecAz[i,j,k,2]*sigma_xy + vecAz[i,j,k,3]*sigma_xz) -
                                      (vecAz[i,j,k,1]*tau_xx   + vecAz[i,j,k,2]*tau_xy   + vecAz[i,j,k,3]*tau_xz)) / volume_av
                G_vis_hat[i,j,k,3] = ((vecAz[i,j,k,1]*sigma_xy + vecAz[i,j,k,2]*sigma_yy + vecAz[i,j,k,3]*sigma_yz) - 
                                      (vecAz[i,j,k,1]*tau_xy   + vecAz[i,j,k,2]*tau_yy   + vecAz[i,j,k,3]*tau_yz)) / volume_av
                G_vis_hat[i,j,k,4] = ((vecAz[i,j,k,1]*sigma_xz + vecAz[i,j,k,2]*sigma_yz + vecAz[i,j,k,3]*sigma_zz) -
                                      (vecAz[i,j,k,1]*tau_xy   + vecAz[i,j,k,2]*tau_yy   + vecAz[i,j,k,3]*tau_zz)) / volume_av
                G_vis_hat[i,j,k,5] = ((vecAz[i,j,k,1]*betax    + vecAz[i,j,k,2]*betay    + vecAz[i,j,k,3]*betaz) -
                                      (vecAz[i,j,k,1]*e_sgs_x  + vecAz[i,j,k,2]*e_sgs_y  + vecAz[i,j,k,3]*e_sgs_z))/ volume_av
            end
        end
    end

    # 乱流粘性の代入
    for k in 1:cellzmax
        for j in 1:cellymax
            for i in 1:cellxmax
                mut[i,j,k] = (mut_bd[i+1,j,k,1] + mut_bd[i,j,k,1] + mut_bd[i,j+1,k,2] 
                            + mut_bd[i,j,k,2] + mut_bd[i,j,k+1,3] + mut_bd[i,j,k,3]) / 6.0
            end
        end
    end

    #=
    println("(u,v ,20)")
    println(Qbase[20,20,2])
    println(Qbase[19,20,2])
    println(Qbase[20,19,2])
    println(Qbase[19,19,2])
    println("(u,v,150)")
    println(Qbase[150,150,2])
    println(Qbase[149,150,2])
    println(Qbase[150,149,2])
    println(Qbase[149,149,2])
    =#
    
    return E_vis_hat, F_vis_hat, G_vis_hat,mut
end

function Smagorinsky_model(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dTdx, dTdy, dTdz,
                            rho_av, u_av, v_av, w_av, volume_av, yplus_av, i, j)
    Cs = 0.18
    Cl = 0.09
    Pr_sgs = 0.6
    k = 0.41

    S11 = 0.5*(dudx + dudx)
    S12 = 0.5*(dudy + dvdx)
    S13 = 0.5*(dudz + dwdx)

    S21 = 0.5*(dvdx + dudy)
    S22 = 0.5*(dvdy + dvdy)
    S23 = 0.5*(dvdz + dwdy)
    
    S31 = 0.5*(dwdx + dudz)
    S32 = 0.5*(dwdy + dvdz)
    S33 = 0.5*(dwdz + dwdz)
    

    absS = ( 2 * (S11^2 + S12^2 + S13^2 + S21^2 + S22^2 + S23^2 + S31^2 + S32^2 + S33^2) )^0.5
    Delta = (volume_av) ^ (1/3) * wallf_Van_Driest(yplus_av)

    #lsgs = 0.5*((Cs*Delta - k*y) - abs(Cs*Delta - k*y)) + k*y # min([ky, Cs*Delta])
    #nu_sgs = (lsgs)^2 * absS
    nu_sgs = (Cs*Delta)^2 * absS
    mu_sgs = nu_sgs * rho_av
    
    #=
    if i == 20 && j == 20
        println("(20,20)")
        println(S11)
        println(S12)
        println(S22)
        println(Delta)
        println(absS)
        println(nu_sgs/volume_av)
    end
    if i == 150 && j == 150
        println("(150,150)")
        println(S11)
        println(S12)
        println(S22)
        println(Delta)
        println(absS)
        println(nu_sgs/volume_av)
    end
    =#



    # SGS stress
    #tau_xx = -2 * rho_av * nu_sgs * 2/3 * dudx +  1/3 * Cl * 2 * rho_av * Delta^2 * absS^2
    #tau_xy = - rho_av * nu_sgs * (dvdx + dudy)
    #tau_yy = -2 * rho_av * nu_sgs * 2/3 * dvdy +  1/3 * Cl * 2 * rho_av * Delta^2 * absS^2
    tau_xx = -2 * rho_av * nu_sgs * S11
    tau_xy = -2 * rho_av * nu_sgs * S12
    tau_xz = -2 * rho_av * nu_sgs * S13

    tau_yy = -2 * rho_av * nu_sgs * S22
    tau_yz = -2 * rho_av * nu_sgs * S23

    tau_zz = -2 * rho_av * nu_sgs * S33

    # SGS energy
    e_sgs_x = -rho_av * nu_sgs / Pr_sgs * dTdx + tau_xx*u_av + tau_xy*v_av + tau_xz*w_av
    e_sgs_y = -rho_av * nu_sgs / Pr_sgs * dTdy + tau_xy*u_av + tau_yy*v_av + tau_yz*w_av
    e_sgs_z = -rho_av * nu_sgs / Pr_sgs * dTdz + tau_xz*u_av + tau_yz*v_av + tau_zz*w_av

    # 
    tau_xx = tau_xx / volume_av   *1.0e0
    tau_xy = tau_xy / volume_av   *1.0e0
    tau_xz = tau_xz / volume_av   *1.0e0
    tau_yy = tau_yy / volume_av   *1.0e0
    tau_yz = tau_yz / volume_av   *1.0e0
    tau_zz = tau_zz / volume_av   *1.0e0
    e_sgs_x = e_sgs_x / volume_av *1.0e0
    e_sgs_y = e_sgs_y / volume_av *1.0e0
    e_sgs_z = e_sgs_z / volume_av *1.0e0
    
    #=
    if i == 51 && j ==3
        println(" iroiro ")
        println(absS)
        println(nu_sgs)
        println(-2 * rho_av * nu_sgs * 2/3 * dudx)
        println(1/3 * Cl * 2 * rho_av * Delta^2 * absS^2)
    end
    =#

    return tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz, e_sgs_x, e_sgs_y, e_sgs_z, mu_sgs
end