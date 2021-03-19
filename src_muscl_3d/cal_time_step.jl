# ------------------------------------
# time integration with explicit methods
# ------------------------------------
function time_integration_explicit(dt, Qcon_hat, RHS, cellxmax, cellymax, cellzmax, nval, icell)
    Threads.@threads for l in 1:nval
        for k in 1+icell:cellzmax-icell
            for j in 1+icell:cellymax-icell
                for i in 1+icell:cellxmax-icell
                    Qcon_hat[i,j,k,l] = Qcon_hat[i,j,k,l] + dt*RHS[i,j,k,l]
                end
            end
        end
    end
    return Qcon_hat
end 

# ------------------------------------
# Calculation of the advection Jacobian matrix by one-wave approximation for lusgs
# ------------------------------------
function one_wave(A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, C_adv_hat_m, C_adv_hat_p, A_beta_shig, B_beta_shig, C_beta_shig, I,
                    Qbase, Qcon, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, specific_heat_ratio, volume, nval)
    beta = 1.1
    kx_av = 0.0
    ky_av = 0.0
    kz_av = 0.0
    
    for k in 2:cellzmax
        for j in 2:cellymax
            for i in 2:cellxmax
                for nn in 1:2   #A, B, C
                    jacob_temp = zeros(nval, nval)
                    if nn == 1
                        kx_av = 0.5*(vecAx[i,j,k,1]+vecAx[i+1,j,k,1]) / volume[i,j,k]
                        ky_av = 0.5*(vecAx[i,j,k,2]+vecAx[i+1,j,k,2]) / volume[i,j,k]
                        kz_av = 0.5*(vecAx[i,j,k,3]+vecAx[i+1,j,k,3]) / volume[i,j,k]
                    elseif nn == 2
                        kx_av = 0.5*(vecAy[i,j,k,1]+vecAy[i,j+1,k,1]) / volume[i,j,k]
                        ky_av = 0.5*(vecAy[i,j,k,2]+vecAy[i,j+1,k,2]) / volume[i,j,k]
                        kz_av = 0.5*(vecAy[i,j,k,3]+vecAy[i,j+1,k,3]) / volume[i,j,k]
                    elseif nn == 3
                        kx_av = 0.5*(vecAz[i,j,k,1]+vecAz[i,j,k+1,1]) / volume[i,j,k]
                        ky_av = 0.5*(vecAz[i,j,k,2]+vecAz[i,j,k+1,2]) / volume[i,j,k]
                        kz_av = 0.5*(vecAz[i,j,k,3]+vecAz[i,j,k+1,3]) / volume[i,j,k]
                    end

                    rho = Qbase[i,j,k,1]
                    u = Qbase[i,j,k,2]
                    v = Qbase[i,j,k,3]
                    w = Qbase[i,j,k,4]
                    e = Qcon[i,j,k,5]
                    p = Qbase[i,j,k,5]

                    g = specific_heat_ratio
                    
                    Z = kx_av*u + ky_av*v + kz_av*w
                    q2   = 0.5*(u^2 + v^2 + w^2)
                    b1c2 = 0.5*q2*(g-1)
                    gebyrho = g*e/rho

                    # advection Jacobian matrix in the general coordinate system
                    jacob_temp[1,1] = 0.0
                    jacob_temp[1,2] = kx_av
                    jacob_temp[1,3] = ky_av
                    jacob_temp[1,4] = kz_av
                    jacob_temp[1,5] = 0.0
                
                    jacob_temp[2,1] = -u*Z + kx_av*b1c2
                    jacob_temp[2,2] = Z - (g-2)*kx_av*u
                    jacob_temp[2,3] = ky_av*u - kx_av*(g-1)*v
                    jacob_temp[2,4] = kz_av*u - kx_av*(g-1)*w
                    jacob_temp[2,5] = (g-1)*kx_av
                
                    jacob_temp[3,1] = -v*Z + ky_av*b1c2
                    jacob_temp[3,2] = kx_av*v - ky_av*(g-1)*u
                    jacob_temp[3,3] = Z - ky_av*(g-2)*v
                    jacob_temp[3,4] = kz_av*v - ky_av*(g-1)*w
                    jacob_temp[3,5] = (g-1)*ky_av

                    jacob_temp[4,1] = -w*Z + kz_av*b1c2
                    jacob_temp[4,2] = kx_av*w - kz_av*(g-1)*u
                    jacob_temp[4,3] = ky_av*w - kz_av*(g-1)*v
                    jacob_temp[4,4] = Z - kz_av*(g-2)*w
                    jacob_temp[4,5] = (g-1)*ky_av
                
                    jacob_temp[5,1] = Z*(-gebyrho + 2*b1c2)
                    jacob_temp[5,2] = kx_av*(gebyrho-b1c2) - (g-1)*u*Z
                    jacob_temp[5,3] = ky_av*(gebyrho-b1c2) - (g-1)*v*Z
                    jacob_temp[5,4] = kz_av*(gebyrho-b1c2) - (g-1)*w*Z
                    jacob_temp[5,5] = g*Z

                    c = (g*rho/p)^0.5
                    shigma = abs(Z) + c*(kx_av^2 + ky_av^2 + kz_av^2)^0.5

                    for l in 1:nval
                        I[l,l] = beta * shigma
                    end

                    if nn == 1
                        for m in 1:nval
                            for l in 1:nval                        
                                A_adv_hat_p[i,j,k,l,m] = 0.5*(jacob_temp[l,m] + I[l,m])
                                A_adv_hat_m[i,j,k,l,m] = 0.5*(jacob_temp[l,m] - I[l,m])
                            end
                        end
                        A_beta_shig[i,j,k] = beta * shigma
                    elseif nn == 2 
                        for m in 1:nval
                            for l in 1:nval
                                B_adv_hat_p[i,j,k,l,m] = 0.5*(jacob_temp[l,m] + I[l,m])
                                B_adv_hat_m[i,j,k,l,m] = 0.5*(jacob_temp[l,m] - I[l,m])
                            end
                        end
                        B_beta_shig[i,j,k] = beta * shigma
                    elseif nn == 3
                        for m in 1:nval
                            for l in 1:nval
                                C_adv_hat_p[i,j,k,l,m] = 0.5*(jacob_temp[l,m] + I[l,m])
                                C_adv_hat_m[i,j,k,l,m] = 0.5*(jacob_temp[l,m] - I[l,m])
                            end
                        end
                        C_beta_shig[i,j,k] = beta * shigma
                    end

                end
            end
        end
    end
    
    # reset
    for l in 1:nval
        I[l,l] = 1.0
    end

    return A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, C_adv_hat_p, C_adv_hat_m, A_beta_shig, B_beta_shig, C_beta_shig
end

# ------------------------------------
# Calculation of the viscosity Jacobian matrix by center differential for lusgs
# ------------------------------------
function central_diff_jacobian(jalphaP, jbetaP, jgammaP, Qbase, Qcon, cellxmax, cellymax, cellzmax, mu, lambda,
                            vecAx, vecAy, vecAz, specific_heat_ratio, volume, nval)
    
    for k in 1:cellzmax
        for j in 1:cellymax
            for i in 1:cellxmax
                rho = Qbase[i,j,k,1]
                
                # alpha
                xi_x_av = 0.5*(vecAx[i,j,k,1] + vecAx[i+1,j,k,1]) / volume[i,j,k]
                xi_y_av = 0.5*(vecAx[i,j,k,2] + vecAx[i+1,j,k,2]) / volume[i,j,k]
                xi_z_av = 0.5*(vecAx[i,j,k,3] + vecAx[i+1,j,k,3]) / volume[i,j,k]
                alpha = (xi_x_av^2 + xi_y_av^2 + xi_z_av^2)^0.5
                jalphaP[i,j,k] = alpha / volume[i,j,k] * (2*mu[i,j,k]/rho)
                
                # beta
                eta_x_av = 0.5*(vecAy[i,j,k,1] + vecAy[i,j+1,k,1]) / volume[i,j,k]
                eta_y_av = 0.5*(vecAy[i,j,k,2] + vecAy[i,j+1,k,2]) / volume[i,j,k]
                eta_z_av = 0.5*(vecAy[i,j,k,3] + vecAy[i,j+1,k,3]) / volume[i,j,k]
                beta = (eta_x_av^2 + eta_y_av^2 + eta_z_av^2)^0.5
                jbetaP[i,j,k] = beta / volume[i,j,k] * (2*mu[i,j,k]/rho)

                # gamma
                zeta_x_av = 0.5*(vecAz[i,j,k,1] + vecAz[i,j,k+1,1]) / volume[i,j,k]
                zeta_y_av = 0.5*(vecAz[i,j,k,2] + vecAz[i,j,k+1,2]) / volume[i,j,k]
                zeta_z_av = 0.5*(vecAz[i,j,k,3] + vecAz[i,j,k+1,3]) / volume[i,j,k]
                gamma = (zeta_x_av^2 + zeta_y_av^2 + zeta_z_av^2)^0.5
                jgammaP[i,j,k] = gamma / volume[i,j,k] * (2*mu[i,j,k]/rho)
            end
        end
    end
    return jalphaP, jbetaP, jgammaP
end

# ------------------------------------
# lusgs
# ------------------------------------
function lusgs(D, Lx, Ly, Lz, Ux, Uy, Uz, LdQ, UdQ, RHS_temp, I, dt, dtau, Qcon_hat, Qconn_hat, delta_Q,
                A_adv_hat_p, A_adv_hat_m, B_adv_hat_p, B_adv_hat_m, C_adv_hat_p, C_adv_hat_m, A_beta_shig, B_beta_shig, C_beta_shig,
                jalphaP, jbetaP, jgammaP, RHS, cellxmax, cellymax, cellzmax, volume, nval, icell)
       
    # calculate L and U
    for m in 1:nval
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:cellymax
                    for i in 1:cellxmax            
                        Lx[i,j,k,l,m] = dt*(A_adv_hat_p[i,j,k,l,m] + jalphaP[i,j,k]*I[l,m])
                        Ly[i,j,k,l,m] = dt*(B_adv_hat_p[i,j,k,l,m] + jbetaP[i,j,k]*I[l,m])
                        Lz[i,j,k,l,m] = dt*(C_adv_hat_p[i,j,k,l,m] + jgammaP[i,j,k]*I[l,m])
                        Ux[i,j,k,l,m] = dt*(A_adv_hat_m[i,j,k,l,m] - jalphaP[i,j,k]*I[l,m])
                        Uy[i,j,k,l,m] = dt*(B_adv_hat_m[i,j,k,l,m] - jbetaP[i,j,k]*I[l,m])
                        Uz[i,j,k,l,m] = dt*(C_adv_hat_m[i,j,k,l,m] - jgammaP[i,j,k]*I[l,m])
                    end
                end
            end
        end
    end

    for k in 1+icell:cellzmax-icell
        for j in 1+icell:cellymax-icell
            for i in 1+icell:cellxmax-icell
                for l in 1:nval
                    for m in 1:nval
                        LdQ[i,j,k,l] += Lx[i-1,j,k,l,m]*delta_Q[i-1,j,k,m] 
                                      + Ly[i,j-1,k,l,m]*delta_Q[i,j-1,k,m]
                                      + Lz[i,j,k-1,l,m]*delta_Q[i,j,k-1,m]
                        UdQ[i,j,k,l] += Ux[i+1,j,k,l,m]*delta_Q[i+1,j,k,m]
                                      + Uy[i,j+1,k,l,m]*delta_Q[i,j+1,k,m]
                                      + Uz[i,j,k+1,l,m]*delta_Q[i,j,k+1,m]
                    end
                end
            end
        end
    end

    # diagonal
    for k in 1+icell:cellzmax-icell
        for j in 1+icell:cellymax-icell
            for i in 1+icell:cellxmax-icell
                D[i,j,k] = dt/dtau[i,j,k] + 1.0 + dt*(A_beta_shig[i,j,k]+2*jalphaP[i,j,k] + B_beta_shig[i,j,k]+2*jbetaP[i,j,k] + C_beta_shig[i,j,k]+2*jgammaP[i,j,k])
            end
        end
    end
    
    # RHS
    for l in 1:nval    
        for k in 1+icell:cellzmax-icell
            for j in 1+icell:cellymax-icell
                for i in 1+icell:cellxmax-icell            
                    RHS_temp[i,j,k,l] = - (Qcon_hat[i,j,k,l]-Qconn_hat[i,j,k,l])*1.0 + dt*RHS[i,j,k,l]
                end
            end
        end
    end

    # lower sweep
    for l in 1:nval
        for k in 1+icell:cellzmax-icell
            for j in 1+icell:cellymax-icell
                for i in 1+icell:cellxmax-icell
                    delta_Q[i,j,k,l] = D[i,j,k]^(-1) * (RHS_temp[i,j,k,l]+LdQ[i,j,k,l])
                end
            end
        end
    end             
    
    # upepr sweep
    for l in 1:nval
        for k in 1+icell:cellzmax-icell
            for j in 1+icell:cellymax-icell
                for i in 1+icell:cellxmax-icell
                    delta_Q[i,j,k,l] = delta_Q[i,j,k,l] - D[i,j,k]^(-1) * UdQ[i,j,k,l]
                end
            end
        end
    end

    # reset
    for l in 1:nval
        for k in 1+icell:cellzmax-icell
            for j in 1+icell:cellymax-icell
                for i in 1+icell:cellxmax-icell
                    LdQ[i,j,k,l] = 0.0
                    UdQ[i,j,k,l] = 0.0                
                end
            end
        end
    end
    
    return delta_Q
end

