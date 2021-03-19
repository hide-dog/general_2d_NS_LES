# ------------------------------------
# advection term by AUSM
# ------------------------------------
function AUSM(E_adv_hat, F_adv_hat, G_adv_hat, QbaseF, QbaseB, QbaseU, QbaseD, QbaseL, QbaseR, 
                QconF, QconB, QconU, QconD, QconL, QconR, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz,
                specific_heat_ratio, volume, nval, Minf, ad_scheme, icell)
    g         = specific_heat_ratio
    temp_vec = zeros(nval)   # vecAx used for the pressure term
    Lpsi = zeros(nval)
    Rpsi = zeros(nval)

    mdot = 0.0
    ph = 0.0
    
    tvA = zeros(3)
    vA  = zeros(3)

    for k in 1+icell:cellzmax -icell
        for j in 1+icell:cellymax -icell
            for i in 1+icell:cellxmax+1 -icell
                volume_av = 0.5*(volume[i-1,j,k] + volume[i,j,k])

                # i-1 cell
                # 規格化
                tvA[1] = 0.5*(vecAx[i-1,j,k,1] + vecAx[i,j,k,1])
                tvA[2] = 0.5*(vecAx[i-1,j,k,2] + vecAx[i,j,k,2])
                tvA[3] = 0.5*(vecAx[i-1,j,k,3] + vecAx[i,j,k,3])
                vA[1] = tvA[1] / (tvA[1]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[2] = tvA[2] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[3] = tvA[3] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                
                rhoL = QbaseL[i,j,k,1]
                UL = QbaseL[i,j,k,2]*vA[1] + QbaseL[i,j,k,3]*vA[2] + QbaseL[i,j,k,4]*vA[3]
                pL = QbaseL[i,j,k,5]
                
                # i cell
                # 規格化
                tvA[1] = 0.5*(vecAx[i,j,k,1] + vecAx[i+1,j,k,1])
                tvA[2] = 0.5*(vecAx[i,j,k,2] + vecAx[i+1,j,k,2])
                tvA[3] = 0.5*(vecAx[i,j,k,3] + vecAx[i+1,j,k,3])
                vA[1] = tvA[1] / (tvA[1]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[2] = tvA[2] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[3] = tvA[3] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                
                rhoR = QbaseR[i,j,k,1]
                UR = QbaseR[i,j,k,2]*vA[1] + QbaseR[i,j,k,3]*vA[2] + QbaseR[i,j,k,4]*vA[3]
                pR = QbaseR[i,j,k,5]              
                
                # scheme
                if  ad_scheme == 1
                    mdot, ph = AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g, i, j)
                elseif ad_scheme == 2
                    mdot, ph = AUSM_plusup_half(rhoL, rhoR, UL, UR, pL, pR, Minf, g, i, j)
                elseif ad_scheme == 4
                    velocity = (0.5 * (QbaseL[i,j,k,2]^2 + QbaseL[i,j,k,3]^2 + QbaseL[i,j,k,4]^2 
                                        + QbaseR[i,j,k,2]^2 + QbaseR[i,j,k,3]^2 + QbaseR[i,j,k,4]^2))^0.5
                    mdot, ph = SLAU_half(rhoL, rhoR, UL, UR, pL, pR, velocity, g)
                end
                
                # flux half
                sqAx = (vecAx[i,j,k,1]^2 + vecAx[i,j,k,2]^2 + vecAx[i,j,k,3]^2)^0.5
                temp_vec[2] = vecAx[i,j,k,1] / sqAx
                temp_vec[3] = vecAx[i,j,k,2] / sqAx
                temp_vec[4] = vecAx[i,j,k,3] / sqAx

                for l in 1:nval
                    Lpsi[l] = QconL[i,j,k,l] / QconL[i,j,k,1]
                    Rpsi[l] = QconR[i,j,k,l] / QconR[i,j,k,1]
                end
                Lpsi[5] = (QconL[i,j,k,5] + QbaseL[i,j,k,5]) / QconL[i,j,k,1]
                Rpsi[5] = (QconR[i,j,k,5] + QbaseR[i,j,k,5]) / QconR[i,j,k,1]

                if mdot > 0
                    for l in 1:nval
                        E_adv_hat[i,j,k,l] = (mdot * Lpsi[l] + ph * temp_vec[l]) * sqAx
                    end
                else
                    for l in 1:nval
                        E_adv_hat[i,j,k,l] = (mdot * Rpsi[l] + ph * temp_vec[l]) * sqAx
                    end
                end            
            end
        end
    end

    for k in 1+icell:cellzmax -icell
        for j in 1+icell:cellymax+1 -icell
            for i in 1+icell:cellxmax -icell
                volume_av = 0.5*(volume[i,j-1,k] + volume[i,j,k])

                # j-1 cell
                # 規格化
                tvA[1] = 0.5*(vecAy[i,j-1,k,1] + vecAy[i,j,k,1])
                tvA[2] = 0.5*(vecAy[i,j-1,k,2] + vecAy[i,j,k,2])
                tvA[3] = 0.5*(vecAy[i,j-1,k,3] + vecAy[i,j,k,3])
                vA[1] = tvA[1] / (tvA[1]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[2] = tvA[2] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[3] = tvA[3] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                
                rhoD = QbaseD[i,j,k,1]
                VD = QbaseD[i,j,k,2]*vA[1] + QbaseD[i,j,k,3]*vA[2] + QbaseD[i,j,k,4]*vA[3]
                pD = QbaseD[i,j,k,5]
                
                # j cell
                tvA[1] = 0.5*(vecAy[i,j,k,1] + vecAy[i,j+1,k,1])
                tvA[2] = 0.5*(vecAy[i,j,k,2] + vecAy[i,j+1,k,2])
                tvA[3] = 0.5*(vecAy[i,j,k,3] + vecAy[i,j+1,k,3])
                vA[1] = tvA[1] / (tvA[1]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[2] = tvA[2] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[3] = tvA[3] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                
                rhoU = QbaseU[i,j,k,1]
                VU = QbaseU[i,j,k,2]*vA[1] + QbaseU[i,j,k,3]*vA[2] + QbaseU[i,j,k,4]*vA[3]
                pU = QbaseU[i,j,k,5]

                # scheme
                if  ad_scheme == 1
                    mdot, ph = AUSM_plus_half(rhoD, rhoU, VD, VU, pD, pU, g, i, j)
                elseif ad_scheme == 2
                    mdot, ph = AUSM_plusup_half(rhoD, rhoU, VD, VU, pD, pU, Minf, g, i, j)
                elseif ad_scheme == 4
                    velocity = (0.5 * (QbaseD[i,j,k,2]^2 + QbaseD[i,j,k,3]^2 + QbaseD[i,j,k,4]^2 + 
                                    QbaseU[i,j,k,2]^2 + QbaseU[i,j,k,3]^2 + QbaseU[i,j,k,4]^2))^0.5
                    mdot, ph = SLAU_half(rhoD, rhoU, VD, VU, pD, pU, velocity, specific_heat_ratio)
                end

                # flux half
                sqAy = (vecAy[i,j,k,1]^2 + vecAy[i,j,k,2]^2 + vecAy[i,j,k,3]^2)^0.5
                temp_vec[2] = vecAy[i,j,k,1] / sqAy
                temp_vec[3] = vecAy[i,j,k,2] / sqAy
                temp_vec[4] = vecAy[i,j,k,3] / sqAy

                for l in 1:nval
                    Lpsi[l] = QconD[i,j,k,l] / QconD[i,j,k,1]
                    Rpsi[l] = QconU[i,j,k,l] / QconU[i,j,k,1]
                end
                Lpsi[5] = (QconD[i,j,k,5] + QbaseD[i,j,k,5]) / QconD[i,j,k,1]
                Rpsi[5] = (QconU[i,j,k,5] + QbaseU[i,j,k,5]) / QconU[i,j,k,1]
                
                if mdot > 0
                    for l in 1:nval
                        F_adv_hat[i,j,k,l] = (mdot * Lpsi[l] + ph * temp_vec[l]) * sqAy
                    end
                else
                    for l in 1:nval
                        F_adv_hat[i,j,k,l] = (mdot * Rpsi[l] + ph * temp_vec[l]) * sqAy
                    end
                end
            end
        end
    end

    for k in 1+icell:cellzmax+1 -icell
        for j in 1+icell:cellymax -icell
            for i in 1+icell:cellxmax -icell
                volume_av = 0.5*(volume[i,j,k-1] + volume[i,j,k])

                # k-1 cell
                # 規格化
                tvA[1] = 0.5*(vecAz[i,j,k-1,1] + vecAz[i,j,k,1])
                tvA[2] = 0.5*(vecAz[i,j,k-1,2] + vecAz[i,j,k,2])
                tvA[3] = 0.5*(vecAz[i,j,k-1,3] + vecAz[i,j,k,3])
                vA[1] = tvA[1] / (tvA[1]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[2] = tvA[2] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[3] = tvA[3] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                
                rhoB = QbaseB[i,j,k,1]
                UB = QbaseB[i,j,k,2]*vA[1] + QbaseB[i,j,k,3]*vA[2] + QbaseB[i,j,k,4]*vA[3]
                pB = QbaseB[i,j,k,5]
                
                # k cell
                # 規格化
                tvA[1] = 0.5*(vecAz[i,j,k,1] + vecAz[i,j,k+1,1])
                tvA[2] = 0.5*(vecAz[i,j,k,2] + vecAz[i,j,k+1,2])
                tvA[3] = 0.5*(vecAz[i,j,k,3] + vecAz[i,j,k+1,3])
                vA[1] = tvA[1] / (tvA[1]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[2] = tvA[2] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                vA[3] = tvA[3] / (tvA[2]^2 + tvA[2]^2 + tvA[3]^2)^0.5
                
                rhoF = QbaseF[i,j,k,1]
                UF = QbaseF[i,j,k,2]*vA[1] + QbaseF[i,j,k,3]*vA[2] + QbaseF[i,j,k,4]*vA[3]
                pF = QbaseF[i,j,k,5]              
                
                # scheme
                if  ad_scheme == 1
                    mdot, ph = AUSM_plus_half(rhoB, rhoF, UB, UF, pB, pF, g, i, j)
                elseif ad_scheme == 2
                    mdot, ph = AUSM_plusup_half(rhoB, rhoF, UB, UF, pB, pF, Minf, g, i, j)
                elseif ad_scheme == 4
                    velocity = (0.5 * (QbaseF[i,j,k,2]^2 + QbaseF[i,j,k,3]^2 + QbaseF[i,j,k,4]^2 
                                        + QbaseB[i,j,k,2]^2 + QbaseB[i,j,k,3]^2 + QbaseB[i,j,k,4]^2))^0.5
                    mdot, ph = SLAU_half(rhoB, rhoF, UB, UF, pB, pF, velocity, g)
                end
                
                # flux half
                sqAz = (vecAz[i,j,k,1]^2 + vecAz[i,j,k,2]^2 + vecAz[i,j,k,3]^2)^0.5
                temp_vec[2] = vecAz[i,j,k,1] / sqAz
                temp_vec[3] = vecAz[i,j,k,2] / sqAz
                temp_vec[4] = vecAz[i,j,k,3] / sqAz

                for l in 1:nval
                    Lpsi[l] = QconB[i,j,k,l] / QconB[i,j,k,1]
                    Rpsi[l] = QconF[i,j,k,l] / QconF[i,j,k,1]
                end
                Lpsi[5] = (QconB[i,j,k,5] + QbaseB[i,j,k,5]) / QconB[i,j,k,1]
                Rpsi[5] = (QconF[i,j,k,5] + QbaseF[i,j,k,5]) / QconF[i,j,k,1]

                if mdot > 0
                    for l in 1:nval
                        G_adv_hat[i,j,k,l] = (mdot * Lpsi[l] + ph * temp_vec[l]) * sqAz
                    end
                else
                    for l in 1:nval
                        G_adv_hat[i,j,k,l] = (mdot * Rpsi[l] + ph * temp_vec[l]) * sqAz
                    end
                end            
            end
        end
    end
    
    return E_adv_hat, F_adv_hat, G_adv_hat
end

function AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g,i,j)
    # param
    beta  = 1/8
    alpha = 3/16
    
    # L, R
    aL = (g * pL / rhoL)^0.5
    aR = (g * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )
    rhoh = 0.5*( rhoL + rhoR )

    ML = UL/ah
    MR = UR/ah

    M_p4 = 0
    M_m4 = 0
    p_p5 = 0
    p_m5 = 0
    if abs(ML) >= 1
        M_p4 = 0.5*(ML + abs(ML))
        p_p5 = 0.5*(ML + abs(ML)) / ML
    else
        M_p4 = 0.25*(ML+1)^2 + beta*(ML^2-1)^2
        p_p5 = 0.25*(ML+1)^2 * (2-ML) + alpha*ML*(ML^2-1)^2
    end
    
    if abs(MR) >= 1
        M_m4 = 0.5*(MR - abs(MR))
        p_m5 = 0.5*(MR - abs(MR)) / MR
    else
        M_m4 = -0.25*(MR-1)^2 - beta*(MR^2-1)^2
        p_m5 = 0.25*(MR-1)^2 * (2+MR) - alpha*MR*(MR^2-1)^2
    end
    
    # Mh = M_p4 + M_m4 + Mp/fa
    Mh = M_p4 + M_m4
#=
    if i==100 && j==200
        println(" MMM_AUSM+ ")
        println(M_p4)
        println(M_m4)
    end
  =#  
    # mdot half
    mdot = ah * Mh
    if Mh > 0
        mdot = mdot * rhoL
    else
        mdot = mdot * rhoR
    end
    
    # p half
    ph = p_p5*pL + p_m5*pR

    return mdot, ph
end

function AUSM_plusup_half(rhoL, rhoR, UL, UR, pL, pR, Minf, g, i,j)
    # param
    beta  = 1/8
    kp    = 0.25
    sigma = 1.0
    ku    = 0.75
    
    # L, R
    aL = (g * pL / rhoL)^0.5
    aR = (g * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )
    rhoh = 0.5*( rhoL + rhoR )

    Mbar = (( UL^2 + UR^2 ) / ( 2 * ah^2 ))^0.5
    Mo   = (min(1,max( Mbar^2, Minf^2 )))^0.5
    fa   = Mo * (2-Mo)

    alpha = 3/16 * ( -4 + 5*fa^2 )

    ML = UL/ah
    MR = UR/ah

    M_p4 = 0
    M_m4 = 0
    p_p5 = 0
    p_m5 = 0
    if abs(ML) >= 1
        M_p4 = 0.5*(ML + abs(ML))
        p_p5 = 0.5*(ML + abs(ML)) / ML
    else
        Mtp  = 0.25*(ML + 1)^2
        Mtm  = -0.25*(ML - 1)^2
        M_p4 = Mtp * (1 - 16*beta*Mtm)
        p_p5 = Mtp * ((2-ML) - 16*alpha*ML*Mtm)
    end
    

    if abs(MR) >= 1
        M_m4 = 0.5*(MR - abs(MR))
        p_m5 = 0.5*(MR - abs(MR)) / MR
    else
        Mtp  = 0.25*(MR + 1)^2
        Mtm  = -0.25*(MR - 1)^2
        M_m4 = Mtm * (1 + 16*beta*Mtp)
        p_m5 = Mtm * ((-2-MR) + 16*alpha*MR*Mtp)
    end
    
    # M half
    Mp = -kp * max( 1 - sigma*Mbar^2, 0.0) * (pR - pL)/(rhoh*ah^2)
    
    Mh = M_p4 + M_m4 + Mp/fa
    
    #=
    if i==100 && j==200
        println(" MMM_AUSM+up ")
        println(M_p4)
        println(M_m4)
        println(Mp/fa)
    end
    =#
        
    # mdot half
    mdot = ah * Mh
    if Mh > 0
        mdot = mdot * rhoL
    else
        mdot = mdot * rhoR
    end
    
    # p half
    pu = -ku * p_p5 * p_m5 * (rhoL + rhoR) * ah *(UR - UL)

    ph = p_p5*pL + p_m5*pR + fa * pu
    return mdot, ph
end

function set_Minf(bdcon, specific_heat_ratio, Rd, nval)
    Minf = 0.0
    check = 0
    for i in 1:6
        if Int(bdcon[i][1]) == 0 || Int(bdcon[i][1]) == 5
            rho = bdcon[i][2]
            u = bdcon[i][3]
            v = bdcon[i][4]
            w = bdcon[i][5]
            p = bdcon[i][6]

            a = (specific_heat_ratio * p / rho)^0.5
            Minf = (u^2 + v^2 + w^2)^0.5 / a

            println(" Minf = " * string(Minf))

            check = 1

        elseif Int(bdcon[i][1]) == 6 || Int(bdcon[i][1]) == 88
            rho = bdcon[i][2]
            u = bdcon[i][3]
            v = bdcon[i][4]
            w = bdcon[i][5]
            T = bdcon[i][nval+2]

            p = (rho*Rd) * T

            a = (specific_heat_ratio * p / rho)^0.5
            Minf = (u^2 + v^2 + w^2)^0.5 / a
            
            println(" Minf = " * string(Minf))

            check = 1
        end
    end

    if check == 0
        println(" Minf error ")
        println("  ")
        println(" if you don't use inlet condition, ")
        println(" don't use AUSM+up ")
        println("  ")
        throw(UndefVarError(:x))
    end

    return Minf
end


function SLAU_half(rhoL, rhoR, UL, UR, pL, pR, velocity, specific_heat_ratio)
    mdot =1.0
    ph =1.0

    # L, R
    aL = (specific_heat_ratio * pL / rhoL)^0.5
    aR = (specific_heat_ratio * pR / rhoR)^0.5

    # half
    ah   = 0.5*( aL + aR )

    ML = UL/ah
    MR = UR/ah

    temp1 = -max( min( ML, 0.0 ), -1.0)
    temp2 =  min( max( MR, 0.0 ), 1.0)
    g = temp1 * temp2

    barUL = (1-g) * (rhoL*abs(UL) + rhoR*abs(UR))/(rhoL + rhoR) + g*abs(UL)
    barUR = (1-g) * (rhoL*abs(UL) + rhoR*abs(UR))/(rhoL + rhoR) + g*abs(UR)

    hatM = min(1.0, 1/ah * velocity)
    chi  = (1-hatM)^2
        
    beta_a = 0
    beta_b = 0
    if abs(MR) >= 1
        beta_a = 0.5 * (1+sign(-MR))
    else
        beta_a = 0.25 * (2+MR) * (MR-1)^2
    end
    if abs(ML) >= 1
        beta_b = 0.5 * (1+sign(ML))
    else
        beta_b = 0.25 * (2-ML) * (ML+1)^2
    end
    
    # mdot half
    mdot = 0.5 * (rhoL*(UL+abs(barUL)) + rhoR*(UR-abs(barUR)) - chi/ah*(pR-pL))
    
    # p half
    ph = 0.5*(pR+pL) + 0.5*(beta_b-beta_a)*(pL-pR) + (1-chi) * (beta_b+beta_a-1) * 0.5*(pR+pL)

    # SLAU2
    #ph = 0.5*(pR+pL) + 0.5*(beta_b-beta_a) + ((UR^2+UL^2)/2)^0.5 * (beta_b+beta_a-1) * 0.5*(rhoR+rhoL)*ah

    return mdot, ph
end
