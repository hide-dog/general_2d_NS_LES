# ------------------------------------
# advection term by AUSM+
# ------------------------------------
function AUSM_plus(E_adv_hat, F_adv_hat, Qbase, Qcon, cellxmax, cellymax, vecAx, vecAy, specific_heat_ratio, volume, nval)
    g         = specific_heat_ratio
    cell_vecAx = zeros(2)     # vecAx at cell center
    cell_vecAy = zeros(2)     # vecAy at cell center
    temp_vecX = zeros(nval)   # vecAx used for the pressure term
    temp_vecY = zeros(nval)   # vecAx used for the pressure term
    Lpsi = zeros(nval)        
    Rpsi = zeros(nval)
    
    for j in 2:cellymax -1
        for i in 2:cellxmax+1 -1
            # i-1 cell
            cell_vecAx[1] = 0.5*(vecAx[i-1,j,1]+vecAx[i,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i-1,j,2]+vecAx[i,j,2])
            
            rhoL = Qbase[i-1,j,1]
            UL = Qcon[i-1,j,2] / Qcon[i-1,j,1]*cell_vecAx[1] + Qcon[i-1,j,3] / Qcon[i-1,j,1]*cell_vecAx[2]
            pL = Qbase[i-1,j,4]
            
            # i cell
            cell_vecAx[1] = 0.5*(vecAx[i,j,1]+vecAx[i+1,j,1])
            cell_vecAx[2] = 0.5*(vecAx[i,j,2]+vecAx[i+1,j,2])
            
            rhoR = Qbase[i,j,1]
            UR = Qcon[i,j,2] / Qcon[i,j,1]*cell_vecAx[1] + Qcon[i,j,3] / Qcon[i,j,1]*cell_vecAx[2]
            pR = Qbase[i,j,4]
            
            # AUSM+
            mdot, ph = AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g)

            # flux half
            temp_vecX[2] = vecAx[i,j,1]
            temp_vecX[3] = vecAx[i,j,2]

            for l in 1:nval
                Lpsi[l] = Qcon[i-1,j,l] / Qcon[i-1,j,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i-1,j,4] + Qbase[i-1,j,4]) / Qcon[i-1,j,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]

            if mdot > 0
                for l in 1:nval
                    E_adv_hat[i,j,l] = mdot * Lpsi[l] + ph * temp_vecX[l]
                end
            else
                for l in 1:nval
                    E_adv_hat[i,j,l] = mdot * Rpsi[l] + ph * temp_vecX[l]
                end
            end            
        end
    end

    for j in 2:cellymax+1 -1
        for i in 2:cellxmax -1
            # j-1 cell
            cell_vecAy[1] = 0.5*(vecAy[i,j-1,1]+vecAy[i,j,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j-1,2]+vecAy[i,j,2])
            
            rhoL = Qbase[i,j-1,1]
            VL = Qbase[i,j-1,2]*cell_vecAy[1] + Qbase[i,j-1,3]*cell_vecAy[2]
            pL = Qbase[i,j-1,4]
            
            # j cell
            cell_vecAy[1] = 0.5*(vecAy[i,j,1]+vecAy[i,j+1,1])
            cell_vecAy[2] = 0.5*(vecAy[i,j,2]+vecAy[i,j+1,2])
            
            rhoR = Qbase[i,j,1]
            VR = Qbase[i,j,2]*cell_vecAy[1] + Qbase[i,j,3]*cell_vecAy[2]
            pR = Qbase[i,j,4]

            # AUSM+up
            mdot, ph = AUSM_plus_half(rhoL, rhoR, VL, VR, pL, pR, g)
            
            # flux half
            temp_vecY[2] = vecAy[i,j,1]
            temp_vecY[3] = vecAy[i,j,2]

            for l in 1:nval
                Lpsi[l] = Qcon[i,j-1,l] / Qcon[i,j-1,1]
                Rpsi[l] = Qcon[i,j,l]   / Qcon[i,j,1]
            end
            Lpsi[4] = (Qcon[i,j-1,4] + Qbase[i,j-1,4]) / Qcon[i,j-1,1]
            Rpsi[4] = (Qcon[i,j,4]   + Qbase[i,j,4])   / Qcon[i,j,1]
            
            if mdot > 0
                for l in 1:nval
                    F_adv_hat[i,j,l] = mdot * Lpsi[l] + ph * temp_vecY[l]
                end
            else
                for l in 1:nval
                    F_adv_hat[i,j,l] = mdot * Rpsi[l] + ph * temp_vecY[l]
                end
            end
        end
    end
    
    return E_adv_hat, F_adv_hat
end

function AUSM_plus_half(rhoL, rhoR, UL, UR, pL, pR, g)
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
        M_p4 = 0.5*(ML+1)^2 + beta*(ML^2-1)^2
        p_p5 = 0.25*(ML+1)^2 * (2-ML) + alpha*ML*(ML^2-1)^2
    end
    
    if abs(MR) >= 1
        M_m4 = 0.5*(MR - abs(MR))
        p_m5 = 0.5*(MR - abs(MR)) / MR
    else
        M_m4 = -0.5*(MR-1)^2 - beta*(MR^2-1)^2
        p_m5 = 0.25*(MR-1)^2 * (2+MR) - alpha*MR*(MR^2-1)^2
    end
    
    # Mh = M_p4 + M_m4 + Mp/fa
    Mh = M_p4 + M_m4
    
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