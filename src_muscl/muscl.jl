# ------------------------------------
# muscl
# ------------------------------------
function muscl(Qbase, QbaseU, QbaseD, QbaseL, QbaseR, QconU, QconD, QconL, QconR, 
                cellxmax, cellymax, vecAx, vecAy, nval, icell, specific_heat_ratio, volume, nodes)
    
    #=
    Q system
    VecA system is sama
    Q is defined on(near) boundary

            nodes[2,3]  QU[2,3]
                |--------------------|
                |       QD[2,3]      |
                |                    |
        QL[2,2] | QR[2,2]    QL[3,2] | QR[3,2]
                |      cell(2,2)     |
                |                    |
                |       QU[2,2]      |
                |--------------------|
            nodes[2,2]  QD[2,2]   nodes[3,2]

    =#
    for l in 1:nval
        for j in 1+icell-1:cellymax-icell+1
            for i in 1+icell:cellxmax-icell+1
                1+icell-1:cellymax+1 -icell +1
                #vecAx_xav = 0.25*( vecAx[i,j,1] + vecAx[i-1,j,1] + vecAx[i,j+1,1] + vecAx[i-1,j+1,1] )
                #vecAx_yav = 0.25*( vecAx[i,j,2] + vecAx[i-1,j,2] + vecAx[i,j+1,2] + vecAx[i-1,j+1,2] )
                #volume_av = 0.5*( volume[i-1,j] + volume[i-1,j+1])

                # ------------------------------
                # L
                # ------------------------------
                #=
                vim = 0.5*( volume[i-2,j] + volume[i-2,j-1])
                vi  = 0.5*( volume[i-1,j] + volume[i-1,j-1])
                vip = 0.5*( volume[i,  j] + volume[i,  j-1])
                
                vim_p = 0.5*( volume[i-2,j+1] + volume[i-2,j])
                vi_p  = 0.5*( volume[i-1,j+1] + volume[i-1,j])
                vip_p = 0.5*( volume[i,  j+1] + volume[i,  j])


                
                lim = 0.5 * ((1/vecAy[i-2,  j,1]^2 + 1/vecAy[i-2,  j,2]^2)^0.5 * vim_p
                            +(1/vecAy[i-2,j-1,1]^2 + 1/vecAy[i-2,j-1,2]^2)^0.5 * vim)
                li  = 0.5 * ((1/vecAy[i-1,  j,1]^2 + 1/vecAy[i-1,  j,2]^2)^0.5 * vi_p
                            +(1/vecAy[i-1,j-1,1]^2 + 1/vecAy[i-1,j-1,2]^2)^0.5 * vi)
                lip = 0.5 * ((1/vecAy[  i,  j,1]^2 + 1/vecAy[  i,  j,2]^2)^0.5 * vip_p
                            +(1/vecAy[  i,j-1,1]^2 + 1/vecAy[  i,j-1,2]^2)^0.5 * vip)
                =#

                lim = 0.5 * (((nodes[i-2,j+1,1] - nodes[i-1,j+1,1])^2 + (nodes[i-2,j+1,2] - nodes[i-1,j+1,2])^2)^0.5 
                            +((nodes[i-2,  j,1] - nodes[i-1,  j,1])^2 + (nodes[i-2,  j,2] - nodes[i-1,  j,2])^2)^0.5)
                li  = 0.5 * (((nodes[i-1,j+1,1] - nodes[  i,j+1,1])^2 + (nodes[i-1,j+1,2] - nodes[  i,j+1,2])^2)^0.5 
                            +((nodes[i-1,  j,1] - nodes[  i,  j,1])^2 + (nodes[i-1,  j,2] - nodes[  i,  j,2])^2)^0.5) 
                lip = 0.5 * (((nodes[  i,j+1,1] - nodes[i+1,j+1,1])^2 + (nodes[  i,j+1,2] - nodes[i+1,j+1,2])^2)^0.5 
                            +((nodes[  i,  j,1] - nodes[i+1,  j,1])^2 + (nodes[  i,  j,2] - nodes[i+1,  j,2])^2)^0.5)

                dip = 2*(Qbase[i,j,l]   - Qbase[i-1,j,l]) / (li  + lip)
                dim = 2*(Qbase[i-1,j,l] - Qbase[i-2,j,l]) / (lim + li)
                QbaseL[i,j,l] = Qbase[i-1,j,l] + 0.5*li*minmod(dim, dip) 
                

                # ------------------------------
                # R
                # ------------------------------
                #=
                vip   = 0.5*( volume[i+1,  j] + volume[i+1,j-1])
                vip_p = 0.5*( volume[i+1,j+2] + volume[i+1,  j])

                lim = li
                li  = lip
                lip = 0.5 * ((1/vecAy[i,j+2,1]^2 + 1/vecAy[i,j+2,2]^2)^0.5 * vip_p
                            +(1/vecAy[i,j+1,1]^2 + 1/vecAy[i,j+1,2]^2)^0.5 * vip)
                =#
                lim = li
                li  = lip
                lip = 0.5 * (((nodes[i+1,j+1,1] - nodes[i+2,j+1,1])^2 + (nodes[i+1,j+1,2] - nodes[i+2,j+1,2])^2)^0.5 
                            +((nodes[i+1,  j,1] - nodes[i+2,  j,1])^2 + (nodes[i+1,  j,2] - nodes[i+2,  j,2])^2)^0.5)
                
                dip = 2*(Qbase[i+1,j,l] - Qbase[i,j,l]) / (li  + lip)
                dim = 2*(Qbase[i,j,l]   - Qbase[i-1,j,l]) / (lim + li)
                QbaseR[i,j,l] = Qbase[i,j,l]  - 0.5*li*minmod(dip, dim) 
            end
        end
    end

    for l in 1:nval
        for j in 1+icell:cellymax-icell + 1
            for i in 1+icell-1:cellxmax-icell+1
                # ------------------------------
                # D
                # ------------------------------
                #=
                vim = 0.5*( volume[i,j-2] + volume[i-1,j-2])
                vi  = 0.5*( volume[i,j-1] + volume[i-1,j-1])
                vip = 0.5*( volume[i,  j] + volume[i-1,  j])
                
                vim_p = 0.5*( volume[i+1,j-2] + volume[i,j-2])
                vi_p  = 0.5*( volume[i+1,j-1] + volume[i,j-1])
                vip_p = 0.5*( volume[i+1,  j] + volume[i,  j])

                lim = 0.5 * ((1/vecAx[  i,j-2,1]^2 + 1/vecAx[  i,j-2,2]^2)^0.5 * vim_p
                            +(1/vecAx[i-1,j-2,1]^2 + 1/vecAx[i-1,j-2,2]^2)^0.5 * vim)
                li  = 0.5 * ((1/vecAx[  i,j-1,1]^2 + 1/vecAx[  i,j-1,2]^2)^0.5 * vi_p
                            +(1/vecAx[i-1,j-1,1]^2 + 1/vecAx[i-1,j-1,2]^2)^0.5 * vi)
                lip = 0.5 * ((1/vecAx[  i,  j,1]^2 + 1/vecAx[  i,  j,2]^2)^0.5 * vip_p 
                            +(1/vecAx[i-1,  j,1]^2 + 1/vecAx[i-1,  j,2]^2)^0.5 * vip)
                =#

                lim = 0.5 * (((nodes[i+1,j-2,1] - nodes[i+1,j-1,1])^2 + (nodes[i+1,j-2,2] - nodes[i+1,j-1,2])^2)^0.5 
                            +((nodes[  i,j-2,1] - nodes[  i,j-1,1])^2 + (nodes[  i,j-2,2] - nodes[  i,j-1,2])^2)^0.5)
                li  = 0.5 * (((nodes[i+1,j-1,1] - nodes[i+1,  j,1])^2 + (nodes[i+1,j-1,2] - nodes[i+1,  j,2])^2)^0.5 
                            +((nodes[  i,j-1,1] - nodes[  i,  j,1])^2 + (nodes[  i,j-1,2] - nodes[  i,  j,2])^2)^0.5) 
                lip = 0.5 * (((nodes[i+1,  j,1] - nodes[i+1,j+1,1])^2 + (nodes[i+1,  j,2] - nodes[i+1,j+1,2])^2)^0.5 
                            +((nodes[  i,  j,1] - nodes[  i,j+1,1])^2 + (nodes[  i,  j,2] - nodes[  i,j+1,2])^2)^0.5)

                dip = 2*(Qbase[i,j,l]   - Qbase[i,j-1,l]) / (li  + lip)
                dim = 2*(Qbase[i,j-1,l] - Qbase[i,j-2,l]) / (lim + li)
                QbaseD[i,j,l] = Qbase[i-1,j,l] + 0.5*li*minmod(dim, dip) 

                # ------------------------------
                # U
                # ------------------------------
                #=
                vip   = 0.5*( volume[  i,j+1] + volume[i-1,j+1])
                vip_p = 0.5*( volume[i+1,j+1] + volume[  i,j+1])

                lim = li
                li  = lip
                lip = 0.5 * ((1/vecAx[i+2,j,1]^2 + 1/vecAx[i+2,j,2]^2)^0.5 * vip_p
                            +(1/vecAx[i+1,j,1]^2 + 1/vecAx[  i,j,2]^2)^0.5 * vip)
                =#
                lim = li
                li  = lip
                lip = 0.5 * (((nodes[i+1,j+1,1] - nodes[i+1,j+2,1])^2 + (nodes[i+1,j+1,2] - nodes[i+1,j+2,2])^2)^0.5 
                            +((nodes[  i,j+1,1] - nodes[  i,j+2,1])^2 + (nodes[  i,j+1,2] - nodes[  i,j+2,2])^2)^0.5)

                dip = 2*(Qbase[i,j+1,l] - Qbase[i,j,l])   / (li  + lip)
                dim = 2*(Qbase[i,j,l]   - Qbase[i,j-1,l]) / (lim + li)
                QbaseU[i,j,l] = Qbase[i,j,l]  - 0.5*li*minmod(dip, dim) 
            end
        end
    end

    QconU = base_to_conservative(QbaseU, QconU, cellxmax, cellymax, specific_heat_ratio)
    QconD = base_to_conservative(QbaseD, QconD, cellxmax, cellymax, specific_heat_ratio)
    QconL = base_to_conservative(QbaseL, QconL, cellxmax, cellymax, specific_heat_ratio)
    QconR = base_to_conservative(QbaseR, QconR, cellxmax, cellymax, specific_heat_ratio)
    
    #println(Qbase[5,5,1])
    #println(QbaseL[:,5,1])
    #println(QbaseL[5,:,1])
    #println(QbaseD[:,5,1])
    #println(QbaseD[5,:,1])
    #println(QbaseR[5,5,1])
    #println(Qbase[6,5,1])
    #=
    println(Qbase[5,1,1])
    println(QbaseL[5,2,1])
    println(Qbase[5,2,1])
    =#
    
    return QbaseU, QbaseD, QbaseL, QbaseR, QconU, QconD, QconL, QconR
end

function minmod(a,b)
    ans = sign(a) * max(0, min(a*sign(b), b*sign(a)))
    return ans
end