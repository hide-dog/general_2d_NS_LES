# ------------------------------------
# muscl
# ------------------------------------
function muscl(Qbase, QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB, QconU, QconD, QconL, QconR, QconF, QconB, 
                cellxmax, cellymax, cellzmax, nval, icell, specific_heat_ratio, volume, nodes)
    
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
        for k in 1+icell-1:cellzmax-icell+1
            for j in 1+icell-1:cellymax-icell+1
                for i in 1+icell:cellxmax-icell+1
                    # ------------------------------
                    # L
                    # ------------------------------
                    
                    lim = 0.25 * (((nodes[i-2,j+1,k+1,1] - nodes[i-1,j+1,k+1,1])^2 
                                 + (nodes[i-2,j+1,k+1,2] - nodes[i-1,j+1,k+1,2])^2 
                                 + (nodes[i-2,j+1,k+1,3] - nodes[i-1,j+1,k+1,3])^2)^0.5 
                                 +((nodes[i-2,  j,k+1,1] - nodes[i-1,  j,k+1,1])^2 
                                 + (nodes[i-2,  j,k+1,2] - nodes[i-1,  j,k+1,2])^2 
                                 + (nodes[i-2,  j,k+1,3] - nodes[i-1,  j,k+1,3])^2)^0.5
                                 +((nodes[i-2,j+1,  k,1] - nodes[i-1,j+1,  k,1])^2 
                                 + (nodes[i-2,j+1,  k,2] - nodes[i-1,j+1,  k,2])^2 
                                 + (nodes[i-2,j+1,  k,3] - nodes[i-1,j+1,  k,3])^2)^0.5 
                                 +((nodes[i-2,  j,  k,1] - nodes[i-1,  j,  k,1])^2 
                                 + (nodes[i-2,  j,  k,2] - nodes[i-1,  j,  k,2])^2 
                                 + (nodes[i-2,  j,  k,3] - nodes[i-1,  j,  k,3])^2)^0.5)
                    
                    li  = 0.25 * (((nodes[i-1,j+1,k+1,1] - nodes[  i,j+1,k+1,1])^2 
                                 + (nodes[i-1,j+1,k+1,2] - nodes[  i,j+1,k+1,2])^2 
                                 + (nodes[i-1,j+1,k+1,3] - nodes[  i,j+1,k+1,3])^2)^0.5 
                                 +((nodes[i-1,  j,k+1,1] - nodes[  i,  j,k+1,1])^2 
                                 + (nodes[i-1,  j,k+1,2] - nodes[  i,  j,k+1,2])^2 
                                 + (nodes[i-1,  j,k+1,3] - nodes[  i,  j,k+1,3])^2)^0.5
                                 +((nodes[i-1,j+1,  k,1] - nodes[  i,j+1,  k,1])^2 
                                 + (nodes[i-1,j+1,  k,2] - nodes[  i,j+1,  k,2])^2 
                                 + (nodes[i-1,j+1,  k,3] - nodes[  i,j+1,  k,3])^2)^0.5 
                                 +((nodes[i-1,  j,  k,1] - nodes[  i,  j,  k,1])^2 
                                 + (nodes[i-1,  j,  k,2] - nodes[  i,  j,  k,2])^2 
                                 + (nodes[i-1,  j,  k,3] - nodes[  i,  j,  k,3])^2)^0.5)
                    
                    lip = 0.25 * (((nodes[  i,j+1,k+1,1] - nodes[i+1,j+1,k+1,1])^2 
                                 + (nodes[  i,j+1,k+1,2] - nodes[i+1,j+1,k+1,2])^2 
                                 + (nodes[  i,j+1,k+1,3] - nodes[i+1,j+1,k+1,3])^2)^0.5 
                                 +((nodes[  i,  j,k+1,1] - nodes[i+1,  j,k+1,1])^2 
                                 + (nodes[  i,  j,k+1,2] - nodes[i+1,  j,k+1,2])^2 
                                 + (nodes[  i,  j,k+1,3] - nodes[i+1,  j,k+1,3])^2)^0.5
                                 +((nodes[  i,j+1,  k,1] - nodes[i+1,j+1,  k,1])^2 
                                 + (nodes[  i,j+1,  k,2] - nodes[i+1,j+1,  k,2])^2 
                                 + (nodes[  i,j+1,  k,3] - nodes[i+1,j+1,  k,3])^2)^0.5 
                                 +((nodes[  i,  j,  k,1] - nodes[i+1,  j,  k,1])^2 
                                 + (nodes[  i,  j,  k,2] - nodes[i+1,  j,  k,2])^2 
                                 + (nodes[  i,  j,  k,3] - nodes[i+1,  j,  k,3])^2)^0.5)
                    

                    dip = 2*(Qbase[i,j,k,l]   - Qbase[i-1,j,k,l]) / (li  + lip)
                    dim = 2*(Qbase[i-1,j,k,l] - Qbase[i-2,j,k,l]) / (lim + li)
                    QbaseL[i,j,k,l] = Qbase[i-1,j,k,l] + 0.5*li*minmod(dim, dip) 
                    

                    # ------------------------------
                    # R
                    # ------------------------------

                    lim = li
                    li  = lip
                    lip = 0.25 * (((nodes[i+1,j+1,k+1,1] - nodes[i+2,j+1,k+1,1])^2 
                                 + (nodes[i+1,j+1,k+1,2] - nodes[i+2,j+1,k+1,2])^2 
                                 + (nodes[i+1,j+1,k+1,3] - nodes[i+2,j+1,k+1,3])^2)^0.5 
                                 +((nodes[i+1,  j,k+1,1] - nodes[i+2,  j,k+1,1])^2 
                                 + (nodes[i+1,  j,k+1,2] - nodes[i+2,  j,k+1,2])^2 
                                 + (nodes[i+1,  j,k+1,3] - nodes[i+2,  j,k+1,3])^2)^0.5
                                 +((nodes[i+1,j+1,  k,1] - nodes[i+2,j+1,  k,1])^2 
                                 + (nodes[i+1,j+1,  k,2] - nodes[i+2,j+1,  k,2])^2 
                                 + (nodes[i+1,j+1,  k,3] - nodes[i+2,j+1,  k,3])^2)^0.5 
                                 +((nodes[i+1,  j,  k,1] - nodes[i+2,  j,  k,1])^2 
                                 + (nodes[i+1,  j,  k,2] - nodes[i+2,  j,  k,2])^2 
                                 + (nodes[i+1,  j,  k,3] - nodes[i+2,  j,  k,3])^2)^0.5)
                    
                    dip = 2*(Qbase[i+1,j,k,l] - Qbase[i,j,k,l]) / (li  + lip)
                    dim = 2*(Qbase[i,j,k,l]   - Qbase[i-1,j,k,l]) / (lim + li)
                    QbaseR[i,j,k,l] = Qbase[i,j,k,l]  - 0.5*li*minmod(dip, dim) 
                end
            end
        end
    end

    for l in 1:nval
        for k in 1+icell-1:cellzmax-icell+1
            for j in 1+icell:cellymax-icell + 1
                for i in 1+icell-1:cellxmax-icell+1
                    # ------------------------------
                    # D
                    # ------------------------------

                    lim = 0.25 * (((nodes[i+1,j-2,k+1,1] - nodes[i+1,j-1,k+1,1])^2 
                                 + (nodes[i+1,j-2,k+1,2] - nodes[i+1,j-1,k+1,2])^2 
                                 + (nodes[i+1,j-2,k+1,3] - nodes[i+1,j-1,k+1,3])^2)^0.5 
                                 +((nodes[  i,j-2,k+1,1] - nodes[  i,j-1,k+1,1])^2 
                                 + (nodes[  i,j-2,k+1,2] - nodes[  i,j-1,k+1,2])^2 
                                 + (nodes[  i,j-2,k+1,3] - nodes[  i,j-1,k+1,3])^2)^0.5
                                 +((nodes[i+1,j-2,  k,1] - nodes[i+1,j-1,  k,1])^2 
                                 + (nodes[i+1,j-2,  k,2] - nodes[i+1,j-1,  k,2])^2 
                                 + (nodes[i+1,j-2,  k,3] - nodes[i+1,j-1,  k,3])^2)^0.5 
                                 +((nodes[  i,j-2,  k,1] - nodes[  i,j-1,  k,1])^2 
                                 + (nodes[  i,j-2,  k,2] - nodes[  i,j-1,  k,2])^2 
                                 + (nodes[  i,j-2,  k,3] - nodes[  i,j-1,  k,3])^2)^0.5)
                    
                    li  = 0.25 * (((nodes[i+1,j-1,k+1,1] - nodes[i+1,  j,k+1,1])^2 
                                 + (nodes[i+1,j-1,k+1,2] - nodes[i+1,  j,k+1,2])^2 
                                 + (nodes[i+1,j-1,k+1,3] - nodes[i+1,  j,k+1,3])^2)^0.5 
                                 +((nodes[  i,j-1,k+1,1] - nodes[  i,  j,k+1,1])^2 
                                 + (nodes[  i,j-1,k+1,2] - nodes[  i,  j,k+1,2])^2 
                                 + (nodes[  i,j-1,k+1,3] - nodes[  i,  j,k+1,3])^2)^0.5
                                 +((nodes[i+1,j-1,  k,1] - nodes[i+1,  j,  k,1])^2 
                                 + (nodes[i+1,j-1,  k,2] - nodes[i+1,  j,  k,2])^2 
                                 + (nodes[i+1,j-1,  k,3] - nodes[i+1,  j,  k,3])^2)^0.5 
                                 +((nodes[  i,j-1,  k,1] - nodes[  i,  j,  k,1])^2 
                                 + (nodes[  i,j-1,  k,2] - nodes[  i,  j,  k,2])^2 
                                 + (nodes[  i,j-1,  k,3] - nodes[  i,  j,  k,3])^2)^0.5)
                    
                    lip = 0.25 * (((nodes[i+1,  j,k+1,1] - nodes[i+1,j+1,k+1,1])^2 
                                 + (nodes[i+1,  j,k+1,2] - nodes[i+1,j+1,k+1,2])^2 
                                 + (nodes[i+1,  j,k+1,3] - nodes[i+1,j+1,k+1,3])^2)^0.5 
                                 +((nodes[  i,  j,k+1,1] - nodes[  i,j+1,k+1,1])^2 
                                 + (nodes[  i,  j,k+1,2] - nodes[  i,j+1,k+1,2])^2 
                                 + (nodes[  i,  j,k+1,3] - nodes[  i,j+1,k+1,3])^2)^0.5
                                 +((nodes[i+1,  j,  k,1] - nodes[i+1,j+1,  k,1])^2 
                                 + (nodes[i+1,  j,  k,2] - nodes[i+1,j+1,  k,2])^2 
                                 + (nodes[i+1,  j,  k,3] - nodes[i+1,j+1,  k,3])^2)^0.5 
                                 +((nodes[  i,  j,  k,1] - nodes[  i,j+1,  k,1])^2 
                                 + (nodes[  i,  j,  k,2] - nodes[  i,j+1,  k,2])^2 
                                 + (nodes[  i,  j,  k,3] - nodes[  i,j+1,  k,3])^2)^0.5)
                    
                    dip = 2*(Qbase[i,j,k,l]   - Qbase[i,j-1,k,l]) / (li  + lip)
                    dim = 2*(Qbase[i,j-1,k,l] - Qbase[i,j-2,k,l]) / (lim + li)
                    QbaseD[i,j,k,l] = Qbase[i,j-1,k,l] + 0.5*li*minmod(dim, dip) 

                    # ------------------------------
                    # U
                    # ------------------------------
                    
                    lim = li
                    li  = lip
                    lip = 0.25 * (((nodes[i+1,j+1,k+1,1] - nodes[i+1,j+2,k+1,1])^2 
                                 + (nodes[i+1,j+1,k+1,2] - nodes[i+1,j+2,k+1,2])^2 
                                 + (nodes[i+1,j+1,k+1,3] - nodes[i+1,j+2,k+1,3])^2)^0.5 
                                 +((nodes[  i,j+1,k+1,1] - nodes[  i,j+2,k+1,1])^2 
                                 + (nodes[  i,j+1,k+1,2] - nodes[  i,j+2,k+1,2])^2 
                                 + (nodes[  i,j+1,k+1,3] - nodes[  i,j+2,k+1,3])^2)^0.5
                                 +((nodes[i+1,j+1,  k,1] - nodes[i+1,j+2,  k,1])^2 
                                 + (nodes[i+1,j+1,  k,2] - nodes[i+1,j+2,  k,2])^2 
                                 + (nodes[i+1,j+1,  k,3] - nodes[i+1,j+2,  k,3])^2)^0.5 
                                 +((nodes[  i,j+1,  k,1] - nodes[  i,j+2,  k,1])^2 
                                 + (nodes[  i,j+1,  k,2] - nodes[  i,j+2,  k,2])^2 
                                 + (nodes[  i,j+1,  k,3] - nodes[  i,j+2,  k,3])^2)^0.5)

                    dip = 2*(Qbase[i,j+1,k,l] - Qbase[i,j,k,l])   / (li  + lip)
                    dim = 2*(Qbase[i,j,k,l]   - Qbase[i,j-1,k,l]) / (lim + li)
                    QbaseU[i,j,k,l] = Qbase[i,j,k,l]  - 0.5*li*minmod(dip, dim) 
                end
            end
        end
    end

    for l in 1:nval
        for k in 1+icell:cellzmax-icell+1
            for j in 1+icell-1:cellymax-icell + 1
                for i in 1+icell-1:cellxmax-icell+1
                    # ------------------------------
                    # B
                    # ------------------------------

                    lim = 0.25 * (((nodes[i+1,j+1,k-2,1] - nodes[i+1,j+1,k-1,1])^2 
                                 + (nodes[i+1,j+1,k-2,2] - nodes[i+1,j+1,k-1,2])^2 
                                 + (nodes[i+1,j+1,k-2,3] - nodes[i+1,j+1,k-1,3])^2)^0.5 
                                 +((nodes[  i,j+1,k-2,1] - nodes[  i,j+1,k-1,1])^2 
                                 + (nodes[  i,j+1,k-2,2] - nodes[  i,j+1,k-1,2])^2 
                                 + (nodes[  i,j+1,k-2,3] - nodes[  i,j+1,k-1,3])^2)^0.5
                                 +((nodes[i+1,  j,k-2,1] - nodes[i+1,  j,k-1,1])^2 
                                 + (nodes[i+1,  j,k-2,2] - nodes[i+1,  j,k-1,2])^2 
                                 + (nodes[i+1,  j,k-2,3] - nodes[i+1,  j,k-1,3])^2)^0.5 
                                 +((nodes[  i,  j,k-2,1] - nodes[  i,  j,k-1,1])^2 
                                 + (nodes[  i,  j,k-2,2] - nodes[  i,  j,k-1,2])^2 
                                 + (nodes[  i,  j,k-2,3] - nodes[  i,  j,k-1,3])^2)^0.5)
                    
                    li  = 0.25 * (((nodes[i+1,j+1,k-1,1] - nodes[i+1,j+1,  k,1])^2 
                                 + (nodes[i+1,j+1,k-1,2] - nodes[i+1,j+1,  k,2])^2 
                                 + (nodes[i+1,j+1,k-1,3] - nodes[i+1,j+1,  k,3])^2)^0.5 
                                 +((nodes[  i,j+1,k-1,1] - nodes[  i,j+1,  k,1])^2 
                                 + (nodes[  i,j+1,k-1,2] - nodes[  i,j+1,  k,2])^2 
                                 + (nodes[  i,j+1,k-1,3] - nodes[  i,j+1,  k,3])^2)^0.5
                                 +((nodes[i+1,  j,k-1,1] - nodes[i+1,  j,  k,1])^2 
                                 + (nodes[i+1,  j,k-1,2] - nodes[i+1,  j,  k,2])^2 
                                 + (nodes[i+1,  j,k-1,3] - nodes[i+1,  j,  k,3])^2)^0.5 
                                 +((nodes[  i,  j,k-1,1] - nodes[  i,  j,  k,1])^2 
                                 + (nodes[  i,  j,k-1,2] - nodes[  i,  j,  k,2])^2 
                                 + (nodes[  i,  j,k-1,3] - nodes[  i,  j,  k,3])^2)^0.5)
                    
                    lip = 0.25 * (((nodes[i+1,j+1,  k,1] - nodes[i+1,j+1,k+1,1])^2 
                                 + (nodes[i+1,j+1,  k,2] - nodes[i+1,j+1,k+1,2])^2 
                                 + (nodes[i+1,j+1,  k,3] - nodes[i+1,j+1,k+1,3])^2)^0.5 
                                 +((nodes[  i,j+1,  k,1] - nodes[  i,j+1,k+1,1])^2 
                                 + (nodes[  i,j+1,  k,2] - nodes[  i,j+1,k+1,2])^2 
                                 + (nodes[  i,j+1,  k,3] - nodes[  i,j+1,k+1,3])^2)^0.5
                                 +((nodes[i+1,  j,  k,1] - nodes[i+1,  j,k+1,1])^2 
                                 + (nodes[i+1,  j,  k,2] - nodes[i+1,  j,k+1,2])^2 
                                 + (nodes[i+1,  j,  k,3] - nodes[i+1,  j,k+1,3])^2)^0.5 
                                 +((nodes[  i,  j,  k,1] - nodes[  i,  j,k+1,1])^2 
                                 + (nodes[  i,  j,  k,2] - nodes[  i,  j,k+1,2])^2 
                                 + (nodes[  i,  j,  k,3] - nodes[  i,  j,k+1,3])^2)^0.5)
                    
                    
                    dip = 2*(Qbase[i,j,k,l]   - Qbase[i,j,k-1,l]) / (li  + lip)
                    dim = 2*(Qbase[i,j,k-1,l] - Qbase[i,j,k-2,l]) / (lim + li)
                    QbaseB[i,j,k,l] = Qbase[i,j,k-1,l] + 0.5*li*minmod(dim, dip) 

                    # ------------------------------
                    # F
                    # ------------------------------
                    
                    lim = li
                    li  = lip
                    lip = 0.25 * (((nodes[i+1,j+1,k+1,1] - nodes[i+1,j+1,k+2,1])^2 
                                 + (nodes[i+1,j+1,k+1,2] - nodes[i+1,j+1,k+2,2])^2 
                                 + (nodes[i+1,j+1,k+1,3] - nodes[i+1,j+1,k+2,3])^2)^0.5 
                                 +((nodes[  i,j+1,k+1,1] - nodes[  i,j+1,k+2,1])^2 
                                 + (nodes[  i,j+1,k+1,2] - nodes[  i,j+1,k+2,2])^2 
                                 + (nodes[  i,j+1,k+1,3] - nodes[  i,j+1,k+2,3])^2)^0.5
                                 +((nodes[i+1,  j,k+1,1] - nodes[i+1,  j,k+2,1])^2 
                                 + (nodes[i+1,  j,k+1,2] - nodes[i+1,  j,k+2,2])^2 
                                 + (nodes[i+1,  j,k+1,3] - nodes[i+1,  j,k+2,3])^2)^0.5 
                                 +((nodes[  i,  j,k+1,1] - nodes[  i,  j,k+2,1])^2 
                                 + (nodes[  i,  j,k+1,2] - nodes[  i,  j,k+2,2])^2 
                                 + (nodes[  i,  j,k+1,3] - nodes[  i,  j,k+2,3])^2)^0.5)

                    dip = 2*(Qbase[i,j,k+1,l] - Qbase[i,j,k,l])   / (li  + lip)
                    dim = 2*(Qbase[i,j,k,l]   - Qbase[i,j,k-1,l]) / (lim + li)
                    QbaseF[i,j,k,l] = Qbase[i,j,k,l]  - 0.5*li*minmod(dip, dim) 
                end
            end
        end
    end

    QconU = base_to_conservative(QbaseU, QconU, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    QconD = base_to_conservative(QbaseD, QconD, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    QconL = base_to_conservative(QbaseL, QconL, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    QconR = base_to_conservative(QbaseR, QconR, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    QconF = base_to_conservative(QbaseF, QconF, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    QconB = base_to_conservative(QbaseB, QconB, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    
    #=
    println("")
    println(Qbase[49,50,2])
    println(QbaseL[50,50,2])
    println(Qbase[50,50,2])
    =#
    
    
    return QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB, QconU, QconD, QconL, QconR, QconF, QconB
end

function minmod(a,b)
    ans = sign(a) * max(0, min(a*sign(b), b*sign(a)))
    return ans
end