# ------------------------------------
# set boundary conditions
# ------------------------------------
function set_boundary(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    """
    bdcon[i][j]
    i : bd number(x-, x+, y-, y+)
    j=1-6 : "bd1_con":"2",
            "bd1_rho":"1.0",
            "bd1_u"  :"300.0",
            "bd1_v"  :"0.0",
            "bd1_p"  :"1.0",
            "bd1_T"  :"300.0",
    """

    Qbase = bd_x_m(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    Qbase = bd_x_p(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    Qbase = bd_y_m(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    Qbase = bd_y_p(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    Qbase = bd_z_m(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    Qbase = bd_z_p(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    return Qbase
end

function bd_x_m(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    # -------------------------------------
    # bd1 = x-
    # -------------------------------------
    ii = icell

    if Int(bdcon[1][1]) == 0
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:cellymax
                    Qbase[ii,j,k,l] = bdcon[1][l+1]
                end
            end
        end
    elseif Int(bdcon[1][1]) == 1
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:cellymax
                    Qbase[ii,j,k,l] = Qbase[ii+1,j,k,l]
                end
            end
        end
    elseif Int(bdcon[1][1]) == 3
        for k in 1:cellzmax
            for j in 1:cellymax
                for l in 1:nval
                    Qbase[ii,j,k,l] = Qbase[ii+1,j,k,l]
                end

                u = Qbase[ii+1,j,k,2]
                v = Qbase[ii+1,j,k,3]
                w = Qbase[ii+1,j,k,4]

                Qbase[ii,j,k,2] = -u
                Qbase[ii,j,k,3] = -v
                Qbase[ii,j,k,4] = -w
            end
        end
    elseif Int(bdcon[1][1]) == 4
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:cellymax
                    for i in 1:icell
                        ii = cellxmax - icell - (icell - i)
                        Qbase[i,j,k,l] = Qbase[ii,j,k,l]
                    end
                end
            end
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    if Int(bdcon[1][1]) != 4
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:cellymax
                    for i in 1:icell-1
                        Qbase[i,j,k,l] = Qbase[i+1,j,k,l]
                    end
                end 
            end
        end
    end

    return Qbase
end

function bd_x_p(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    # -------------------------------------
    # bd2 = x+
    # -------------------------------------
    ii = cellxmax - (icell-1)

    if Int(bdcon[2][1]) == 0
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:cellymax
                    Qbase[ii,j,k,l] = bdcon[2][l+1]
                end
            end
        end
    elseif Int(bdcon[2][1]) == 1
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:cellymax
                    Qbase[ii,j,k,l] = Qbase[ii-1,j,k,l]
                end
            end
        end
    elseif Int(bdcon[2][1]) == 3
        for k in 1:cellzmax
            for j in 1:cellymax
                for l in 1:nval
                    Qbase[ii,j,k,l] = Qbase[ii-1,j,k,l]
                end
                u = Qbase[ii-1,j,k,2]
                v = Qbase[ii-1,j,k,3]
                w = Qbase[ii-1,j,k,4]

                Qbase[ii,j,k,2] = -u
                Qbase[ii,j,k,3] = -v
                Qbase[ii,j,k,4] = -w
            end
        end
    elseif Int(bdcon[2][1]) == 4
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:cellymax
                    for i in 1:icell
                        ii = cellxmax - (i-1)
                        iii = icell + (icell-i+1)
                        Qbase[ii,j,k,l] = Qbase[iii,j,k,l]
                    end
                end
            end
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    if Int(bdcon[2][1]) != 4
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:cellymax
                    for i in 1:icell-1
                        ii = cellxmax - (icell-1) +i
                        Qbase[ii,j,k,l] = Qbase[ii-1,j,k,l]
                    end
                end 
            end
        end
    end

    return Qbase
end

function bd_y_m(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    # -------------------------------------
    # bd3 = y-
    # -------------------------------------
    jj = icell

    if Int(bdcon[3][1]) == 0
        for l in 1:nval
            for k in 1:cellzmax
                for i in 1:cellxmax
                    Qbase[i,jj,k,l] = bdcon[3][l+1]
                end
            end
        end
    elseif Int(bdcon[3][1]) == 1
        for k in 1:cellzmax
            for l in 1:nval
                for i in 1:cellxmax
                    Qbase[i,jj,k,l] = Qbase[i,jj+1,k,l]
                end
            end
        end
    elseif Int(bdcon[3][1]) == 3
        for k in 1:cellzmax
            for i in 1:cellxmax
                for l in 1:nval
                    Qbase[i,jj,k,l] = Qbase[i,jj+1,k,l]
                end

                u = Qbase[i,jj+1,k,2]
                v = Qbase[i,jj+1,k,3]
                w = Qbase[i,jj+1,k,4]

                Qbase[i,jj,k,2] = -u
                Qbase[i,jj,k,3] = -v
                Qbase[i,jj,k,4] = -w
            end
        end
    elseif Int(bdcon[3][1]) == 4
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:icell
                    for i in 1:cellxmax
                        jj = cellymax - icell - (icell - j)
                        Qbase[i,j,k,l] = Qbase[i,jj,k,l]
                    end
                end 
            end
        end
    elseif Int(bdcon[3][1]) == 7
        for k in 1:cellzmax
            for i in 1:cellxmax
                for l in 1:nval
                    Qbase[i,jj,k,l] = Qbase[i,jj+1,k,l]
                end

                u = Qbase[i,jj+1,k,2]
                v = Qbase[i,jj+1,k,3]
                w = Qbase[i,jj+1,k,4]

                Qbase[i,jj,k,2] = -u
                Qbase[i,jj,k,3] = -v
                Qbase[i,jj,k,4] = -w

                Tw = bdcon[3][nval+2]
                p  = Qbase[i,jj,k,5]
                Qbase[i,jj,k,1] = p/(Tw*Rd)
            end
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    if Int(bdcon[3][1]) != 4
        for k in 1:cellzmax
            for l in 1:nval
                for j in 1:icell-1
                    for i in 1:cellxmax
                        Qbase[i,j,k,l] = Qbase[i,j+1,k,l]
                    end
                end 
            end
        end
    end

    return Qbase
end

function bd_y_p(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    # -------------------------------------
    # bd4 = y+
    # -------------------------------------
    jj = cellymax - (icell-1)

    if Int(bdcon[4][1]) == 0
        for l in 1:nval
            for k in 1:cellzmax
                for i in 1:cellxmax
                    Qbase[i,jj,k,l] = bdcon[4][l+1]
                end
            end
        end
    elseif Int(bdcon[4][1]) == 1
        for l in 1:nval
            for k in 1:cellzmax
                for i in 1:cellxmax 
                    Qbase[i,jj,k,l] = Qbase[i,jj-1,k,l]
                end
            end
        end
    elseif Int(bdcon[4][1]) == 3
        for k in 1:cellzmax
            for i in 1:cellxmax
                for l in 1:nval
                    Qbase[i,jj,k,l] = Qbase[i,jj-1,k,l]
                end

                u = Qbase[i,jj-1,k,2]
                v = Qbase[i,jj-1,k,3]
                w = Qbase[i,jj-1,k,4]

                Qbase[i,jj,k,2] = -u
                Qbase[i,jj,k,3] = -v
                Qbase[i,jj,k,4] = -w
            end
        end
    elseif Int(bdcon[4][1]) == 4
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:icell
                    for i in 1:cellxmax
                        jj  = cellymax - (j-1)
                        jjj = icell + (icell-j+1)
                        Qbase[i,jj,k,l] = Qbase[i,jjj,k,l]
                    end
                end
            end 
        end
    elseif Int(bdcon[4][1]) == 5
        temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
        temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
        for l in 1:nval
            for k in 1:cellzmax
                for i in 1:temp1
                    Qbase[i,jj,k,l] = Qbase[i,jj-1,k,l]
                end
            end
        end
        for l in 1:nval
            for k in 1:cellzmax
                for i in temp1+1:temp2
                    Qbase[i,jj,k,l] = bdcon[4][l+1]
                end
            end
        end
        for l in 1:nval
            for k in 1:cellzmax
                for i in temp2+1:cellxmax
                    Qbase[i,jj,k,l] = Qbase[i,jj-1,k,l]
                end
            end
        end
    elseif Int(bdcon[4][1]) == 6
        temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
        temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
        for l in 1:nval
            for k in 1:cellzmax
                for i in 1:temp1
                    Qbase[i,jj,k,l] = Qbase[i,jj-1,k,l]
                end
            end
        end
        for k in 1:cellzmax
            for i in temp1+1:temp2
                for l in 1:nval
                    Qbase[i,jj,k,l] = bdcon[4][l+1]
                end
                T = bdcon[4][nval+2]
                p = (Qbase[i,jj,k,1]*Rd) * T
                Qbase[i,jj,k,5] = p
            end
        end
        for l in 1:nval
            for k in 1:cellzmax
                for i in temp2+1:cellxmax
                    Qbase[i,jj,k,l] = Qbase[i,jj-1,k,l]
                end
            end
        end
    elseif Int(bdcon[4][1]) == 88
        for k in 1:cellzmax
            for i in 1:cellxmax
                for l in 1:nval
                    Qbase[i,jj,k,l] = bdcon[4][l+1]
                end
                T = bdcon[4][nval+2]
                p = (Qbase[i,jj,k,1]*Rd) * T
                Qbase[i,jj,k,5] = p
            end
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    if Int(bdcon[4][1]) != 4
        for l in 1:nval
            for k in 1:cellzmax
                for j in 1:icell-1
                    for i in 1:cellymax
                        jj = cellymax - (icell-1)+j
                        Qbase[i,jj,k,l] = Qbase[i,jj-1,k,l]
                    end
                end 
            end
        end
    end

    return Qbase
end

function bd_z_m(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    # -------------------------------------
    # bd5 = z-
    # -------------------------------------
    kk = icell

    if Int(bdcon[5][1]) == 0
        for l in 1:nval
            for j in 1:cellymax
                for i in 1:cellxmax
                    Qbase[i,j,kk,l] = bdcon[5][l+1]
                end
            end
        end
    elseif Int(bdcon[5][1]) == 1
        for l in 1:nval
            for j in 1:cellymax
                for i in 1:cellxmax
                    Qbase[i,j,kk,l] = Qbase[i,j,kk+1,l]
                end
            end
        end
    elseif Int(bdcon[5][1]) == 3
        for j in 1:cellymax
            for i in 1:cellxmax
                for l in 1:nval
                    Qbase[i,j,kk,l] = Qbase[i,j,kk+1,l]
                end

                u = Qbase[i,j,kk+1,2]
                v = Qbase[i,j,kk+1,3]
                w = Qbase[i,j,kk+1,4]

                Qbase[i,j,kk,2] = -u
                Qbase[i,j,kk,3] = -v
                Qbase[i,j,kk,4] = -w
            end
        end
    elseif Int(bdcon[5][1]) == 4
        for l in 1:nval
            for k in 1:icell
                for j in 1:cellymax
                    for i in 1:cellxmax
                        kk = cellzmax - icell - (icell - k)
                        Qbase[i,j,k,l] = Qbase[i,j,kk,l]
                    end
                end
            end
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    if Int(bdcon[1][1]) != 4
        for l in 1:nval
            for k in 1:icell-1
                for j in 1:cellymax
                    for i in 1:cellxmax
                        Qbase[i,j,k,l] = Qbase[i,j,k+1,l]
                    end
                end 
            end
        end
    end

    return Qbase
end

function bd_z_p(Qbase, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, bdcon, Rd, g, nval, icell)
    # -------------------------------------
    # bd6 = z+
    # -------------------------------------
    kk = cellxmax - (icell-1)

    if Int(bdcon[6][1]) == 0
        for l in 1:nval
            for j in 1:cellymax
                for i in 1:cellxmax
                    Qbase[i,j,kk,l] = bdcon[6][l+1]
                end
            end
        end
    elseif Int(bdcon[6][1]) == 1
        for l in 1:nval
            for j in 1:cellymax
                for i in 1:cellxmax
                    Qbase[i,j,kk,l] = Qbase[i,j,kk-1,l]
                end
            end
        end
    elseif Int(bdcon[6][1]) == 3
        for j in 1:cellymax
            for i in 1:cellxmax
                for l in 1:nval
                    Qbase[i,j,kk,l] = Qbase[i,j,kk-1,l]
                end
                u = Qbase[i,j,kk-1,2]
                v = Qbase[i,j,kk-1,3]
                w = Qbase[i,j,kk-1,4]

                Qbase[i,j,kk,2] = -u
                Qbase[i,j,kk,3] = -v
                Qbase[i,j,kk,4] = -w
            end
        end
    elseif Int(bdcon[6][1]) == 4
        for l in 1:nval
            for k in 1:icell
                for j in 1:cellymax
                    for i in 1:cellxmax
                        kk = cellzmax - (k-1)
                        kkk = icell + (icell-k+1)
                        Qbase[i,j,kk,l] = Qbase[i,j,kkk,l]
                    end
                end
            end
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    if Int(bdcon[6][1]) != 4
        for l in 1:nval
            for k in 1:icell-1
                for j in 1:cellymax
                    for i in 1:cellxmax
                        kk = cellzmax - (icell-1) +k
                        Qbase[i,j,kk,l] = Qbase[i,j,kk-1,l]
                    end
                end 
            end
        end
    end

    return Qbase
end

function check_bd(bdcon)
    for l in 1:6
        if Int(bdcon[l][1]) == 0
        elseif Int(bdcon[l][1]) == 1
        #elseif Int(bdcon[l][1]) == 2
        elseif Int(bdcon[l][1]) == 3
        elseif Int(bdcon[l][1]) == 4
        elseif Int(bdcon[l][1]) == 5
        elseif Int(bdcon[l][1]) == 6
        elseif Int(bdcon[l][1]) == 7
        elseif Int(bdcon[l][1]) == 88
        else
            println("\n check boundary condition ! \n")
            throw(UndefVarError(:x))
        end
    end
end