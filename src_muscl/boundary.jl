# ------------------------------------
# set boundary conditions
# ------------------------------------
function set_boundary(Qbase, cellxmax, cellymax, vecAx, vecAy, bdcon, Rd, g, nval, icell)
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
    
    # -------------------------------------
    # bd1 = x-
    # -------------------------------------
    ii = icell

    if Int(bdcon[1][1]) == 0
        for j in 1:cellymax
            for l in 1:nval
                Qbase[ii,j,l] = bdcon[1][l+1]
            end
        end
    elseif Int(bdcon[1][1]) == 1
        for j in 1:cellymax
            for l in 1:nval
                Qbase[ii,j,l] = Qbase[ii+1,j,l]
            end
        end
    elseif Int(bdcon[1][1]) == 2
        for j in 1:cellymax
            for l in 1:nval
                Qbase[ii,j,l] = Qbase[ii+1,j,l]
            end

            xvec = vecAx[ii+1,j,1]
            yvec = vecAx[ii+1,j,2]
            u = Qbase[ii+1,j,2]
            v = Qbase[ii+1,j,3]

            Qbase[ii,j,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[ii,j,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)            
        end
    elseif Int(bdcon[1][1]) == 3
        for j in 1:cellymax
            for l in 1:nval
                Qbase[ii,j,l] = Qbase[ii+1,j,l]
            end

            u = Qbase[ii+1,j,2]
            v = Qbase[ii+1,j,3]

            Qbase[ii,j,2] = -u
            Qbase[ii,j,3] = -v
        end
    elseif Int(bdcon[1][1]) == 4
        for j in 1:cellymax
            for l in 1:nval
                for i in 1:icell
                    ii = cellxmax - icell - (icell - i)
                    Qbase[i,j,l] = Qbase[ii,j,l]
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
            for j in 1:cellymax
                for i in 1:icell-1
                    Qbase[i,j,l] = Qbase[i+1,j,l]
                end
            end 
        end
    end
    

    # -------------------------------------
    # bd2 = x+
    # -------------------------------------
    ii = cellxmax - (icell-1)

    if Int(bdcon[2][1]) == 0
        for j in 1:cellymax
            for l in 1:nval
                Qbase[ii,j,l] = bdcon[2][l+1]
            end
        end
    elseif Int(bdcon[2][1]) == 1
        for j in 1:cellymax
            for l in 1:nval
                Qbase[ii,j,l] = Qbase[ii-1,j,l]
            end
        end
    elseif Int(bdcon[2][1]) == 2
        for j in 1:cellymax
            for l in 1:nval
                Qbase[ii,j,l] = Qbase[ii-1,j,l]
            end

            xvec = vecAx[cellymax,j,1]
            yvec = vecAx[cellymax,j,2]
            u = Qbase[ii-1,j,2]
            v = Qbase[ii-1,j,3]

            Qbase[ii,j,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[ii,j,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
        end
    elseif Int(bdcon[2][1]) == 3
        for j in 1:cellymax
            for l in 1:nval
                Qbase[ii,j,l] = Qbase[ii-1,j,l]
            end
            u = Qbase[ii-1,j,2]
            v = Qbase[ii-1,j,3]

            Qbase[ii,j,2] = -u
            Qbase[ii,j,3] = -v
        end
    elseif Int(bdcon[2][1]) == 4
        for j in 1:cellymax
            for l in 1:nval
                for i in 1:icell
                    ii = cellxmax - (i-1)
                    iii = icell + (icell-i+1)
                    Qbase[ii,j,l] = Qbase[iii,j,l]
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
            for j in 1:cellymax
                for i in 1:icell-1
                    ii = cellxmax - (icell-1) +i
                    Qbase[ii,j,l] = Qbase[ii-1,j,l]
                end
            end 
        end
    end
    
    # -------------------------------------
    # bd3 = y-
    # -------------------------------------
    jj = icell

    if Int(bdcon[3][1]) == 0
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = bdcon[3][l+1]
            end
        end
    elseif Int(bdcon[3][1]) == 1
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj+1,l]
            end
        end
    elseif Int(bdcon[3][1]) == 2
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj+1,l]
            end

            xvec = vecAy[i,jj+1,1]
            yvec = vecAy[i,jj+1,2]
            u = Qbase[i,jj+1,2]
            v = Qbase[i,jj+1,3]
            
            Qbase[i,jj,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,jj,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
        end
    elseif Int(bdcon[3][1]) == 3
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj+1,l]
            end

            u = Qbase[i,jj+1,2]
            v = Qbase[i,jj+1,3]

            Qbase[i,jj,2] = -u
            Qbase[i,jj,3] = -v
        end
    elseif Int(bdcon[3][1]) == 4
        for i in 1:cellxmax
            for l in 1:nval
                for j in 1:icell
                    jj = cellymax - icell - (icell - j)
                    Qbase[i,j,l] = Qbase[i,jj,l]
                end
            end 
        end
    elseif Int(bdcon[3][1]) == 7
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj+1,l]
            end

            u = Qbase[i,jj+1,2]
            v = Qbase[i,jj+1,3]

            Qbase[i,jj,2] = -u
            Qbase[i,jj,3] = -v

            Tw = bdcon[3][nval+2]
            p  = Qbase[i,jj,4]
            Qbase[i,jj,1] = p/(Tw*Rd)
        end
    else
        println("------------------------")
        println(" boundary error ")
        println(" please check boundary number ")
        println("------------------------")
        throw(UndefVarError(:x))
    end

    if Int(bdcon[3][1]) != 4
        for l in 1:nval
            for i in 1:cellxmax
                for j in 1:icell-1
                    Qbase[i,j,l] = Qbase[i,j+1,l]
                end
            end 
        end
    end
    
    # -------------------------------------
    # bd4 = y+
    # -------------------------------------
    jj = cellymax - (icell-1)

    if Int(bdcon[4][1]) == 0
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = bdcon[4][l+1]
            end
        end
    elseif Int(bdcon[4][1]) == 1
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj-1,l]
            end
        end
    elseif Int(bdcon[4][1]) == 2
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj-1,l]
            end

            xvec = vecAy[i,jj,1]
            yvec = vecAy[i,jj,2]
            u = Qbase[i,jj-1,2]
            v = Qbase[i,jj-1,3]

            Qbase[i,jj,2] = ((-xvec^2+yvec^2)*u-2*xvec*yvec*v)/(xvec^2+yvec^2)
            Qbase[i,jj,3] = (-2*xvec*yvec*u+(xvec^2-yvec^2)*v)/(xvec^2+yvec^2)
        end
    elseif Int(bdcon[4][1]) == 3
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj-1,l]
            end

            u = Qbase[i,jj-1,2]
            v = Qbase[i,jj-1,3]

            Qbase[i,jj,2] = -u
            Qbase[i,jj,3] = -v
        end
    elseif Int(bdcon[4][1]) == 4
        for i in 1:cellxmax
            for l in 1:nval
                for j in 1:icell
                    jj  = cellymax - (j-1)
                    jjj = icell + (icell-j+1)
                    Qbase[i,jj,l] = Qbase[i,jjj,l]
                end
            end 
        end
    elseif Int(bdcon[4][1]) == 5
        temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
        temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
        for i in 1:temp1
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj-1,l]
            end
        end
        for i in temp1+1:temp2
            for l in 1:nval
                Qbase[i,jj,l] = bdcon[4][l+1]
            end
        end
        for i in temp2+1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj-1,l]
            end
        end
    elseif Int(bdcon[4][1]) == 6
        temp1 = Int64(round(cellxmax/4)+1) # 1/4流出
        temp2 = Int64((temp1-1)*3+1)       # 1/4~3/4流入
        for i in 1:temp1
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj-1,l]
            end
        end
        for i in temp1+1:temp2
            for l in 1:nval
                Qbase[i,jj,l] = bdcon[4][l+1]
            end
            T = bdcon[4][nval+2]
            p = (Qbase[i,jj,1]*Rd) * T
            Qbase[i,jj,4] = p
        end
        for i in temp2+1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = Qbase[i,jj-1,l]
            end
        end
    elseif Int(bdcon[4][1]) == 88
        for i in 1:cellxmax
            for l in 1:nval
                Qbase[i,jj,l] = bdcon[4][l+1]
            end
            T = bdcon[4][nval+2]
            p = (Qbase[i,jj,1]*Rd) * T
            Qbase[i,jj,4] = p
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
            for i in 1:cellxmax
                for j in 1:icell-1
                    jj = cellymax - (icell-1)+j
                    Qbase[i,jj,l] = Qbase[i,jj-1,l]
                end
            end 
        end
    end

    return Qbase
end

function check_bd(bdcon)
    for l in 1:4
        if Int(bdcon[l][1]) == 0
        elseif Int(bdcon[l][1]) == 1
        elseif Int(bdcon[l][1]) == 2
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