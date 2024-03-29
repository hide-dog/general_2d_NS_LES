# ------------------------------------
# set volume
# ------------------------------------
function set_volume(nodes, cellxmax, cellymax, cellzmax, vecAx, vecAy, vecAz, volume)
    for k in 1:cellzmax
        for j in 1:cellymax
            for i in 1:cellxmax
                #----------------------------
                #----------------------------
                #----------------------------
                rx = nodes[i+1,j+1,k+1,1] - nodes[i,j,k,1]
                ry = nodes[i+1,j+1,k+1,2] - nodes[i,j,k,2]
                rz = nodes[i+1,j+1,k+1,3] - nodes[i,j,k,3]

                Ax = 0.5 * (vecAx[i,j,k,1] + vecAx[i+1,j,k,1]) +
                     0.5 * (vecAy[i,j,k,1] + vecAy[i,j+1,k,1]) +
                     0.5 * (vecAz[i,j,k,1] + vecAz[i,j,k+1,1])
                Ay = 0.5 * (vecAx[i,j,k,2] + vecAx[i+1,j,k,2]) +
                     0.5 * (vecAy[i,j,k,2] + vecAy[i,j+1,k,2]) +
                     0.5 * (vecAz[i,j,k,2] + vecAz[i,j,k+1,2])
                Az = 0.5 * (vecAx[i,j,k,3] + vecAx[i+1,j,k,3]) +
                     0.5 * (vecAy[i,j,k,3] + vecAy[i,j+1,k,3]) +
                     0.5 * (vecAz[i,j,k,3] + vecAz[i,j,k+1,3])

                volume[i,j,k] = (Ax*rx + Ay*ry + Az*rz) /3
            end
        end
    end
    return volume
end 

# ------------------------------------
# set cellcenter
# ------------------------------------
function set_cellcenter(cellcenter, nodes, cellxmax, cellymax, cellzmax)
    for l in 1:3
        for k in 1:cellzmax
            for j in 1:cellymax
                for i in 1:cellxmax
                    xl = 0.125*(nodes[i+1,j+1,  k,l] +
                                nodes[i+1,  j,  k,l] + 
                                nodes[  i,j+1,  k,l] + 
                                nodes[  i,  j,  k,l] +
                                nodes[i+1,j+1,k+1,l] + 
                                nodes[i+1,  j,k+1,l] +
                                nodes[  i,j+1,k+1,l] +
                                nodes[  i,  j,k+1,l])
                    cellcenter[i,j,k,l] = xl
                end
            end
        end
    end
    return cellcenter
end

# ------------------------------------
# set dx and dy for local time stepping
# ------------------------------------
function set_dx_lts(dx, dy, dz, nodes, cellxmax, cellymax, cellzmax, icell)
    # define at cell boundaries
    # dx = zeros(cellxmax+1, cellymax)
    # dy = zeros(cellxmax, cellymax+1)

    point = zeros(3,3)
    
    for k in 1+icell:cellzmax -icell
        for j in 1+icell:cellymax -icell
            for i in 1+icell:cellxmax+1 -icell
                #= 点と直線の距離
                   ax+by+cz+d = 0
                   A(x0, y0, z0)
                   => D = abs(ax0 + by0 + cz0 + d) / abs(a^2 + b^2 + c^2)
                =#
                for l in 1:3
                    point[1,l] = nodes[i,  j,  k,l]
                    point[2,l] = nodes[i,j+1,  k,l]
                    point[3,l] = nodes[i,  j,k+1,l]
                end

                a,b,c,d = dis_point_plane(point)

                ccx = 0.125 * (nodes[i,  j,  k,1] + nodes[i+1,  j,  k,1] +
                               nodes[i,j+1,  k,1] + nodes[i+1,j+1,  k,1] +
                               nodes[i,  j,k+1,1] + nodes[i+1,  j,k+1,1] )
                ccy = 0.125 * (nodes[i,  j,  k,2] + nodes[i+1,  j,  k,2] +
                               nodes[i,j+1,  k,2] + nodes[i+1,j+1,  k,2] +
                               nodes[i,  j,k+1,2] + nodes[i+1,  j,k+1,2] )
                ccz = 0.125 * (nodes[i,  j,  k,3] + nodes[i+1,  j,  k,3] +
                               nodes[i,j+1,  k,3] + nodes[i+1,j+1,  k,3] +
                               nodes[i,  j,k+1,3] + nodes[i+1,  j,k+1,3] )

                dx1 = abs(a*ccx + b*ccy + c*ccz + d) / (a^2 + b^2 + c^2)^0.5

                ccx = 0.125 * (nodes[i,  j,  k,1] + nodes[i-1,  j,  k,1] +
                               nodes[i,j+1,  k,1] + nodes[i-1,j+1,  k,1] +
                               nodes[i,  j,k+1,1] + nodes[i-1,  j,k+1,1] )
                ccy = 0.125 * (nodes[i,  j,  k,2] + nodes[i-1,  j,  k,2] +
                               nodes[i,j+1,  k,2] + nodes[i-1,j+1,  k,2] +
                               nodes[i,  j,k+1,2] + nodes[i-1,  j,k+1,2] )
                ccz = 0.125 * (nodes[i,  j,  k,3] + nodes[i-1,  j,  k,3] +
                               nodes[i,j+1,  k,3] + nodes[i-1,j+1,  k,3] +
                               nodes[i,  j,k+1,3] + nodes[i-1,  j,k+1,3] )
                
                dx2 = abs(a*ccx + b*ccy + c*ccz + d) / (a^2 + b^2 + c^2)^0.5
                
                dx[i,j,k] = dx1 + dx2
            end
        end
    end

    for k in 1+icell:cellzmax -icell
        for j in 1+icell:cellymax+1 -icell
            for i in 1+icell:cellxmax -icell
                for l in 1:3
                    point[1,l] = nodes[  i,  j,  k,l]
                    point[2,l] = nodes[i+1,  j,  k,l]
                    point[3,l] = nodes[  i,  j,k+1,l]
                end

                a,b,c,d = dis_point_plane(point)

                ccx = 0.125 * (nodes[i,  j,  k,1] + nodes[  i,j+1,  k,1] +
                               nodes[i+1,j,  k,1] + nodes[i+1,j+1,  k,1] +
                               nodes[i,  j,k+1,1] + nodes[  i,j+1,k+1,1] )
                ccy = 0.125 * (nodes[i,  j,  k,2] + nodes[  i,j+1,  k,2] +
                               nodes[i+1,j,  k,2] + nodes[i+1,j+1,  k,2] +
                               nodes[i,  j,k+1,2] + nodes[  i,j+1,k+1,2] )
                ccz = 0.125 * (nodes[i,  j,  k,3] + nodes[  i,j+1,  k,3] +
                               nodes[i+1,j,  k,3] + nodes[i+1,j+1,  k,3] +
                               nodes[i,  j,k+1,3] + nodes[  i,j+1,k+1,3] )

                dx1 = abs(a*ccx + b*ccy + c*ccz + d) / (a^2 + b^2 + c^2)^0.5

                ccx = 0.125 * (nodes[i,  j,  k,1] + nodes[  i,j-1,  k,1] +
                               nodes[i+1,j,  k,1] + nodes[i+1,j-1,  k,1] +
                               nodes[i,  j,k+1,1] + nodes[  i,j-1,k+1,1] )
                ccy = 0.125 * (nodes[i,  j,  k,2] + nodes[  i,j-1,  k,2] +
                               nodes[i+1,j,  k,2] + nodes[i+1,j-1,  k,2] +
                               nodes[i,  j,k+1,2] + nodes[  i,j-1,k+1,2] )
                ccz = 0.125 * (nodes[i,  j,  k,3] + nodes[  i,j-1,  k,3] +
                               nodes[i+1,j,  k,3] + nodes[i+1,j-1,  k,3] +
                               nodes[i,  j,k+1,3] + nodes[  i,j-1,k+1,3] )
                
                dx2 = abs(a*ccx + b*ccy + c*ccz + d) / (a^2 + b^2 + c^2)^0.5
                
                dy[i,j,k] = dx1 + dx2
            end
        end
    end

    for k in 1+icell:cellzmax+1 -icell
        for j in 1+icell:cellymax -icell
            for i in 1+icell:cellxmax -icell
                for l in 1:3
                    point[1,l] = nodes[  i,  j,  k,l]
                    point[2,l] = nodes[i+1,  j,  k,l]
                    point[3,l] = nodes[  i,j+1,  k,l]
                end

                a,b,c,d = dis_point_plane(point)

                ccx = 0.125 * (nodes[i,    j,  k,1] + nodes[  i,  j,k+1,1] +
                               nodes[i+1,  j,  k,1] + nodes[i+1,  j,k+1,1] +
                               nodes[  i,j+1,  k,1] + nodes[  i,j+1,k+1,1] )
                ccy = 0.125 * (nodes[i,    j,  k,2] + nodes[  i,  j,k+1,2] +
                               nodes[i+1,  j,  k,2] + nodes[i+1,  j,k+1,2] +
                               nodes[  i,j+1,  k,2] + nodes[  i,j+1,k+1,2] )
                ccz = 0.125 * (nodes[i,    j,  k,3] + nodes[  i,  j,k+1,3] +
                               nodes[i+1,  j,  k,3] + nodes[i+1,  j,k+1,3] +
                               nodes[  i,j+1,  k,3] + nodes[  i,j+1,k+1,3] )

                dx1 = abs(a*ccx + b*ccy + c*ccz + d) / (a^2 + b^2 + c^2)^0.5

                ccx = 0.125 * (nodes[i,    j,  k,1] + nodes[  i,  j,k-1,1] +
                               nodes[i+1,  j,  k,1] + nodes[i+1,  j,k-1,1] +
                               nodes[  i,j+1,  k,1] + nodes[  i,j+1,k-1,1] )
                ccy = 0.125 * (nodes[i,    j,  k,2] + nodes[  i,  j,k-1,2] +
                               nodes[i+1,  j,  k,2] + nodes[i+1,  j,k-1,2] +
                               nodes[  i,j+1,  k,2] + nodes[  i,j+1,k-1,2] )
                ccz = 0.125 * (nodes[i,    j,  k,3] + nodes[  i,  j,k-1,3] +
                               nodes[i+1,  j,  k,3] + nodes[i+1,  j,k-1,3] +
                               nodes[  i,j+1,  k,3] + nodes[  i,j+1,k-1,3] )
                
                dx2 = abs(a*ccx + b*ccy + c*ccz + d) / (a^2 + b^2 + c^2)^0.5
                
                dz[i,j,k] = dx1 + dx2
            end
        end
    end

    return dx, dy, dz
end


# ------------------------------------
# Distance between point and plane
# ------------------------------------
function dis_point_plane(point)
    #= 点と直線の距離
        ax+by+cz+d = 0
        A(x0, y0, z0)
        => D = abs(ax0 + by0 + cz0 + d) / abs(a^2 + b^2 + c^2)

        平面に対する法線ベクトル:vec{n} = (p, q, r)
        p(x−x0) + q(y−y0) + r(z−z0) = 0

        法線ベクトルは外積から求める。
    =#
    a1 = point[2,1] - point[1,1]
    a2 = point[2,2] - point[1,2]
    a3 = point[2,3] - point[1,3]
    b1 = point[3,1] - point[1,1]
    b2 = point[3,2] - point[1,2]
    b3 = point[3,3] - point[1,3]

    p = a2*b3 - a3*b2
    q = a3*b1 - a1*b3
    r = a1*b2 - a2*b1

    a = p
    b = q
    c = r
    d = -p*point[1,1] - q*point[1,2] - r*point[1,3]

    return a, b, c, d
end

# ------------------------------------
# set dtau for local time stepping
# ------------------------------------
function set_lts(dtau, lambda_facex, lambda_facey, lambda_facez, Qbase, cellxmax, cellymax, cellzmax, mu, dx, dy, dz,
                vecAx, vecAy, vecAz, volume, specific_heat_ratio, cfl, icell)
    g = specific_heat_ratio
    
    for k in 1+icell:cellzmax -icell
        for j in 1+icell:cellymax -icell
            for i in 1+icell:cellxmax+1 -icell
                rho_av = 0.5 * (Qbase[i,j,k,1] + Qbase[i-1,j,k,1])
                u_av   = 0.5 * (Qbase[i,j,k,2] + Qbase[i-1,j,k,2])
                v_av   = 0.5 * (Qbase[i,j,k,3] + Qbase[i-1,j,k,3])
                w_av   = 0.5 * (Qbase[i,j,k,4] + Qbase[i-1,j,k,4])
                mu_av  = 0.5 * (   mu[i,j,k]   +    mu[i-1,j,k]  )
                
                ap = (g * Qbase[  i,j,k,5] / Qbase[  i,j,k,1])^0.5
                am = (g * Qbase[i-1,j,k,5] / Qbase[i-1,j,k,1])^0.5
                a_av = 0.5 * (ap + am)

                U   = u_av*vecAx[i,j,k,1] + v_av*vecAx[i,j,k,2] + w_av*vecAx[i,j,k,3]
                lambda_facex[i,j,k] = abs(U) + a_av + 2*mu_av/(rho_av*dx[i,j,k])
            end
        end
    end
    
    for k in 1+icell:cellzmax -icell
        for j in 1+icell:cellymax+1 -icell
            for i in 1+icell:cellxmax -icell
                rho_av = 0.5 * (Qbase[i,j,k,1] + Qbase[i,j-1,k,1])
                u_av   = 0.5 * (Qbase[i,j,k,2] + Qbase[i,j-1,k,2])
                v_av   = 0.5 * (Qbase[i,j,k,3] + Qbase[i,j-1,k,3])
                w_av   = 0.5 * (Qbase[i,j,k,4] + Qbase[i,j-1,k,4])
                mu_av  = 0.5 * (   mu[i,j,k]   +    mu[i,j-1,k]  )
                
                ap = (g * Qbase[  i,j,k,5] / Qbase[  i,j,k,1])^0.5
                am = (g * Qbase[i,j-1,k,5] / Qbase[i,j-1,k,1])^0.5
                a_av = 0.5 * (ap + am)

                V   = u_av*vecAy[i,j,k,1] + v_av*vecAy[i,j,k,2] + w_av*vecAy[i,j,k,3]
                lambda_facey[i,j,k] = abs(V) + a_av + 2*mu_av/(rho_av*dy[i,j,k])
            end
        end
    end

    for k in 1+icell:cellzmax+1 -icell
        for j in 1+icell:cellymax -icell
            for i in 1+icell:cellxmax -icell
                rho_av = 0.5 * (Qbase[i,j,k,1] + Qbase[i,j,k-1,1])
                u_av   = 0.5 * (Qbase[i,j,k,2] + Qbase[i,j,k-1,2])
                v_av   = 0.5 * (Qbase[i,j,k,3] + Qbase[i,j,k-1,3])
                w_av   = 0.5 * (Qbase[i,j,k,4] + Qbase[i,j,k-1,4])
                mu_av  = 0.5 * (   mu[i,j,k]   +    mu[i,j,k-1]  )
                
                ap = (g * Qbase[  i,j,k,5] / Qbase[  i,j,k,1])^0.5
                am = (g * Qbase[i,j,k-1,5] / Qbase[i,j,k-1,1])^0.5
                a_av = 0.5 * (ap + am)

                W   = u_av*vecAz[i,j,k,1] + v_av*vecAz[i,j,k,2] + w_av*vecAz[i,j,k,3]
                lambda_facez[i,j,k] = abs(W) + a_av + 2*mu_av/(rho_av*dz[i,j,k])
            end
        end
    end
    
    for k in 1+icell:cellzmax-icell
        for j in 1+icell:cellymax-icell
            for i in 1+icell:cellxmax-icell
                a1 = lambda_facex[  i,  j,  k]
                a2 = lambda_facex[i+1,  j,  k]
                a3 = lambda_facey[  i,  j,  k]
                a4 = lambda_facey[  i,j+1,  k]
                a5 = lambda_facez[  i,  j,  k]
                a6 = lambda_facez[  i,  j,k+1]
                lmax = maximum([a1, a2, a3, a4, a5, a6])

                dtau[i,j,k] = cfl * volume[i,j,k] / lmax
            end
        end
    end
    return dtau
end

# ------------------------------------
# set viscosity by Sutherland's formula
# https://cattech-lab.com/science-tools/sutherland/
# ------------------------------------
function set_mu(mu, Qbase, cellxmax, cellymax, cellzmax,specific_heat_ratio, Rd)
    mu0 = 1.82e-5     # Reference Viscosity, Pa s
    T0  = 293.15      # Reference Temperature, K
    C   = 117         # Sutherland's constant, K
    
    for k in 1:cellzmax
        for j in 1:cellymax
            for i in 1:cellxmax
                # T = p/(rho Rd )
                T = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd)
                mu[i,j,k] = mu0 * (T/T0)^1.5 * (T0+C)/(T+C)
            end
        end
    end
    return mu
end

# ------------------------------------
# set thermal Conductivity by Sutherland's formula
# https://doi.org/10.11357/jsam1937.37.694
# ------------------------------------
function set_lambda(lambda, Qbase, cellxmax, cellymax, cellzmax, mu, specific_heat_ratio, Rd)
    lam0 = 22.3*10^(-3)  # Reference thermal Conductivity, W/(m K)
    T0   = 273.15        # Reference Temperature, K
    C    = 125           # Sutherland's constant, K

    for k in 1:cellzmax
        for j in 1:cellymax
            for i in 1:cellxmax
                T = Qbase[i,j,k,5]/(Qbase[i,j,k,1]*Rd)
                lambda[i,j,k] = lam0*((T0+C)/(T+C))*(T/T0)^1.5
            end
        end
    end

    return lambda
end

# ------------------------------------
# Conversion from primitive variables to conserved variables
# ------------------------------------
function base_to_conservative(Qbase, Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """ 
    for k in 1:cellzmax
        for j in 1:cellymax
            for i in 1:cellxmax
                Qcon[i,j,k,1] = Qbase[i,j,k,1]
                Qcon[i,j,k,2] = Qbase[i,j,k,1]*Qbase[i,j,k,2]
                Qcon[i,j,k,3] = Qbase[i,j,k,1]*Qbase[i,j,k,3]
                Qcon[i,j,k,4] = Qbase[i,j,k,1]*Qbase[i,j,k,4]
                Qcon[i,j,k,5] = Qbase[i,j,k,5]/(specific_heat_ratio-1) + Qbase[i,j,k,1]*(Qbase[i,j,k,2]^2+Qbase[i,j,k,3]^2+Qbase[i,j,k,3]^2)/2
            end
        end
    end
    return Qcon
end

# ------------------------------------
# Conversion from conserved variables to primitive variables
# ------------------------------------
function conservative_to_base(Qbase, Qcon, cellxmax, cellymax, cellzmax, specific_heat_ratio)
    """
    Qbase=[rho,u,v,p]
    Qcon=[rho,rhou,rhov,e]
    """
    for k in 1:cellzmax
        for j in 1:cellymax
            for i in 1:cellxmax
                Qbase[i,j,k,1] = Qcon[i,j,k,1]
                Qbase[i,j,k,2] = Qcon[i,j,k,2]/Qcon[i,j,k,1]
                Qbase[i,j,k,3] = Qcon[i,j,k,3]/Qcon[i,j,k,1]
                Qbase[i,j,k,4] = Qcon[i,j,k,4]/Qcon[i,j,k,1]
                Qbase[i,j,k,5] = (Qcon[i,j,k,5]-Qcon[i,j,k,1]*(Qbase[i,j,k,2]^2+Qbase[i,j,k,3]^2+Qbase[i,j,k,4]^2)/2)*(specific_heat_ratio-1)
            end
        end
    end
    return Qbase
end

# ------------------------------------
# Conversion from Cartesian coordinate system to general coordinate system
# ------------------------------------
function setup_Qcon_hat(Qcon, Qcon_hat, cellxmax, cellymax, cellzmax, volume, nval)
    for l in 1:nval
        for k in 1:cellzmax
            for j in 1:cellymax
                for i in 1:cellxmax
                    Qcon_hat[i,j,k,l] = Qcon[i,j,k,l] * volume[i,j,k]
                end
            end
        end
    end
    return Qcon_hat
end

# ------------------------------------
# Conversion from general coordinate system to Cartesian coordinate system
# ------------------------------------
function Qhat_to_Q(Qcon, Qcon_hat, cellxmax, cellymax, cellzmax, volume, nval)
    for l in 1:nval
        for k in 1:cellzmax
            for j in 1:cellymax
                for i in 1:cellxmax
                    Qcon[i,j,k,l] = Qcon_hat[i,j,k,l] / volume[i,j,k]
                end
            end
        end
    end
    return Qcon
end

# ------------------------------------
# cal average convective value
# ------------------------------------
function cal_Qave(Qbase_ave, Qbase, cellxmax, cellymax, cellzmax, nval, stepnum)
    for l in 1:nval
        for k in 1:cellzmax
            for j in 1:cellymax
                for i in 1:cellxmax
                    temp = (Qbase_ave[i,j,k,l]*(stepnum-1) + Qbase[i,j,k,l]) / stepnum
                    Qbase_ave[i,j,k,l] = temp
                end
            end
        end
    end
    return Qbase_ave
end