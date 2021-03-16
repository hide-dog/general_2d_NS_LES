# ------------------------------------
# calculate Residuals
# ------------------------------------
function set_res(res, Delta_Qcon_hat, Delta_Qcon_hat_temp, cellxmax, cellymax, cellzmax, nval, icell)
    for l in 1:nval
        for k in 1+icell:cellzmax-icell
            for j in 1+icell:cellymax-icell
                for i in 1+icell:cellxmax-icell
                    res[i,j,k,l] = abs(Delta_Qcon_hat[i,j,k,l]-Delta_Qcon_hat_temp[i,j,k,l])
                end
            end
        end
    end 
    return res
end

# ------------------------------------
# calculate Residuals by norm-2
# ------------------------------------
function check_converge(res, RHS, cellxmax, cellymax, cellzmax, init_small, nval, icell)
    norm2 = zeros(nval)

    tempAxb = zeros(nval)
    tempb = zeros(nval)
    for l in 1:nval
        for k in 1+icell:cellzmax-icell
            for j in 1+icell:cellymax-icell
                for i in 1+icell:cellxmax-icell            
                    tempAxb[l] = tempAxb[l] + res[i,j,k,l]^2
                    tempb[l] = tempb[l] + RHS[i,j,k,l]^2
                end
            end
        end
    end
    for l in 1:nval
        norm2[l] = (tempAxb[l]/(tempb[l]+init_small)) ^ 0.5
    end
    return norm2
end