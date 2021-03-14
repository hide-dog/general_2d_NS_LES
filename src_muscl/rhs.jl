# ------------------------------------
# calculate right hand side
# ------------------------------------
function setup_RHS(RHS, cellxmax, cellymax, E_adv_hat, F_adv_hat, E_vis_hat, F_vis_hat, nval, volume, icell)
    Threads.@threads for l in 1:nval
        for j in 1+icell:cellymax-icell
            for i in 1+icell:cellxmax-icell
                RHS[i,j,l] = - (E_adv_hat[i+1,  j, l] - E_vis_hat[i+1,  j, l]) +
                             (E_adv_hat[  i,  j, l] - E_vis_hat[  i,  j, l]) -
                             (F_adv_hat[  i,j+1, l] - F_vis_hat[  i,j+1, l]) +
                             (F_adv_hat[  i,  j, l] - F_vis_hat[  i,  j, l])
            end
        end
    end

    #=
    println("RHS[20,50,:]")
    i = 50
    j = 50
    println(RHS[i,j,2])
    println((E_adv_hat[i+1,j,2] - E_adv_hat[i,j,2])/RHS[i,j,2])
    println((E_vis_hat[i+1,j,2] - E_vis_hat[i,j,2])/RHS[i,j,2])
    println((F_adv_hat[i,j+1,2] - F_adv_hat[i,j,2])/RHS[i,j,2])
    println((F_vis_hat[i,j+1,2] - F_vis_hat[i,j,2])/RHS[i,j,2])
    =#

    return RHS
end