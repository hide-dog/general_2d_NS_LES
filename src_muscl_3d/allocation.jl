# ------------------------------------
# Secure the array for explicit and implicit
# ------------------------------------
function common_allocation(cellxmax, cellymax, cellzmax, nval)
    # define at cell center
    Qbase     = zeros(cellxmax, cellymax, cellzmax, nval)        # primitive variables
    Qbase_ave = zeros(cellxmax, cellymax, cellzmax, nval)        # average primitive variables 
    volume    = zeros(cellxmax, cellymax, cellzmax)              # volume

    cellcenter = zeros(cellxmax, cellymax, cellzmax, 3)          # cell center
    
    wally = zeros(cellxmax, cellymax, cellzmax)                  # dis from wall
    yplus =  ones(cellxmax, cellymax, cellzmax) * 1e5            # yplus

    Qcon     = zeros(cellxmax, cellymax, cellzmax, nval)         # conserved variables
    Qcon_hat = zeros(cellxmax, cellymax, cellzmax, nval)         # conserved variables in the general coordinate system
    
    QbaseU = zeros(cellxmax+1, cellymax, cellzmax, nval)        # Q up
    QbaseD = zeros(cellxmax+1, cellymax, cellzmax, nval)        # Q down
    QbaseL = zeros(cellxmax, cellymax+1, cellzmax, nval)        # Q left
    QbaseR = zeros(cellxmax, cellymax+1, cellzmax, nval)        # Q right
    QbaseF = zeros(cellxmax, cellymax, cellzmax+1, nval)        # Q front
    QbaseB = zeros(cellxmax, cellymax, cellzmax+1, nval)        # Q back

    QconU = zeros(cellxmax+1, cellymax, cellzmax, nval)        # Q up
    QconD = zeros(cellxmax+1, cellymax, cellzmax, nval)        # Q down
    QconL = zeros(cellxmax, cellymax+1, cellzmax, nval)        # Q left
    QconR = zeros(cellxmax, cellymax+1, cellzmax, nval)        # Q right
    QconF = zeros(cellxmax, cellymax, cellzmax+1, nval)        # Q front
    QconB = zeros(cellxmax, cellymax, cellzmax+1, nval)        # Q back

    mu     = zeros(cellxmax, cellymax, cellzmax)                 # viscosity
    lambda = zeros(cellxmax, cellymax, cellzmax)                 # thermal Conductivity
    mut    = zeros(cellxmax, cellymax, cellzmax)                 # turb viscosity
    mut_bd = zeros(cellxmax+1, cellymax+1, cellzmax+1, 2)        # turb viscosity on bd

    RHS = zeros(cellxmax, cellymax, cellzmax, nval)              # right hand side
    
    # define at cell boundaries
    dx = zeros(cellxmax+1, cellymax, cellzmax)                   # distance from cell center
    dy = zeros(cellxmax, cellymax+1, cellzmax)                   # distance from cell center
    dz = zeros(cellxmax, cellymax, cellzmax+1)                   # distance from cell center

    E_adv_hat = zeros(cellxmax+1,   cellymax,   cellzmax, nval)    # flux of advection in the x-direction
    F_adv_hat = zeros(  cellxmax, cellymax+1,   cellzmax, nval)    # flux of advection in the y-direction
    G_adv_hat = zeros(  cellxmax,   cellymax, cellzmax+1, nval)    # flux of advection in the z-direction

    E_vis_hat = zeros(cellxmax+1,   cellymax,   cellzmax, nval)    # flux of viscosity in the x-direction
    F_vis_hat = zeros(  cellxmax, cellymax+1,   cellzmax, nval)    # flux of viscosity in the y-direction
    G_vis_hat = zeros(  cellxmax,   cellymax, cellzmax+1, nval)    # flux of viscosity in the z-direction

    return Qbase, Qbase_ave, volume, cellcenter, wally, yplus, dx, dy, dz, Qcon, Qcon_hat, mu, mut, mut_bd, lambda, 
            E_adv_hat, F_adv_hat, G_adv_hat, E_vis_hat, F_vis_hat, G_vis_hat, RHS, QbaseU, QbaseD, QbaseL, QbaseR, QbaseF, QbaseB,
            QconU, QconD, QconL, QconR, QconF, QconB
            
end 

function allocation_implicit(cellxmax, cellymax, cellzmax, nval)
    # define at cell center
    Qbasen     = zeros(cellxmax, cellymax, cellzmax, nval)        # primitive variables for inner iteration
    Qconn      = zeros(cellxmax, cellymax, cellzmax, nval)        # conserved variables for inner iteration
    Qconn_hat  = zeros(cellxmax, cellymax, cellzmax, nval)        # conserved variables in general coordinate system for inner iteration
    Qbasem     = zeros(cellxmax, cellymax, cellzmax, nval)        # primitive variables for inner iteration

    dtau         = zeros(cellxmax, cellymax, cellzmax)            # computational time steps

    A_adv_hat_p = zeros(cellxmax, cellymax, cellzmax, nval, nval) # Jacobian matrix A+ for one-wave approximation
    A_adv_hat_m = zeros(cellxmax, cellymax, cellzmax, nval, nval) # Jacobian matrix A- for one-wave approximation
    B_adv_hat_p = zeros(cellxmax, cellymax, cellzmax, nval, nval) # Jacobian matrix B+ for one-wave approximation
    B_adv_hat_m = zeros(cellxmax, cellymax, cellzmax, nval, nval) # Jacobian matrix B- for one-wave approximation
    C_adv_hat_p = zeros(cellxmax, cellymax, cellzmax, nval, nval) # Jacobian matrix C+ for one-wave approximation
    C_adv_hat_m = zeros(cellxmax, cellymax, cellzmax, nval, nval) # Jacobian matrix C- for one-wave approximation
    A_beta_shig = zeros(cellxmax, cellymax, cellzmax)             # A sigma for one-wave approximation
    B_beta_shig = zeros(cellxmax, cellymax, cellzmax)             # B sigma for one-wave approximation
    C_beta_shig = zeros(cellxmax, cellymax, cellzmax)             # C sigma for one-wave approximation

    jalphaP = zeros(cellxmax, cellymax, cellzmax)                 # Jacobian matrix for viscosity
    jbetaP  = zeros(cellxmax, cellymax, cellzmax)                 # Jacobian matrix for viscosity
    jgammaP = zeros(cellxmax, cellymax, cellzmax)                 # Jacobian matrix for viscosity

    delta_Q      = zeros(cellxmax, cellymax, cellzmax, nval)      # delta Q for lusgs
    delta_Q_temp = zeros(cellxmax, cellymax, cellzmax, nval)      # temporary delta Q for lusgs

    D  = zeros(cellxmax, cellymax, cellzmax)                      # Diagonal of LHS
    Lx = zeros(cellxmax, cellymax, cellzmax, nval, nval)          # lower of LHS
    Ly = zeros(cellxmax, cellymax, cellzmax, nval, nval)          # lower of LHS
    Ux = zeros(cellxmax, cellymax, cellzmax, nval, nval)          # upper of LHS
    Uy = zeros(cellxmax, cellymax, cellzmax, nval, nval)          # upper of LHS

    LdQ = zeros(cellxmax, cellymax, cellzmax, nval)               # lower of LHS
    UdQ = zeros(cellxmax, cellymax, cellzmax, nval)               # upper of LHS

    RHS_temp = zeros(cellxmax, cellymax, cellzmax, nval)          # temporary RHS
    res = zeros(cellxmax, cellymax, cellzmax, nval)               # residual

    # define at cell boundaries
    lambda_facex = zeros(cellxmax+1, cellymax, cellzmax)          # lambda for computational time steps
    lambda_facey = zeros(cellxmax, cellymax+1, cellzmax)          # lambda for computational time steps
    lambda_facez = zeros(cellxmax, cellymax, cellzmax+1)          # lambda for computational time steps

    # misc
    norm2 = zeros(nval)                                 # Residuals by norm-2
    I = zeros(nval, nval)                               # identity matrix
    for l in 1:nval
        I[l,l] = 1.0
    end

    return Qbasen, Qconn, Qconn_hat, Qbasem, dtau, lambda_facex, lambda_facey, lambda_facez,
            A_adv_hat_m, A_adv_hat_p, B_adv_hat_m, B_adv_hat_p, C_adv_hat_p, C_adv_hat_m, A_beta_shig, B_beta_shig, C_beta_shig,
            jalphaP, jbetaP, jgammaP, delta_Q, delta_Q_temp, D, Lx, Ly, Ux, Uy, LdQ, UdQ, RHS_temp, res,
            norm2, I
end