listk <- function(...){
    # get variable names from input expressions
    cols <- as.list(substitute(list(...)))[-1]
    vars <- names(cols)
    Id_noname <- if (is.null(vars)) seq_along(cols) else which(vars == "")

    if (length(Id_noname) > 0)
        vars[Id_noname] <- sapply(cols[Id_noname], deparse)
    # ifelse(is.null(vars), Id_noname <- seq_along(cols), Id_noname <- which(vars == ""))
    x <- setNames(list(...), vars)
    return(x)
}


ET_main <- function(
    ae, Eec, m,
    L_0, E_0, H_0,
    uecr, k, z, d, zom, zov,
    Tr_K,
    VPD,  W, E_Energy, Qne,
    maxiters=100, show_progress = FALSE, ...)
{
    n <- 0
    L_upd <- L_0+1

    while (n == 0 || (abs(L_upd-L_0)>0.01 && n <= maxiters)) {
        n <- n + 1
        if (show_progress) sprintf('running i=%3d, diff = %.2f\n', n, abs(L_upd-L_0)) %>% cat()

        if (n > 1){
            E_0 <- Epa_upd
            H_0 <- H_upd
            L_0 <- L_upd
            ustar_0 <-ustar_upd
        }
        L_0_abs <- abs(L_0)

        if (L_0_abs<100) {
            y <-(z-d)/L_0
            if (L_0<0) {
                # "unstable"
                x <-(1-16*y)^0.25
                fsv <-2*log((1+x^2)/2)
                fsm <-2*log((1+x)/2)+log((1+x^2)/2)-2*atan(x)+pi/2
            } else {
                # "stable"
                if(y<=1){
                    fsv <-5*(zom/L_0)-y
                    fsm <-5*(zom/L_0)-y
                }else{
                    fsv <- -(5*log((z-d)/zom))
                    fsm <- -(5*log((z-d)/zom))
                }
            }
        } else {
            # "neutral"
            fsv <- 0
            fsm <- 0
        }
        ustar_upd = uecr*k/(log((z-d)/zom)-fsm*((z-d)/L_0))
        L_upd     = -ustar_upd^3/k/10/(H_0+0.61*Tr_K*1004*E_0/2450000)*1004*Tr_K
        fu_neutral= 0.622*k^2*uecr/287.04/Tr_K/(log((z-d)/zov)-fsv)/(log((z-d)/zom)-fsm)
        EA_upd    = (1-W)*fu_neutral*1000*VPD*30*60/2450000
        Epa_upd   = EA_upd + E_Energy
        H_upd     = Qne-Epa_upd*2450000/30/60
        E_upd     = (2*Epa_upd-m*ae)*(m*ae/Epa_upd)^2
        # browser()
    }

    listk(ustar_upd, L_upd, fu_neutral, EA_upd, Epa_upd, H_upd, E_upd)
}

optim_func <- function(ae, dm){
    params <- dm[, .(Eec, m, L_0, E_0, H_0, uecr, k, z, d, zom, zov, Tr_K, VPD,  W, E_Energy, Qne)] %>% as.list()
    dp  <- do.call(ET_main, c(list(ae = ae, maxiters=100), params))
    dp
}

ae_opt <- function(ae, dm){
    dp <- optim_func(ae, dm)

    x1   <- ae*dm$m
    x2   <- (x1/dp$Epa_upd)^2*(2*dp$Epa_upd-x1)
    dif2 <- (x2-dm$Eec)^2
    RMSE <- sqrt(mean(dif2, na.rm = T))

    return(RMSE)
}
