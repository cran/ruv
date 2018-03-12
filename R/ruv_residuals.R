ruv_residuals <-
function (fit, type = c("residuals", "adjusted.Y"), subset_and_sort = TRUE) 
{
    if (fit$misc$method == "RUVinv") 
        return(FALSE)
    WZX = cbind(fit$misc$W, fit$misc$Z, fit$misc$X)
    kq = ncol(WZX) - ncol(fit$misc$X)
    Y1 = RUV1(fit$Y, fit$misc$eta, fit$misc$ctl, include.intercept = fit$misc$include.intercept)
    agb = solve(t(WZX) %*% WZX) %*% t(WZX) %*% Y1
    if (type[1] == "residuals") 
        Y1 = Y1 - WZX %*% agb
    if (type[1] == "adjusted.Y") 
        Y1 = Y1 - WZX[, 1:kq, drop = FALSE] %*% agb[1:kq, , drop = FALSE]
    if (subset_and_sort) 
        Y1 = (Y1[, fit$misc$colsubset])[, fit$misc$colorder]
    return(Y1)
}
