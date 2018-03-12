RUV1 <-
function (Y, eta, ctl, include.intercept = TRUE) 
{
    if (is.null(eta)) 
        return(Y)
    m = nrow(Y)
    n = ncol(Y)
    ctl = tological(ctl, n)
    if (is.numeric(eta)) 
        if (length(eta) == 1) 
            eta = matrix(1, n, 1)
    if (is.matrix(eta)) 
        if (nrow(eta) != n) 
            eta = t(eta)
    eta = design.matrix(eta, name = "eta", include.intercept = include.intercept)
    eta = t(eta)
    Yc = Y[, ctl, drop = FALSE]
    etac = eta[, ctl, drop = FALSE]
    if (sum(is.na(Yc)) == 0) 
        return(Y - Yc %*% t(etac) %*% solve(etac %*% t(etac)) %*% 
            eta)
    else {
        for (i in 1:m) {
            keep = !is.na(Yc[i, ])
            Yci = Yc[i, keep, drop = FALSE]
            etaci = etac[, keep, drop = FALSE]
            Y[i, ] = Y[i, ] - Yci %*% t(etaci) %*% solve(etaci %*% 
                t(etaci)) %*% eta
        }
        return(Y)
    }
}
