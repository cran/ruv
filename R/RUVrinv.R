RUVrinv <-
function (Y, X, ctl, Z = 1, eta = NULL, include.intercept = TRUE, 
    fullW0 = NULL, invsvd = NULL, lambda = NULL, k = NULL, l = NULL, 
    randomization = FALSE, iterN = 1e+05, inputcheck = TRUE) 
{
    if (is.data.frame(Y)) 
        Y = data.matrix(Y)
    m = nrow(Y)
    n = ncol(Y)
    X = design.matrix(X, include.intercept = FALSE)
    p = ncol(X)
    if (is.numeric(Z)) 
        if (length(Z) == 1) 
            Z = matrix(1, m, 1)
    if (!is.null(Z)) {
        Z = design.matrix(Z, name = "Z", include.intercept = include.intercept)
        q = ncol(Z)
    }
    else q = 0
    ctl = tological(ctl, n)
    if (inputcheck) 
        inputcheck1(Y, X, Z, ctl)
    Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
    if (!is.null(lambda)) 
        return(RUVinv(Y, X, ctl, Z = Z, fullW0 = fullW0, invsvd = invsvd, 
            lambda = lambda, randomization = randomization, iterN = iterN))
    if (is.null(k)) {
        if (is.null(l) & (ncol(X) > 1)) 
            warning("Neither lambda nor k are specified, so a call to getK will be made.  But p > 1 and l is not specified.  Arbitrarily setting l = 1.")
        temp = getK(Y, X, ctl, Z = Z, fullW0 = fullW0, inputcheck = FALSE)
        k = temp$k
        fullW0 = temp$fullW0
    }
    ruv4fit = RUV4(Y, X, ctl, k, Z = Z, fullW0 = fullW0)
    lambda = sum(ruv4fit$sigma2[ctl])
    return(RUVinv(Y, X, ctl, Z = Z, fullW0 = fullW0, invsvd = invsvd, 
        lambda = lambda, randomization = randomization, iterN = iterN))
}
