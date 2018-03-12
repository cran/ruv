getK <-
function (Y, X, ctl, Z = 1, eta = NULL, include.intercept = TRUE, 
    fullW0 = NULL, cutoff = NULL, method = "select", l = 1, inputcheck = TRUE) 
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
    nc = sum(ctl)
    if (inputcheck) 
        inputcheck1(Y, X, Z, ctl)
    if (p > 1) 
        return(getK(Y, X[, l, drop = FALSE], ctl, Z = cbind(X[, 
            -l, drop = FALSE], Z), eta = eta, include.intercept = include.intercept, 
            fullW0 = fullW0, cutoff = cutoff, method = method, 
            inputcheck = FALSE))
    Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
    if (q > 0) {
        Y = residop(Y, Z)
        X = residop(X, Z)
    }
    if (method == "select") {
        if (nc <= 1000 & m <= 300) 
            method = "leave1out"
        else method = "fast"
        if (method == "fast" & (nc < 3 * m + 100 | m > 500)) 
            warning("Either m or n_c (or both) are too large.  This algorithm (\"getK\") will be slow, provide poor results, or both.  It is recommended to choose K in some other manner.  To override this warning, set \"method\" to either \"fast\" or \"leave1out.\"  Currently defaulting to method \"fast.\"", 
                immediate. = TRUE)
    }
    K1 = min(m - p - q - 1, nc - 2)
    X = X/sqrt(sum(X^2))
    Yc = Y[, ctl]
    Y0 = residop(Y, X)
    if (is.null(fullW0)) {
        fullW0 = svd(Y0 %*% t(Y0))$u[, 1:(m - p - q), drop = FALSE]
    }
    W0 = fullW0[, 1:K1, drop = FALSE]
    alpha = solve(t(W0) %*% W0) %*% t(W0) %*% Y0
    bycx = solve(t(X) %*% X) %*% t(X) %*% Yc
    ac = alpha[, ctl, drop = FALSE]
    sizeratios = getsizeratios(bycx, ac, method = method)
    if (is.null(cutoff)) 
        cutoff = getKcutoff(m, n)
    keep = sizeratios > cutoff
    K = sum(keep)
    return(list(k = K, cutoff = cutoff, sizeratios = sizeratios, 
        fullW0 = fullW0))
}
