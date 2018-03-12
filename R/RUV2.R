RUV2 <-
function (Y, X, ctl, k, Z = 1, eta = NULL, include.intercept = TRUE, 
    fullW = NULL, svdyc = NULL, do_projectionplot = TRUE, inputcheck = TRUE) 
{
    if (is.data.frame(Y)) 
        Y = data.matrix(Y)
    m = nrow(Y)
    n = ncol(Y)
    X = rX = design.matrix(X, include.intercept = FALSE)
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
    if (k > sum(ctl)) 
        stop("k must not be larger than the number of controls")
    Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
    if (q > 0) {
        Y = residop(Y, Z)
        X = residop(X, Z)
    }
    if (is.null(fullW)) {
        if (is.null(svdyc)) 
            svdyc = lsvd(Y[, ctl, drop = FALSE])
        fullW = svdyc$u[, 1:min((m - p - q), sum(ctl)), drop = FALSE]
    }
    W = alpha = byx = bwx = projectionplotW = projectionplotalpha = NULL
    if (k > 0) {
        W = fullW[, 1:k, drop = FALSE]
        XZW = cbind(X, Z, W)
        if (do_projectionplot) {
            ppvars = projectionplotvariables(Y, X, W)
            byx = ppvars$byx
            bwx = ppvars$bwx
            projectionplotalpha = ppvars$projectionplotalpha
            projectionplotW = ppvars$projectionplotW
        }
    }
    else {
        XZW = cbind(X, Z)
    }
    A = solve(t(XZW) %*% XZW)
    AXZW = A %*% t(XZW)
    betagammaalphahat = AXZW %*% Y
    resids = Y - XZW %*% betagammaalphahat
    betahat = betagammaalphahat[1:p, , drop = FALSE]
    if (k > 0) 
        alpha = betagammaalphahat[(p + q + 1):(p + q + k), , 
            drop = FALSE]
    multiplier = as.matrix(diag(A)[1:p])
    df = m - p - q - k
    sigma2 = apply(resids^2, 2, sum)/df
    sigma2 = as.vector(sigma2)
    se = sqrt(multiplier %*% t(sigma2))
    tvals = betahat/se
    pvals = tvals
    for (i in 1:nrow(pvals)) pvals[i, ] = 2 * pt(-abs(tvals[i, 
        ]), df)
    Fstats = apply(betahat * (solve(AXZW[1:p, , drop = FALSE] %*% 
        t(AXZW[1:p, , drop = FALSE])) %*% betahat), 2, sum)/p/sigma2
    Fpvals = pf(Fstats, p, df, lower.tail = FALSE)
    return(list(betahat = betahat, sigma2 = sigma2, t = tvals, 
        p = pvals, Fstats = Fstats, Fpvals = Fpvals, multiplier = multiplier, 
        df = df, W = W, alpha = alpha, byx = byx, bwx = bwx, 
        X = rX, k = k, ctl = ctl, Z = Z, eta = eta, fullW = fullW, 
        projectionplotW = projectionplotW, projectionplotalpha = projectionplotalpha, 
        include.intercept = include.intercept, method = "RUV2"))
}
