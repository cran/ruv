ruv_summary <-
function (Y, fit, rowinfo = NULL, colinfo = NULL, colsubset = NULL, 
    sort.by = "F.p", var.type = c("ebayes", "standard", "pooled"), 
    p.type = c("standard", "rsvar", "evar"), min.p.cutoff = 1e-24) 
{
    fit = variance_adjust(fit)
    m.orig = data.frame(mean.orig = apply(Y, 2, mean))
    m = data.frame(mean = apply(RUV1(Y, fit$eta, fit$ctl, include.intercept = fit$include.intercept), 
        2, mean))
    R = data.frame(fit$X)
    if (!is.null(rowinfo)) 
        R = cbind(R, data.frame(rowinfo))
    if (!is.null(fit$Z)) 
        R = cbind(R, data.frame(fit$Z))
    rownames(R) = rownames(Y)
    if (var.type[1] == "ebayes") {
        F.p = data.frame(fit$Fpvals.ebayes)
        F.p.BH = data.frame(fit$Fpvals.ebayes.BH)
        sigma2 = data.frame(fit$sigma2.ebayes)
        if (p.type[1] == "standard") {
            p = data.frame(t(fit$p.ebayes))
            p.BH = data.frame(t(fit$p.ebayes.BH))
            var.b = data.frame(t(fit$varbetahat.ebayes))
        }
        if (p.type[1] == "rsvar") {
            p = data.frame(t(fit$p.rsvar.ebayes))
            p.BH = data.frame(t(fit$p.rsvar.ebayes.BH))
            var.b = data.frame(t(fit$varbetahat.rsvar.ebayes))
        }
        if (p.type[1] == "evar") {
            p = data.frame(t(fit$p.evar))
            p.BH = data.frame(t(fit$p.evar.BH))
            var.b = data.frame(t(fit$varbetahat.evar))
        }
    }
    else if (var.type[1] == "pooled") {
        F.p = data.frame(fit$Fpvals.pooled)
        F.p.BH = data.frame(fit$Fpvals.pooled.BH)
        sigma2 = data.frame(fit$sigma2.pooled)
        if (p.type[1] == "standard") {
            p = data.frame(t(fit$p.pooled))
            p.BH = data.frame(t(fit$p.pooled.BH))
            var.b = data.frame(t(fit$varbetahat.pooled))
        }
        if (p.type[1] == "rsvar") {
            p = data.frame(t(fit$p.rsvar.pooled))
            p.BH = data.frame(t(fit$p.rsvar.pooled.BH))
            var.b = data.frame(t(fit$varbetahat.rsvar.pooled))
        }
        if (p.type[1] == "evar") {
            p = data.frame(t(fit$p.evar.pooled))
            p.BH = data.frame(t(fit$p.evar.pooled.BH))
            var.b = data.frame(t(fit$varbetahat.evar.pooled))
        }
    }
    else {
        F.p = data.frame(fit$Fpvals)
        F.p.BH = data.frame(fit$Fpvals.BH)
        sigma2 = data.frame(fit$sigma2)
        if (p.type[1] == "standard") {
            p = data.frame(t(fit$p))
            p.BH = data.frame(t(fit$p.BH))
            var.b = data.frame(t(fit$varbetahat))
        }
        if (p.type[1] == "rsvar") {
            p = data.frame(t(fit$p.rsvar))
            p.BH = data.frame(t(fit$p.rsvar.BH))
            var.b = data.frame(t(fit$varbetahat.rsvar))
        }
        if (p.type[1] == "evar") {
            p = data.frame(t(fit$p.evar))
            p.BH = data.frame(t(fit$p.evar.BH))
            var.b = data.frame(t(fit$varbetahat.evar))
        }
    }
    names(F.p) = "F.p"
    names(F.p.BH) = "F.p.BH"
    names(sigma2) = "sigma2"
    names(p) = paste("p", colnames(fit$X), sep = "_")
    names(p.BH) = paste("p.BH", colnames(fit$X), sep = "_")
    names(var.b) = paste("var.b", colnames(fit$X), sep = "_")
    b = data.frame(t(fit$betahat))
    names(b) = paste("b", colnames(fit$X), sep = "_")
    fit.ctl = data.frame(fit.ctl = fit$ctl)
    C = cbind(F.p, F.p.BH, p, p.BH, b, sigma2, var.b, fit.ctl, 
        m)
    if (!is.null(fit$eta)) 
        C = cbind(C, m.orig)
    if (!is.null(colinfo)) {
        if (ncol(colinfo) == ncol(Y)) 
            colinfo = t(colinfo)
        C = cbind(C, data.frame(colinfo))
    }
    rownames(C) = colnames(Y)
    misc = list()
    misc$method = fit$method
    misc$X = fit$X
    misc$Z = fit$Z
    misc$W = fit$W
    misc$eta = fit$eta
    misc$ctl = fit$ctl
    misc$k = fit$k
    misc$lambda = fit$lambda
    misc$include.intercept = fit$include.intercept
    misc$multiplier = fit$multiplier
    misc$byx = fit$byx
    misc$bwx = fit$bwx
    if (!is.null(fit$projectionplotalpha)) 
        misc$ppalpha = fit$projectionplotalpha
    else misc$ppalpha = fit$alpha
    if (var.type[1] == "ebayes") 
        misc$df = fit$df.ebayes
    else if (var.type[1] == "pooled") 
        misc$df = fit$df.pooled
    else misc$df = fit$df
    misc$colsubset = rep(TRUE, ncol(Y))
    if (!is.null(colsubset)) {
        C = C[colsubset, , drop = FALSE]
        misc$ppalpha = misc$ppalpha[, colsubset, drop = FALSE]
        misc$byx = fit$byx[, colsubset, drop = FALSE]
        misc$colsubset = colsubset
    }
    misc$colorder = 1:nrow(C)
    if (!is.null(sort.by)) {
        colorder = order(C[[sort.by]])
        C = C[colorder, , drop = FALSE]
        misc$ppalpha = misc$ppalpha[, colorder, drop = FALSE]
        misc$byx = fit$byx[, colorder, drop = FALSE]
        misc$colorder = colorder
    }
    pvalcolnames = c("F.p", "F.p.BH", paste("p", colnames(fit$X), 
        sep = "_"), paste("p.BH", colnames(fit$X), sep = "_"))
    for (pvalcol in pvalcolnames) {
        C[C[, pvalcol] < min.p.cutoff, pvalcol] = min.p.cutoff
    }
    return(list(Y = Y, R = R, C = C, misc = misc))
}
