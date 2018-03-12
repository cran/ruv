variance_adjust <-
function (fit, ctl.idx = NULL, ebayes = TRUE, pooled = TRUE, 
    evar = TRUE, rsvar = TRUE, bin = 10, rescaleconst = NULL) 
{
    n = ncol(fit$betahat)
    p = nrow(fit$betahat)
    if (is.null(ctl.idx)) {
        if (is.list(fit$ctl)) 
            rsvar = FALSE
        else ctl.idx = fit$ctl
    }
    if (length(fit$multiplier) == p) 
        fit$multiplier = matrix(fit$multiplier, p, n)
    if (length(fit$df) == 1) 
        fit$df = rep(fit$df, n)
    if (TRUE) {
        varbetahat = p.BH = matrix(NA, p, n)
        for (l in 1:p) {
            varbetahat[l, ] = fit$sigma2 * fit$multiplier[l, 
                ]
            p.BH[l, ] = p.adjust(fit$p[l, ], method = "BH")
        }
        fit$p.BH = p.BH
        fit$varbetahat = varbetahat
        fit$Fpvals.BH = p.adjust(fit$Fpvals, method = "BH")
    }
    if (rsvar) {
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            multiplier = mean(fit$betahat[l, ctl.idx]^2/fit$sigma2[ctl.idx])
            varbetahat[l, ] = fit$sigma2 * multiplier
            tvals[l, ] = fit$betahat[l, ]/sqrt(varbetahat[l, 
                ])
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), fit$df)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.rsvar = pvals
        fit$p.rsvar.BH = p.BH
        fit$varbetahat.rsvar = varbetahat
    }
    if (evar) {
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            varbetahat[l, ] = get_empirical_variances(fit$sigma2, 
                fit$betahat[l, ], bin = bin, rescaleconst = rescaleconst)
            tvals[l, ] = fit$betahat[l, ]/sqrt(varbetahat[l, 
                ])
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), Inf)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.evar = pvals
        fit$p.evar.BH = p.BH
        fit$varbetahat.evar = varbetahat
    }
    if (ebayes) {
        temp = sigmashrink(fit$sigma2, fit$df)
        fit$sigma2.ebayes = temp$sigma2
        fit$df.ebayes = temp$df
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            varbetahat[l, ] = fit$sigma2.ebayes * fit$multiplier[l, 
                ]
            tvals[l, ] = fit$betahat[l, ]/sqrt(varbetahat[l, 
                ])
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), fit$df.ebayes)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.ebayes = pvals
        fit$p.ebayes.BH = p.BH
        fit$varbetahat.ebayes = varbetahat
        fit$Fpvals.ebayes = pf(fit$Fstats * (fit$sigma2/fit$sigma2.ebayes), 
            p, fit$df.ebayes, lower.tail = FALSE)
        fit$Fpvals.ebayes.BH = p.adjust(fit$Fpvals.ebayes, method = "BH")
    }
    if (rsvar & ebayes) {
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            multiplier = mean(fit$betahat[l, ctl.idx]^2/fit$sigma2.ebayes[ctl.idx])
            varbetahat[l, ] = fit$sigma2.ebayes * multiplier
            tvals[l, ] = fit$betahat[l, ]/sqrt(fit$sigma2.ebayes * 
                multiplier)
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), fit$df.ebayes)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.rsvar.ebayes = pvals
        fit$p.rsvar.ebayes.BH = p.BH
        fit$varbetahat.rsvar.ebayes = varbetahat
    }
    if (pooled) {
        fit$sigma2.pooled = rep(mean(fit$sigma2, na.rm = TRUE), 
            n)
        fit$df.pooled = rep(sum(fit$df, na.rm = TRUE), n)
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            varbetahat[l, ] = fit$sigma2.pooled * fit$multiplier[l, 
                ]
            tvals[l, ] = fit$betahat[l, ]/sqrt(varbetahat[l, 
                ])
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), fit$df.pooled)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.pooled = pvals
        fit$p.pooled.BH = p.BH
        fit$varbetahat.pooled = varbetahat
        fit$Fpvals.pooled = pf(fit$Fstats * (fit$sigma2/fit$sigma2.pooled), 
            p, fit$df.pooled, lower.tail = FALSE)
        fit$Fpvals.pooled.BH = p.adjust(fit$Fpvals.pooled, method = "BH")
    }
    if (rsvar & pooled) {
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            multiplier = mean(fit$betahat[l, ctl.idx]^2/fit$sigma2.pooled[ctl.idx])
            varbetahat[l, ] = fit$sigma2.pooled * multiplier
            tvals[l, ] = fit$betahat[l, ]/sqrt(fit$sigma2.pooled * 
                multiplier)
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), fit$df.pooled)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.rsvar.pooled = pvals
        fit$p.rsvar.pooled.BH = p.BH
        fit$varbetahat.rsvar.pooled = varbetahat
    }
    if (evar & pooled) {
        varbetahat = pvals = p.BH = tvals = matrix(0, p, n)
        for (l in 1:p) {
            multiplier = mean(fit$betahat[l, ]^2/fit$sigma2.pooled)
            varbetahat[l, ] = fit$sigma2.pooled * multiplier
            tvals[l, ] = fit$betahat[l, ]/sqrt(fit$sigma2.pooled * 
                multiplier)
            pvals[l, ] = 2 * pt(-abs(tvals[l, ]), fit$df.pooled)
            p.BH[l, ] = p.adjust(pvals[l, ], method = "BH")
        }
        fit$p.evar.pooled = pvals
        fit$p.evar.pooled.BH = p.BH
        fit$varbetahat.evar.pooled = varbetahat
    }
    return(fit)
}
