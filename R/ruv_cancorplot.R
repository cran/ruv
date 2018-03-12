ruv_cancorplot <-
function (Y, X, ctl, W1 = NULL, W2 = NULL) 
{
    checks = check.ggplot()
    if (checks) {
        X = design.matrix(X, include.intercept = FALSE)
        if (is.null(W1)) {
            W1 = svd(Y)$u
            K1 = min(nrow(Y), ncol(Y))
            W1 = W1[, 1:K1]
        }
        else K1 = ncol(W1)
        if (is.null(W2)) {
            Yc = Y[, ctl]
            W2 = svd(Yc)$u
            K2 = min(nrow(Yc), ncol(Yc))
            W2 = W2[, 1:K2]
        }
        else K2 = ncol(W2)
        K = nrow(Y)
        cc1 = cc2 = rep(NA, nrow(Y))
        for (k in 1:K1) cc1[k] = cancor(X, W1[, 1:k, drop = FALSE])$cor[1]^2
        for (k in 1:K2) cc2[k] = cancor(X, W2[, 1:k, drop = FALSE])$cor[1]^2
        df = data.frame(featureset = as.factor(rep(c("All", "Control"), 
            each = K)), K = rep(1:K, 2), cc = c(cc1, cc2))
        ccplot = ggplot(df) + coord_cartesian(ylim = c(0, 1)) + 
            theme_bw() + geom_point(aes_string(x = "K", y = "cc", 
            color = "featureset"), na.rm = TRUE) + geom_line(aes_string(x = "K", 
            y = "cc", color = "featureset"), na.rm = TRUE) + 
            labs(x = "k", y = bquote(cancor^2)) + scale_color_manual(name = "", 
            values = c("black", "green"))
        return(ccplot)
    }
    else return(FALSE)
}
