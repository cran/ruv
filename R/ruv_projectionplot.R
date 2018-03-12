ruv_projectionplot <-
function (fit, X.col = 1, factor = "gradient", adjusted = TRUE) 
{
    checks = check.ggplot() & check.X.col(fit, X.col) & check.factor(fit, 
        factor)
    if (checks) {
        if (is.character(X.col)) 
            X.col = which(colnames(fit$misc$X) == X.col)[1]
        y = as.vector(fit$misc$byx[X.col, ])
        bwx = fit$misc$bwx
        alpha = fit$misc$ppalpha
        if (factor == "gradient") {
            graddir = bwx[X.col, ]/sqrt(sum(bwx[X.col, ]^2))
            slope = sum(bwx[X.col, ] * graddir)
            x = as.vector(t(matrix(graddir)) %*% alpha)
            xlabel = "Gradient Factor"
            ylabel = bquote(b[YX] ~ ~(.(colnames(fit$misc$X)[X.col])))
        }
        else {
            factor = as.numeric(factor)
            x = alpha[factor, ]
            slope = bwx[X.col, factor]
            xlabel = paste0("Factor ", factor)
            ylabel = bquote(b[YX] ~ ~(.(colnames(fit$misc$X)[X.col])))
            if (adjusted) {
                if (nrow(alpha) > 1) 
                  y = y - as.vector(bwx[X.col, -factor, drop = FALSE] %*% 
                    alpha[-factor, , drop = FALSE])
                ylabel = bquote(Adjusted ~ b[YX] ~ ~(.(colnames(fit$misc$X)[X.col])))
            }
        }
        df = cbind(fit$C, data.frame(pplot.x = x, pplot.y = y, 
            pplot.yhat = x * slope))
        pplot = ggplot(df, aes_string(x = "pplot.x", y = "pplot.y")) + 
            theme_bw() + xlab(xlabel) + ylab(ylabel) + geom_line(aes_string(x = "pplot.x", 
            y = "pplot.yhat"), color = "black") + geom_point(data = df, 
            aes_string(x = "pplot.x", y = "pplot.y"))
        return(pplot)
    }
    else return(FALSE)
}
