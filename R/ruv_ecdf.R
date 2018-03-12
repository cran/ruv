ruv_ecdf <-
function (fit, X.col = "all", power = 1, uniform.lines = 0) 
{
    checks = check.ggplot() & check.X.col_2(fit, X.col) & check.power(power)
    if (checks) {
        if (is.numeric(X.col)) 
            X.col = colnames(fit$misc$X)[X.col]
        if (X.col == "all") {
            pvar = "F.p"
            xlab = "P-value (F test)"
        }
        else {
            pvar = paste0("p_", X.col)
            xlab = paste0("P-value (", X.col, ")")
        }
        pvalplot = ggplot(fit$C, aes_string(paste0("`", pvar, 
            "`"))) + theme_classic() + xlab(xlab) + ylab("ECDF") + 
            coord_trans(x = trans_new("power", function(x) x^power, 
                function(x) x^(1/power)), y = trans_new("power", 
                function(x) x^power, function(x) x^(1/power)), 
                limx = c(0, 1), limy = c(0, 1))
        if (!is.null(uniform.lines)) 
            for (a in uniform.lines) pvalplot = pvalplot + geom_line(data = data.frame(x = seq(0, 
                1, length.out = 100), y = seq(a, 1, length.out = 100)), 
                aes_string(x = "x", y = "y"), color = "gray", 
                alpha = 0.25) + geom_line(data = data.frame(x = seq(a, 
                1, length.out = 100), y = seq(0, 1, length.out = 100)), 
                aes_string(x = "x", y = "y"), color = "gray", 
                alpha = 0.25)
        pvalplot = pvalplot + stat_ecdf(size = 1.5, alpha = 0.5)
        return(pvalplot)
    }
    else return(FALSE)
}
