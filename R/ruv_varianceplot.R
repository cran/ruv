ruv_varianceplot <-
function (fit, X.col = 1, power = 1/4) 
{
    checks = check.ggplot() & check.X.col(fit, X.col) & check.power(power)
    if (checks) {
        if (is.numeric(X.col)) 
            X.col = colnames(fit$misc$X)[X.col]
        sigma2var = "sigma2"
        bvar = paste0("b_", X.col)
        var.bvar = paste0("var.b_", X.col)
        vplot = ggplot(fit$C, aes(get(sigma2var), get(bvar)^2)) + 
            theme_bw() + xlab(bquote(sigma^2)) + ylab(bquote(hat(beta)^2 ~ 
            ~(.(X.col))))
        if (sd(fit$C$sigma2) != 0) 
            vplot = vplot + coord_trans(x = trans_new("power", 
                function(x) x^power, function(x) x^(1/power)), 
                y = trans_new("power", function(x) x^power, function(x) x^(1/power)))
        else vplot = vplot + coord_trans(y = trans_new("power", 
            function(x) x^power, function(x) x^(1/power)))
        vplot = vplot + geom_point() + geom_line(aes(get(sigma2var), 
            get(var.bvar)), color = "black")
        return(vplot)
    }
    else return(FALSE)
}
