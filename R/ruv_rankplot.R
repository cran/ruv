ruv_rankplot <-
function (fit, pctl, X.col = "all", uniform.lines = 0) 
{
    checks = check.ggplot() & check.X.col_2(fit, X.col)
    if (checks) {
        if (is.numeric(X.col)) 
            X.col = colnames(fit$misc$X)[X.col]
        if (X.col == "all") {
            pvar = "F.p"
            xlab = "Rank (F test p-value)"
        }
        else {
            pvar = paste0("p_", X.col)
            xlab = paste0("Rank (", X.col, " p-value)")
        }
        if (length(pctl) == 1) 
            pctl = fit$C[[pctl]]
        ranks = rank(fit$C[[pvar]], ties.method = "first")[pctl]
        counts = ecdf(ranks)(1:nrow(fit$C)) * length(ranks)
        df = data.frame(x = 1:nrow(fit$C), y = counts)
        rankplot = ggplot(df, aes_string(x = "x", y = "y")) + 
            theme_bw() + xlab(xlab) + ylab("# of Positive Controls") + 
            coord_cartesian(xlim = c(0, nrow(fit$C)))
        if (!is.null(uniform.lines)) 
            for (a in uniform.lines) rankplot = rankplot + geom_line(data = data.frame(x = seq(0, 
                1, length.out = 100) * nrow(fit$C), y = seq(a, 
                1, length.out = 100) * length(ranks)), aes_string(x = "x", 
                y = "y"), color = "gray", alpha = 0.25) + geom_line(data = data.frame(x = seq(a, 
                1, length.out = 100) * nrow(fit$C), y = seq(0, 
                1, length.out = 100) * length(ranks)), aes_string(x = "x", 
                y = "y"), color = "gray", alpha = 0.25)
        rankplot = rankplot + geom_point(size = 1.5, shape = 3)
        return(rankplot)
    }
    else return(FALSE)
}
