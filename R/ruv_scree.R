ruv_scree <-
function (Y = NULL, Z = 1, Y.svd = NULL) 
{
    checks = check.ggplot()
    if (checks) {
        if (is.null(Y.svd)) {
            if (is.numeric(Z)) 
                if (length(Z) == 1) 
                  Z = matrix(1, nrow(Y), 1)
            if (!is.null(Z)) {
                Y = residop(Y, Z)
                q = ncol(Z)
            }
            else q = 0
            svd.Y = svd(Y)
        }
        plotval = svd.Y$d[1:(length(svd.Y$d) - q)]
        zero = plotval < plotval[floor(length(plotval)/2)] * 
            1e-08
        plotval[zero] = min(plotval[!zero])
        screeplot = ggplot(data.frame(plotval, K = 1:length(plotval), 
            zero = zero), aes_string(x = "K", y = "plotval")) + 
            theme_bw() + xlab("K") + ylab("Singular Value") + 
            theme(axis.title.x = element_text(size = rel(1.3)), 
                axis.text.x = element_text(size = rel(1.3))) + 
            theme(axis.title.y = element_text(size = rel(1.3)), 
                axis.text.y = element_blank()) + coord_trans(y = "log", 
            x = "log") + geom_point()
        if (sum(zero) > 0) 
            screeplot = screeplot + aes(color = zero) + scale_color_manual(name = "Zero", 
                values = c("black", "red"))
        return(screeplot)
    }
    else return(FALSE)
}
