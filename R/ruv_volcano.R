ruv_volcano <-
function (fit, X.col = 1) 
{
    checks = check.ggplot() & check.X.col(fit, X.col)
    reverselog_trans <- function(base = exp(1)) {
        trans <- function(x) -log(x, base)
        inv <- function(x) base^(-x)
        trans_new(paste0("reverselog-", format(base)), trans, 
            inv, log_breaks(base = base), domain = c(1e-100, 
                Inf))
    }
    if (checks) {
        if (is.numeric(X.col)) 
            X.col = colnames(fit$misc$X)[X.col]
        bvar = paste0("b_", X.col)
        pvar = paste0("p_", X.col)
        vplot = ggplot(fit$C, aes_string(x = paste0("`", bvar, 
            "`"), y = paste0("`", pvar, "`"))) + theme_bw() + 
            xlab(bquote(hat(beta) ~ ~(.(X.col)))) + ylab("P-Value (-log scale)") + 
            scale_y_continuous(trans = reverselog_trans(10)) + 
            geom_point()
        return(vplot)
    }
    else return(FALSE)
}
