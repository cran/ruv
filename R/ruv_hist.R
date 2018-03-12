ruv_hist <-
function (fit, X.col = "all", breaks = c(0, 0.001, 0.01, 0.05, 
    seq(0.1, 1, by = 0.1))) 
{
    checks = check.ggplot() & check.X.col_2(fit, X.col)
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
        pvalhist = ggplot(fit$C, aes_string(paste0("`", pvar, 
            "`"))) + theme_classic() + theme(axis.title.y = element_blank()) + 
            xlab(xlab) + geom_histogram(aes_string(y = "..density.."), 
            breaks = breaks, position = "identity", color = "black") + 
            geom_histogram(aes_string(y = "..density.."), breaks = breaks, 
                position = "identity")
        return(pvalhist)
    }
    else return(FALSE)
}
