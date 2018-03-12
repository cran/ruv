check.X.col <-
function (fit, X.col) 
{
    if (is.null(X.col)) {
        warning("X.col must be specified.")
        return(FALSE)
    }
    if (length(X.col) > 1) {
        warning("Only one column of X may be chosen (X.col must be length 1).")
        return(FALSE)
    }
    if (X.col == "all") {
        warning("X.col cannot be \"all\".  A specific column of X must be chosen.")
        return(FALSE)
    }
    ok = colnames(fit$misc$X)
    if (X.col %in% ok | X.col %in% 1:length(ok)) 
        return(TRUE)
    warning(paste("X.col must be an integer between 1 and", length(ok), 
        "or one of the following character strings:", paste(ok, 
            collapse = " ")))
    return(FALSE)
}
