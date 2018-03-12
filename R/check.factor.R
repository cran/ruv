check.factor <-
function (fit, factor) 
{
    if (is.null(factor)) {
        warning("factor must be specified.")
        return(FALSE)
    }
    if (length(factor) > 1) {
        warning("Only one factor may be chosen.")
        return(FALSE)
    }
    if (factor == "gradient" | factor %in% 1:nrow(fit$misc$ppalpha)) 
        return(TRUE)
    warning(paste("factor must be either \"gradient\" or an integer between 1 and", 
        nrow(fit$misc$ppalpha)))
    return(FALSE)
}
