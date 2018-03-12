inputcheck1 <-
function (Y, X, Z, ctl, check.na = TRUE) 
{
    m = nrow(Y)
    n = ncol(Y)
    if (m > n) 
        warning("m is greater than n!  This is not a problem itself, but may indicate that you need to transpose your data matrix.  Please ensure that rows correspond to observations (e.g. microarrays) and columns correspond to features (e.g. genes).")
    if (check.na & sum(is.na(Y)) > 0) 
        warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) > 
        0) 
        warning("Y contains infinities.  This is not supported.")
    XZ = cbind(X, Z)
    d = svd(t(XZ) %*% XZ, nu = 0, nv = 0)$d
    if (d[1]/d[length(d)] > 10^10) 
        warning("There appears to be linear dependence between the columns of X and Z.")
    if (sum(ctl) == 0) 
        warning("No genes are defined as control genes.  This is not supported.")
    return(NULL)
}
