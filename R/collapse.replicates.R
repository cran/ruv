collapse.replicates <-
function (df, M) 
{
    M = replicate.matrix(M)
    keepcol = apply(abs(residop(data.matrix(data.frame(lapply(df, 
        as.factor))), M)), 2, sum) < 1e-11
    keeprow = apply(M, 2, function(x) match(1, x))
    df = df[keeprow, keepcol, drop = FALSE]
    rownames(df) = colnames(M)
    return(df)
}
