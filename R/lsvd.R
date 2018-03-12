lsvd <-
function (Y) 
{
    temp = svd(Y %*% t(Y))
    return(list(u = temp$u, d = sqrt(temp$d)))
}
