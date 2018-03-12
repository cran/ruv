projectionplotvariables <-
function (Y, X, W) 
{
    W0 = residop(W, X)
    bwx = solve(t(X) %*% X) %*% t(X) %*% W
    svdw0 = svd(W0)
    vd = t((1/svdw0$d) * t(svdw0$v))
    W0 = W0 %*% vd
    bwx = bwx %*% vd
    W0Y = t(W0) %*% Y
    u = svd(W0Y %*% t(W0Y))$u
    W0 = W0 %*% u
    bwx = bwx %*% u
    projectionplotalpha = t(W0) %*% Y
    byx = solve(t(X) %*% X) %*% t(X) %*% Y
    projectionplotW = W0 + X %*% bwx
    return(list(byx = byx, bwx = bwx, projectionplotalpha = projectionplotalpha, 
        projectionplotW = projectionplotW))
}
