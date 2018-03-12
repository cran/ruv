getinvsvd <-
function (Y, XZ, ctl) 
{
    m = nrow(Y)
    if (is.null(XZ)) {
        pq = 0
    }
    else if (length(XZ) == 1) {
        if (XZ == 1) {
            XZ = matrix(1, m, 1)
            pq = 1
        }
    }
    else {
        pq = ncol(XZ)
    }
    k = m - pq
    if (pq > 0) 
        temp = residop(Y, XZ)
    else temp = Y
    temp2 = svd(temp[, ctl] %*% t(temp[, ctl]))
    if ((k > sum(ctl)) & (pq > 0)) {
        U1 = temp2$u[, 1:sum(ctl), drop = FALSE]
        U2 = svd(residop(diag(m), cbind(XZ, U1)))$u[, 1:(m - 
            pq - sum(ctl)), drop = FALSE]
        d1 = temp2$d[1:sum(ctl)]
        d2 = rep(0, ncol(U2))
        invsvd = list(u = cbind(U1, U2), d = c(d1, d2))
    }
    else invsvd = temp2
    return(invsvd)
}
