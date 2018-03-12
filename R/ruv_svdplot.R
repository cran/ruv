ruv_svdplot <-
function (Y.data, Y.space = NULL, info = NULL, k = c(1, 2), Z = 1, 
    left = TRUE) 
{
    checks = check.ggplot()
    if (checks) {
        if (is.data.frame(Y.space)) 
            Y.space = data.matrix(Y.space)
        if (is.numeric(Z)) 
            if (length(Z) == 1) 
                Z = matrix(1, nrow(Y.data), 1)
        if (!is.null(Z)) 
            Y.data = residop(Y.data, Z)
        if (is.null(Y.space)) 
            Y.space = Y.data
        if (!is.list(Y.space)) {
            if (!is.null(Z)) 
                Y.space = residop(Y.space, Z)
            Y.space = svd(Y.space)
        }
        if (left) {
            Ysuv = Y.space$v
            Ysuv2 = Y.space$u
            N = nrow(Y.data)
        }
        else {
            Ysuv = Y.space$u
            Ysuv2 = Y.space$v
            N = ncol(Y.data)
        }
        UV = matrix(0, N, 2)
        uvlim = matrix(0, 2, 2)
        for (i in 1:2) {
            k1 = floor(abs(k[i]))
            k2 = ceiling(abs(k[i]))
            a = c(1 - (abs(k[i]) - k1), abs(k[i]) - k1)
            a = a/sqrt(sum(a^2))
            a[1] = a[1] * sign(k[i])^k1
            a[2] = a[2] * sign(k[i])^k2
            uv = a[1] * Ysuv[, k1] + a[2] * Ysuv[, k2]
            uvlimvect = a[1] * Y.space$d[k1] * Ysuv2[, k1] + 
                a[2] * Y.space$d[k2] * Ysuv2[, k2]
            uvlim[i, ] = c(min(uvlimvect), max(uvlimvect))
            if (left) 
                UV[, i] = as.vector(Y.data %*% uv)
            else UV[, i] = as.vector(t(uv) %*% Y.data)
        }
        if (left) 
            side = "Left"
        else side = "Right"
        df = data.frame(x = UV[, 1], y = UV[, 2])
        if (!is.null(info)) {
            info = data.frame(info)
            df = cbind(df, info)
        }
        svdplot = ggplot(df, aes_string(x = "x", y = "y")) + 
            theme_bw() + xlab(paste(side, "singular vector", 
            k[1])) + ylab(paste(side, "singular vector", k[2])) + 
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
            coord_cartesian(xlim = uvlim[1, ], ylim = uvlim[2, 
                ]) + geom_point()
        if (!is.null(info)) {
            if (ncol(info) == 1) 
                svdplot = svdplot + aes(color = info) + labs(color = "")
            if (ncol(info) == 2) 
                svdplot = svdplot + aes(color = info[, 1], shape = info[, 
                  2]) + labs(color = "", shape = "")
        }
        return(svdplot)
    }
    else return(FALSE)
}
