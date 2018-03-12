ruv_svdgridplot <-
function (Y.data, Y.space = NULL, rowinfo = NULL, colinfo = NULL, 
    k = 1:3, Z = 1, left.additions = NULL, right.additions = NULL, 
    factor.labels = paste("S.V.", k)) 
{
    checks = check.ggplot() & check.gridExtra()
    if (checks) {
        get_legend = function(myggplot) {
            tmp = ggplot_gtable(ggplot_build(myggplot))
            leg = which(sapply(tmp$grobs, function(x) x$name) == 
                "guide-box")
            if (length(leg) != 0) 
                return(tmp$grobs[[leg]])
            else return(ggplot() + theme_void())
        }
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
        d.scaled = Y.space$d/sqrt(sum(Y.space$d^2))
        K = length(k)
        U = matrix(0, nrow(Y.data), K)
        V = matrix(0, ncol(Y.data), K)
        D = rep(0, K)
        ulim = matrix(0, K, 2)
        vlim = matrix(0, K, 2)
        for (i in 1:K) {
            k1 = floor(abs(k[i]))
            k2 = ceiling(abs(k[i]))
            a = c(1 - (abs(k[i]) - k1), abs(k[i]) - k1)
            a = a/sqrt(sum(a^2))
            a[1] = a[1] * sign(k[i])^k1
            a[2] = a[2] * sign(k[i])^k2
            u = a[1] * Y.space$u[, k1] + a[2] * Y.space$u[, k2]
            v = a[1] * Y.space$v[, k1] + a[2] * Y.space$v[, k2]
            ulimvect = a[1] * Y.space$d[k1] * Y.space$u[, k1] + 
                a[2] * Y.space$d[k2] * Y.space$u[, k2]
            vlimvect = a[1] * Y.space$d[k1] * Y.space$v[, k1] + 
                a[2] * Y.space$d[k2] * Y.space$v[, k2]
            ulim[i, ] = c(min(ulimvect), max(ulimvect))
            vlim[i, ] = c(min(vlimvect), max(vlimvect))
            U[, i] = as.vector(Y.data %*% v)
            V[, i] = as.vector(t(u) %*% Y.data)
            D[i] = sqrt(a[1]^2 * d.scaled[k1]^2 + a[2]^2 * d.scaled[k2]^2)
        }
        if (!is.null(rowinfo)) 
            rowinfo = data.frame(rowinfo)
        if (!is.null(colinfo)) 
            colinfo = data.frame(colinfo)
        plots = rep(list(NA), K^2)
        for (i in 1:K) for (j in 1:K) {
            if (i == j) {
                df.rect = data.frame(x = 0, y = 0, xmin = -D[i], 
                  ymin = -D[i], xmax = D[i], ymax = D[i])
                thisplot = ggplot() + theme_classic() + theme(axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank()) + theme(axis.text.y = element_blank(), 
                  axis.ticks.y = element_blank()) + theme(axis.line = element_blank()) + 
                  xlab(factor.labels[i]) + ylab(factor.labels[i]) + 
                  scale_x_continuous(position = "top") + geom_rect(data = df.rect, 
                  aes_string(xmin = "xmin", ymin = "ymin", xmax = "xmax", 
                    ymax = "ymax"), alpha = 0.4) + coord_cartesian(xlim = c(-1, 
                  1), ylim = c(-1, 1))
                plots[[(i - 1) * K + j]] = thisplot
            }
            else if (i > j) {
                df = data.frame(x = U[, j], y = U[, i])
                if (!is.null(rowinfo)) 
                  df = cbind(df, rowinfo)
                thisplot = ggplot(data = df, aes_string(x = "x", 
                  y = "y")) + theme_bw() + theme(axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
                  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                    axis.title.y = element_blank()) + theme(legend.position = "none") + 
                  coord_cartesian(xlim = ulim[j, ], ylim = ulim[i, 
                    ]) + geom_point()
                if (!is.null(rowinfo) & is.null(left.additions)) {
                  if (ncol(rowinfo) == 1) 
                    left.additions = list(aes(color = rowinfo), 
                      labs(color = ""))
                  if (ncol(rowinfo) == 2) 
                    left.additions = list(aes(color = rowinfo[[1]], 
                      shape = rowinfo[[2]]), labs(color = "", 
                      shape = ""))
                }
                if (!is.null(left.additions)) 
                  thisplot = thisplot + left.additions
                plots[[(i - 1) * K + j]] = thisplot
            }
            else {
                df = data.frame(x = V[, j], y = V[, i])
                if (!is.null(colinfo)) 
                  df = cbind(df, colinfo)
                thisplot = ggplot(df, aes_string(x = "x", y = "y")) + 
                  theme_bw() + theme(axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
                  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                    axis.title.y = element_blank()) + theme(legend.position = "none") + 
                  coord_cartesian(xlim = vlim[j, ], ylim = vlim[i, 
                    ]) + geom_point()
                if (!is.null(colinfo) & is.null(right.additions)) {
                  if (ncol(colinfo) == 1) 
                    right.additions = list(aes(color = colinfo, 
                      alpha = 0.2), scale_alpha_identity(), labs(color = ""))
                  if (ncol(colinfo) == 2) 
                    right.additions = list(aes(color = colinfo[[1]], 
                      shape = colinfo[[2]], alpha = 0.2), scale_alpha_identity(), 
                      labs(color = "", shape = ""))
                }
                if (!is.null(right.additions)) 
                  thisplot = thisplot + right.additions
                plots[[(i - 1) * K + j]] = thisplot
            }
        }
        layout_matrix = kronecker(t(matrix(1:K^2, K, K)), matrix(1, 
            3, 3))
        if (!is.null(left.additions) | !is.null(right.additions)) {
            layout_matrix = rbind(layout_matrix, K^2 + 1)
            layout_matrix = cbind(layout_matrix, K^2 + 2)
            layout_matrix = cbind(layout_matrix, K^2 + 2)
            layout_matrix[K * 3 + 1, ((K) * 3 + 1):ncol(layout_matrix)] = K^2 + 
                3
            layout_matrix[((K) * 3 + 1):nrow(layout_matrix), 
                K * 3 + (1:2)] = K^2 + 4
            plots[[K^2 + 1]] = plots[[K^2 + 2]] = ggplot() + 
                theme_void()
            if (!is.null(left.additions)) 
                plots[[K^2 + 1]] = get_legend(plots[[K + 1]] + 
                  theme(legend.position = "bottom"))
            if (!is.null(right.additions)) 
                plots[[K^2 + 2]] = get_legend(plots[[2]] + theme(legend.position = "right"))
        }
        return(grid.arrange(grobs = plots, layout_matrix = layout_matrix))
    }
    else return(FALSE)
}
