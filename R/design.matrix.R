design.matrix <-
function (a, name = "X", remove.collinear = TRUE, include.intercept = TRUE) 
{
    if (is.vector(a)) 
        a = matrix(a)
    if (is.matrix(a)) {
        if (is.null(colnames(a))) 
            colnames(a) = paste(name, 1:ncol(a), sep = "")
        for (i in 1:ncol(a)) {
            if (is.na(colnames(a)[i])) 
                colnames(a)[i] = paste(name, i, sep = "")
            if (colnames(a)[i] == "") 
                colnames(a)[i] = paste(name, i, sep = "")
        }
    }
    if (is.numeric(a)) 
        A = a
    else {
        if (is.factor(a)) {
            a = data.frame(a)
            names(a)[1] = paste(name, 1, sep = "")
        }
        a = data.frame(a)
        varnames = colnames(a)
        for (i in 1:length(varnames)) if (varnames[i] == paste("X", 
            i, sep = "")) 
            varnames[i] = paste(name, i, sep = "")
        A = design.column(a[, 1], name = varnames[1])
        if (ncol(a) > 1) 
            for (i in 2:ncol(a)) A = cbind(A, design.column(a[, 
                i], name = varnames[i]))
    }
    if (remove.collinear) 
        if (ncol(A) > 1) {
            if (ncol(A) > nrow(A)) 
                A = A[, 1:nrow(A)]
            A0 = scale(A, center = FALSE, scale = TRUE)
            d = svd(A)$d
            if (d[1]/d[length(d)] > 10^9) {
                warning("Collinearity detected.  Removing some columns.")
                toremove = NULL
                for (i in 2:ncol(A0)) {
                  A1 = A0[, 1:(i - 1)]
                  if (!is.null(toremove)) 
                    A1 = A1[, -toremove]
                  if (mean(residop(A0[, i, drop = FALSE], A1)^2) < 
                    10^(-9)) 
                    toremove = c(toremove, i)
                }
                A = A[, -toremove, drop = FALSE]
            }
        }
    if (include.intercept) {
        if (sum(residop(matrix(1, nrow(A), 1), A)^2) > 10^(-8)) {
            A = cbind(1, A)
            colnames(A)[1] = paste(name, "0", sep = "")
        }
    }
    return(A)
}
