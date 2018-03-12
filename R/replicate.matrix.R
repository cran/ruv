replicate.matrix <-
function (a, burst = NULL, return.factor = FALSE, name = "M", 
    sep = "_", burstsep = "_") 
{
    if (is.matrix(a) & is.numeric(a)) {
        m = nrow(a)
        b = abs(a) > 10^(-8)
        if ((sum(abs(apply(a, 1, sum) - rep(1, m))) < 10^(-8)) & 
            (sum(abs(apply(b, 1, sum) - rep(1, m))) < 10^(-8))) {
            M = a
            if (is.null(colnames(M))) 
                colnames(M) = paste(name, 1:ncol(M), sep = "")
            for (i in 1:ncol(M)) {
                if (is.na(colnames(M)[i])) 
                  colnames(M)[i] = paste(name, i, sep = "")
                if (colnames(M)[i] == "") 
                  colnames(M)[i] = paste(name, i, sep = "")
            }
            colnames(M) = make.names(colnames(M), unique = TRUE)
            if (is.null(burst) & (!return.factor)) 
                return(M)
            a = rep("", m)
            for (i in 1:ncol(M)) a[M[, i] > 0.5] = colnames(M)[i]
            a = as.factor(a)
        }
    }
    a = data.frame(a)
    A = matrix("", nrow(a), ncol(a))
    for (i in 1:ncol(a)) A[, i] = as.character(a[, i])
    a = make.names(apply(A, 1, paste, collapse = sep))
    if (!is.null(burst)) {
        for (b in burst) {
            bn = sum(a == b)
            if (bn == 0) 
                warning(paste("Unable to burst factor level: ", 
                  b, "  -- not found.", sep = ""))
            a[a == b] = paste(b, 1:bn, sep = burstsep)
        }
    }
    a = as.factor(a)
    if (return.factor) 
        return(a)
    M = model.matrix(~a - 1)
    colnames(M) = levels(a)
    return(M)
}
