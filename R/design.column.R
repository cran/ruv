design.column <-
function (a, name = "X") 
{
    if (is.numeric(a)) {
        A = matrix(a)
        colnames(A) = name
        return(A)
    }
    else {
        a = as.factor(a)
        if (length(levels(a)) == 1) 
            return(NULL)
        A = model.matrix(~a)[, -1, drop = FALSE]
        colnames(A) = paste(name, levels(a)[-1], sep = ".")
        return(A)
    }
}
