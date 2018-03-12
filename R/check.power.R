check.power <-
function (power) 
{
    if (power <= 0) {
        warning("The \"power\" variable must be positive.")
        return(FALSE)
    }
    if (power > 1) 
        warning("The \"power\" variable is normally less than 1.  Are you sure this is what you want?")
    return(TRUE)
}
