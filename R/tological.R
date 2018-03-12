tological <-
function (ctl, n) 
{
    ctl2 = rep(FALSE, n)
    ctl2[ctl] = TRUE
    return(ctl2)
}
