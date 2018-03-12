google_search <-
function (a) 
{
    a = as.character(a)
    b = gsub(" ", "+", as.character(a))
    return(paste0("<a href='http://www.google.com/#hl=en&q=", 
        b, "' target='_blank'>", a, "</a>"))
}
