ruv_rle <-
function (Y, rowinfo = NULL, probs = c(0.05, 0.25, 0.5, 0.75, 
    0.95), ylim = c(-0.5, 0.5)) 
{
    checks = check.ggplot()
    if (checks) {
        rle = t(apply(t(Y) - apply(Y, 2, median), 2, quantile, 
            probs = probs))
        colnames(rle) = c("min", "lower", "middle", "upper", 
            "max")
        df = cbind(data.frame(rle.x.factor = 1:nrow(rle)), data.frame(rle))
        if (!is.null(rowinfo)) {
            rowinfo = data.frame(rowinfo)
            df = cbind(df, rowinfo)
        }
        rleplot = ggplot(df, aes_string(x = "rle.x.factor")) + 
            geom_boxplot(aes_string(lower = "lower", middle = "middle", 
                upper = "upper", max = "max", min = "min", group = "rle.x.factor"), 
                stat = "identity") + theme_bw() + theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
            theme(axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.5))) + 
            geom_hline(yintercept = 0) + coord_cartesian(ylim = ylim)
        if (!is.null(rowinfo)) 
            if (ncol(rowinfo) == 1) 
                rleplot = rleplot + aes(fill = rowinfo) + labs(fill = "")
        return(rleplot)
    }
    else return(FALSE)
}
