rocplot.multiple <- function(test.data.list, groupName = "grp", predName = "res", title = "ROC Plot", p.value = TRUE){
    require(plyr)
    require(ggplot2)
    plotdata <- llply(test.data.list, function(x) with(x, rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
    plotdata <- list(roc = ldply(plotdata, function(x) x$roc),
                     stats = ldply(plotdata, function(x) x$stats)
    )
    
    if (p.value == TRUE){
        annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
    } else {
        annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
    }
    
    p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
        geom_line(aes(colour = .id)) +
        geom_abline (intercept = 0, slope = 1) +
        theme_bw() +
        scale_x_continuous("False Positive Rate (1-Specificity)") +
        scale_y_continuous("True Positive Rate (Sensitivity)") +
        scale_colour_brewer(palette="Set1", breaks = names(test.data.list), labels = paste(names(test.data.list), ": ", annotation, sep = "")) +
        opts(title = title,
             plot.title = theme_text(face="bold", size=14), 
             axis.title.x = theme_text(face="bold", size=12),
             axis.title.y = theme_text(face="bold", size=12, angle=90),
             panel.grid.major = theme_blank(),
             panel.grid.minor = theme_blank(),
             legend.justification=c(1,0), 
             legend.position=c(1,0),
             legend.title=theme_blank(),
             legend.key = theme_blank()
        )
    return(p)
}