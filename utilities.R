###############
## UTILITIES ##
###############

## PLOTTING FUNCTIONS

## http://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
grid_arrange_shared_legend <- function(..., ncol = 1) {
  require(grid)
  plots <- list(...)
  if (length(plots) == 1) { ## already a list of plots
    plots <- plots[[1]]
  }
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, c(lapply(plots, function(x)
      x + theme(legend.position="none")), ncol = ncol)),
    legend,
    nrow = 2, 
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

#' Title Generate a list of ggplot barplots (one per factor level) using a template p
#'
#' @param p  ggplot template
#' @param data data.frame passed on to ggplot
#' @param vbl factor variable (character string) to splice data in smaller datasets
#' @param ylab y-axis label
#'
#' @return a list of ggplot objects
#' @export
#'
#' @examples
generate_barplot <- function(p, data, vbl, ylab = "share of communities") {

  levs <- levels(data[[vbl]])
  result <- lapply(levs, 
                   function(lev) { p + geom_bar(data = filter_(data, 
                                                               lazyeval::interp(~variable == value, variable = as.name(vbl), value = lev))) + 
                       labs(y = ylab) + ggtitle(lev)})
  return(result)
}

#' Title Generate a list of ggplot scatterplots (one per factor level) using a template p
#'
#' @param p  ggplot template
#' @param data data.frame passed on to ggplot
#' @param vbl factor variable (character string) to splice data in smaller datasets
#' @param ylab y-axis label
#'
#' @return a list of ggplot objects
#' @export
#'
generate_pointplot <- function(p, data, vbl, ylab = "Divergence competitor") {
  levs <- levels(data[[vbl]])
  result <- lapply(levs, 
                   function(lev) { p + geom_point(data = filter_(data, 
                                                                 lazyeval::interp(~variable == value, variable = as.name(vbl), value = lev))) + 
                       labs(y = ylab) + ggtitle(lev)})
  return(result)
}

#' Title Generate a list of ggplot boxplots (one per factor level) using a template p
#'
#' @param p  ggplot template
#' @param data data.frame passed on to ggplot
#' @param vbl factor variable (character string) to splice data in smaller datasets
#' @param ylab y-axis label
#'
#' @return a list of ggplot objects
#' @export
#'
generate_boxplot <- function(p, data, vbl, ylab = "Divergence competitor") {
  levs <- levels(data[[vbl]])
  result <- lapply(levs, 
                   function(lev) { p + geom_boxplot(data = filter_(data, 
                                                                 lazyeval::interp(~variable == value, variable = as.name(vbl), value = lev)), 
                                                    outlier.size = 0.8) + 
                       geom_boxplot(data = filter_(data, lazyeval::interp(~variable == value, variable = as.name(vbl), value = lev)), 
                                    aes(color = NULL), outlier.color = "transparent") +
                       labs(y = ylab) + ggtitle(lev)})
  return(result)
}

## COMPARISONS FUNCTIONS
## Determine best of frogs (x) or competitor (y) in a Dunnett test. 
determine_best <- function(x, y) {
  local.df <- data.frame(method = x, divergence = y)
  model <- lm(divergence ~ method, data = local.df)
  results <- fortify(summary(glht(model, linfct = mcp(method = "Dunnett"))))
  if (all(results$estimate >= 0 & results$p < 0.05)) {
    return("frogs")
  } else {
    better.method <- subset(results, p < 0.05 & estimate < 0)
    if (nrow(better.method)) {
      return("competitor")
    } else {
      return("tied")
    }
  }
}

## Determine best of frogs (x) or competitor (y) in a Dunn test with BH p-values correction
determine_best_np <- function(x, g, method = "BH") {
  ## set frogs as reference factor
  g <- relevel(g, ref = "frogs")
  dunn.results <- dunn.test.control(x, g, p.adjust = method)
  results <- data.frame(method = rownames(dunn.results$p.value), 
                        p = as.numeric(dunn.results$p.value), 
                        estimate = as.numeric(dunn.results$statistic))
  if (all(results$estimate >= 0 & results$p < 0.05)) {
    return("frogs")
  } else {
    better.method <- subset(results, p < 0.05 & estimate < 0)
    if (nrow(better.method)) {
      return("competitor")
    } else {
      return("tied")
    }
  }
}

## Robust t-test
my.paired.t.test <- function(x, y, ...) {
  if (sd(x - y)) {
    t.test(x, y, paired = TRUE, ...)
  } else {
    return(list(p.value = 0, estimate = mean(x -y)))
  }
}

## Robust Wilcoxon test
my.paired.wilcox.test <- function(x, y, ...) {
  if (sd(x - y)) {
    wilcox.test(x, y, paired = TRUE, ...)
  } else {
    return(list(p.value = 0))
  }
}

###################
## COLOR PALETTE ##
###################
manual.palette <- c("frogs"      = rgb(0, 1, 0, alpha = 0.6, maxColorValue = 1), 
                    "tied"       = "grey60", 
                    ## "competitor" = rgb(1, 0, 0, alpha = 0.6, maxColorValue = 1), 
                    "competitor" = rgb(0, 0.5, 1, alpha = 0.6, maxColorValue = 1), 
                    "mothur (MA)"     = rgb(1, 0.2, 0, alpha = 0.6, maxColorValue = 1), 
                    "mothur (SOP)" = rgb(1, 0.5, 0, alpha = 0.6, maxColorValue = 1),                   
                    "uparse (MA)"     = rgb(0, 0.2, 1, alpha = 0.6, maxColorValue = 1), 
                    "uparse (SOP)" = rgb(0, 0.5, 1, alpha = 0.6, maxColorValue = 1), 
                    "qiime (MA)"      = rgb(1, 0.2, 1, alpha = 0.6, maxColorValue = 1), 
                    "qiime (SOP)"  = rgb(1, 0.5, 1, alpha = 0.6, maxColorValue = 1))