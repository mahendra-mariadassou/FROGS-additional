## Figures Bioinformatics submission

## Load Packages ----------------------------------------------------------------
library(reshape2)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
source("utilities.R")

## Import and format data -------------------------------------------------------
data <- read.table("utax.tsv", sep = "\t", header = TRUE)

data <- data %>% mutate(amplicon = ifelse(amplicon == "V4V4", "V4", as.character(amplicon)))

data$nb_OTU <- as.numeric(sub("sp", "", data$nb_OTU))
if ("OTU.lost" %in% names(data)) {
  data <- data %>% mutate(OTU.lost = -OTU.lost) %>%
    melt(id.vars = c("databank", "nb_OTU", "dataset", 
                     "set_number", "amplicon", "abundance_law", "method"), 
         value.name = "divergence", 
         variable.name = "rank")
  levels(data$rank) <- c(levels(data$rank)[1:(length(levels(data$rank))-2)], "FN", "FP")
} else {
  data <- data %>% 
    melt(id.vars = c("databank", "nb_OTU", "dataset", 
                     "set_number", "amplicon", "abundance_law", "method"), 
         value.name = "divergence", 
         variable.name = "rank")
}
data.2 <- data %>% dcast(databank + nb_OTU + dataset + set_number + amplicon + abundance_law + rank ~ method, value.var = "divergence", fun.aggregate = mean) %>% 
  melt(id.vars = c("databank", "nb_OTU", "dataset", 
                   "set_number", "amplicon", "abundance_law", "rank", "frogs"), 
       value.name = "divergence", 
       variable.name = "method")


## Order methods and keep only frogs, uparse, mothur and qiime
data.otus <- filter(data, (rank %in% c("FP", "FN"))) %>% mutate(error_type = rank) %>% 
  mutate(method = factor(method, levels = c("frogs", "uparse_sop", "uparse", "mothur_sop", "mothur", "qiime_sop", "qiime"), 
                         labels = c("frogs", "uparse", "uparse (MA)", "mothur", "mothur (MA)", "qiime", "qiime (MA)"))) %>% 
  filter(method %in% c("frogs", "uparse", "mothur", "qiime"))
data.2 <- filter(data.2, !(rank %in% c("FP", "FN"))) %>% mutate(rank = factor(rank)) %>% 
  mutate(method = factor(method, levels = c("uparse_sop", "uparse", "mothur_sop", "mothur", "qiime_sop", "qiime"),
                         labels = c("uparse (SOP)", "uparse (MA)", "mothur (SOP)", "mothur (MA)", "qiime (SOP)", "qiime (MA)")))

## Custom palette --------------------------------------------------------------------------
manual.palette <- c("frogs"      = rgb(0, 1, 0, alpha = 0.6, maxColorValue = 1), 
                    "tied"       = "grey60", 
                    "competitor" = rgb(0, 0.5, 1, alpha = 0.6, maxColorValue = 1), 
                    "mothur (MA)"     = rgb(1, 0.2, 0, alpha = 0.6, maxColorValue = 1), 
                    "mothur" = rgb(1, 0.5, 0, alpha = 0.6, maxColorValue = 1),                   
                    "uparse (MA)"     = rgb(0, 0.2, 1, alpha = 0.6, maxColorValue = 1), 
                    "uparse" = rgb(0, 0.5, 1, alpha = 0.6, maxColorValue = 1), 
                    "qiime (MA)"      = rgb(1, 0.2, 1, alpha = 0.6, maxColorValue = 1), 
                    "qiime"  = rgb(1, 0.5, 1, alpha = 0.6, maxColorValue = 1))

## Compute statistics -----------------------------------------------------------------------
data.wilcox.test <- data.2 %>% group_by(databank, nb_OTU, dataset, amplicon, abundance_law, rank, method) %>% 
  summarize(pval = wilcox.test(frogs, divergence, paired = TRUE)$p.value, 
            measure = wilcox.test(frogs, divergence, paired = TRUE, conf.int = TRUE)$estimate) %>% 
  mutate(best = ifelse(measure >0, "competitor", "frogs")) %>% 
  mutate(best = ifelse(pval < 0.05, best, "tied")) %>%
  mutate(best = factor(best, levels = c("frogs", "tied", "competitor")))

## Plot graphic ----------------------------------------------------------------------------

## Figure 3 (main text)
p <- ggplot(mapping = aes(x = rank, fill = best)) + facet_grid(abundance_law + amplicon ~ nb_OTU) + theme_bw() + 
  scale_fill_manual(values = manual.palette) + theme(axis.text.x = element_text(angle = 90), 
                                                     axis.text = element_text(size = 6), 
                                                     axis.title = element_text(size = 8), 
                                                     strip.text = element_text(size = 6), 
                                                     plot.title = element_text(size = 10, hjust = 0.5))
plot.list <- generate_barplot(p, 
                              data = data.wilcox.test %>% 
                                filter(method %in% c("uparse (SOP)", "mothur (SOP)", "qiime (SOP)")) %>% 
                                mutate(method = droplevels(method)), 
                              vbl = "method", ylab = "Number of communities")
g <- grid_arrange_shared_legend(plot.list, ncol = 1, plot = FALSE)
## export to tiff
ggsave(plot = g, file = "Escudie_FROGS_Fig3.tiff", width = 86, height = 231, units = "mm", dpi = 350)

## Figure 4 (main text)
p <- ggplot(data.otus, mapping = aes(x = method, fill = method, y = divergence, alpha = amplicon)) + 
  facet_grid(nb_OTU ~ abundance_law, scales = "free_y") + 
  scale_y_continuous(limits = c(0, NA)) + 
  scale_fill_manual(values = manual.palette) + 
  scale_color_manual(values = manual.palette) + 
  scale_alpha_discrete(guide = guide_legend(override.aes = list(fill = "grey10")), range = c(0.2, 0.6)) + 
  theme_bw() + 
  theme(axis.text = element_text(color = "black"), 
        axis.text.x = element_text(angle = 90)) + 
  ylab(NULL)
p.fn <- p + geom_boxplot(data = filter(data.otus, error_type == "FN")) + ggtitle("False Negative OTUs")
p.fp <- p + geom_boxplot(data = filter(data.otus, error_type == "FP")) + ggtitle("False Positive OTUs")
## Export to tiff
g <- grid_arrange_shared_legend(p.fp, p.fn, ncol = 2, plot = FALSE)
ggsave(plot = g, file = "Escudie_FROGS_Fig4.tiff", width = 178, height = 127, units = "mm", dpi = 350)

## Figure 
