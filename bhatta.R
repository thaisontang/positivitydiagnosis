# Input Matrix - psmat
exposure.lab <- "Drug"
exposure <- dataset[,exposure.lab]
c.exp <- as.character(exposure)
dataset <- dataset
labels <- unique(c.exp)
# labels <- paste(c("Inst"), 1:nhosp)
psmat <- ps


# Check if Categorical Exposure / Bin in Case of Continuous Exposures
# bin.width


# Resolution
grid.res <- 100
grid <- seq(0, 1, by = 1/grid.res)

# Labels
if(!(is.null(labels))) {
  if(length(labels) == ncol(psmat)) {
    levels <- labels
  } else {
    stop("Length of labels and psmat do not match")
  }
} else if(!(is.null(colnames(psmat)))) {
  levels <- colnames(psmat)
} else { levels <- 1:ncol(psmat) }

rownames(psmat) <- colnames(psmat) <- levels

# Bhattacharrya Coefficients
bhatta <- matrix(NA, nrow = length(levels), ncol = length(levels))
rownames(bhatta) <- colnames(bhatta) <- levels
for (i in 1:length(levels)) {
  for (j in 1:length(levels)) {
    index <- psmat[c.exp == levels[i], levels[i]]
    ref <- psmat[c.exp == levels[j], levels[i]]
    pindex <- pref <- rep(0, length(grid))
    tab <- table(findInterval(index, grid))/length(index)
    pindex[as.numeric(names(tab))] <- tab
    tab <- table(findInterval(ref, grid))/length(ref)
    pref[as.numeric(names(tab))] <- tab
    bhatta[i,j] <- sum(sqrt(pindex * pref))
  }
}

# Matrix to Long-Format Data
mat_pairs <- function(mat) {
  xvals <- rownames(mat)
  yvals <- colnames(mat)
  grid.vals <- expand.grid(xvals, yvals)
  data.fit2 <- data.frame(`RowVar` = grid.vals$Var1, `ColVar` = grid.vals$Var2, `Value` = as.vector(mat))
  data.fit2[complete.cases(data.fit2),]
}

bhatta_values <- as_tibble(mat_pairs(bhatta)) 

plot.sig <- ggplot(bhatta_values, aes(`RowVar`, `ColVar`)) +
  geom_tile(aes(fill = `Value`), colour = "white") +
  scale_fill_gradient(low = "blue", high = "white", limits=c(0,1)) + theme_bw() +
  ggtitle("Bhattacharyya Coefficients between Distributions of\nIndex Drugs & Reference Drugs") +
  labs(x = "Reference Drug", y = "Index Drug") +
  theme(legend.title = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 11, margin = margin(t = 18)),
        axis.title.y = element_text(size = 11, margin = margin(r = 18)),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11,
                                  margin = ggplot2::margin(b = 12)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_colourbar(
    frame.colour = "black", label.vjust = 1, title = "Bhattacharyya\nCoefficient\n")) 
print(plot.sig)