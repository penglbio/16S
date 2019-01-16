################################################################################
#
# Schlaeppi et al., PNAS, 2013
#
# Functions used in the CAP analysis and to generate Figures 4B, 5B, S15 & S22 
#
# Ruben Garrido Oter
#
# garridoo@uni-duesseldorf.de
#
################################################################################

library(vegan)
library(calibrate)
library(Biostrings)
library(RColorBrewer)

# variability_table generates a table with the relative variability explained 
# by each transformation given a CCA object

variability_table <- function(cca){

    chi <- c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
    variability_table <- cbind(chi, chi / chi[1])
    colnames(variability_table) <- c("inertia", "proportion")
    rownames(variability_table) <- c("total", "constrained", "unconstrained") 
    return(variability_table)

}


# cap_var_props returns the variance explained by the CPCoA transformation 
# for a given CCA object

cap_var_props <- function(cca){

        eig_tot <- sum(cca$CCA$eig)
        var_propdf <- cca$CCA$eig / eig_tot
        return(var_propdf)
}


# cap_var_props returns the variance explained by the PCoA transformation 
# for a given CCA object

pca_var_props <- function(cca) {

        eig_tot <- sum(cca$CA$eig)
        var_propdf <- cca$CA$eig / eig_tot
        return(var_propdf)

}


# cca_ci returns the confidence interval for the variance explained by a 
# CPCoA transformation for a given CCA object using a ANOVA-based permutation 
# test with 5000 permutations by default

cca_ci <- function(cca, permutations = 5000) {

    var_tbl <- variability_table(cca)
    p <- permutest(cca, permutations=permutations)
    ci <- quantile(p$F.perm, 
                   c(.05, .95)) * p$chi[1] / var_tbl["total", "inertia"]
    return(ci)

}

tr <- function(m, threshold = 2, scale = 1000) {
    
    m <- m / rowSums(m, na=T) * scale
    m <- log2(m + 1)
    return(m)

}

tra<-function(m, threshold = 2, scale = 1000) {
  m <- m / rowSums(m, na=T) * scale
  return(m)
}
# plot_cap plots the two principal components of a CPCoA transformation 
# subsetting different variables by color and / or by shape as well as plotting 
# the respective percentaje of the variance explained by each component as well 
# as by the transformation overall
# 
# p: coordinates of the samples in the CAP space
# d: matrix containing the environmental variables of the dataset
# col_var: vector of variables to subset each different color in the figure
# pch_var: vector of variables to subset each different shape
# col_comp: vector of values to subset each different color
# col_pch: vector of values each different shape
# shapes: vector of shapes to assign to each subset
# colors: vector of colors to assign to each subset
# file_name: name of the file to be used as output (.pdf by default)
# svg: boolean indicating if the output file should be written in svg format
# constraint: name of the variable constrained by the CPCoA
# ci: confindence interval for the perc. of variance explained by the transform.
# var_tbl: table containing the relative variabilities of the transformation
# cap_var: vector of variances explained by each component
# perm_anova: output of appliying the perm_anova test to the CCA object
# cex.main: cex of the main text of the figure
# cex: cex of the points plotted in the figure

plot_cap <- function(p, d, col_var, pch_var, col_comp, pch_comp, shapes, 
                     colors, file_name, svg = F, constraint, ci, var_tbl, 
                     cap_var, perm_anova, cex.main = 0.8, cex = 20,mar_var,leg_x,leg_y,leg_var,leg_point_type) {
	
	col <- rep("black", dim(p)[1])
	pch <- rep(19, dim(p)[1])
	for (i in 1:length(col_var)) {
		 index <- grep(col_var[i], d[rownames(p), col_comp[i]]) 
		 col[index] <- colors[i]
	}
	for (i in 1:length(pch_var)) {
		index <- grep(pch_var[i], d[rownames(p), pch_comp[i]]) 
		pch[index] <- shapes[i]
	}
		xlab <- paste("Constrained PCoA 1 (", 
                      format(cap_var[1] * 100, digits = 4), " %)", sep = "")
		ylab <- paste("Constrained PCoA 2 (", 
                      format(cap_var[2] * 100, digits = 4), " %)", sep = "")
		main <- paste(constraint, ": [", 
                      format(var_tbl["constrained", "proportion"] * 100, 
                             digits=2), 
                      "% of variance; P < ", 
                      format(perm_anova[1, 4], digits = 2),
				      "; 95% CI = ", format(ci[1] * 100, digits = 2), 
				      "%, ", format(ci[2] * 100, digits = 2), "%]", sep = "")
		if(svg) svg(file = file_name) else pdf(file = file_name)
		par(mar=mar_var)
		plot(p, col = col, pch = pch, xlab = xlab, 
             ylab = ylab, main = main, cex.main = cex.main, cex = cex)
		abline(v = 0, h = 0, lty = "dotted")
		legend(leg_x,leg_y,legend =leg_var,col =unique(col),pch=leg_point_type,bty ="n",xpd = TRUE,cex=0.9,pt.cex = 1.2)
		dev.off()
}


# plot_cap plots the two principal components of a PCoA transformation 
# subsetting different variables by color and / or by shape as well as 
# plotting the respective percentaje of the variance explained by each 
# component as well as by the transformation overall

plot_pcoa <- function(p, d, variables, col_var, pch_var, col_comp, 
                      pch_comp, shapes, colors, file_name, svg=F, 
                      pcoa_var, cex.main = 0.8, cex = 1.5) {
	
	col <- rep("black", dim(p)[1])
	pch <- rep(19, dim(p)[1])
	for (i in 1:length(col_var)) {
		index <- grep(col_var[i], d[rownames(p), col_comp[i]]) 
		col[index] <- colors[i]
	}
	for (i in 1:length(pch_var)) {
		index <- grep(pch_var[i], d[rownames(p), pch_comp[i]]) 
		pch[index] <- shapes[i]
	}
		xlab <- paste("PCo 1 (", 
                      format(pcoa_var[1] * 100, digits = 4), " %)", sep = "")
		ylab <- paste("PCo 2 (", 
                      format(pcoa_var[2] * 100, digits = 4), " %)", sep = "")
		main <- "Unconstrained PCoA"	
		if(svg) 
            svg(file=file_name) 
        else 
            pdf(file=file_name)
		plot(p, col = col, pch = pch, xlab = xlab, ylab = ylab, 
             main = main, cex.main = cex.main, cex = cex)
		abline(v = 0, h = 0, lty = "dotted")
  	legend(0.2,0.4,legend = col_var ,col =col,pch=16,bty ="n",xpd = TRUE,cex=0.9,pt.cex = 1.2,ncol=6)
dev.off()
}


# plot_pcoa_spp plots the species scores obtained from a CAP analysis 
# contained in the vector of coordinates p. Allows to subset species 
# scores by taxa and to plot different sizes by relative abundances

plot_pcoa_spp <- function(p, centroids = NULL, pch, file_name, svg = F, 
                          pcoa_var, abundances = NULL, otu_subset = NULL, phy,
                          taxonomy = NULL, cex.main = 0.8, cex = 0.8) {

	col <- rep("black", dim(p)[1])

	if (!is.null(taxonomy)) {
	  for (i in 1:length(phy)) {
	    colors[taxonomy[rownames(p),]$Phylum == phy[i]] <- ful_cols[i]
	  }
	  
	}
	

	if(!is.null(abundances)) {
		cex <- colMeans(abundances[, rownames(p)] * 0.5)
	}


		xlab <- paste("PCo 1 (", 
                      format(pcoa_var[1] * 100, digits=4), " %)", sep="")
		ylab <- paste("PCo 2 (", 
                      format(pcoa_var[2] * 100, digits=4), " %)", sep="")
		main <- "Unconstrained PCoA"	
		if(svg) svg(file = file_name) else pdf(file = file_name)



		plot(p, col = col, pch = pch, xlab = xlab, ylab = ylab, main = main, 
             cex.main = cex.main, cex = cex)
		abline(v = 0, h = 0, lty = "dotted")	

		if(!is.null(otu_subset)) {
			idx <- which(rownames(p) %in% otu_subset)
			points(p[idx, 1:2], pch = 16, col = col[idx], cex = cex[idx])
		}


		dev.off()
	
}

# plot_biplots plots the species scores as well as the biplots obtained from 
# a CAP analysis. The vector p should contain the coordinates of the species 
# scores and the centroids of the clusters from which the biplots are obtained 
# should be given in the variable centroids. Allows to subset by color and shape 
# as well as to plot different sizes by relative abundances

plot_biplots <- function(p, centroids, pch, col, file_name, svg = F, constraint, 
                         ci, var_tbl, cap_var, perm_anova, abundances = NULL, phy=NULL,
                         otu_subset = NULL, taxonomy = NULL, cex.main = 0.8, 
                         cex = 10,labels_var,mar_var,leg_x,leg_y,leg_var) {
	
	colors <- rep("grey22", dim(p)[1]) 
	ful_cols <-brewer.pal(9,"Set1")
		xlab <- paste("Constrained PCoA 1 (", 
                  format(cap_var[1] * 100, digits = 4), " %)", sep = "")
	ylab <- paste("Constrained PCoA 2 (", 
                  format(cap_var[2] * 100, digits = 4), " %)", sep = "")
	main <- paste(constraint, ": [", 
                  format(var_tbl["constrained", "proportion"] * 100, digits = 2),
      			  "% of variance; P < ", format(perm_anova[1, 4], digits = 2), 
			      "; 95% CI = ", format(ci[1] * 100, digits = 2), 
			      "%, ", format(ci[2] * 100, digits = 2), "%]", sep="")	
	if (svg) svg(file = file_name) else pdf(file = file_name)

	if (!is.null(taxonomy)) {
	  for (i in 1:length(phy)) {
	    colors[taxonomy[rownames(p),]$Phylum == phy[i]] <- ful_cols[i]
	  }
	
	}

	if (!is.null(abundances)) {
		cex <- colMeans(abundances[, rownames(p)] * 0.5)
	}

	par(mar=mar_var)
	plot(rbind(centroids, p), col = "white", pch = 19, xlab = xlab, ylab = ylab, 
            main = main, cex.main = cex.main)
	points(p, col = colors, cex = cex, pch = pch)
	arrows(x0 = c(0, 0, 0, 0), y0 = c(0, 0, 0, 0), 
           x1 = centroids[, 1], y1 = centroids[, 2],
		   col = col, pch = 17, length = 0.1, cex = 1.5, lwd = 2)
	abline(v=0, h=0, lty="dotted")
  text(centroids[,1],centroids[,2],labels = labels_var)
	if(!is.null(otu_subset)){
		idx <- which(rownames(p) %in% otu_subset)
		points(p[idx, 1:2], pch = 16, col = colors[idx], cex = cex[idx])
	}
	legend(leg_x,leg_y,legend =leg_var,col =ful_cols[1:length(phy)],pch=16,bty ="n",xpd = TRUE,cex=0.9,pt.cex = 1.2)
	dev.off()

}

