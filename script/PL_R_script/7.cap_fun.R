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


# variability_table generates a table with the relative variability explained 
# by each transformation given a CCA object

variability_table <- function(cca){
  
  chi <- c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi / chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained") 
  return(variability_table)
  
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
# cap_var_props returns the variance explained by the CPCoA transformation 
# for a given CCA object

cap_var_props <- function(cca){
  
  eig_tot <- sum(cca$CCA$eig)
  var_propdf <- cca$CCA$eig / eig_tot
  return(var_propdf)
}
##plot
plot_cap <- function(p, d,centroids, col_var, pch_var, col_comp, pch_comp, shapes, 
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
		arrows(x0 = c(0, 0, 0, 0), y0 = c(0, 0, 0, 0), 
		       x1 = centroids[, 1], y1 = centroids[, 2],
		       col = unique(col), pch = 17, length = 0.1, cex = 1.5, lwd = 2)
		text(centroids[, 1], centroids[, 2],leg_var)
		legend(leg_x,leg_y,legend =leg_var,col =unique(col),pch=leg_point_type,bty ="n",xpd = TRUE,cex=0.9,pt.cex = 1.2)
		dev.off()
	
}


