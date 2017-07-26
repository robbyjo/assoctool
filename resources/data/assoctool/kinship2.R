# Association analysis tool
# Version: 0.1
# By: Roby Joehanes
#
# Copyright 2016-2017 Roby Joehanes
# This file is distributed under the GNU General Public License version 3.0.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

suppressMessages(library(kinship2));
suppressMessages(library(coxme));
eval(parse(text=param_cmd));
# If method is not specified by user, default to REML
if (is.null(param_list[['method']])) {
	param_list$method <- "REML";
}
param_list$y <- TRUE; # We need this for Rsq computation
# I added the option of use.pt (TRUE/FALSE). If TRUE, then the p-value will be recomputed using
# T statistics with conservative degrees of freedom for errors (DFE). This is the default for this tool.
# The default in kinship2 uses pchisq, which is equivalent to discarding the degrees of freedom (more liberal).
# Please use the use.pt=TRUE.
opt$use_pt <- TRUE;
if (!is.null(param_list[['use.pt']])) {
	if (is.logical(param_list$use.pt)) {
		opt$use_pt <- param_list$use.pt;
	}
	param_list$use.pt <- NULL;
}
if (is(ped, "pedigreeList")) ped <- kinship(ped); # Pedigree data is in pedigreeList, convert to matrix
if (is.null(param_list[['varlist']]) & isDefined(ped)) {
	# Users did not specify varlist, but he/she loaded the pedigree
	cat("varlist is missing, but pedigree is loaded. Fill in automatically.\n");
	ped <- ped[rownames(ped) %in% pdata[, opt$pedigree_id_col], colnames(ped) %in% pdata[, opt$pedigree_id_col]];
	param_list$varlist <- list();
	param_list$varlist[[opt$pedigree_id_col]] <- ped;
	param_list$formula <- combineFormulas(param_list$formula, as.formula(paste("y~(1|", opt$pedigree_id_col, ")", sep="")));
	cat("Adding pedigree into formula:\n");
	print(param_list$formula);
}
doOne <- function(i) {
	param_list$data[, opt$omics_var_name] <- get(mdata, i);
	result <- do.call(lmekin, param_list);
	# The author does NOT provide table output, so compute ourselves
	fixedef <- result$coefficients$fixed;
	stderr <- sqrt(diag(result$var));
	zval <- fixedef/stderr;
	df_est2 <- result$n - length(result$coefficients$fixed) - length(result$vcoef) - 1;
	if (opt$use_pt) {
		pval <- 2*pt(abs(zval), df_est2, lower.tail=FALSE);
	} else {
		pval <- pchisq(zval^2, df=1, lower.tail=FALSE);
	}
	
	#reduced_y <- result$model[, attr(result$terms, "response")];
	# attr(result$terms, "response") seems to be always 1
	sst <- var(result$y) * (length(result$y) - 1);
	ssr <- zval^2 * (sum(result$residuals^2) / df_est2);
	rsq <- ssr / sst;
	result <- as.matrix(cbind(P=pval, Fx = fixedef, SE=stderr, T=zval, RSq=rsq));
	rownames(result) <- names(fixedef);
	
	# Returns P-value, Effect size, Standard Error, T-statistics, R^2
	if (rownames(result)[1] == "(Intercept)") result <- result[-1,];
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "T", "RSq"), length(..patterns)), rep(..patterns, each=5), sep="_");
	
	return (result);
}
