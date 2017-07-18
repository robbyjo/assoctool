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

suppressMessages(library(kinship));
# NOTE: Kinship1 is always fitted using ML, and NOT REML!
eval(parse(text=param_cmd));
param_list$fixed <- param_list$formula;
param_list$formula <- NULL;
if (is.null(param_list[['random']]) & is.null(param_list[['varlist']])) {
	cat("varlist and random options are missing. Fill in automatically.\n");
	ped <- ped[colnames(ped) %in% pdata[, opt$pedigree_id_col], colnames(ped) %in% pdata[, opt$pedigree_id_col]];
	param_list$varlist <- list();
	param_list$varlist[[opt$pedigree_id_col]] <- ped;
	param_list$random <- as.formula(paste("~1|", opt$pedigree_id_col, sep=""));
} else {
	if (is.null(param_list[['random']]) | is.null(param_list[['varlist']])) stop("You need to fill in both random and varlist options, not just one of them!");
}
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
doOne <- function(i) {
	param_list$data[, opt$omics_var_name] <- txFun(get(mdata, i));
	result <- do.call(lmekin, param_list);
	tbl <- result$ctable;

	if (opt$use_pt) {
		tbl[,4] <- 2*pt(abs(tbl[,3]), result$df.residual, lower.tail=FALSE);
	}
	#reduced_y <- result$model[, attr(result$terms, "response")];
	# attr(result$terms, "response") seems to be always 1
	reduced_y <- result$fitted.values + result$residuals;
	sst <- var(reduced_y) * (length(reduced_y) - 1);
	ssr <- tbl[,3]^2 * (sum(result$residuals^2) / result$df.residual);
	rsq <- ssr / sst;
	
	# Returns P-value, Effect size, Standard Error, T-statistics, R^2
	result <- cbind(tbl[,4], tbl[,1], tbl[,2], tbl[,3], rsq);
	rownames(result) <- rownames(tbl);
	if (rownames(result)[1] == "(Intercept)") result <- result[-1,];
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "T", "RSq"), length(..patterns)), rep(..patterns, each=5), sep="_");
	return (result);
}
