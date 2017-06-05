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

suppressMessages(library(pscl));
eval(parse(text=param_cmd));
param_list$y <- TRUE; # We need this for Rsq computation

doOne <- function(i) {
	param_list$data[, opt$omics_var_name] <- txFun(as.numeric(mdata[i,]));
	result <- do.call(zeroinfl, param_list);
	tbl1 <- summary(result)$coef$count;
	tbl2 <- summary(result)$coef$zero;
	rownames(tbl2) <- paste("Zero", rownames(tbl2), sep="_");
	rownames(tbl2)[1] <- paste("ZeroIntercept");
	
	sst <- var(result$y) * (length(result$y) - 1);
	ssr <- tbl1[,3]^2 * (sum(result$residuals^2) / result$df.residual);
	rsq <- ssr / sst;
	
	# Returns P-value, Effect size, Standard Error, Z-statistics, R^2 for the count part
	result <- cbind(tbl1[,4], tbl1[,1], tbl1[,2], tbl1[,3], rsq);
	if (rownames(result)[1] == "(Intercept)") result <- result[-1,];
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "Z", "RSq"), length(..patterns)), rep(..patterns, each=5), sep="_");

	# Return the zero part P-value, Effect size, Standard Error, Z-statistics
	result2 <- cbind(tbl2[,4], tbl2[,1], tbl2[,2], tbl2[,3]);
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result2)[grep(opt$result_var_pattern, rownames(result2))];
		result2 <- result2[c(1, ..patterns),]; # Always include zero intercept because it is usually important
	}
	if (NROW(result2) > 1) result2 <- as.vector(t(result2));
	names(result2) <- paste(rep(c("P", "Fx", "SE", "Z"), length(..patterns)), rep(..patterns, each=4), sep="_");
	result <- c(result, result2);

	return (result);
}
