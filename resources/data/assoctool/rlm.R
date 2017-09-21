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

# RLM
suppressMessages(library(robust));
eval(parse(text=param_cmd));
# Set iteration number to 1000 to overcome convergence problem (if not set yet).
if (is.null(param_list[['control']])) {
	param_list$control <- lmRob.control(mxr=1000,mxf=1000,mxs=1000);
}
if (is.null(param_list[['x']])) {
	param_list$x <- TRUE;
}
if (is.null(param_list[['y']])) {
	param_list$y <- TRUE;
}
doOne <- function(i, mdata) {
	param_list$data[, opt$omics_var_name] <- get(mdata, i);
	result <- do.call(lmRob, param_list);
	sst <- var(result$y) * (length(result$y) - 1);
	tbl <- summary(result)$coef;
	if ("weights" %in% names(result)) {
		ssr <- unlist(lapply(1:NROW(tbl), function (x) sum(result$weights * as.vector(result$x[,x])) * tbl[x, 1]^2));
	} else {
		ssr <- unlist(lapply(1:NROW(tbl), function (x) sum(as.vector(result$x[,x])) * tbl[x, 1]^2));
	}
	rsq <- ssr / sst;
	# Returns P-value, Effect size, Standard Error, T-statistics, R^2
	result <- cbind(tbl[,4], tbl[,1], tbl[,2], tbl[,3], rsq);
	if (rownames(result)[1] == "(Intercept)") result <- result[-1,];
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "T", "RSq"), length(..patterns)), rep(..patterns, each=5), sep="_");

	return (result);
}
