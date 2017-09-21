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

# GEE
suppressMessages(library(gee));
eval(parse(text=param_cmd));
if (is.null(param_list[['id']])) {
	param_list$id <- ifelse (is.numeric(pdata[, opt$id_col]), pdata[, opt$id_col], as.numeric(factor(pdata[, opt$id_col])));
}
if (!is(param_list[['id']], "numeric")) {
	cat("Function gee does not like non-numerical ID. Converting to numerical ID.\n");
	param_list$id <- as.numeric(factor(param_list[['id']]));
}

doOne <- function(i, mdata) {
	param_list$data[, opt$omics_var_name] <- get(mdata, i);
	suppressMessages(result <- do.call(gee, param_list));
	tbl <- summary(result)$coef[, c(1,4,5)]; # Always the robust estimate

	# Returns P-value, Effect size, Standard Error, Z-statistics
	result <- cbind(2*pnorm(abs(tbl[,3]), lower.tail=FALSE), tbl[,1], tbl[,2], tbl[,3]);
	if (rownames(result)[1] == "(Intercept)") result <- result[-1,];
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "Z"), length(..patterns)), rep(..patterns, each=4), sep="_");
	return (result);
}
