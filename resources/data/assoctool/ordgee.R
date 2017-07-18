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

# ORDGEE
suppressMessages(library(geepack));
eval(parse(text=param_cmd));
cat("Note: ordgee only works with complete data. assoctool does NOT provide any prefiltering for you!\n");
if (is.null(param_list[['id']])) {
	param_list$data[, opt$id_col] <- as.numeric(factor(param_list$data[, opt$id_col]));
	param_list$id <- as.name(opt$id_col);
} else if (is(param_list$id, "character")) {
	if (!(param_list$id %in% colnames(param_list$data))) stop("Column name specified as the ID is not in your phenotype data.");
	param_list$data[, param_list$id] <- as.numeric(factor(param_list$data[, param_list$id]));
	param_list$id <- as.name(param_list$id);
} else stop("You need to specify which column name in your phenotype data that will serve as cluster ID for ordgee.");

doOne <- function(i) {
	param_list$data[, opt$omics_var_name] <- txFun(get(mdata, i));
	result <- do.call(ordgee, param_list);
	tbl <- summary(result)$mean;

	# No RSq here since it is an ordered value
	# Returns P-value, Effect size, Standard Error, Z-statistics
	result <- cbind(tbl[,4], tbl[,1], tbl[,2], tbl[,3]);
	rownames(result) <- rownames(tbl);
	# We do not remove intercept here since the intercept may be useful
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "Z"), length(..patterns)), rep(..patterns, each=4), sep="_");
	return (result);
}
