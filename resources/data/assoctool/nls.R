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

eval(parse(text=param_cmd));
..lmFun <- match.fun(tolower(opt$method));
doOne <- function(i, mdata) {
	param_list$data[, opt$omics_var_name] <- get(mdata, i);
	result <- do.call(..lmFun, param_list);
	tbl <- summary(result)$coef;
	
	# Returns P-value, Effect size, Standard Error, T-statistics
	result <- cbind(tbl[,4], tbl[,1], tbl[,2], tbl[,3]);
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "T"), length(..patterns)), rep(..patterns, each=4), sep="_");
	
	return (result);
}
