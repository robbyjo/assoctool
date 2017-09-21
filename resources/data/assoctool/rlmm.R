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


suppressMessages(library(robustlmm));
eval(parse(text=param_cmd));

# RLMM
# Compute DFE (useful for p-value computation) using R's lm
# Since we only cares about the DFE, we can use any dummy values for y
# Fixed effect version of the formula above
pdata[, opt$omics_var_name] <- rnorm(NROW(pdata));
..DFE <- lm(convertToLMFormula(opt$fm), data=pdata)$df.residual;
if (is.null(param_list[['control']])) {
	param_list$control <- lmerControl(optimizer="bobyqa", calc.derivs=FALSE, check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev="ignore", check.nobs.vs.nRE="ignore");
}

doOne <- function(i, mdata) {
	param_list$data[, opt$omics_var_name] <- get(mdata, i);
	result <- do.call(rlmer, param_list);
	tbl <- summary(result)$coef;

	sst <- var(result@resp$y) * (length(result@resp$y) - 1);
	ssr <- tbl[,3]^2 * (sum(esult@resp$weights * resid(result)^2) / ..DFE);
	rsq <- ssr / sst;
	
	# Returns P-value, Effect size, Standard Error, T-statistics, R^2
	result <- cbind(2*pt(abs(tbl[,3]), df=..DFE, lower.tail=FALSE), tbl[,1], tbl[,2], tbl[,3], rsq);
	if (rownames(result)[1] == "(Intercept)") result <- result[-1,];
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "T", "RSq"), length(..patterns)), rep(..patterns, each=5), sep="_");
	
	return (result);
}
