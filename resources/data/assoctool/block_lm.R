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

if (!isFixedEffectFormula(opt$fm)) {
	stop("***ERROR: BLOCKLM cannot accept non-fixed effect formula!");
}


# If omics is the ..response variable
if (..is_lhs <- isLeftTerm(opt$fm, opt$omics_var_name)) {
	pdata[, opt$omics_var_name] <- rnorm(NROW(pdata));
	..cur_fm <- opt$fm; 
	cat("Omics variable is on the left hand side of the formula (i.e., dependent variable).\n");
} else {
	..cur_fm <- removeTerm(opt$fm, opt$omics_var_name);
	cat("Omics variable is on the right hand side of the formula (i.e., independent variable). Formula without the omics variable is:\n");
	print(..cur_fm);
	if (!is.null(opt$result_var_name)) {
		cat("Note: BLOCKLM ignores request for variable output when the omics variable is one of the independent variables.\n");
	}
}

..model_fram <- model.frame(..cur_fm, data=pdata, ..weights=param_list$..weights, filter=1:NROW(pdata));
..filter <- model.extract(..model_fram, "filter");
# filter will contain the indices that do not have missing data

# Any ..weights defined?
..weights <- model.extract(..model_fram, "..weights");

if (is.null(..weights)) {
	..qr_pdat <- qr(model_mat);
} else {
	..weights <- sqrt(..weights);
	..qr_pdat <- qr(model_mat * ..weights);
}

..diag_r <- diag(chol2inv(qr.R(..qr_pdat)));
..DFE <- NROW(a) - ..qr_pdat$rank;

..response <- model.response(..model_fram);
if (!is.null(..response) & !is.null(..weights)) ..response <- ..response * ..weights;
if (!..is_lhs) {
	..response <- qr.resid(..qr_pdat, ..response);
	..res_resid <- sum(..response^2);
}

doBlock <- function(i, mdata) {
	# mdata is m x n at this point; m = num genes; n = num samples
	if (!is.null(..weights)) {
		mdata <- t(mdata[, ..filter]);
	} else {
		mdata <- t(sweep(mdata[, ..filter], MARGIN=2, ..weights, `*`));
	}
	# mdata is now n x m
	if (..is_lhs) {
		# This one is slower, but it works with ..weights and/or zero intercept (if needed)
		beta <- qr.coef(..qr_pdat, mdata); # p x m ; p = num of covariates
		resid <- qr.resid(..qr_pdat, mdata); # n x m
		se <- sqrt(..diag_r %o% (colSums(resid^2) / ..DFE)); # p x m ; %o% is outer product

		# If omics at LHS, we can pick and choose which covariates to output
		# This step assumes beta and se have row names (double check).
		if (rownames(beta)[1] == "(Intercept)") {
			beta <- beta[-1,];
			se <- se[-1,];
		}
		if (!is.null(opt$result_var_name)) {
			if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
			beta <- beta[..patterns,];
			se <- se[..patterns,];
		}
	} else {
		mdata <- qr.resid(..qr_pdat, mdata); # n x m # Residualize the covariants
		var_resid <- colSums(mdata^2); # 1 x m
		beta <- (..response %*% mdata) / var_resid; # 1 x m # Beta contains only the coef for omics variable
		se <- sqrt((..res_resid/var_resid - beta^2)/ DFE);
	}

	tstat <- beta / se;
	pval <- 2*pt(abs(tstat), df=..DFE, lower.tail=FALSE);

	param_list$data[, opt$omics_var_name] <- get(mdata, i);
	result <- do.call(..lmFun, param_list);
	tbl <- summary(result)$coef;
	#reduced_y <- result$model[, attr(result$terms, "..response")];
	# attr(result$terms, "..response") seems to be always 1
	reduced_y <- result$model[, 1];
	sst <- var(reduced_y) * (length(reduced_y) - 1);
	ssr <- tbl[,3]^2 * (sum(resid(result)^2) / result$df.residual);
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
