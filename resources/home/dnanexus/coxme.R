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

if (is.na(pedigree_type)) {
	# For pedigree-related calls, coxme function is in kinship library for KINSHIP1, but in coxme library for KINSHIP2
	# If there is no pedigree involved, assume coxme library
	suppressMessages(library(coxme));
}
eval(parse(text=param_cmd));
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
	param_list$data[, opt$omics_var_name] <- txFun(as.numeric(mdata[i,]));
	x <- do.call(coxme, param_list);
	# coxme does not provide coefficient tables like the others, unfortunately. So, we need to copy and paste from print.coxme.R
	beta <- x$coefficients;
	nvar <- length(beta);
	nfrail<- nrow(x$var) - nvar;

	# Log likelihood computation, not really necessary
	#loglik <- x$loglik + c(0,0, x$penalty);
	#chi1 <- 2*diff(loglik[1:2]);
	#chi2 <- 2*diff(loglik[c(1,3)]);
	#ll <- rbind(c(chi1, x$df[1], pchisq(chi1, x$df[1], lower.tail=FALSE),
	#	chi1 - 2*x$df[1], chi1 - log(x$n[1])*x$df[1]),
	#	c(chi2, x$df[2], pchisq(chi2,x$df[2], lower.tail=FALSE),
	#	chi2 - 2*x$df[2], chi2 - log(x$n[1])*x$df[2]));
	#dimnames(ll) <- list(c("Integrated loglik", " Penalized loglik"), c("Chisq", "df", "p", "AIC", "BIC"));

	se <- sqrt(diag(x$var)[nfrail+1:nvar]);
	result <- cbind(pchisq((beta/ se)^2, 1, lower.tail=FALSE), beta, se, beta/se);

	# Returns P-value, Effect size, Standard Error, Z-statistics
	# We do not have R^2 here
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "Z"), length(..patterns)), rep(..patterns, each=4), sep="_");

	return (result);
}
