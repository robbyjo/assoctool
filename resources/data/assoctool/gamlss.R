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

# GAM
# NOTE: You can NOT specify the method like in vanilla gamlss call (e.g., method=mixed(1,20)).
# This will throw an error because of the R's different treatment on parameter execution.
# In essence, outside of function (and especially more so using eval), the parameters are
# evaluated prior to function execution. Thus, the specification of mixed() or RS() or CG()
# functions in assoctool parameters will throw an error because they are only defined
# *inside* of the gamlss function. To work around this, use: method=substitute(mixed(1,20))
# to let R know that this is a deferred execution.

suppressMessages(library(gamlss));
eval(parse(text=param_cmd));
..robust <- FALSE;
if (!is.null(param_list[['robust']])) {
	if (is.logical(param_list$robust)) {
		..robust <- param_list$robust;
	}
	param_list$robust <- NULL;
}
..hessian.fun <- "R";
if (!is.null(param_list[['hessian.fun']])) {
	if (param_list$hessian.fun %in% c("R", "PB")) {
		..hessian.fun <- param_list$hessian.fun;
	}
	param_list$hessian.fun <- NULL;
}

# Tell gamlss to be quiet
if (is.null(param_list[['control']])) {
	param_list$control <- gamlss.control(trace=FALSE);
}
if (is.null(param_list[['i.control']])) {
	param_list$i.control <- glim.control(bf.trace=FALSE);
}

# TODO: Check this function
doOne <- function(i, mdata) {
	param_list$data[, opt$omics_var_name] <- get(mdata, i);
	suppressMessages(result <- do.call(gamlss, param_list));

	# Much of the code below is taken from SUMMARY.R of gamlss package. The reason is that there is no way to extract
	# coefficients quietly with summary() function. The suppressMessages does not work. In addition, there is a bug below.
	# This is why I cannot use summary() function as-is.
	type <- "vcov";
	# This is a workaround so that function gamlss:::gen.likelihood (gen-Likelihood.R) does not throw an error on basically an unused branch of code.
	# If it throws an error, then the code will switch to QR approximation
	# This workaround is true as of version 5.0-1 (Nov 2016) release of gamlss.
	result$call["data"] <- NULL;
	covmat <- try(suppressWarnings(vcov(result, type="all", robust=..robust,  hessian.fun = ..hessian.fun)), silent = TRUE);
	if (any(class(covmat)%in%"try-error"||any(is.na(covmat$se)))) {
		type <- "qr";
	}
	tbl_mu <- tbl_sigma <- tbl_nu <- tbl_tau <- NULL;

	if (type=="vcov") {
		tbl <- cbind(2*pt(-abs(covmat$coef/covmat$se), result$df.res), covmat$coef, covmat$se, covmat$coef/covmat$se);
		idx <- 0;
		if ("mu"%in%result$parameters & length(eval(parse(text="result$mu.fix==TRUE")))==0 & result$mu.df != 0) {
			tbl_mu <- tbl[1:result$mu.qr$rank, , drop=FALSE];
			idx <- result$mu.qr$rank;
		}
		if ("sigma"%in%result$parameters & length(eval(parse(text="result$sigma.fix==TRUE")))==0 & result$sigma.df != 0) {
			tbl_sigma <- tbl[(idx+1):(idx+result$sigma.qr$rank), , drop=FALSE];
			idx <- idx + result$sigma.qr$rank;
		}
		if ("nu"%in%result$parameters & length(eval(parse(text="result$nu.fix==TRUE")))==0 & result$nu.df != 0) {
			tbl_nu <- tbl[(idx+1):(idx+result$nu.qr$rank), , drop=FALSE];
			idx <- idx + result$nu.qr$rank;
		}
		if ("tau"%in%result$parameters & length(eval(parse(text="result$tau.fix==TRUE")))==0 & result$tau.df != 0) {
			tbl_tau <- tbl[(idx+1):(idx+result$tau.qr$rank), , drop=FALSE];
		}
	} else {
		est.disp <- FALSE;
		estimatesgamlss<-function (coef.p, s.err, est.disp, df.r) {
			tvalue <- coef.p/s.err;
			if (!est.disp)  {
				coef.table <- cbind(2*pnorm(-abs(tvalue)), coef.p, s.err, tvalue);
			} else if (df.r > 0) {
				coef.table <- cbind(2*pt(-abs(tvalue), df.r), coef.p, s.err, tvalue);
			} else {
				coef.table <- cbind(Inf, coef.p, NA, NA);
			}
			rownames(coef.table) <- names(coef.p);
			return(coef.table)
		}
		if ("mu"%in%result$parameters & result$mu.df != 0) {
			Qr <- result$mu.qr;
			df.r <- result$noObs - result$mu.df;
			if (df.r > 0) est.disp <- TRUE;
			p1 <- 1:(result$mu.df-result$mu.nl.df);
			covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE]);
			tbl_mu <- estimatesgamlss(coef.p=result$mu.coefficients[Qr$pivot[p1]], s.err=sqrt(diag(covmat.unscaled)), est.disp=est.disp, df.r=df.r);
		}
		if ("sigma"%in%result$parameters & result$sigma.df != 0) {
			Qr <- result$sigma.qr;
			df.r <- result$noObs - result$sigma.df;
			p1 <- 1:(result$sigma.df-result$sigma.nl.df);
			covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE]);
			tbl_sigma <- estimatesgamlss(coef.p=result$sigma.coefficients[Qr$pivot[p1]], s.err=sqrt(diag(covmat.unscaled)), est.disp=est.disp, df.r=df.r);
		}
		if ("nu"%in%result$parameters & result$nu.df != 0) {
			Qr <- result$nu.qr;
			df.r <- result$noObs - result$nu.df;
			p1 <- 1:(result$nu.df-result$nu.nl.df);
			covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE]);
			tbl_nu <- estimatesgamlss(coef.p=result$nu.coefficients[Qr$pivot[p1]], s.err=sqrt(diag(covmat.unscaled)), est.disp=est.disp, df.r=df.r);
		}
		if ("tau"%in%result$parameters & result$tau.df != 0) {
			Qr <- result$tau.qr;
			df.r <- result$noObs - result$tau.df;
			p1 <- 1:(result$tau.df-result$tau.nl.df);
			covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE]);
			tbl_tau <- estimatesgamlss(coef.p=result$tau.coefficients[Qr$pivot[p1]], s.err=sqrt(diag(covmat.unscaled)), est.disp=est.disp, df.r=df.r);
		}
	}

	new_tbl <- c();
	if (!is.null(tbl_mu)) {
		if (rownames(tbl_mu)[1] == "(Intercept)") rownames(tbl_mu)[1] <- "Intercept";
		rownames(tbl_mu) <- paste("Mu-", rownames(tbl_mu), sep="");
		new_tbl <- rbind(new_tbl, tbl_mu);
	}
	if (!is.null(tbl_sigma)) {
		if (rownames(tbl_sigma)[1] == "(Intercept)") rownames(tbl_sigma)[1] <- "Intercept";
		rownames(tbl_sigma) <- paste("Sigma-", rownames(tbl_sigma), sep="");
		new_tbl <- rbind(new_tbl, tbl_sigma);
	}
	if (!is.null(tbl_nu)) {
		if (rownames(tbl_nu)[1] == "(Intercept)") rownames(tbl_nu)[1] <- "Intercept";
		rownames(tbl_nu) <- paste("Nu-", rownames(tbl_nu), sep="");
		new_tbl <- rbind(new_tbl, tbl_nu);
	}
	if (!is.null(tbl_tau)) {
		if (rownames(tbl_tau)[1] == "(Intercept)") rownames(tbl_tau)[1] <- "Intercept";
		rownames(tbl_tau) <- paste("Tau-", rownames(tbl_tau), sep="");
		new_tbl <- rbind(new_tbl, tbl_tau);
	}
	tbl <- new_tbl;
	rm(new_tbl);


	# Returns P-value, Effect size, Standard Error, Z-statistics
	result <- tbl;
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "Z"), length(..patterns)), rep(..patterns, each=4), sep="_");
	return (result);
}
