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

suppressMessages(library(logistf));

# Fix the logistf p-value calculation
logistf <- function(formula = attr(data, "formula"), data = sys.parent(), pl = TRUE, alpha = 0.05,
		control, plcontrol, firth = TRUE, init, weights, plconf=NULL, dataout=TRUE, ...)
{
	#n <- nrow(data)
#    if (is.null(weights)) weights<-rep(1,nrow(data))
	call <- match.call()
	if(missing(control)) control<-logistf.control()
	if(pl==TRUE & missing(plcontrol)) plcontrol<-logistpl.control()
	
	mf <- match.call(expand.dots =FALSE)
	m <- match(c("formula", "data","weights", "na.action", "offset"), names(mf), 0L)
	#   mf<-model.frame(formula, data=data, weights=weights)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	y <- model.response(mf)
	n <- length(y)
	x <- model.matrix(formula, data = data) ## Model-Matrix 
	k <- ncol(x)    ## Anzahl Effekte
	cov.name <- labels(x)[[2]]
	weight <- as.vector(model.weights(mf)  )
	offset <- as.vector(model.offset(mf)   )
	if (is.null(offset)) offset<-rep(0,n)
	if (is.null(weight)) weight<-rep(1,n)
	
	if (missing(init)) init<-rep(0,k)
	if (is.null(plconf) & pl==TRUE) plconf<-1:k
	
	if (dimnames(x)[[2]][1] == "(Intercept)")  {
		int <- 1
		coltotest <- 2:k
	} else {
		int <- 0
		coltotest <-1:k
	}
	# FIXME: The logistf.fit has a very nasty bug. It will hang if x has some collinear columns
	fit.full<-logistf:::logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=1:k, init, control=control)
	fit.null<-logistf:::logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=int, init, control=control)
	fit <- list(coefficients = fit.full$beta, alpha = alpha, terms=colnames(x), var = fit.full$var, df = (k-int), loglik =c(fit.null$loglik, fit.full$loglik),
			iter = fit.full$iter, n = sum(weight), y = y, formula = formula(formula), call=match.call(), conv=fit.full$conv)
	names(fit$conv)<-c("LL change","max abs score","beta change")
	beta<-fit.full$beta
	covs<-fit.full$var
	pi<-fit.full$pi
	fit$firth<-firth
	fit$linear.predictors <- as.vector(x %*% beta + offset)
	fit$predict <- fit.full$pi
	fit$hat.diag <- fit.full$Hdiag
	if(firth)
		fit$method <- "Penalized ML"
	else fit$method <- "Standard ML"
	vars <- diag(covs)
	fit$prob <- 1 - pchisq((beta^2/vars), 1)
	fit$method.ci <- rep("Wald",k)
	fit$ci.lower <- as.vector(beta + qnorm(alpha/2) * vars^0.5)
	fit$ci.upper <- as.vector(beta + qnorm(1 - alpha/2) * vars^0.5)
	fit$alpha<-alpha
	fit$conflev<-1-alpha
	if(pl) {
		betahist.lo<-vector(length(plconf),mode="list")
		betahist.up<-vector(length(plconf),mode="list")
		pl.conv<-matrix(0,length(plconf),4)
		dimnames(pl.conv)[[1]]<-as.list(plconf)
		dimnames(pl.conv)[[2]]<-as.list(c("lower, loglik","lower, beta", "upper, loglik", "upper, beta"))
		LL.0 <- fit.full$loglik - qchisq(1 - alpha, 1)/2
		pl.iter<-matrix(0,k,2)
#        fit$ci.lower <- fit$ci.upper <- rep(0, k)
		icount<-0
		for(i in plconf) {
			icount<-icount+1
			inter<-logistf:::logistpl(x, y, beta, i, LL.0, firth, -1, offset, weight, plcontrol)
			fit$ci.lower[i] <- inter$beta
			pl.iter[i,1]<-inter$iter
			betahist.lo[[icount]]<-inter$betahist
			pl.conv.lower<-t(inter$conv)
			inter<-logistf:::logistpl(x, y, beta, i, LL.0, firth, 1, offset, weight, plcontrol)
			fit$ci.upper[i] <- inter$beta
			pl.iter[i,2]<-inter$iter
			betahist.up[[icount]]<-inter$betahist
			pl.conv.upper<-t(inter$conv)
			pl.conv[icount,]<-cbind(pl.conv.lower,pl.conv.upper)
			fit.i<-logistf:::logistf.fit(x,y, weight=weight, offset=offset, firth, col.fit=(1:k)[-i], control=control)
			#fit$prob[i] <- 1-pchisq(2*(fit.full$loglik-fit.i$loglik),1) # RJ modification
			fit$prob[i] <- pchisq(2*(fit.full$loglik-fit.i$loglik),1, lower.tail=FALSE)
			fit$method.ci[i] <- "Profile Likelihood"
		}
		fit$pl.iter<-pl.iter
		fit$betahist<-list(lower=betahist.lo, upper=betahist.up)
		fit$pl.conv<-pl.conv
	}
	names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- names(fit$coefficients) <- dimnames(x)[[2]]
	if(dataout) {
		fit$data<-data
		fit$weights<-weight
	}
	attr(fit, "class") <- c("logistf")
	fit
}

eval(parse(text=param_cmd));

doOne <- function(i) {
	param_list$data[, opt$omics_var_name] <- txFun(get(mdata, i));
	result <- do.call(logistf, param_list);
	tbl <- cbind(result$coefficients, sqrt(diag(result$var)), qchisq(result$prob, 1, lower.tail=FALSE), result$prob);
	
	# Returns P-value, Effect size, Standard Error, Chi-square
	result <- cbind(tbl[,4], tbl[,1], tbl[,2], tbl[,3]);
	if (rownames(result)[1] == "(Intercept)") result <- result[-1,];
	if (!is.null(opt$result_var_name)) {
		if (is.null(..patterns)) ..patterns <<- rownames(result)[grep(opt$result_var_pattern, rownames(result))];
		result <- result[..patterns,];
	}
	if (NROW(result) > 1) result <- as.vector(t(result));
	names(result) <- paste(rep(c("P", "Fx", "SE", "ChiSq"), length(..patterns)), rep(..patterns, each=4), sep="_");
	
	return (result);
}
