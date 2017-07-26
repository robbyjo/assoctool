# Modified pedigreemm source
# Based on pedigreemm 0.3-3
# By Roby Joehanes

# What is modified:
# 1. Allowing optimizer options to pass through
# 2. Allow direct specification of pedigree matrix
# 3. Shortcut to reconstruct necessary matrices instead of calling lmer/glmer

pedigreemm <-
		function(formula, data, family = NULL, REML = TRUE, pedigree = list(),
				control = list(), start = NULL, verbose = FALSE, 
				subset, weights, na.action, offset, contrasts = NULL,
				model = TRUE, x = TRUE, ...)
{
	gaus <- FALSE
	if (is.null(family)) {
		gaus <- TRUE
	} else {
		## copied from glm()
		if (is.character(family)) 
			family <- get(family, mode = "function", envir = parent.frame())
		if (is.function(family)) 
			family <- family()
		if (!inherits(family, "family")) stop("unknown family type")
		gaus <- family$family == "gaussian" && family$link == "identity"
	}
	mc <- match.call()
	lmerc <- mc                         # create a call to lmer
	lmerc[[1]] <- if (gaus) as.name("lmer") else as.name("glmer")
	lmerc$pedigree <- NULL
	if (!gaus) lmerc$REML <- NULL
	
	if (!length(pedigree))              # call [g]lmer instead
		return(eval.parent(lmerc))
	
	stopifnot(is.list(pedigree),        # check the pedigree argument
			length(names(pedigree)) == length(pedigree),
			all(sapply(pedigree, function(x) is(x, "pedigree") | is(x, "matrix") | is(x, "Matrix")))) # RJ modification
	
	lmerc[[1]] <- if (gaus) quote(lme4::lFormula) else quote(lme4::glFormula) # RJ modification
	lmf <- eval(lmerc, parent.frame())

	relfac <- pedigree          # copy the pedigree list for relfactor
	pnms <- names(pedigree)
	#pp <- lmf@pp # RJ modification
	#resp <- lmf@resp # RJ modification
	fl <- lmf$reTrms$flist # RJ modification
	stopifnot(all(pnms %in% names(fl)))
	asgn <- attr(fl, "assign")
	#Zt <- pp$Zt # RJ modification
	for (i in seq_along(pedigree)) {
		tn <- which(match(pnms[i], names(fl)) == asgn)
		if (length(tn) > 1)
			stop("a pedigree factor must be associated with only one r.e. term")
		ind <- (lmf$reTrms$Gp)[tn:(tn+1L)] # RJ modification
		rowsi <- (ind[1]+1L):ind[2]
		if (is(relfac[[i]], "pedigree")) relfac[[i]] <- relfactor(pedigree[[i]], rownames(lmf$reTrms$Zt)[rowsi]) # RJ modification
		lmf$reTrms$Zt[rowsi,] <- relfac[[i]] %*% lmf$reTrms$Zt[rowsi,] # RJ modification
	};
	#reTrms <- list(Zt=Zt,theta=lmf@theta,Lambdat=pp$Lambdat,Lind=pp$Lind,  # RJ modification starts
	#		lower=lmf@lower,flist=lmf@flist,cnms=lmf@cnms, Gp=lmf@Gp)
	dfl <- list(fr=lmf$fr, X=lmf$X, reTrms=lmf$reTrms, start=lmf$reTrms$theta)
	optimizer_opt <- ifelse(gaus, control$optimizer, control$optimizer[[2]]);
	if (is.null(optimizer_opt)) optimizer_opt <- "Nelder_Mead"; # RJ modification ends
	if (gaus) {
		dfl$REML = lmf$REML > 0L # RJ modification
		devfun <- do.call(mkLmerDevfun,dfl)
		opt <- optimizeLmer(devfun, optimizer=optimizer_opt,  # RJ modification starts
			restart_edge = control$restart_edge,
			boundary.tol = control$boundary.tol,
			control = control$optCtrl,
			calc.derivs=control$calc.derivs,
			use.last.params=control$use.last.params,...) # RJ modification ends
	} else {
		dfl$family <- family
		devfun <- do.call(mkGlmerDevfun,dfl)
		opt <- optimizeGlmer(devfun, optimizer=optimizer_opt, # RJ modification starts
			restart_edge=control$restart_edge,
			boundary.tol=control$boundary.tol,
			control = control$optCtrl,
			stage=2,
			calc.derivs=control$calc.derivs,
			use.last.params=control$use.last.params,...) # RJ modification ends
	}
	mm <- mkMerMod(environment(devfun), opt, lmf$reTrms, lmf$fr, mc)  # RJ modification
	cls <- if (gaus) "lmerpedigreemm" else "glmerpedigreemm"
	ans <- do.call(new, list(Class=cls, relfac=relfac,
					frame=mm@frame, flist=mm@flist, cnms=mm@cnms, Gp=mm@Gp,
					theta=mm@theta, beta=mm@beta,u=mm@u,lower=mm@lower,
					devcomp=mm@devcomp, pp=mm@pp,resp=mm@resp,optinfo=mm@optinfo))
	ans@call <- evalq(mc)
	ans
}

suppressMessages(library(lme4));
suppressMessages(library(pedigreemm));
eval(parse(text=param_cmd));
if (is.null(param_list[['family']])) {
	param_list$control <- lmerControl(optimizer="bobyqa", calc.derivs=FALSE, check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev="ignore", check.nobs.vs.nRE="ignore");
} else {
	param_list$control <- glmerControl(optimizer="bobyqa", calc.derivs=FALSE, check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev="ignore", check.nobs.vs.nRE="ignore");
}

pdata[, opt$omics_var_name] <- rnorm(NROW(pdata));
..DFE <- lm(convertToLMFormula(opt$fm), data=pdata)$df.residual;
param_list$pedigree <- ped;
..xx <- grep("\\(1\\|new_ids)", gsub(" ", "", opt$formula_str))
if (length(..xx) == 0) {
	opt$fm <- combineFormulas(opt$fm, as.formula("y ~ (1|new_ids)"));
	cat("Adding pedigree into formula:\n");
	print(opt$fm);
	param_list$formula <- opt$fm;
}
rm(..xx);
# If ped is a list of matrices, make sure the diagonal is 1
# Matrix from kinship MUST be scaled by two, i.e., kmat <- 2*kmat, where kmat is the matrix you get from kinship call.
# Pedigreemm always uses the Cholesky form of the pedigree matrix. We assume raw matrix is not yet decomposed.
if (is(ped, "bdsmatrix")) {
	ped <- as.matrix(ped);
}
if (is(ped, "matrix") | is(ped, "Matrix")) {
	if (opt$pedigree_type %in% c("kinship1", "kinship2")) {
		diag(ped) <- 0.5;
		cat("Pedigree matrix is from kinship. Multiply by 2 and perform Cholesky decomposition.\n");
		ped <- list(new_ids = chol(2 * ped));
	} else if (opt$pedigree_type %in% c("sparse_matrix", "dense_matrix")) {
		cat("Pedigree matrix is loaded through binary. Perform Cholesky decomposition.\n");
		ped <- list(new_ids = chol(ped));
	}
}

doOne <- function(i) {
	param_list$data[, opt$omics_var_name] <- get(mdata, i);
	result <- do.call(pedigreemm, param_list);
	tbl <- summary(result)$coef;
	
	sst <- var(result@resp$y) * (length(result@resp$y) - 1);
	ssr <- tbl[,3]^2 * (sum(resid(result)^2) / ..DFE);
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
