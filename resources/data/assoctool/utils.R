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


# Because there are hooks in prologue / epilogue and custom analysis code, it is VERY important
# to establish a good variable naming convention. 
# Variable conventions:
# Temporary variables start with double dots
# mdata holds main data
# pdata holds phenotype data
# ped_data holds pedigree data (if any)
# ped holds the pedigree object (if any). The type depends on which package is requested.
# annot_data holds annotation data (if any)
# All function names must be in camelCase
# All matrices must have data suffix

# args holds the raw options passed into this program and will be NOT deleted once known options have been parsed.
# This will allow custom programs to pass parameters through the command line.

# Utility functions
trim <- function (x) unlist(lapply(as.list(x), function(v) gsub("^\\s+|\\s+$", "", v)));
isDefined <- function(object) exists(as.character(substitute(object)));
isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol;
isFixedEffectFormula <- function(formula) {
	vars <- grep("\\|", attr(terms(formula), "term.labels"));
	return(length(vars) == 0);
}
paramToList <- function(param_str) {
	ll <- eval(parse(text = paste("list(", param_str, ")")));
	return(ll);
}
processBooleanArg <- function(b, var_name, default=FALSE) {
	if (is.na(b)) return (default);
	if (!is.logical(b)) {
		if (is.numeric(b)) {
			return (b != 0);
		} else if (is.character(b)) {
			if (b == "TRUE" | b == "T") {
				return (TRUE);
			} else if (b == "FALSE" | b == "F") {
				return(FALSE);
			}
			if (tolower(gsub("^-*", "", b)) == tolower(var_name)) return(TRUE);
		}
	}
	cat(paste("Unknown characters in ", var_name, "!", sep=""), "\n");
	return (default);
}
processIntegerArg <- function(v, var_name, default=NA) {
	if (is.numeric(v)) {
		if (isWholeNumber(v)) return (v);
		return (default);
	}
	if (is.na(v)) return(v);
	cat(paste("Unknown characters in ", var_name, "!", sep=""), "\n");
	return (default);
}

# Very simple command line argument parsing
processArgs <- function(args, show=TRUE) {
	stopifnot (length(args) > 0);
	args <- do.call(rbind, regmatches(args, regexpr("=", args), invert = TRUE)); # split by the first equal sign only
	stopifnot (NCOL(args) == 2);
	cat("This program is called with the following parameters:\n");
	colnames(args) <- c("Parameter", "Value");
	args[, 1] <- gsub("^-*", "", args[, 1]); # If option is prefixed by dashes, remove the dashes
	if (show) {
		x <- data.frame(args);
		ff <- x[x[, 1] == "formula", 2];
		x <- x[x[, 1] != "formula", ];
		print(data.frame(x), row.names=FALSE);
		cat(paste("formula =", ff), "\n");
		rm(x, ff);
	}
	rownames(args) <- args[, 1];
	args <- args[, 2];
	return (args);
}
lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE), na.rm=TRUE) / qchisq(0.5, df=1);

convertToLMFormula <- function(ff) {
	form <- ff;
	# The rest of the code in this function is taken from lme4
	nobars <- function(term) {
		if (!('|' %in% all.names(term))) return(term);
		#if (is.call(term) && term[[1]] == as.name('|')) return(NULL);
		if (is.call(term) && term[[1]] == as.name('|')) return(nobars(term[[3]]));
		if (length(term) == 2) {
			nb <- nobars(term[[2]]);
			if (is.null(nb)) return(NULL);
			term[[2]] <- nb;
			return(term);
		}
		nb2 <- nobars(term[[2]]);
		nb3 <- nobars(term[[3]]);
		if (is.null(nb2)) return(nb3);
		if (is.null(nb3)) return(nb2);
		term[[2]] <- nb2;
		term[[3]] <- nb3;
		term;
	}
	form[[3]] <- if (is.null(nb <- nobars(form[[3]]))) 1 else nb;
	return(form);
}
combineFormulas <- function(f1, f2, op="+") {
	ff = as.character(f1);
	return(formula(paste(ff[[2]], ff[[1]], ff[[3]], op, paste(deparse(f2[[3]]), collapse=""))));
}

# Pedigree construction; for pedigreemm only
constructPedigree <- function(all_ids, dad_ids, mom_ids) {
	library(pedigreemm);
	n <- length(all_ids);
	stopifnot(n == length(dad_ids), n == length(mom_ids));
	
	# Start from root. Roots are also Dad / Mom Ids that are not in all_ids.
	root <- unique(setdiff(c(dad_ids, mom_ids), all_ids));
	
	# Also, root has no parents.
	unique_ids <- unique(all_ids[is.na(dad_ids) & is.na(mom_ids)]);
	unique_ids <- unique(c(unique_ids, root));
	unique_ids <- sort(unique_ids[!is.na(unique_ids)]);
	
	# Successively find and add the IDs of children of the roots in a loop.
	# No incestuous relationship is assumed. Violation will break this algorithm.
	cur_gen_ids <- unique_ids;
	
	while (TRUE) {
		cur_kid_ids <- unique(all_ids[(dad_ids %in% cur_gen_ids) | (mom_ids %in% cur_gen_ids)]);
		if (length(cur_kid_ids) == 0) break;
		# If there's a skip in a generation, put the younger generation later in the sequence
		# Example:
		# Grandpa + Grandma
		#         |
		#        Dad + Mom
		#            |
		#           Kid
		# At the first iteration, grandpa, grandma, and mom are selected and
		# assumed to be in the same generation.
		# At the second iteration, dad and kid are selected as the second generation
		# At the third iteration, only kid is selected.
		# Since kid has been selected as the second generation, it has to be removed from
		# the queue and has to be assigned as the third generation
		unique_ids <- setdiff(unique_ids, cur_kid_ids);
		cur_kid_ids <- sort(cur_kid_ids[!is.na(cur_kid_ids)]);
		cur_gen_ids <- cur_kid_ids;
		unique_ids <- c(unique_ids, cur_kid_ids);
	}
	
	dads <- match(dad_ids, unique_ids);
	moms <- match(mom_ids, unique_ids);
	kids <- match(all_ids, unique_ids);
	origtbl <- cbind(kids, dads, moms);
	colnames(origtbl) <- c("new_ids", "fathers_ids", "mothers_ids");
	
	non_kids <- setdiff(1:length(unique_ids), kids);
	kids <- c(non_kids, kids);
	dads <- c(rep(NA, length(non_kids)), dads);
	moms <- c(rep(NA, length(non_kids)), moms);
	tbl <- cbind(kids, dads, moms);
	tbl <- tbl[order(tbl[,"kids"]),];
	rm(n, root, unique_ids, cur_gen_ids, cur_kid_ids, dads, moms, kids, non_kids);
	
	return (list(ped = pedigree(sire=tbl[,"dads"], dam=tbl[,"moms"], label=tbl[,"kids"]), tbl=origtbl));
}

autoDetectTar <- function() {
	if (Sys.getenv("tar") != "") return(Sys.getenv("tar"));
	tar_fn <- "";
	if (.Platform$OS.type == "unix") {
		if (file.exists("/bin/tar")) {
			tar_fn <- "/bin/tar";
		} else if (file.exists("/usr/bin/tar")) {
			tar_fn <- "/usr/bin/tar";
		}
	} else if (.Platform$OS.type == "windows") {
		if (file.exists("C:/cygwin64/bin/tar.exe")) {
			tar_fn <- "C:/cygwin64/bin/tar.exe";
		} else if (file.exists("C:/cygwin/bin/tar.exe")) {
			tar_fn <- "C:/cygwin/bin/tar.exe";
		}
	}
	if (tar_fn == "") stop("Automatic detection of tar failed. You need to set Sys.setenv(\"tar\"). For Windows users, install Cygwin.");
	Sys.setenv(tar = tar_fn);
	return(tar_fn);
}

fixCygwinFilename <- function(fn) {
	if (grepl("\\:\\\\", fn)) fn <- paste("/cygdrive/", gsub("\\:\\\\", "/", fn), sep="");
	if (grepl("\\:/", fn)) fn <- paste("/cygdrive/", gsub("\\:/", "/", fn), sep="");
	if (grepl("\\\\", fn)) fn <- gsub("\\\\", "/", fn);
	return (fn);
}

convertTextFileToBMat <- function(filename, prefix, delimiter="\t", sliceSize=10000, ...) {
	require(filematrix);
	fm <- fm.create.from.text.file(filename, prefix, delimiter=delimiter, sliceSize=sliceSize, ...);
	return(fm);
}

convertMatrixToBMat <- function(filename, prefix, delimiter="\t", sliceSize=10000, ...) {
	require(filematrix);
	fm <- fm.create.from.text.file(filename, prefix, delimiter=delimiter, sliceSize=sliceSize, ...);
	return(fm);
}

convertToBMatTar <- function(filename, prefix, delimiter="\t", sliceSize=10000, ...) {
	require(filematrix);
	tar_fn <- autoDetectTar();
	fm <- fm.create.from.text.file(filename, prefix, delimiter=delimiter, sliceSize=sliceSize, ...);
	files <- c(fm$info.filename, fm$data.filename, fm$rnamefile, fm$cnamefile);
	b <- unlist(lapply(files, function(x) file.exists(x)));
	files <- files[b];
	fn <- paste(fm$data.filename, ".tar", sep="");
	dd <- dirname(prefix);
	if (.Platform$OS.type == "windows") {
		fn <- fixCygwinFilename(fn);
		files <- unlist(lapply(files, fixCygwinFilename));
		dd <- fixCygwinFilename(dd);
	}
	#tar(tarfile = fn, files = files, compression="none", tar=tar_fn);
	files <- unlist(lapply(files, function(x) gsub(dd, "", x)));
	out <- system(command = paste(tar_fn, "cf", fn, "-C", dd, paste(files, collapse=" ")));
	if(out == 0) closeAndDeleteFiles(fm);
	return(out);
}

untar <- function(filename) {
	tar_fn <- autoDetectTar();
	filename <- fixCygwinFilename(filename);
	out <- system(command = paste(tar_fn, "xf", filename));
	return(out);
}

