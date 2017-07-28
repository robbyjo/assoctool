# QQ plot tool
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


# Variable conventions:
# Temporary variables start with double dots
# mdata holds main data
# input_cols holds all the P-value columns one wants to plot
# All function names must be in camelCase
# All matrices must have data suffix

# args holds the raw options passed into this program and will be NOT deleted once known options have been parsed.
# This will allow custom programs to pass parameters through the command line.

default_code_path <- Sys.getenv("ASSOCTOOL_DIR");
if (is.null(default_code_path) | default_code_path == "") default_code_path <- "/data/assoctool/";
source(paste(default_code_path, "utils.R", sep=""));
args <- processArgs(commandArgs(trailingOnly=TRUE));

{
	opt <- list();
	opt$input_file <- args["input_file"];
	opt$input_cols <- args["input_cols"];
	opt$qq_file <- args["output_file"];
	opt$qq_fmt <- args["output_format"];
	opt$qq_param_list <- args["param_list"];
	opt$qq_fmt_param_list <- args["format_param_list"];
	opt$qq_p_threshold <- processFloatArg(args["qq_p_threshold"], "qq_p_threshold", log10(5e-8));

	opt$recognized_formats <- c("png", "pdf", "tiff", "bmp", "jpeg");
	if (is.na(opt$input_file)) stop("Input file is missing!");
	if (is.na(opt$input_cols)) stop("Input column specification is missing!");
	if (is.na(opt$qq_file)) stop("No Q-Q file is specified!");
	if (is.na(opt$qq_fmt)) {
		cat("Q-Q format missing, defaulting to png.\n");
		opt$qq_fmt <- "png";
	} else {
		opt$qq_fmt <- tolower(trim(opt$qq_fmt));
		stopifnot(opt$qq_fmt %in% opt$recognized_formats);
	}
}

# TODO option for multiple page format
opt$input_cols <- trim(unlist(strsplit(opt$input_cols, ",")));

library(data.table);

# Input file loading
# Assumption: rows = num samples, cols = num phenotypes
..fn <- tolower(opt$input_file);
if (endsWith(..fn, ".rds")) {
	cat("Loading", opt$input_file, "as RDS...\n");
	mdata <- readRDS(opt$input_file);
	if (class(mdata) != "data.frame") mdata <- data.frame(mdata, check.names=FALSE, stringsAsFactors=FALSE);
} else if (endsWith(..fn, ".rda") | endsWith(..fn, ".rdata")) {
	cat("Loading", ..fn, "as RDa...\n");
	..vv <- load(opt$input_file);
	..ii <- 1;
	..lv <- length(..vv);
	if (..lv == 0) stop(paste(opt$input_file, "contains no data!"));
	if (..lv > 1) {
		cat("NOTE: File", opt$input_file, "contains multiple objects\n");
		..ss <- rep(0, ..lv);
		for (i in 1:..lv) ..ss[i] <- object.size(eval(parse(text=..vv[i])));
		..ii <- which.max(..ss);
		rm(..ss);
	}
	cat("Taking the largest object as mdata:", ..vv[..ii], "\n");
	mdata <- eval(parse(text=..vv[..ii]));
	rm(list=..vv[..ii]); # We will not delete the other objects
	rm(..vv, ..ii, ..lv);
	if (class(mdata) != "data.frame") mdata <- data.frame(mdata, check.names=FALSE, stringsAsFactors=FALSE);
} else if (endsWith(..fn, ".bz2")) {
	cat("Loading", opt$input_file, "as bzipped text...\n");
	mdata <- read.csv(bzfile(opt$input_file), check.names=FALSE, stringsAsFactors=FALSE);
} else if (endsWith(..fn, ".gz")) {
	cat("Loading", opt$input_file, "as gzipped text...\n");
	mdata <- read.csv(gzfile(opt$input_file), check.names=FALSE, stringsAsFactors=FALSE);
} else if (endsWith(..fn, ".xz")) {
	cat("Loading", opt$input_file, "as xzipped text...\n");
	mdata <- read.csv(xzfile(opt$input_file), check.names=FALSE, stringsAsFactors=FALSE);
} else if (endsWith(..fn, ".zip")) {
	cat("Loading", opt$input_file, "as zipped text...\n");
	mdata <- read.csv(unz(opt$input_file), check.names=FALSE, stringsAsFactors=FALSE);
} else {
	# Assume text
	cat("Loading", opt$input_file, "as text...\n");
	mdata <- data.frame(fread(opt$input_file), check.names=FALSE, stringsAsFactors=FALSE);
}
rm(..fn);
cat("Input file has been loaded. Dimension:", dim(mdata), "\n");

if (!all(opt$input_cols %in% colnames(mdata))) {
	cat("Some of the input columns are missing": opt$input_cols[!(opt$input_cols %in% colnames(mdata))], "\n");
}

## Q-Q plot

qqPlot <- function(pvector, p0 = opt$qq_p_threshold, col=c("#A0A0A0", "#000000"), ...) {
	p_order <- order(pvector,decreasing=FALSE);
	if (any(pvector == 0)) {
		pvector[pvector == 0] <- .Machine$double.xmin;
	}
	o <- -log10(pvector[p_order]);
	n <- length(o);
	e <- -log10( 1:n/n );
	b <- o >= p0;
	
	# Make sure that the dots are not too crowded
	ob <- duplicated(round(o, digits=2)) & !b;
	o <- o[!ob];
	e <- e[!ob];
	b <- b[!ob];
	
	plot(e[!b],o[!b],pch=19,cex=0.7, ...,
			xlab=expression(Expected~~-log[10](italic(p))),
			ylab=expression(Observed~~-log[10](italic(p))),
			xlim=c(0,e[1]), ylim=c(0,o[1]),col=col[1]);
	points(e[b],o[b],pch=19,cex=0.7,col=col[2]);
	abline(a=0,b=1,col=rgb(1,0.65,0),lty=1);
	#abline(lm(o~e), col=rgb(0.5,0,0),lty=2);
}

for (colname in opt$input_cols) {
	cat("Column", colname, "\n");
	if (!(colname %in% colnames(mdata))) {
		cat("not found in the data, skipping...\n");
		next;
	}
	if (class(mdata[, colname]) != "numeric") stop("The column MUST be numeric!");
	if (any(mdata[, colname] < 0 | mdata[, colname] > 1)) stop("Not a valid p-value column. All values must be between 0 and 1!");

	# Compute lambda:
	..ll <- lambda(mdata[,colname]);
	cat("Lambda is:", ..ll, "\n");

	qq_fmt <- match.fun(opt$qq_fmt);
	if (!is.na(opt$qq_fmt_param_list)) {
		..param_list <- paramToList(opt$qq_fmt_param_list);
		..param_list$filename <- paste(opt$qq_file, "-", colname, sep="");
		do.call(qq_fmt, ..param_list);
		rm(..param_list);
	} else {
		qq_fmt(paste(opt$qq_file, "-", colname, sep=""));
	}
	if (!is.na(opt$qq_param_list)) {
		..param_list <- paramToList(opt$qq_param_list);
		..param_list$pvector <- mdata[, colname];
		do.call(qqPlot, ..param_list);
		rm(..param_list);
	} else {
		qqPlot(mdata[,colname]);
	}
	legend("topleft", legend=paste("lambda =", round(..ll, digits=4)));
	dev.off()
}

