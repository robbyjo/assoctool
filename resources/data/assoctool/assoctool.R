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

# TODO see Econometrics / Psychometrics
# GAM (mgcv / vgam)
# Nonparametric regression (np / loess)
# Generalized method of moment (gmm)
# Principal component regression (pls package, mvr)
# Extreme bounds analysis? (eba)

# Utility functions
#args <- c("--omics_file=try-gene.txt", "--pheno_file=masterpheno-3rd.txt", "--id_col=cel_files", "--output_file=output.txt",
#	"--pedigree_file=sabre_ped_0407_v1_rj.csv", "--pedigree_id_col=SabreID", "--pedigree_id=sabreid", "--pedigree_father=fid", "--pedigree_mother=mid", "--pedigree_type=kinship1",
#	"--result_var_name=Sex,Age,BMI,Glucose", "--factors_list=Batch_Lump",
#	"--annot_file=HuEx-1_0-st-v2.na33.1.hg19.transcript-core.txt", "--annot_marker_id=transcript_cluster_id", "--annot_cols=transcript_cluster_id,seqname,strand,start,stop,GeneSymbol",
#	"--method=lmm", "--formula='y ~ Sex + Age + BMI + Glucose + WBC_Pred + PLT_Pred + LY_PER_Pred + MO_PER_Pred + EO_PER_Pred + BA_PER_Pred + PC1_Gene + all_probeset_mean + all_probeset_stdev + neg_control_mean + neg_control_stdev + pos_control_mean + pos_control_stdev + all_probeset_rle_mean + all_probeset_mad_residual_mean + mm_mean + Residual_all_probeset_mean_Gene + (1|Batch_Lump)'")
#source("utils.R");
#args <- processArgs(args);

default_code_path <- Sys.getenv("ASSOCTOOL_DIR");
if (is.null(default_code_path) | default_code_path == "") default_code_path <- "/data/assoctool/";
default_block_size <- 10000;

source(paste(default_code_path, "utils.R", sep=""));
args <- processArgs(commandArgs(trailingOnly=TRUE));

# Parameter checking and sanitizing (LONG)
{
	opt <- list();
	# Input and output options
	opt$met_file <- args["omics_file"];
	opt$out_file <- args["output_file"];
	opt$save_as_binary <- processBooleanArg(args["save_as_binary"], "save_as_binary");
	opt$block_size <- processIntegerArg(args["block_size"], "block_size", default_block_size);
	opt$load_all <- processBooleanArg(args["load_all"], "load_all");
	
	# Phenotypic file related options
	opt$pheno_file <- args["pheno_file"];
	opt$id_col <- args["id_col"];
	opt$pheno_filter_criteria <- args["pheno_filter_criteria"];
	opt$factors_list <- args["factors_list"];

	# Method and analysis related options
	opt$analysis_type <- tolower(trim(args["analysis_type"]));
	if (is.na(opt$analysis_type) | opt$analysis_type == "") opt$analysis_type <- "gwas";
	opt$method <- args["method"];
	opt$formula_str <- args["formula"];
	opt$tx_fun_str <- args["tx_fun"];
	opt$fn_param_list <- args["fn_param_list"];
	opt$preload_code <- args["preload_code"];
	opt$prologue_code <- args["prologue_code"];
	opt$epilogue_code <- args["epilogue_code"];
	opt$analysis_code <- args["analysis_code"];
	opt$omics_var_name <- args["omics_var_name"];
	opt$result_var_name <- args["result_var_name"];
	opt$compute_fdr <- processBooleanArg(args["compute_fdr"], "compute_fdr");

	# Annotation related options
	opt$annot_file <- args["annot_file"];
	opt$annot_marker_id <- args["annot_marker_id"];
	opt$annot_filter_criteria <- args["annot_filter_criteria"];
	opt$annot_cols <- args["annot_cols"];

	#GWAS-oriented options
	opt$maf_threshold <- processFloatArg(args["maf_threshold"], "maf_threshold");
	opt$mac_threshold <- processIntegerArg(args["mac_threshold"], "mac_threshold", 1);
	opt$chromosome <- tolower(args["chromosome"]);
	if (!(opt$chromosome %in% c("x", "y", "mt"))) opt$chromosome <- "autosome";
	opt$sex <- args["sex"];
	if (is.na(opt$sex) | opt$sex == "") opt$sex <- "sex";

	# GDS related options
	opt$gds_var_id <- args["gds_var_id"];
	if (is.na(opt$gds_var_id) | opt$gds_var_id == "") opt$gds_var_id <- "variant.id";
	opt$gds_sample_id <- args["gds_sample_id"];
	if (is.na(opt$gds_sample_id) | opt$gds_sample_id == "") opt$gds_sample_id <- "sample.id";
	
	# Pedigree related options
	opt$pedigree_file <- args["pedigree_file"];
	opt$pedigree_id_col <- args["pedigree_id_col"];
	opt$pedigree_type <- args["pedigree_type"];
	opt$pedigree_id <- args["pedigree_id"];
	opt$pedigree_father <- args["pedigree_father"];
	opt$pedigree_mother <- args["pedigree_mother"];
	opt$pedigree_sex <- args["pedigree_sex"];
	
	# Miscellaneous options
	opt$num_cores <- processIntegerArg(args["num_cores"], "num_cores");
	opt$progress_bar <- processBooleanArg(args["progress_bar"], "progress_bar");
	opt$debug <- processBooleanArg(args["debug"], "progress_bar");
	#opt$from <- processIntegerArg(args["from"], "from");
	#opt$to <- processIntegerArg(args["to"], "to");
	
	# Parameter validation / sanity check before any loading takes place
	opt$recognized_analysis <- c("gwas", "epigenome", "transcriptome", "protein", "metabolite", "microbiome", "generic");
	opt$recognized_methods <- c("lm", "lmm", "glm", "glmm", "pedigreemm", "kinship1", "kinship2", "nls", "nlmm", "logistf", "rlm", "polr", "survival", "coxme", "gee", "geeglm", "ordgee", "zeroinfl", "glmnb", "censreg", "truncreg", "betareg", "quantreg", "mlogit", "relogit", "gamlss", "zelig", "custom");
	opt$recognized_pedigree <- c("pedigreemm", "kinship1", "kinship2", "sparse_matrix", "dense_matrix");

	if (is.na(opt$met_file)) stop("Omics file is missing!");
	if (!file.exists(opt$met_file)) stop("Omics file does not exist!");
	if (is.na(opt$pheno_file)) stop("Phenotype file is missing!");
	if (!file.exists(opt$pheno_file)) stop("Phenotype file does not exist!");
	if (is.na(opt$out_file)) stop("Output file is missing!");
	if (!is.na(opt$num_cores) & opt$num_cores < 1) {
		cat ("Detected invalid num_cores. Defaulting to using all available cores!\n");
		opt$num_cores <- NA;
	}
	if (!(opt$analysis_type %in% opt$recognized_analysis)) stop(paste("The only analysis methods are:", paste(opt$recognized_analysis, collapse=", ")));
	if (is.na(opt$id_col)) stop("ID column specification is missing!");
	if (!is.na(opt$pedigree_file)) {
		if (is.na(opt$pedigree_id_col)) opt$pedigree_id_col <- opt$id_col;
		if (is.na(opt$pedigree_type)) {
			cat("Pedigree file is specified, but pedigree type is missing. Assuming pedigreemm.\n");
			opt$pedigree_type <- "pedigreemm";
		} else {
			opt$pedigree_type <- tolower(trim(opt$pedigree_type));
			stopifnot(opt$pedigree_type %in% opt$recognized_pedigree);
		}
		if (is.na(opt$pedigree_id)) stop("You MUST specify the ID column in the pedigree file");
		if (is.na(opt$pedigree_father)) stop("You MUST specify the Father ID column in the pedigree file");
		if (is.na(opt$pedigree_mother)) stop("You MUST specify the Mother ID column in the pedigree file");
	}
	if (!is.na(opt$annot_file)) {
		if (is.na(opt$annot_marker_id)) stop("If you specify annotation file, you MUST mention which column name contains the ID.");
	} else {
		cat("No annotation file is specified.\n")
	}
	if (is.na(opt$omics_var_name)) opt$omics_var_name <- "y";
	if (!is.na(opt$result_var_name)) {
		opt$result_var_name <- trim(unlist(strsplit(opt$result_var_name, ",")));
		opt$result_var_pattern <- paste(opt$result_var_name, collapse="|");
	} else {
		opt$result_var_pattern <- NULL;
	}

	if (is.na(opt$formula_str)) stop("Formula needs to be specified. R formula string.");

	#Formula
	# Strip off first and last quotes
	opt$formula_str <- trim(gsub("^\\\"|\\\"$|^'|'$", "", opt$formula_str));
	opt$fm <- as.formula(opt$formula_str);
	..xx <- grep("Surv\\(", opt$formula_str);
	..is_survival <- ifelse (length(..xx) == 0, FALSE, ..xx > 0);
	rm(..xx);

	if (is.na(opt$method)) { # Method not specified, auto detect
		if (isFixedEffectFormula(opt$fm)) {
			opt$method <- ifelse(..is_survival, "survival", "lm");
		} else {
			opt$method <- ifelse(..is_survival, "coxme", ifelse(is.na(opt$pedigree_type), "lmm", opt$pedigree_type));
		}
	} else {
		opt$method <- tolower(trim(opt$method));
		if (!(opt$method %in% opt$recognized_methods)) stop(paste("The only supported methods are:", paste(opt$recognized_methods, collapse=", ")));
		if (opt$method == "custom") {
			if (is.na(opt$analysis_code)) stop("For custom method, you MUST enter in your custom analysis code!");
			checkAnalysisCode <- function(fn) {
				# Load the analysis code locally (meaning: it will be purged after this function exists)
				# All we want to do at this point is to check the sanity of the analysis code.
				# It will be reloaded globally below
				source(fn, local=TRUE);
				# Make sure function doOne is defined
				if (!isDefined(doOne)) stop("Function doOne is not defined!");
			}
			checkAnalysisCode(opt$analysis_code);
		}
	}
	rm(..is_survival);
	
	# Check if transform function is valid (if any)
	if (!is.na(opt$tx_fun_str)) {
		txFun <- eval(parse(text=opt$tx_fun_str));
	} else {
		txFun <- function(x) x;
	}
} # End of parameter check

suppressMessages(library(data.table));
suppressMessages(library(filematrix));

opt$skip_loading_phenotype <- FALSE;
opt$skip_loading_omics <- FALSE;
opt$skip_loading_pedigree <- FALSE;
opt$skip_loading_annotation <- FALSE;

# Preloading code can override the options above
if (!is.na(opt$preload_code)) {
	cat ("Loading pre-loading code...", opt$preload_code, "\n");
	source(opt$preload_code);
}

if (!opt$skip_loading_phenotype) {
	# Phenotype data loading
	# Assumption: rows = num samples, cols = num phenotypes
	pdata <- data.frame(fread(opt$pheno_file), check.names=FALSE, stringsAsFactors=FALSE);
	cat("Phenotype file has been loaded. Dimension:", dim(pdata), "\n");
	# Double check the IDs in the phenotype data and the methylation data!
	if (!(opt$id_col %in% colnames(pdata))) stop(paste(opt$id_col, "is not in phenotype file!"));
	if (!is.na(opt$pheno_filter_criteria)) {
		pdata <- data.table(pdata);
		pdata <- pdata[eval(parse(text=opt$pheno_filter_criteria))]
		pdata <- data.frame(pdata, check.names=FALSE, stringsAsFactors=FALSE);
		cat("Phenotype file has been filtered. Dimension:", dim(pdata), "\n");
	}

	#Make sure chip, row, and column effects are factors
	if (!is.na(opt$factors_list)) {
		cat("Converting to factors:", opt$factors_list, "\n");
		for (..factor in unlist(strsplit(opt$factors_list, ","))) {
			..factor <- trim(..factor);
			pdata[, ..factor] <- as.factor(pdata[, ..factor]);
		}
		rm(..factor);
	}

	if (opt$analysis_type == "gwas") {
		if (opt$sex %in% colnames(pdata)) {
			..avail_sex <- unique(pdata[, opt$sex]);
			if (any(..avail_sex == "")) {
				..avail_sex[..avail_sex == ""] <- "<empty>";
			}
			cat("Sex is encoded as:", ..avail_sex, "\n");
			..b <- ..avail_sex %in% c("M", "F");
			if (!all(..b)) {
				stop(paste("Some of the sex in column '", opt$sex, "' are not encoded in M or F.", sep=""))
			}
			rm(..avail_sex, ..b);
		} else {
			stop(paste("This is a GWAS analysis and the column '", opt$sex, "' to indicate sex is not in the phenotype file!", sep=""))
		}
	}
}

if (!opt$skip_loading_pedigree) {
	# Load pedigree if any
	if (!is.na(opt$pedigree_file)) {
		cat("Loading pedigree file:", opt$pedigree_file, "\n");
		..fn <- tolower(opt$pedigree_file);
		if (endsWith(..fn, ".rds")) {
			cat("Loading", opt$pedigree_file, "as RDS...\n");
			ped <- readRDS(opt$pedigree_file);
			#if (class(ped) != "matrix") ped <- as.matrix(ped);
		} else if (endsWith(..fn, ".rda") | endsWith(..fn, ".rdata")) {
			cat("Loading", opt$pedigree_file, "as RDa...\n");
			..vv <- load(opt$pedigree_file);
			..ii <- 1;
			..lv <- length(..vv);
			if (..lv == 0) stop(paste(opt$pedigree_file, "contains no data!"));
			if (..lv > 1) {
				cat("NOTE: File", opt$pedigree_file, "contains multiple objects\n");
				..ss <- rep(0, ..lv);
				for (i in 1:..lv) ..ss[i] <- object.size(eval(parse(text=..vv[i])));
				..ii <- which.max(..ss);
				rm(..ss);
			}
			cat("Taking the largest object as the pedigree:", ..vv[..ii], "\n");
			ped <- eval(parse(text=..vv[..ii]));
			rm(list=..vv[..ii]); # We will not delete the other objects
			rm(..vv, ..ii, ..lv);
			#if (class(ped) != "matrix") ped <- as.matrix(ped);
		} else { # Text file
			if (opt$pedigree_type %in% c("sparse_matrix", "dense_matrix")) stop("Sparse / dense matrix pedigrees MUST be entered through rds / rda / RData files!");

			cat("Loading", opt$pedigree_file, "as text. Constructing sibship.\n");
			ped_data <- data.frame(fread(opt$pedigree_file), check.names=FALSE, stringsAsFactors=FALSE);
			cat("Pedigree file has been loaded. Dimension:", dim(ped_data), "\n");
			..kid_ids <- ped_data[, opt$pedigree_id];
			..dad_ids <- ped_data[, opt$pedigree_father];
			..mom_ids <- ped_data[, opt$pedigree_mother];
			
			# Is there anyone with missing pedigree?
			..has_missing_ped <- !(as.character(pdata[, opt$pedigree_id_col]) %in% as.character(..kid_ids));
			if (sum(!..has_missing_ped) == 0) stop("The IDs in the pedigree file does NOT match at all with the IDs in the phenotype file!");
			if (sum(..has_missing_ped) > 0) {
				cat("WARNING: Not all IDs in phenotype file has pedigree data. Assuming singletons.\n");
				..missing_ids <- pdata[..has_missing_ped, opt$pedigree_id_col];
				..kid_ids <- c(..kid_ids, ..missing_ids);
				..dad_ids <- c(..dad_ids, rep(NA, length(..missing_ids)));
				..mom_ids <- c(..mom_ids, rep(NA, length(..missing_ids)));
				rm(..missing_ids);
			}
			rm(..has_missing_ped);
			if (opt$pedigree_type == "kinship1") {
				if (!(opt$method %in% c("lmm", "kinship1", "coxme"))) stop("Pedigree type 'kinship1' is only supported in LMM or COXME!");
				suppressMessages(library(kinship));
				..fam_id <- makefamid(factor(..kid_ids), ..dad_ids, ..mom_ids);
				ped <- makekinship(..fam_id, factor(..kid_ids), ..dad_ids, ..mom_ids);
				opt$method <- "kinship1";
			} else if (opt$pedigree_type == "kinship2") {
				if (!(opt$method %in% c("lmm", "kinship2", "coxme"))) stop("Pedigree type 'kinship2' is only supported in LMM or COXME!");
				suppressMessages(library(kinship2));
				suppressMessages(library(coxme));
				if (is.na(opt$pedigree_sex)) {
					cat("WARNING: pedigree_sex is missing, trying to reconstruct the matrix using makekinship. May be significantly slower.")
					..fam_id <- makefamid(factor(..kid_ids), ..dad_ids, ..mom_ids);
					ped <- makekinship(..fam_id, factor(..kid_ids), ..dad_ids, ..mom_ids);
				} else {
					ped <- kinship(pedigree(factor(..kid_ids), ..dad_ids, ..mom_ids, pdata[, opt$pedigree_sex]));
				}
				ped <- makekinship(factor(..kid_ids), ..dad_ids, ..mom_ids);
				diag(ped) <- 0.5;
				opt$method <- "kinship2";
			} else if (opt$pedigree_type == "pedigreemm") {
				if (!(opt$method %in% c("lmm", "glmm", "pedigreemm"))) stop("Pedigree type 'pedigreemm' is only supported in LMM or GLMM!");
				..ped_struct <- constructPedigree(..kid_ids, ..dad_ids, ..mom_ids);
				ped <- list(new_ids = ..ped_struct[["ped"]]);
				opt$method <- "pedigreemm";
				ped_data <- cbind(ped_data, ..ped_struct[["tbl"]]);
				..temp_ped <- ped_data[as.character(ped_data[,opt$pedigree_id]) %in% as.character(pdata[, opt$pedigree_id_col]), ];
				..temp_ped <- ..temp_ped[match(as.character(pdata[, opt$pedigree_id_col]), as.character(..temp_ped[, opt$pedigree_id])), ];
				pdata <- cbind(pdata, ..temp_ped[, c("new_ids", "fathers_ids", "mothers_ids")]);
				rm(..temp_ped, ..ped_struct);
			}
			rm(..kid_ids, ..dad_ids, ..mom_ids);
		}
	}
}

# Data loading
# Assumption: rows = num markers, cols = num samples
# The file is assumed to be big, so let's use fread to speed up the loading

if (!opt$skip_loading_omics) {
	..fn <- tolower(opt$met_file);
	if (endsWith(..fn, ".rds")) {
		cat("Loading", opt$met_file, "as RDS...\n");
		mdata <- readRDS(opt$met_file);
		if (class(mdata) != "matrix") mdata <- as.matrix(mdata);
	} else if (endsWith(..fn, ".rda") | endsWith(..fn, ".rdata")) {
		cat("Loading", opt$met_file, "as RDa...\n");
		..vv <- load(opt$met_file);
		..ii <- 1;
		..lv <- length(..vv);
		if (..lv == 0) stop(paste(opt$met_file, "contains no data!"));
		if (..lv > 1) {
			cat("NOTE: File", opt$met_file, "contains multiple objects\n");
			..ss <- rep(0, ..lv);
			for (i in 1:..lv) ..ss[i] <- object.size(eval(parse(text=..vv[i])));
			..ii <- which.max(..ss);
			rm(..ss);
		}
		cat("Taking the largest object as mdata:", ..vv[..ii], "\n");
		mdata <- eval(parse(text=..vv[..ii]));
		rm(list=..vv[..ii]); # We will not delete the other objects
		rm(..vv, ..ii, ..lv);
		if (class(mdata) != "matrix") mdata <- as.matrix(mdata);
	} else if (endsWith(..fn, ".gds")) {
		cat("Loading", opt$met_file, "as GDS...\n");
		suppressMessages(library(SeqArray));
		suppressMessages(library(SeqVarTools));
		suppressMessages(library(gdsfmt));
		mdata <- seqOpen(opt$met_file);
	} else if (endsWith(..fn, ".bmat.tar")) {
		cat("Loading", opt$met_file, "as .bmat.tar...\n");
		..out <- untar(opt$met_file);
		if (..out != 0) stop("Opening file failed!");
		rm(..out);
		mdata <- fm.open(gsub(".bmat.tar", "", opt$met_file, ignore.case=TRUE));
	} else {
		# Assume text
		cat("Loading", opt$met_file, "as text...\n");
		mdata <- data.frame(fread(opt$met_file), check.names=FALSE, stringsAsFactors=FALSE);
		rownames(mdata) <- mdata[,1];
		mdata <- mdata[, -1];
		mdata <- as.matrix(mdata);
	}

	# Loading the entire dataset to the memory is not wise
	if (is(mdata, "matrix") & !opt$load_all & NROW(mdata) > opt$block_size) {
		cat("Converting data matrix into filematrix .bmat format...\n");
		mdata <- fm.create.from.matrix("mdata-temp", mdata);
		# Force R to relinquish any unused memory
		while(gc(reset=TRUE)[2,3] != gc(reset=TRUE)[2,3]) {}
	}
	rm(..fn);
	cat("Main file has been loaded. Dimension:", dim(mdata), "\n");
}

if (!opt$skip_loading_annotation) {
	# Annotation loading and filtering (if any)
	..included_marker_ids <- ifelse(is(mdata, "SeqVarGDSClass"), seqGetData(mdata, opt$gds_var_id), rownames(mdata));
	if (!is.na(opt$annot_file)) {
		cat("Loading annotation", opt$annot_file, "\n");
		annot_data <- fread(opt$annot_file);
		cat("Annotation file has been loaded. Dimension:", dim(annot_data), "\n");
		if (!is.na(opt$annot_filter_criteria)) {
			annot_data <- annot_data[eval(parse(text=opt$annot_filter_criteria))]
			cat("Annotation file has been filtered. Dimension:", dim(annot_data), "\n");
		}
		annot_data <- data.frame(annot_data, check.names=FALSE, stringsAsFactors=FALSE);
		..included_marker_ids <- annot_data[, opt$annot_marker_id];
		if (!is.na(opt$annot_cols) & opt$annot_cols != "") {
			..annot_cols <- trim(unlist(strsplit(gsub("^\\\"|\\\"$|^'|'$", "", opt$annot_cols), ",")));
			..annot_cols <- unique(c(..annot_cols, opt$annot_marker_id));
			if (!all(..annot_cols %in% colnames(annot_data))) stop("Some of the annot_cols are not in the annotation columns!");
			annot_data <- annot_data[, ..annot_cols];
			rm(..annot_cols);
		}
	}
}

opt$skip_matching <- FALSE;

get <- function(mat, i) {
	return (txFun(as.numeric(mat[i, ])));
}


# Prologue code can skip matching between phenotype and omics dataset.
# (or use a custom matching rule)
if (!is.na(opt$prologue_code)) {
	cat ("Loading pre-processing code...", opt$prologue_code, "\n");
	#eval(parse(text=prologue_code));
	source(opt$prologue_code);
}

if (!opt$skip_matching) {
	if (is(mdata, "SeqVarGDSClass")) {
		..ids <- intersect(seqGetData(mdata, opt$gds_sample_id), pdata[, opt$id_col]);
		if (length(..ids) == 0) stop("There is no common IDs between those specified in the phenotype data vs. the main data!");
		cat("There are", length(..ids), "IDs in common\n");
		pdata <- pdata[match(..ids, pdata[,opt$id_col]), ];
		cat("There are", length(seqGetData(mdata, opt$gds_var_id)), "number of SNPs in the file.\n");
	} else {
		..ids <- intersect(as.character(colnames(mdata)), as.character(pdata[, opt$id_col]));
		if (length(..ids) == 0) stop("There is no common IDs between those specified in the phenotype data vs. the main data!");
		cat("There are", length(..ids), "IDs in common\n");
		pdata <- pdata[match(..ids, pdata[,opt$id_col]), ];
	}
	rm(..ids);
	cat("Phenotype file has been matched by ID. Dimension:", dim(pdata), "\n");
}

# Force R to relinquish any unused memory
while(gc(reset=TRUE)[2,3] != gc(reset=TRUE)[2,3]) {}

#######################
# MAIN ANALYSIS CODE
#######################

####################
# METHOD SELECTION

# Unfortunately, we cannot execute the following commands since the parameters
# may call functions specific to the model (and thereby the library). So, we store
# the command as strings, which will be executed AFTER the loading of the library
param_cmd <- "param_list <- ifelse(is.na(opt$fn_param_list), list(), paramToList(opt$fn_param_list)); param_list$formula <- opt$fm; param_list$data <-pdata; if (NA %in% names(param_list)) { param_list[[which(NA %in% names(param_list))]] <- NULL; }";
..patterns <- NULL;

# NOTE: If you have a big machine with big RAM, use the parallel version.
# If not, use the serial version.

# Show how many cores your machine has:
library(parallel);
if (is.na(opt$num_cores)) opt$num_cores <- detectCores(logical=TRUE);
cat("The number of cores will be used:", opt$num_cores, "\n");
options(mc.cores = opt$num_cores);

# Show the size of your dataset in MB
cat("Used memory (MB):\n");
(used_mem <- gc(reset=TRUE)[2,2]);

# Load the relevant plugins
if (opt$method == "custom") {
	try(eval(parse(text=param_cmd)));
	source(opt$analysis_code);
	if (!isDefined(doOne)) stop("Function doOne is not defined!");
} else if (opt$method == "lmm" | opt$method == "glmm" | opt$method == "nlmm") {
	source(paste(default_code_path, "lmm.R", sep=""));
} else if (opt$method == "lm" | opt$method == "glm") {
	source(paste(default_code_path, "lm.R", sep=""));
} else {
	source(paste(default_code_path, opt$method, ".R", sep=""));
}

result_all <- c();
..ids <- pdata[,opt$id_col];

if (is(mdata, "SeqVarGDSClass")) {
	..gds <- mdata;
	..num_markers <- length(seqGetData(..gds, opt$gds_var_id));
	..num_blocks <- ceiling(..num_markers / opt$block_size);
	..block_start <- ((1:..num_blocks) - 1) * opt$block_size + 1;
	..block_end   <- ((1:..num_blocks) * opt$block_size);
	..block_end[..num_blocks] <- ..num_markers;
	..males <- pdata[, opt$sex] == "M";
	..females <- !..males;
	if (opt$progress_bar) ..pb <- txtProgressBar(max=..num_blocks, style=3);
	for (..block_no in 1:..num_blocks) {
		if (opt$debug) cat("Block #", ..block_no, "\n");
		seqSetFilter(..gds, variant.sel = ..block_start[..block_no]:..block_end[..block_no], sample.id = ..ids, verbose = FALSE);
		mdata <- t(altDosage(..gds));
		mdata <- mdata[rownames(mdata) %in% ..included_marker_ids, ..ids, drop=FALSE];
		if (opt$analysis_type == "gwas") {
			if (opt$chromosome == "X") {
				..count_f <- rowSums(mdata[, ..females, drop = FALSE], na.rm = TRUE);
				..n_f <- rowSums(!is.na(mdata[, ..females, drop = FALSE]));
				..count_m <- rowSums(mdata[, ..males, drop = FALSE], na.rm = TRUE);
				..n_m <- rowSums(!is.na(mdata[, ..males, drop = FALSE]));
				mac <- pmin(..count_f, 2*..n_f - ..count_f) + pmin(..count_m, ..n_m - ..count_m);
				maf <- mac / (2*..n_f + ..n_m);
				..n <- ..n_f + ..n_m;
				rm(..count_f, ..count_m, ..n_f, ..n_m);
			} else if (opt$chromosome == "Y") {
				..count <- rowSums(mdata[, ..males, drop = FALSE], na.rm = TRUE);
				..n <- rowSums(!is.na(mdata[, ..males, drop = FALSE]));
				mac <- pmin(..count, ..n - ..count);
				maf <- mac / ..n;
				rm(..count);
			} else {
				..count <- rowSums(mdata, na.rm = TRUE);
				..n <- rowSums(!is.na(mdata));
				mac <- pmin(..count, 2*..n - ..count);
				maf <- mac / (2 * ..n);
				rm(..count);
			}
			..b <- mac >= opt$mac_threshold;
			if (!is.na(opt$maf_threshold)) {
				..b <- b & (maf >= opt$maf_threshold);
			}
			..sb <- sum(..b);
			if (..sb < NROW(mdata)) {
				mdata <- mdata[..b, , drop=FALSE];
			}
			rm(..sb);
		}
		# If some SNPs still survive filtering, do the analysis. If not, skip the block altogether.
		if (NROW(mdata) > 0) {
			cur_result <- do.call(rbind, mclapply(1:NROW(mdata), doOne, mc.cores=opt$num_cores));
			if (opt$analysis_type == "gwas") {
				cur_result <- cbind(N=..n, MAC=mac, MAF=maf, cur_result);
			}
			rownames(cur_result) <- rownames(mdata);
			result_all <- rbind(result_all, cur_result);
		}
		if (opt$progress_bar) setTxtProgressBar(..pb, ..block_no);
	}
} else if (is(mdata, "filematrix") | is(mdata, "matrix")) {
	result_all <- do.call(rbind, mclapply(1:NROW(mdata), doOne, mc.cores=opt$num_cores));
	#if (is(..fm, "filematrix")) closeAndDeleteFiles(..fm);
	# TODO: Partial loading for text
}

# Dispose progress bar
if (opt$progress_bar) {
	..pb$kill();
	rm(..pb);
}

# Remove ALL temporary variables
rm(list=ls(all.names=TRUE)[grep("^\\.\\.", ls(all.names=TRUE))]);

if (opt$debug) {
	cat("class(result_all):\n");
	print(class(result_all));
	cat("dim(result_all):\n");
	print(dim(result_all));
	cat("colnames(result_all):\n");
	print(colnames(result_all));
}

for (..col in grep("^P_", colnames(result_all))) {
	cat("Lambda of ", ..col, "=",lambda(result_all[,..col]), "\n");
	# Compute FDR (Benjamini & Hochberg)
	if (opt$compute_fdr) {
		result_all[, paste("FDR", colnames(result_all)[..col],sep="_")] <- p.adjust(result_all[, ..col], "BH");
	}
}

if (!is.na(opt$epilogue_code)) {
	cat ("Loading post-processing code...", opt$prologue_code, "\n");
	#eval(parse(text=epilogue_code));
	source(opt$epilogue_code);
}

# Annotation has been requested
if (!is.na(opt$annot_cols)) {
	annot_data <- annot_data[match(rownames(result_all), annot_data[, opt$annot_marker_id]), ];
	result_all <- data.frame(result_all, check.names=FALSE, stringsAsFactors=FALSE);
	result_all <- cbind(annot_data, result_all);
}

if (!opt$save_as_binary) {
	cat("Saving output as text to", opt$out_file, "\n");
	fwrite(result_all, bzfile(opt$out_file, "w"), row.names=TRUE);
} else {
	cat("Saving output as binary to", opt$out_file, "\n");
	saveRDS(result_all, file=opt$out_file, compress="bz2");
}
