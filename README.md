## **Assoctool**

Maintained by Roby Joehanes ([robyjoehanes@hsl.harvard.edu](mailto:robyjoehanes@hsl.harvard.edu); [robbyjo@gmail.com](mailto:robbyjo@gmail.com)) 

Location: Commons_Development:/assoctool_test/ (will be moved once this app is done)

Tool URL: [https://github.com/robbyjo/assoctool](https://github.com/robbyjo/assoctool)

Tags: Statistics

Version: 0.0.1


## **What does this app do?**

Assoctool is a very versatile single marker association analysis app. Current supported models:

1. Linear fixed or mixed effects model with or without pedigree

2. Generalized linear fixed or mixed effects model with or without pedigree with any link function

3. Non-linear fixed or mixed effect model

4. Robust linear fixed or mixed model

5. Survival model

6. Coxme (with or without pedigree)

7. Generalized Estimating Equation (GEE)

8. GEE-GLM

9. Ordinal GEE

10. Zero Inflated models (Poisson or Negative Binomial)

11. Generalized Additive Model (GAM) through GAMLSS

12. Multinomial regression

13. Truncated regression

14. Censored regression

15. Firth's regression

16. Rare event logistic regression

17. Ordinal regression

18. Beta regression

19. Quantile regression

20. Proportional Odds regression

Highly customized formula can be specified using R format. Anything that R can do, it can be done through assoctool.

Once the article for this tool is published, I will put it in.

## **What are typical use cases for this app?**

Single marker association analysis (GWAS) or for transcriptomic (TWAS), epigenomic (EWAS), or microbiome.

## **What data are required for this app to run?**

This app requires genotypic file in GDS format or txt / csv and phenotypic file in csv format.

Optional pedigree matrix can be inputted through csv or RDa format to account for family relatedness or cryptic relatedness. Pedigree matrix can also be constructed through kinship data.

Optional annotation file can be inputted for filtering markers

## **What does this app output?**

This app outputs the effect, standard error, T or Z statistics, p-values, partial R^2 (when applicable) for the marker and all the fixed effect covariates. That is, anything that you see in summary(result) will be collected. This app allows specification which covariates to be output.

## **How does this app work?**

This app runs single marker analysis in parallel.

**Inputs**

<table>
  <tr>
    <td>COMMON (Required)</td>
    <td></td>
  </tr>
  <tr>
    <td><tt>omics_file</tt>  file</td>
    <td>Omics file. For GWAS, this should be a GDS file. For any others (TWAS / EWAS / microbiome), you can use txt / csv / RDa / RDs format</td>
  </tr>
  <tr>
    <td><tt>pheno_file</tt>  file</td>
    <td>Phenotypic file in txt / csv format containing all covariates to use in the model. Must contain an ID column (which is by default assumed to be named "ID" --- this default can be changed; see advanced options) to match with the omics file. Note that any columns not named according to R convention will be automatically renamed (e.g., column <tt>age-orig</tt> will be renamed <tt>age.orig</tt>, as per R convention) </td>
  </tr>
  <tr>
    <td><tt>formula</tt>  string</td>
    <td>The formula for your analysis. By default the marker is called "y" (again, default can be changed in advanced options). Example:
    <ol>
<li><tt>fnbmd ~ y + y*age + age + sex + height</tt><br>
Here is a typical GWAS on femoral neck BMD. The SNP is called y here and this model can account for interaction with age.</li>
<li><tt>y ~ smoking + snp_chrna3 + age + sex + (1|batch)</tt><br>
Here is a typical TWAS/EWAS analysis, where <tt>y</tt> is the expression / methylation levels. The batch term is specifically cast as random effect</ol> </td>
  </tr>
  <tr>
    <td><tt>analysis_type</tt>  combo</td>
    <td>Specify analysis type (one of GWAS, Epigenome, Transcriptome, Protein, Metabolite, Microbiome, Generic). For GWAS analysis, the app will also compute Minimum Allele Frequency information per variant. In version 0.0.1, analyses other than GWAS have no difference.</td>
  </tr>
  <tr>
    <td><tt>method</tt>  combo</td>
    <td>Specify which method to use:<ol>
<li>LM for linear fixed effect model</li>
<li>LMM for linear mixed effect model</li>
<li>GLM for generalized linear fixed effect model (any family except negative binomial)</li>
<li>GLMM for generalized linear mixed effect model (any family)</li>
<li>NLS for non-linear fixed effect model</li>
<li>NLMM for non-linear mixed effect model</li>
<li>PEDIGREEMM for linear mixed effect model with pedigree structure using pedigreemm package</li>
<li>KINSHIP1 for linear mixed effect model with pedigree structure using kinship1 package</li>
<li>KINSHIP2 for linear mixed effect model with pedigree structure using kinship2 package</li>
<li>LOGISTF for Firth's logistic regression</li>
<li>RLM for robust linear model</li>
<li>POLR for proportional odds logistic regression</li>
<li>SURVIVAL for fixed effect survival model</li>
<li>COXME for fixed or mixed effect survival model with or without pedigree</li>
<li>GEE for generalized estimating equation</li>
<li>GEEGLM for GEE using GEEpack</li>
<li>ORDGEE for Ordinal GEE</li>
<li>ZEROINFL for zero inflated models (Poisson / Negative Binomial)</li>
<li>GLMNB for generalized linear fixed effect model on Negative Binomial family</li>
<li>CENSREG for censored regression</li>
<li>TRUNCREG for truncated regression</li>
<li>BETAREG for regression on beta-distributed responses</li>
<li>QUANTREG for regression on quantile-based responses</li>
<li>MLOGIT for multinomial logit regression model</li>
<li>RELOGIT for rare event logistic model</li>
<li>GAMLSS for generalized additive model</li>
<li>ZELIG for any model that you can specify using Zelig package</li>
<li>CUSTOM for custom model. In this case, you MUST specify your analysis code at minimum. See advanced options</li>
</ol></td>
  </tr>
  <tr>
    <td>Common (Optional)</td>
    <td></td>
  </tr>
  <tr>
    <td><tt>pedigree_file</tt>  file</td>
    <td>File containing pedigree information</td>
  </tr>
  <tr>
    <td><tt>pedigree_type</tt>  combo</td>
    <td>Specify what the pedigree format is: whether it is kinship1, kinship2, or pedigree mm format or sparse / dense matrix. In case of matrix format, the pedigree file is expected to be in rds format. In case of the others (kinship1 or 2 or pedigreemm), they are expected to be in text files (txt / csv) and then you need to specify the ID column of the text file for participant, father, and mother. See advanced options.</td>
  </tr>
  <tr>
    <td><tt>annot_file</tt>  file</td>
    <td>File containing annotation information. If you want to filter the markers by annotation (e.g., by position), you can do so by loading annotation file</td>
  </tr>
  <tr>
    <td><tt>annot_marker_id</tt>  string</td>
    <td>The column name in the annotation file containing the ID for the marker (that you can link to your main omics file).</td>
  </tr>
  <tr>
    <td><tt>annot_filter_criteria</tt>  string</td>
    <td>A string containing an R expression to select which marker to include. Example: <tt>Position >= 1000000 && Position <= 2000000</tt>. This means that markers that lie within 1MB to 2MB position will be included in the analysis. Obviously, you will need to have a column called <tt>Position</tt> in your annotation file to make it work.</td>
  </tr>
  <tr>
    <td><tt>annot_cols</tt>  string</td>
    <td>If you wish your result file come pre-annotated using the annotation file, you can specify which columns to include. This is a comma-separated list specifying different columns in the annotation file.</td>
  </tr>
  <tr>
    <td>Common for GWAS analysis only</td>
    <td></td>
  </tr>
  <tr>
    <td><tt>maf_threshold</tt> float</td>
    <td>Minimum Allele Frequency (MAF) threshold under which the marker will not be analyzed. By default it is not specified. Applicable only for GWAS analysis. If you want to analyze only markers with MAF>=5%, then you put 0.05 here.</td>
  </tr>
  <tr>
    <td><tt>mac_threshold</tt> integer</td>
    <td>Minimum Allele Count (MAC) threshold under which the marker will not be analyzed. By default it is set to 1. Applicable only for GWAS analysis. If you want to analyze only markers with MAC>=5, then you put 5 here.</td>
  </tr>
  <tr>
    <td><tt>chromosome</tt> string</td>
    <td>Specify which chromosome the markers are in. Applicable only for GWAS analysis. This information is used to correctly compute MAF. Enter X, Y, or MT for X, Y, or mitochondrial chromosome. Anything else will be treated as autosome. This is not needed if you have a GDS file format because GDS file contains chromosome field built-in (see GDS format section).</td>
  </tr>
  <tr>
    <td><tt>sex</tt> string</td>
    <td>Column name in phenotype file that indicates sex. Only applicable for GWAS analyses; ignored for other analysis types. MUST be encoded as M or F to indicate males and females, respectively (again, only for GWAS case). This is only used to compute MAF / MAC.</td>
  </tr>
  <tr>
    <td><tt>impute_genotype</tt> boolean</td>
    <td>Whether or not missing genotype data should be imputed or not. If set to true, then missing genotype will be set to MAF</td>
  </tr>
  <tr>
    <td>ADVANCED (Essential for pedigree analysis)</td>
    <td></td>
  </tr>
  <tr>
    <td><tt>pedigree_id</tt> string</td>
    <td>The name of column in pedigree file indicating participant ID. By default this is also called "<tt>ID</tt>". Applicable only for kinship / pedigreemm pedigree format. In kinship / pedigreemm format, pedigree file contains at least three columns, the individual ID, father ID, and mother ID, followed optionally by sex information.</td>
  </tr>
  <tr>
    <td><tt>pedigree_father</tt> string</td>
    <td>The name of column in pedigree file indicating the ID of the father of the subject. By default this is called "<tt>father</tt>". Applicable only for kinship / pedigreemm pedigree format.</td>
  </tr>
  <tr>
    <td><tt>pedigree_mother</tt> string</td>
    <td>The name of column in pedigree file indicating the ID of the mother of the subject. By default this is called "<tt>mother</tt>". Applicable only for kinship / pedigreemm pedigree format.</td>
  </tr>
  <tr>
    <td><tt>pedigree_sex</tt> string</td>
    <td>The name of column in pedigree file indicating the sex of the. By default this is called "<tt>sex</tt>". Applicable only for kinship / pedigreemm pedigree format.</td>
  </tr>
  <tr>
    <td><tt>pedigree_id_col</tt> string</td>
    <td>If for some reason the pedigree file is specified using a different set of ID numbers than the one specified in the phenotypic data (a common situation in clinical setting), then you need to provide a linker ID column name in here. In case of Framingham Heart Study cohort, the pedigree is often provided using SabreID, while the genotype data is available in ShareID. The pedigree file contains the <tt>pedigree_id</tt>, <tt>pedigree_father</tt>, and <tt>pedigree_mother</tt> in SabreID, then we add another column ShareID to link to the phenotype file. In this case, we fill in ShareID here.</td>
  </tr>
  <tr>
    <td>ADVANCED (For some analysis / models)</td>
    <td></td>
  </tr>
  <tr>
    <td><tt>fn_param_list</tt> string</td>
    <td>A string containing a list of parameters that you wish to pass to the main analysis function. Example use of this would be to specify the family and link for GLM / GLMM type analysis. If you want to do weighted binomial GLM with logit link, you would put in the following text:<br>
<tt>family=binomial(link="logit"), weights=pdata[, "Weight"]</tt><br>
The pdata is the object in which phenotype data is stored. Certainly, for this example, you will need to make sure that the phenotype file contains a column called <tt>Weight</tt>.
This string will be passed to the glm function (along with formula and the pdata).</td>
  </tr>
  <tr>
    <td><tt>tx_fun</tt> string</td>
    <td>Transformation function, specified in R format. In EWAS / TWAS / Microbiome analysis, sometimes it is necessary to transform the marker due to non conformance to normality. One example of such transformation is the inverse normal transform (INT), which can be specified as:<br>
<tt>function(x) qqnorm(as.numeric(x), plot.it=FALSE)$x</tt><br>
Or, Z-transform, which can be specified as:<br>
<tt>function(x) scale(as.numeric(x), center=TRUE, scale=TRUE)</tt></td>
  </tr>
  <tr>
    <td><tt>id_col</tt> string</td>
    <td>The name of column in phenotype file indicating subject ID. By default this is called "<tt>ID</tt>". If you name it something else (e.g., SabreID), then you should specify it here</td>
  </tr>
  <tr>
    <td><tt>omics_var_name</tt> string</td>
    <td>By default, the omics will be called "<tt>y</tt>" in formula. If you wish to call the omics something else in the formula (for clarity), then you specify it here</td>
  </tr>
  <tr>
    <td><tt>compute_fdr</tt> boolean</td>
    <td>If set to true, then Benjamini & Hochberg's False Discovery Rate (FDR) will be computed for each p-value column.</td>
  </tr>
  <tr>
    <td><tt>result_var_name</tt> string</td>
    <td>A string containing a list of terms you wish to output into the result file. For example, you are running the following model:<br>
<tt>fnbmd ~ y + y*age + age + sex + height</tt><br>
If <tt>result_var_name</tt> is left blank, then assoctool will output the statistics for y (i.e., the marker), y*age, age, sex, and height (i.e., everything). Suppose you are interested in the statistics of the marker and its interaction with age. Then, you would fill in "<tt>y,y:age</tt>" here.</td>
  </tr>
  <tr>
    <td><tt>factors_list</tt> string</td>
    <td>A comma separated list of column names in phenotype file that should be treated as factors, rather than integers. By default, phenotype file will be loaded as data.frame. R usually guesses the data type for each column. If a column, say batch number, is entirely integer, then R will treat it as numeric. In regression equation, treating batch as numeric variable will be wrong. To override that behavior, you specify which columns should be treated as factors. In terms of R code, basically assoctool will iterate the list and apply <tt>as.factor</tt> function to each column.</td>
  </tr>
  <tr>
    <td>ADVANCED (For custom GDS format)</td>
    <td></td>
  </tr>
  <tr>
    <td><tt>gds_var_id</tt> string</td>
    <td>Field name indicating variant ID. By default it is "<tt>variant.id</tt>". If you have a GDS file not prepared by seqVCF2GDS (of seqArray R package), then the resulting GDS file may be structured differently and thus you need to specify which field contains variant ID name.</td>
  </tr>
  <tr>
    <td><tt>gds_sample_id</tt> string</td>
    <td>Field name indicating sample ID. By default it is "<tt>sample.id</tt>".</td>
  </tr>
  <tr>
    <td><tt>gds_chrom_id</tt> string</td>
    <td>Field name indicating chromosome name for each variant. Default to "<tt>chromosome</tt>"</td>
  </tr>
  <tr>
    <td><tt>gds_pos_id</tt> string</td>
    <td>Field name indicating position of each variant. Default to "<tt>position</tt>"</td>
  </tr>
  <tr>
    <td><tt>gds_allele_id</tt> string</td>
    <td>Field name indicating the alleles for each variant. Default to "<tt>allele</tt>"</td>
  </tr>
  <tr>
    <td>ADVANCED (For custom analysis model not covered above)</td>
    <td></td>
  </tr>
  <tr>
    <td><tt>analysis_code</tt> file</td>
    <td>If you choose a CUSTOM analysis method, you will need to upload your own R code in here. This app expects the code to contain a function called <tt>doOne</tt>, with two parameters: <tt>i</tt> and <tt>mdata</tt>. The variable <tt>i</tt> indicates the index of the marker in <tt>mdata</tt> matrix. This function must return the statistics for marker <tt>i</tt>. For an example, you can view assoctool's <tt>lm.R</tt> file.</td>
  </tr>
  <tr>
    <td>preload_code file</td>
    <td>If you have a custom data format or require a custom library, you can put it here. This code will be loaded right BEFORE the any of the data (main, phenotype, pedigree, annotation) are loaded and prepared. If you do choose to load the data yourself, here is the convention for assoctool:<ol>
<li>Omics data is loaded into mdata.</li>
<li>Phenotype data is loaded into pdata, assumed a data frame.</li>
<li>Annotation data is loaded into annot_data, assumed a data frame.</li>
<li>Pedigree data is loaded into ped object.</li>
</ol>
Variables mdata and ped can be of any object AS LONG AS the analysis code understands them. Assoctool's default for mdata is matrix object, and ped can be either matrix or kinship/pedigreemm object.

If you custom-load the data, you can then skip assoctool's default loading data routine by setting appropriate flag:<ol>
<li><tt>opt$skip_loading_omics</tt> --> If set to <tt>TRUE</tt>, assoctool will skip loading omics file.</li>
<li><tt>opt$skip_loading_phenotype</tt> --> If set to <tt>TRUE</tt>, assoctool will skip loading phenotype file.</li>
<li><tt>opt$skip_loading_pedigree</tt> --> If set to <tt>TRUE</tt>, assoctool will skip loading pedigree file.</li>
<li><tt>opt$skip_loading_annotation</tt> --> If set to <tt>TRUE</tt>, assoctool will skip loading annotation file.</li></ol></td>
  </tr>
  <tr>
    <td><tt>prologue_code</tt> file</td>
    <td>If you have a custom library to load or prepare the dataset differently, this is the place. The R code specified here will be run AFTER the loading of all datasets, but BEFORE any analysis run. An example use is to define the get function in assoctool, which is defined as follows:

</tt>get <- function(mat, i) txFun(as.numeric(mat[i, ]));</tt>

This function extracts marker i from omics data mat, assuming that mat is a matrix and that the markers are stored row wise. If you have different arrangement, this is where to override it.

Another use of a prologue code is to define ID matching between omics data and phenotype data. By default, assoctool will attempt to match by ID as specified by <tt>id_col</tt>. However, this process is not foolproof and that you may have duplicate IDs that have to be resolved differently (such as ID joining of sort in repeat measures data, for example). If you define the data matching differently, then you will need to set <tt>opt$skip_matching</tt> to TRUE.</td>
  </tr>
  <tr>
    <td><tt>epilogue_code</tt> file</td>
    <td>If filled, this is an R code that will be executed AFTER analysis is run. A good use case for this is to pare down unused columns, or to use your favorite False Discovery Rate routine on certain p-values.</td>
  </tr>
  <tr>
    <td>ADVANCED </td>
    <td></td>
  </tr>
  <tr>
    <td><tt>compress_output</tt> string</td>
    <td>Set to none if you do not wish to compress the output. You can choose GZIP, BZ2, or XZ compression if you do.</td>
  </tr>
  <tr>
    <td><tt>save_as_binary</tt> boolean</td>
    <td>Save the result file as a binary file (RDS format) that you can later load in R using <tt>readRDS</tt> function. If you choose to save as binary, please disable the compression option above because R automatically compress the binary file.</td>
  </tr>
  <tr>
    <td><tt>from</tt> integer; 
</tt>to</tt> integer</td>
    <td>If you decide to only analyze a subset of markers, you can specify that here. This is useful to distribute the computational load of one file to several servers. Example: Server 1 compute markers 1 to 100000 --> <tt>from=1, to=100000</tt>; Server 2 compute markers 100001 to 200000, etc.</td>
  </tr>
  <tr>
    <td><tt>block_size</tt> integer</td>
    <td>Assoctool tries to load the omics data file per block because the file tends to be very big and do not fit in RAM simultaneously. This parameter specifies how many markers assoctool may load at a time. By default it is set to 2500. If your run returns with an out of memory error, try reducing this block size.</td>
  </tr>
  <tr>
    <td><tt>num_cores</tt> integer</td>
    <td>The number of cores that you want to use for computation. Defaults to the maximum number available in that server.</td>
  </tr>
  <tr>
    <td><tt>load_all</tt> boolean</td>
    <td>If set to true, then assoctool will attempt to load the entirety of omics data to memory. Useful if the data file is relatively small. If this option is set to true, then <tt>block_size</tt> will be ignored.</td>
  </tr>
  <tr>
    <td><tt>debug</tt> integer</td>
    <td>Set 0 to disable debug. Higher number increases verbosity of debug message.</td>
  </tr>
</table>


**Outputs**

<table>
  <tr>
    <td><tt>output_file</tt>  string</td>
    <td>The desired name for your output file. This output will contain all the statistics of all markers desired.</td>
  </tr>
</table>


### **To run this app**

Click on Projects and select a project in which you wish to run this app. Then click the Start Analysis button and select this app.

### **To run this app from the command line**

Example use case:

```
$ dx run Commons_Development:/assoctool_test/assoctool --watch --yes \
--destination Commons_Development:/assoctool_test/ --instance-type mem3_ssd1_x32 \
-iomics_file=cdbg:/Genotypes/freeze.4/GDS/freeze4.chr21.pass.gtonly.minDP10.genotypes.gds \
-ipheno_file=Commons_Development:/assoctool_test/dummy-pheno.txt \
-isex=Sex -ichromosome=autosome -iblock_size=2000 \
-ioutput_file=chr21-dummy-results.txt.bz2 -icompress_output=BZ2 -idebug=1 \
-iid_col="ID" -iformula="Pheno1~SNP+Sex+Age+Pheno2" \
-imethod="LM" -iomics_var_name="SNP" -iresult_var_name="SNP,Age"
```
