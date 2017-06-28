#!/bin/bash
# assoctool 0.0.1
# Association analysis tool
# Version: 0.0.1
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

# Do not forget to do this:
# sudo chmod 666 /var/run/docker.sock
# dx-docker add-to-applet robbyjo/r-mkl-bioconductor:3.4.0 assoctool

main() {

    echo "Value of omics_file: '$omics_file'"
    echo "Value of output_file: '$output_file'"
    echo "Value of save_as_binary: '$save_as_binary'"
    echo "Value of compress_output: '$compress_output'"
    echo "Value of load_all: '$load_all'"

    echo "Value of pheno_file: '$pheno_file'"
    echo "Value of id_col: '$id_col"
    echo "Value of pheno_filter_criteria: '$pheno_filter_criteria'"
    echo "Value of factors_list: '$factors_list'"

    echo "Value of method: '$method'"
    echo "Value of formula: '$formula'"
    echo "Value of tx_fun: '$tx_fun'"
    echo "Value of fn_param_list: '$fn_param_list'"
    echo "Value of preload_code: '$preload_code'"
    echo "Value of prologue_code: '$prologue_code'"
    echo "Value of epilogue_code: '$epilogue_code'"
    echo "Value of analysis_code: '$analysis_code'"
    echo "Value of omics_var_name: '$omics_var_name'"
    echo "Value of result_var_name: '$result_var_name'"
    echo "Value of compute_fdr: '$compute_fdr'"

    echo "Value of annot_file: '$annot_file'"
    echo "Value of annot_marker_id: '$annot_marker_id'"
    echo "Value of annot_filter_criteria: '$annot_filter_criteria'"
    echo "Value of annot_cols: '$annot_cols'"

    echo "Value of pedigree_file: '$pedigree_file'"
    echo "Value of pedigree_type: '$pedigree_type'"
    echo "Value of pedigree_id: '$pedigree_id'"
    echo "Value of pedigree_id_col: '$pedigree_id_col'"
    echo "Value of pedigree_father: '$pedigree_father'"
    echo "Value of pedigree_mother: '$pedigree_mother'"
    echo "Value of pedigree_sex: '$pedigree_sex'"

    echo "Value of num_cores: '$num_cores'"
    echo "Value of block_size: '$block_size'"
    echo "Value of debug: '$debug'"

    mkdir -p /data/
    omics_file_name=`dx ls "$omics_file"`
    echo "Real name of omics file: $omics_file_name"
    pheno_file_name=`dx ls "$pheno_file"`
    echo "Real name of phenotype file: $pheno_file_name"

    PARMS=(--omics_file="/data/${omics_file_name}")
    dx download "$omics_file" -o "/data/${omics_file_name}" &
    PARMS+=(--pheno_file="/data/${pheno_file_name}")
    dx download "$pheno_file" -o "/data/${pheno_file_name}" &

# Required parameters
    PARMS+=(--output_file=results)
    PARMS+=(--id_col="$id_col")
    PARMS+=(--method="$method")
    PARMS+=(--formula="$formula")

# Optional parameters
    if [[ "$save_as_binary" != "" ]] ; then
        PARMS+=(--save_as_binary="$save_as_binary")
    fi
 
    if [[ "$pheno_filter_criteria" != "" ]] ; then
        PARMS+=(--pheno_filter_criteria="$pheno_filter_criteria")
    fi
    if [[ "$factors_list" != "" ]] ; then
        PARMS+=(--factors_list="$factors_list")
    fi

    if [[ "$tx_fun" != "" ]] ; then
        PARMS+=(--tx_fun="$tx_fun")
    fi
    if [[ "$fn_param_list" != "" ]] ; then
        PARMS+=(--fn_param_list="$fn_param_list")
    fi
    if [[ "$preload_code" != "" ]] ; then
        PARMS+=(--preload_code="/data/preload_code.R")
        dx download "$preload_code" -o "/data/preload_code.R" &
    fi
    if [[ "$prologue_code" != "" ]] ; then
        PARMS+=(--prologue_code="/data/prologue_code")
        dx download "$prologue_code" -o "/data/prologue_code.R" &
    fi
    if [[ "$epilogue_code" != "" ]] ; then
        PARMS+=(--epilogue_code="/data/epilogue_code.R")
        dx download "$epilogue_code" -o "/data/epilogue_code.R" &
    fi
    if [[ "$analysis_code" != "" ]] ; then
        PARMS+=(--analysis_code="/data/$analysis_code.R")
        dx download "$analysis_code" -o "/data/analysis_code.R" &
    fi
    if [[ "$omics_var_name" != "" ]] ; then
        PARMS+=(--omics_var_name="$omics_var_name")
    fi
    if [[ "$result_var_name" != "" ]] ; then
        PARMS+=(--result_var_name="$result_var_name")
    fi
    if [[ "$compute_fdr" != "" ]] ; then
        PARMS+=(--compute_fdr="$compute_fdr")
    fi

    if [[ "$annot_file" != "" ]] ; then
    	annot_file_name=`dx ls "$annot_file"`
    	echo "Real name of annotation file: $annot_file_name"
        PARMS+=(--annot_file="/data/${annot_file_name}")
        dx download "$annot_file" -o "/data/${annot_file_name}" &
    fi
    if [[ "$annot_marker_id" != "" ]] ; then
        PARMS+=(--annot_marker_id="$annot_marker_id")
    fi
    if [[ "$annot_filter_criteria" != "" ]] ; then
        PARMS+=(--annot_filter_criteria="$annot_filter_criteria")
    fi
    if [[ "$annot_cols" != "" ]] ; then
        PARMS+=(--annot_cols="$annot_cols")
    fi

    if [[ "$pedigree_file" != "" ]] ; then
    	pedigree_file_name=`dx ls "$pedigree_file"`
    	echo "Real name of pedigree file: $pedigree_file_name"
        PARMS+=(--pedigree_file="/data/${pedigree_file_name}")
        dx download "$pedigree_file" -o "/data/${pedigree_file_name}" &
    fi
    if [[ "$pedigree_id_col" != "" ]] ; then
        PARMS+=(--pedigree_id_col="$pedigree_id_col")
    fi
    if [[ "$pedigree_type" != "" ]] ; then
        PARMS+=(--pedigree_type="$pedigree_type")
    fi
    if [[ "$pedigree_id" != "" ]] ; then
        PARMS+=(--pedigree_id="$pedigree_id")
    fi
    if [[ "$pedigree_father" != "" ]] ; then
        PARMS+=(--pedigree_father="$pedigree_father")
    fi
    if [[ "$pedigree_mother" != "" ]] ; then
        PARMS+=(--pedigree_mother="$pedigree_mother")
    fi
    if [[ "$pedigree_sex" != "" ]] ; then
        PARMS+=(--pedigree_sex="$pedigree_sex")
    fi

    if [[ "$num_cores" != "" ]] ; then
        PARMS+=(--num_cores="$num_cores")
    fi

    if [[ "$block_size" != "" ]] ; then
        PARMS+=(--block_size="$block_size")
    fi

    if [[ "$load_all" != "" ]] ; then
        PARMS+=(--load_all="$load_all")
    fi

    if [[ "$debug" != "" ]] ; then
        PARMS+=(--debug="$debug")
    fi

    echo "Waiting for all file transfers to complete..."
    wait
    sudo chmod o+rw /tmp
    echo "===== Listing /data/ ====="
    ls -R /data/
    echo "===== Listing /data/ ends ====="
    echo "Working directory is:"
    pwd
    echo "Checking phenotype file" 
    if [ -e "/data/${pheno_file_name}" ]
    then
       head -n1 /data/${pheno_file_name}
    else
       echo "The phenofile is not ready"
    fi
 
    echo "Rscript assoctool.R ${PARMS[@]}"
    echo "Running code"
    dx-docker run -v /data/:/data/ robbyjo/r-mkl-bioconductor:3.4.0 /usr/bin/Rscript --vanilla /data/assoctool/assoctool.R "${PARMS[@]}"
    echo "Finished running code"
    results=$(dx upload results --brief)
    dx-jobutil-add-output results "$results" --class=file
    if [[ $compress_output == "GZIP" ]] ; then
       gzip -9 $results
       results="${results}.gz"
    elif [[ $compress_output == "BZ2" ]] ; then
       bzip2 -9 $results
       results="${results}.bz2"
    elif [[ $compress_output == "XZ" ]] ; then
       xz -9 $results
       results="${results}.xz"
    fi
    dx mv ${results} ${output_file}
}
