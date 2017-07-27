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
# dx-docker add-to-applet robbyjo/r-mkl-bioconductor:3.4.1 assoctool

main() {
    if [ $debug == "true" ] ; then
       echo "Setting up profiling"
       DEBIAN_FRONTEND=noninteractive apt-get update
       DEBIAN_FRONTEND=noninteractive apt-get -y install sysstat
       sed -e "s/false/true/g" /etc/default/sysstat > /etc/default/sysstat.bak
       mv /etc/default/sysstat.bak /etc/default/sysstat
       /etc/init.d/sysstat start
    fi
    echo "Value of omics_file: '$omics_file'"
    echo "Value of output_file: '$output_file'"
    echo "Value of save_as_binary: '$save_as_binary'"
    echo "Value of compress_output: '$compress_output'"
    echo "Value of load_all: '$load_all'"

    echo "Value of pheno_file: '$pheno_file'"
    echo "Value of id_col: '$id_col"
    echo "Value of pheno_filter_criteria: '$pheno_filter_criteria'"
    echo "Value of factors_list: '$factors_list'"

    echo "Value of analysis_type: '$analysis_type'"
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

    echo "Value of maf_thershold: '$maf_thershold'"
    echo "Value of mac_thershold: '$mac_thershold'"
    echo "Value of chromosome: '$chromosome'"
    echo "Value of sex: '$sex'"
    echo "Value of gds_var_id: '$gds_var_id'"
    echo "Value of gds_sample_id: '$gds_sample_id'"
    echo "Value of gds_chrom_id: '$gds_chrom_id'"
    echo "Value of gds_pos_id: '$gds_pos_id'"
    echo "Value of gds_allele_id: '$gds_allele_id'"

    echo "Value of pedigree_file: '$pedigree_file'"
    echo "Value of pedigree_type: '$pedigree_type'"
    echo "Value of pedigree_id: '$pedigree_id'"
    echo "Value of pedigree_id_col: '$pedigree_id_col'"
    echo "Value of pedigree_father: '$pedigree_father'"
    echo "Value of pedigree_mother: '$pedigree_mother'"
    echo "Value of pedigree_sex: '$pedigree_sex'"

    #echo "Value of from: '$from'"
    #echo "Value of to: '$to'"
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
    PARMS+=(--output_file=/data/results)
    PARMS+=(--id_col="$id_col")
    PARMS+=(--method="$method")
    PARMS+=(--analysis_type="$analysis_type")
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
        preload_code_file_name=`dx ls "$preload_code"`
        echo "Real name of preload code file: $preload_code_file_name"
        PARMS+=(--preload_code="/data/${preload_code_file_name}")
        dx download "$preload_code" -o "/data/${preload_code_file_name}" &
    fi
    if [[ "$prologue_code" != "" ]] ; then
        prologue_code_file_name=`dx ls "$prologue_code"`
        echo "Real name of prologue code file: $prologue_code_file_name"
        PARMS+=(--prologue_code="/data/${prologue_code_file_name}")
        dx download "$prologue_code" -o "/data/${prologue_code_file_name}" &
    fi
    if [[ "$epilogue_code" != "" ]] ; then
        epilogue_code_file_name=`dx ls "$epilogue_code"`
        echo "Real name of epilogue code file: $epilogue_code_file_name"
        PARMS+=(--epilogue_code="/data/${epilogue_code_file_name}")
        dx download "$epilogue_code" -o "/data/${epilogue_code_file_name}" &
    fi
    if [[ "$analysis_code" != "" ]] ; then
        analysis_code_file_name=`dx ls "$analysis_code"`
        echo "Real name of analysis code file: $analysis_code_file_name"
        PARMS+=(--analysis_code="/data/${analysis_code_file_name}")
        dx download "$analysis_code" -o "/data/${analysis_code_file_name}" &
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
    if [[ "$maf_threshold" != "" ]] ; then
        PARMS+=(--maf_threshold="$maf_threshold")
    fi
    if [[ "$mac_threshold" != "" ]] ; then
        PARMS+=(--mac_threshold="$mac_threshold")
    fi

    if [[ "$gds_var_id" != "" ]] ; then
        PARMS+=(--gds_var_id="$gds_var_id")
    fi
    if [[ "$gds_sample_id" != "" ]] ; then
        PARMS+=(--gds_sample_id="$gds_sample_id")
    fi
    if [[ "$gds_chrom_id" != "" ]] ; then
        PARMS+=(--gds_chrom_id="$gds_chrom_id")
    fi
    if [[ "$gds_pos_id" != "" ]] ; then
        PARMS+=(--gds_pos_id="$gds_pos_id")
    fi
    if [[ "$gds_allele_id" != "" ]] ; then
        PARMS+=(--gds_allele_id="$gds_allele_id")
    fi
    if [[ "$chromosome" != "" ]] ; then
        PARMS+=(--chromosome="$chromosome")
    fi
    if [[ "$sex" != "" ]] ; then
        PARMS+=(--sex="$sex")
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

    #if [[ "$from" != "" ]] ; then
    #    PARMS+=(--from="$from")
    #fi
    #if [[ "$to" != "" ]] ; then
    #    PARMS+=(--to="$to")
    #fi

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

    if [[ "$progress_bar" != "" ]] ; then
        PARMS+=(--progress_bar="$progress_bar")
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

    if [ $debug == "true" ] ; then
       echo "Running CPU benchmark code... (sar)"
       sar -u 60 > /data/cpu_benchmark.prof &
       echo "Running RAM benchmark code... (vmstat)"
       vmstat -Sm 60 /data/ram_benchmark.prof & 
    fi
    echo '#!/bin/bash' > /data/runme.sh
    echo 'export MKL_NUM_THREADS=1' >> /data/runme.sh
    x="Rscript assoctool.R ${PARMS[@]}"
    echo "echo \"$x\"" >> /data/runme.sh
    echo 'echo "Running code"' >> /data/runme.sh
    x="/usr/bin/Rscript --vanilla /data/assoctool/assoctool.R ${PARMS[@]}"
    echo "$x" >> /data/runme.sh
    echo 'echo "Finished running code"' >> /data/runme.sh
    chmod 700 /data/runme.sh

    echo "============================ Script begin ============================ "
    cat /data/runme.sh
    echo "============================= Script end ============================= "

    dx-docker run -v /data/:/data/ robbyjo/r-mkl-bioconductor:3.4.1 /bin/bash /data/runme.sh

    if [ $debug == "true" ] ; then
       echo "CPU benchmark results:"
       echo "=========================== BEGIN ==========================="
       cat /data/cpu_benchmark.prof
       echo "============================ END ============================"
       pkill sar
       echo "RAM benchmark results:"
       echo "=========================== BEGIN ==========================="
       cat /data/ram_benchmark.prof
       echo "============================ END ============================"
       pkill vmstat
    fi
    results="/data/results"
    if [ -e $results ] ; then
       if [[ $compress_output == "GZIP" ]] ; then
          echo "Gzipping result file..."
          gzip -9 $results
          results="/data/results.gz"
       elif [[ $compress_output == "BZ2" ]] ; then
          echo "Bzipping result file..."
          bzip2 -9 $results
          results="/data/results.bz2"
       elif [[ $compress_output == "XZ" ]] ; then
          echo "Xzipping result file..."
          xz -9 $results
          results="/data/results.xz"
       fi

       results=$(dx upload ${results} --brief)
       echo "Uploaded results: '$results'"
       dx-jobutil-add-output results "$results" --class=file
       if [ $debug == "true" ] ; then
          echo "Working directory is:"
          pwd
          echo "===== Listing pwd ====="
          ls -R
          echo "===== Listing pwd ends ====="
       fi
       echo "Moving results to ${output_file}"
       dx mv ${results} ${output_file}
    else
       echo "No result file is detected"
    fi
}
