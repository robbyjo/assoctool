R_archive=/tmp/R-3.3.0_MKL_all.tar.gz
mkl_archive=/tmp/mkl_redist.tar.gz


.PHONY : all
all: extract_r setup_mkl install_missing create_compile_links

.PHONY : extract_r
extract_r :  ${R_archive}
	cd /opt && sudo tar xzf ${R_archive} 
	

.PHONY : extract_mkl
extract_mkl: ${mkl_archive}  
	cd / && sudo tar xzf ${mkl_archive}


.PHONY : setup_mkl
setup_mkl: extract_mkl
	sudo sh -c "echo /opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64 > /etc/ld.so.conf.d/mkl_share.conf" 
	sudo ldconfig


.PHONY : install_missing
install_missing:
	sudo apt-get update --yes
	sudo apt-get install --no-install-recommends --yes libnetcdfc7


.PHONY : create_compile_links
create_compile_links:
	sudo ln -s /opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64/libmkl_intel_lp64.so /usr/lib/libmkl_intel_lp64.so
	sudo ln -s /opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64/libmkl_core.so /usr/lib/libmkl_core.so
	sudo ln -s /opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64/libmkl_gf_lp64.so /usr/lib/libmkl_gf_lp64.so
	sudo ln -s /opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64/libmkl_gnu_thread.so /usr/lib/libmkl_gnu_thread.so
