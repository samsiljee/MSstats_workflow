########################################################################################
## Install required packages for Advanced R course section from CRAN and Bioconductor
packages <- c('reshape2','MSstats')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(packages)

########################################################################################
## Tests whether all required packages 
## Press the 'Source' button in the top right corner of this pane and check 
## whether the output in the Console pane confirms that all packages are installed

required_packages <- c(packages)
installed_packages <- required_packages %in% installed.packages()[,"Package"]
missing_packages <- required_packages[!installed_packages]

if ( length(missing_packages) > 0 ) {
    warning(sprintf('FOLLOWING PACKAGES NEED TO BE INSTALLED STILL:\n\t%s\n',
                    paste(missing_packages, collapse=', ')))
} else{
    message('ALL PACKAGES ARE INSTALLED, WE\'RE GOOD TO GO!\n')
}
