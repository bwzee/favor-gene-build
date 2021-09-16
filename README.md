# favor-gene-build
Simple automation of FAVOR GENE


## Installation

Currently this code is implemented in the R language

```
R version 3.5.0 (2018-04-23) -- "Joy in Playing"

library(devtools)

devtools::install_github("bwzee/favor-gene-build")

```


## Running the Build Process

A number of input files are used in this process. These need to be downloaded first.


```

cd src/


#Run the build process testing
build_main_genedb(testing=T,samplesize=10000)


```
