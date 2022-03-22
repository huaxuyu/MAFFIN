# MAFFIN

## Overview

MAFFIN is an integrated R package for metabolomics sample normalization. MAFFIN main algorithm consists of three modules, **high-quality feature selection, MS signal intensity correction, and maximal density fold change normalization**. This package also implements commonly used normalization methods and normalization evaluation functions.

## Installation

**Check if R package "devtools" is installed.**
```
# If not, install R package "devtools" first
install.packages("devtools")
```

**Install MAFFIN**
```
devtools::install_github("Waddlessss/MAFFIN")
```

## Usage

MAFFIN takes feature intensity table (dataframe) as input, with features in row and samples in column by default. 

Specifically, the column names are sample names. The first row includes sample group names. The input feature intensity table can be prepared in Excel and read into R using

```
# Read the input feature intensity table.
inputTable = read.csv(filename)
```

Your input data table should be similar to this:

<img src='man/figures/ExampleDataTable.PNG' align="left" height="200"/></a>


An example input data table can also be viewed in R:
```
# View the example input data table
View(TestingData)
```





