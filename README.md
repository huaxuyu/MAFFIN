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

MAFFIN takes feature intensity table (dataframe) as input with features in row and samples in column by default. Specifically, the column names are recognized as sample names, and the first row is recognized as sample group names. 
