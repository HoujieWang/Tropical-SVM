# Tropical-SVM
This is R code for computing a tropical support vector machine SVM) for classification over the tropical projective torus.  We can also apply this code for classifying gene trees into two classes.  A tropical SVM is, like a classical SVM, a tropical hyperplane which maximizes the minimum distance of the data points to the tropical hyperplane to separate the data points.  
## Reference
https://arxiv.org/abs/2003.00677
## Installation
To run this code we need the [ape](https://cran.r-project.org/web/packages/ape/index.html), [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html), [parallel](https://CRAN.R-project.org/view=HighPerformanceComputing), [gtools](https://cran.r-project.org/web/packages/gtools/index.html), [ggplot2](https://ggplot2.tidyverse.org) and [e1071](https://cran.r-project.org/web/packages/e1071/index.html) packages.
## Contents
### 15% test data
This folder contains the data sets with training and 15% testing data pre-separated(data_matrix_list_15%.RData) and all assignments (asgn_1~4_15%.RData) of each algorithm when the maximum accuracy is reached. Folder of 20% and 25% contain silimar files.
#### data_matrix_list_15%.RData
When loaded in R, this file will be named "data_matrix_list", which is a list of length 12 with each sublist containing 10 data sets with each created by randomly choosing 15% data points as tesing set from one overall data set created by [Mesquite](http://www.mesquiteproject.org) at one of 12 different C values. 
#### asgn_1_15%.Rdata
When loaded in R, this file will be named "assignment1", which is a list of length 12 with each element a matrix of 10 rows. Each row is an assignment at which algorithm1 outputs the best classification accuracy for a corresponding testing data set. Other files named starting by "asgn" correspond to the maximum assignment table of other algorithms.
### Algorithm1~4.R
These algorithms are Algorithms 1~4, respectively, in Section 5 in the reference.
### graph reproducer.R
This R script allows users to reprocude the graphs in section 6 by simplying chaning the size of tesing data size.

