# Tropical-SVM
This is R code for computing a tropical support vector machine for classification over the tropical projective torus.  We can also apply this code for classifying gene trees into two classes.  A tropical SVM is, like a classical SVM, a tropical hyperplane which maximizes the minimum distance of the data points to the tropical hyperplane to separate the data points.  
## Reference
https://arxiv.org/abs/2003.00677
## Installation
To run this code please unzip the [Data.zip] and put all files inside in your R directory. Also, please make sure that [ape](https://cran.r-project.org/web/packages/ape/index.html), [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html), [parallel](https://CRAN.R-project.org/view=HighPerformanceComputing), [gtools](https://cran.r-project.org/web/packages/gtools/index.html), [ggplot2](https://ggplot2.tidyverse.org) and [e1071](https://cran.r-project.org/web/packages/e1071/index.html) packages are installed.
## Contents
### Data.zip
This folder contains all the data sets with training and testing data set pre-separated(data_matrix_list_15~25%.RData) and all assignments (asgn_1\~4_15\~25%.RData) of each algorithm and each percentage of test data when the maximum accuracy is reached. 
#### Genetree Data/data_matrix_list_15%.RData
When loaded in R, this file will be named "data_matrix_list", which is a list of length 12 with each sublist containing 10 data sets with each created by randomly choosing 15% data points as tesing data set from one overall data set created by [Mesquite](http://www.mesquiteproject.org) at one of 12 different C values. 

#### asgn_1_15%.Rdata
When loaded in R, this file will be named "assignment1", which is a list of length 12 with each element a matrix of 10 rows. Each row is an assignment at which algorithm1 outputs the best classification accuracy for a corresponding testing data set. Other files named starting by "asgn" correspond to the maximum assignment table of other algorithms.
### Algorithm1~4.R
These algorithms are Algorithms 1~4, respectively, in Section 5 in the reference.
### graph producer.R
This R script allows users to reprocude the graphs in section 6 by simplying changing the size of tesing data size. To run this script, please make sure all required packages installed and Rdata put in the directory.
### unbounded.RData
When loaded in R, this file will be named "data", which allows users to verify that the objective function is not an upper bound when C is smaller than 0.00625.
