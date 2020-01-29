# Tropical-SVM
This is R code for computing a tropical support vector machine SVM) for classification over the tropical projective torus.  We can also apply this code for classifying gene trees into two classes.  A tropical SVM is, like a classical SVM, a tropical hyperplane which maximizes the minimum distance of the data points to the tropical hyperplane to separate the data points.  
## Reference
Our arxiv link will be posted here (TBA)
## Installation
To run this code we need the ape package.... 
## Contents
### Algorithm1~5.R
These algorithm are Algorithms 1~5, respectively, in Section 5 in the reference.
### Sample.R
This codes provide sample codes to demonstrate the accuracy of our tropical svm under four different theorems at all possible assignments when large enough (i.e. more than 10 points from each category) feasible data sets exist.
### lp_solver.R
This algorithm, which is much slower and needs feasible data set input, generates constraint sets of the linear programming of section 4 in the reference and inputs them to the function [lp](https://www.rdocumentation.org/packages/lpSolve/versions/5.6.15/topics/lp). A detailed instruction of the output of this algorithm is available here: [lp.object](https://www.rdocumentation.org/packages/lpSolve/versions/5.6.13.3/topics/lp.object).
