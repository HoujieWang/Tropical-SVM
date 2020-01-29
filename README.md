# Tropical-SVM
This is R code for computing a tropical support vector machine SVM) for classification over the tropical projective torus.  We can also apply this code for classifying gene trees into two classes.  A tropical SVM is, like a classical SVM, a tropical hyperplane which maximizes the minimum distance of the data points to the tropical hyperplane to separate the data points.  
## Reference
Our arxiv link will be posted here (TBA)
## Installation
To run this code we need the [ape](https://cran.r-project.org/web/packages/ape/index.html), [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html) and [parallel](https://CRAN.R-project.org/view=HighPerformanceComputing) packages. 
## Contents
### Algorithm1~5.R
These algorithms are Algorithms 1~5, respectively, in Section 5 in the reference.
### Sample.R
This R script displays the accuracy of our tropical svm by applying algorithm 1~4 to all possible assignments of the four theorems respectively and using algorithm 5 to calculate the accuracy. A user could directly run the codes as long as the packages are installed and the data (LittleTree6A.txt and LittleTree6B.txt, each contains 100 trees simulated at C=6) are put in the R directory. If a different data set is preferred, a user would have to maually set the number of leaves and load the data set by typing in the data name. 
### lp_solver.R
This algorithm solves the linear programming (equation (4)~(7)) in section 4 in the reference. It makes use of the functionality of the [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html) and its output is the [lp.object](https://www.rdocumentation.org/packages/lpSolve/versions/5.6.13.3/topics/lp.object).
