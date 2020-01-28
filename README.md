# Tropical-SVM
Support Vector Machine applicable to phylogenetic tree classification.
## Reference
Our arxiv link will be posted here
## Installation
We recommend downloading this repository to your local folder.
## Contents
### Algorithms 1~5
Algorithms 1~5 are all listed on section 5 of the reference.
### Sample
Here we provided some sample codes to demonstrate the accuracy of our tropical svm under four different theorems at all possible assignments when large enough (i.e. more than 10 points from each category) feasible data sets exist.
### lp_solver
This algorithm solves the linear programming on section 4 of the reference, which is much slower needs feasible data set input. It depends on the package lpSolve and its output is the [lp.object](https://www.rdocumentation.org/packages/lpSolve/versions/5.6.13.3/topics/lp.object)
