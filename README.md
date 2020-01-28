# Tropical-SVM
Support Vector Machine applicable to phylogenetic tree classification.
## Reference
Our arxiv link will be posted here
## Installation
We recommend downloading this repository to your local folder.
## Contents
### Algorithms 1~5
Algorithms 1~5 are all listed at section 5 in the reference.
### Sample
Here we provided some sample codes to demonstrate the accuracy of our tropical svm under four different theorems at all possible assignments when large enough (i.e. more than 10 points from each category) feasible data sets exist.
### lp_solver
This algorithm, which is much slower and needs feasible data set input, generates constraint sets of the linear programming of section 4 in the reference and inputs them to the function [lp](https://www.rdocumentation.org/packages/lpSolve/versions/5.6.15/topics/lp). A detailed instruction of the output of this algorithm is available here: [lp.object](https://www.rdocumentation.org/packages/lpSolve/versions/5.6.13.3/topics/lp.object).
