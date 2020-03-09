# Tropical-SVM
This is R code for computing a tropical support vector machine SVM) for classification over the tropical projective torus.  We can also apply this code for classifying gene trees into two classes.  A tropical SVM is, like a classical SVM, a tropical hyperplane which maximizes the minimum distance of the data points to the tropical hyperplane to separate the data points.  
## Reference
https://arxiv.org/abs/2003.00677
## Installation
To run this code we need the [ape](https://cran.r-project.org/web/packages/ape/index.html), [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html), [parallel](https://CRAN.R-project.org/view=HighPerformanceComputing), [gtools](https://cran.r-project.org/web/packages/gtools/index.html), [ggplot2](https://ggplot2.tidyverse.org) and [e1071](https://cran.r-project.org/web/packages/e1071/index.html) packages.
## Contents
### Trees
Each file contains 100 trees. The file names indicate the C value when it was generated and its categories (A or B). For example:
```{r}
LittleTree6A.txt
```
It means that the 100 trees in this file are produced at C=6 and labeled as category A.
### Algorithm1~4.R
These algorithms are Algorithms 1~4, respectively, in Section 5 in the reference.
### Sample.R
This R script respectively applys algorithm 1~4 to all possible assignments of the four theorems and calculates the accuracy at these assignments with algorithm 5. A user could directly run all the codes as long as the packages are installed and the data (We used LittleTree6A.txt and LittleTree6B.txt for this sample) are put in the R directory. If a different data set is preferred (e.g. miniTree0.8A.txt and miniTree0.8B.txt), a user would have to maually set the number of leaves and load the data set by typing in the data name. 
```{r}
N = 5
P  = matrix(unlist(lapply(read.tree("miniTree0.8A.txt")...
Q  = matrix(unlist(lapply(read.tree("miniTree0.8B.txt")...
```

