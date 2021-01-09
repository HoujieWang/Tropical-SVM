# Sample Execution for algorithm 1
library(gtools); library(lpSolve)

# Load the data. For details of data set, please see the README on Github
load("data_15%.RData")

data <- rbind(data_matrix_list[[1]][[1]][[1]], data_matrix_list[[1]][[1]][[2]])
n <- nrow(data_matrix_list[[1]][[1]][[1]])
tst_data <- rbind(data_matrix_list[[1]][[1]][[3]], data_matrix_list[[1]][[1]][[4]])
ntst <- nrow(data_matrix_list[[1]][[1]][[3]])

# Create all possible assignments in the case of algorithm 1. For algorithm 2 and 4, please set r=3 and for algorithm 3 please set r=2
all_assignment <- permutations(n = ncol(data), r = 4)

# User can run all possible assignments and methods to find the best, which could be time consuming.
# We've performed experiment for the data set on Github. For a quick graphical evaluation, please run graph_producer.R
# beta can be set 1 as default and varying beta does not make obvious changes in performance 
# method_ind takes value from 1 to 70 (same for algorithm 2~4), where each represents a classification method
# User can set method_ind = 1: 70 to automatically compute the prediction accuracy in terms of 70 different methods
# So far, this algorithm only supports the case where two categories contain the same number of points
algorithm1(assignment = 1: 4, data, tst_data, n, ntst, beta = 1, method_ind = 1: 70)
