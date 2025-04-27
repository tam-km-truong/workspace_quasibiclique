\* row_knapsack_problem *\
Minimize
knapsack: 4 row_19 + 4 row_23 + 3 row_30 + 4 row_43 + 4 row_49
Subject To
sensitivity: row_19 + row_23 + 2 row_30 + row_43 + row_49 >= 4.2
Bounds
 row_19 <= 1
 row_23 <= 1
 row_30 <= 1
 row_43 <= 1
 row_49 <= 1
End
