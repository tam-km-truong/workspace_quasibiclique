\* row_knapsack_problem *\
Minimize
knapsack: 2 row_0 + 2 row_1 + row_4
Subject To
sensitivity: row_0 + row_1 + 2 row_4 >= 0.004
Bounds
 row_0 <= 1
 row_1 <= 1
 row_4 <= 1
End
