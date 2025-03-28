\* column_knapsack_problem *\
Minimize
knapsack: 4 col_63 + 6 col_64 + 4 col_79 + 4 col_83 + 5 col_95
Subject To
sensitivity: 2 col_63 + 2 col_79 + 2 col_83 + col_95 >= 0.7
Bounds
 col_63 <= 1
 col_64 <= 1
 col_79 <= 1
 col_83 <= 1
 col_95 <= 1
End
