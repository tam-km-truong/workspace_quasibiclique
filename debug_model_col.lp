\* column_knapsack_problem *\
Minimize
knapsack: col_1 + 2 col_4 + col_6
Subject To
sensitivity: 2 col_1 + col_4 + 2 col_6 >= 3.5
Bounds
 col_1 <= 1
 col_4 <= 1
 col_6 <= 1
End
