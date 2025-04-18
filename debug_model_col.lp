\* column_knapsack_problem *\
Minimize
knapsack: 2 col_0 + 2 col_1
Subject To
sensitivity: __dummy >= 0
Bounds
 __dummy = 0
 col_0 <= 1
 col_1 <= 1
End
