\* column_knapsack_problem *\
Minimize
knapsack: col_10 + col_18 + col_27 + col_42
Subject To
sensitivity: __dummy >= 0
Bounds
 __dummy = 0
 col_10 <= 1
 col_18 <= 1
 col_27 <= 1
 col_42 <= 1
End
