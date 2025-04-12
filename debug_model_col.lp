\* column_knapsack_problem *\
Minimize
knapsack: 5 col_103 + 5 col_169 + 6 col_177 + 5 col_192 + 4 col_38
Subject To
sensitivity: col_103 + col_169 + col_192 + 2 col_38 >= 0.005
Bounds
 col_103 <= 1
 col_169 <= 1
 col_177 <= 1
 col_192 <= 1
 col_38 <= 1
End
