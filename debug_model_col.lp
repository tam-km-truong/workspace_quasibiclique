\* column_knapsack_problem *\
Minimize
knapsack: 5 col_1 + 4 col_2 + 3 col_3 + 6 col_4 + 4 col_6 + 5 col_9
Subject To
sensitivity: 2 col_1 + 3 col_2 + 4 col_3 + col_4 + 3 col_6 + 2 col_9 >= 1.5
Bounds
 col_1 <= 1
 col_2 <= 1
 col_3 <= 1
 col_4 <= 1
 col_6 <= 1
 col_9 <= 1
End
