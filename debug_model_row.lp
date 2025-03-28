\* row_knapsack_problem *\
Minimize
knapsack: 4 row_22 + 4 row_3 + 3 row_4 + 4 row_7 + 3 row_75 + 3 row_83
 + 5 row_89
Subject To
sensitivity: row_22 + row_3 + 2 row_4 + row_7 + 2 row_75 + 2 row_83 >= 0.9
Bounds
 row_22 <= 1
 row_3 <= 1
 row_4 <= 1
 row_7 <= 1
 row_75 <= 1
 row_83 <= 1
 row_89 <= 1
End
