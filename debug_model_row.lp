\* row_knapsack_problem *\
Minimize
knapsack: row_1 + 5 row_3 + 2 row_5 + 3 row_9
Subject To
sensitivity: 5 row_1 + row_3 + 6 row_4 + 4 row_5 + 6 row_6 + 6 row_7 + 3 row_9
 >= 6.4
Bounds
 row_1 <= 1
 row_3 <= 1
 row_4 <= 1
 row_5 <= 1
 row_6 <= 1
 row_7 <= 1
 row_9 <= 1
End
