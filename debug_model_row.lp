\* row_knapsack_problem *\
Minimize
knapsack: 4 row_23 + 4 row_24 + 4 row_39 + 3 row_40 + 5 row_42 + 4 row_5
 + 4 row_8
Subject To
sensitivity: row_23 + row_24 + row_39 + 2 row_40 + row_5 + row_8 >= 0.007
Bounds
 row_23 <= 1
 row_24 <= 1
 row_39 <= 1
 row_40 <= 1
 row_42 <= 1
 row_5 <= 1
 row_8 <= 1
End
