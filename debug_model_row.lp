\* row_knapsack_problem *\
Minimize
knapsack: 20 row_12 + 19 row_29 + 20 row_3 + 21 row_30 + 19 row_37 + 19 row_40
 + 19 row_41 + 19 row_42 + 19 row_47 + 19 row_6
Subject To
sensitivity: 2 row_12 + 3 row_29 + 2 row_3 + row_30 + 3 row_37 + 3 row_40
 + 3 row_41 + 3 row_42 + 3 row_47 + 3 row_6 >= 2.6
Bounds
 row_12 <= 1
 row_29 <= 1
 row_3 <= 1
 row_30 <= 1
 row_37 <= 1
 row_40 <= 1
 row_41 <= 1
 row_42 <= 1
 row_47 <= 1
 row_6 <= 1
End
