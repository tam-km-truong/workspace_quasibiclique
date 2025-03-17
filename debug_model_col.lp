\* column_knapsack_problem *\
Minimize
knapsack: 25 col_11 + 24 col_13 + 24 col_15 + 2 col_16 + 25 col_18 + 24 col_19
 + 24 col_21 + 25 col_22 + col_27 + 3 col_28 + 25 col_3 + col_38 + 25 col_39
 + 26 col_47 + 25 col_6
Subject To
sensitivity: col_11 + 2 col_13 + 2 col_15 + 24 col_16 + col_18 + 2 col_19
 + 2 col_21 + col_22 + 25 col_27 + 23 col_28 + col_3 + 25 col_38 + col_39
 + col_6 >= 11.1
Bounds
 col_11 <= 1
 col_13 <= 1
 col_15 <= 1
 col_16 <= 1
 col_18 <= 1
 col_19 <= 1
 col_21 <= 1
 col_22 <= 1
 col_27 <= 1
 col_28 <= 1
 col_3 <= 1
 col_38 <= 1
 col_39 <= 1
 col_47 <= 1
 col_6 <= 1
End
