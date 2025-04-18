\* König_E *\
Minimize
min_degrees: 2 col_0 + 2 col_1 + 3 col_2 + col_3 + col_4 + 2 row_0 + 2 row_1
 + row_2 + row_3 + 3 row_4
Subject To
edge_0_2: col_2 + row_0 >= 1
edge_0_3: col_3 + row_0 >= 1
edge_0_4: col_4 + row_0 >= 1
edge_1_2: col_2 + row_1 >= 1
edge_1_3: col_3 + row_1 >= 1
edge_1_4: col_4 + row_1 >= 1
edge_2_0: col_0 + row_2 >= 1
edge_2_1: col_1 + row_2 >= 1
edge_2_3: col_3 + row_2 >= 1
edge_2_4: col_4 + row_2 >= 1
edge_3_0: col_0 + row_3 >= 1
edge_3_1: col_1 + row_3 >= 1
edge_3_3: col_3 + row_3 >= 1
edge_3_4: col_4 + row_3 >= 1
edge_4_0: col_0 + row_4 >= 1
edge_4_1: col_1 + row_4 >= 1
Bounds
 0 <= col_0 <= 1
 0 <= col_1 <= 1
 0 <= col_2 <= 1
 0 <= col_3 <= 1
 0 <= col_4 <= 1
 0 <= row_0 <= 1
 0 <= row_1 <= 1
 0 <= row_2 <= 1
 0 <= row_3 <= 1
 0 <= row_4 <= 1
End
