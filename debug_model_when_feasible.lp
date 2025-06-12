\* König_E *\
Minimize
min_degrees: 10 col_102 + 10 col_110 + 10 col_116 + 10 col_125 + 10 col_129
 + 10 col_186 + 10 col_189 + 10 col_191 + 10 col_22 + 10 col_54 + 10 col_58
 + 10 col_6 + 10 col_61 + 10 col_63 + 10 col_97 + 15 row_10 + 14 row_12
 + 12 row_13 + 13 row_16 + 14 row_19 + 11 row_21 + 15 row_29 + 15 row_30
 + 14 row_31 + 13 row_42 + 14 row_5
Subject To
edge_12_189: col_189 + row_12 >= 1
edge_13_54: col_54 + row_13 >= 1
edge_13_58: col_58 + row_13 >= 1
edge_13_63: col_63 + row_13 >= 1
edge_16_6: col_6 + row_16 >= 1
edge_16_61: col_61 + row_16 >= 1
edge_19_97: col_97 + row_19 >= 1
edge_21_102: col_102 + row_21 >= 1
edge_21_110: col_110 + row_21 >= 1
edge_21_116: col_116 + row_21 >= 1
edge_21_186: col_186 + row_21 >= 1
edge_31_129: col_129 + row_31 >= 1
edge_42_125: col_125 + row_42 >= 1
edge_42_191: col_191 + row_42 >= 1
edge_5_22: col_22 + row_5 >= 1
Bounds
 0 <= col_102 <= 1
 0 <= col_110 <= 1
 0 <= col_116 <= 1
 0 <= col_125 <= 1
 0 <= col_129 <= 1
 0 <= col_186 <= 1
 0 <= col_189 <= 1
 0 <= col_191 <= 1
 0 <= col_22 <= 1
 0 <= col_54 <= 1
 0 <= col_58 <= 1
 0 <= col_6 <= 1
 0 <= col_61 <= 1
 0 <= col_63 <= 1
 0 <= col_97 <= 1
 0 <= row_10 <= 1
 0 <= row_12 <= 1
 0 <= row_13 <= 1
 0 <= row_16 <= 1
 0 <= row_19 <= 1
 0 <= row_21 <= 1
 0 <= row_29 <= 1
 0 <= row_30 <= 1
 0 <= row_31 <= 1
 0 <= row_42 <= 1
 0 <= row_5 <= 1
End
