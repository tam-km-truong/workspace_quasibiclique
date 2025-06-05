\* König_E *\
Minimize
min_degrees: 3 col_110 + 4 col_113 + 3 col_121 + 3 col_122 + 4 col_123
 + 4 col_138 + 4 col_140 + 4 col_152 + 4 col_160 + 4 col_161 + 3 col_167
 + 3 col_177 + 3 col_181 + 4 col_19 + 3 col_217 + 3 col_220 + 3 col_234
 + 4 col_251 + 3 col_281 + 4 col_284 + 3 col_289 + 3 col_381 + 3 col_402
 + 4 col_407 + 3 col_414 + 3 col_418 + 4 col_423 + 3 col_430 + 4 col_55
 + 27 row_2 + 27 row_23 + 23 row_24 + 23 row_27
Subject To
edge_23_110: col_110 + row_23 >= 1
edge_23_430: col_430 + row_23 >= 1
edge_24_121: col_121 + row_24 >= 1
edge_24_181: col_181 + row_24 >= 1
edge_24_234: col_234 + row_24 >= 1
edge_24_281: col_281 + row_24 >= 1
edge_24_381: col_381 + row_24 >= 1
edge_24_418: col_418 + row_24 >= 1
edge_27_167: col_167 + row_27 >= 1
edge_27_177: col_177 + row_27 >= 1
edge_27_217: col_217 + row_27 >= 1
edge_27_220: col_220 + row_27 >= 1
edge_27_289: col_289 + row_27 >= 1
edge_27_402: col_402 + row_27 >= 1
edge_2_122: col_122 + row_2 >= 1
edge_2_414: col_414 + row_2 >= 1
Bounds
 0 <= col_110 <= 1
 0 <= col_113 <= 1
 0 <= col_121 <= 1
 0 <= col_122 <= 1
 0 <= col_123 <= 1
 0 <= col_138 <= 1
 0 <= col_140 <= 1
 0 <= col_152 <= 1
 0 <= col_160 <= 1
 0 <= col_161 <= 1
 0 <= col_167 <= 1
 0 <= col_177 <= 1
 0 <= col_181 <= 1
 0 <= col_19 <= 1
 0 <= col_217 <= 1
 0 <= col_220 <= 1
 0 <= col_234 <= 1
 0 <= col_251 <= 1
 0 <= col_281 <= 1
 0 <= col_284 <= 1
 0 <= col_289 <= 1
 0 <= col_381 <= 1
 0 <= col_402 <= 1
 0 <= col_407 <= 1
 0 <= col_414 <= 1
 0 <= col_418 <= 1
 0 <= col_423 <= 1
 0 <= col_430 <= 1
 0 <= col_55 <= 1
 0 <= row_2 <= 1
 0 <= row_23 <= 1
 0 <= row_24 <= 1
 0 <= row_27 <= 1
End
