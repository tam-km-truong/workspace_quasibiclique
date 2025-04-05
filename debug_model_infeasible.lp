\* max_e_h *\
Maximize
max_sum_rows_columns: sumCol_7 + sumRow_39
Subject To
col_slack_7: - 4 col_7 - row_39 + sumCol_7 >= -4
col_sum_7: - row_39 + sumCol_7 = 0
col_threshold: col_7 >= 2
row_slack_39: - col_7 - 4 row_39 + sumRow_39 >= -4
row_sum_39: - col_7 + sumRow_39 = 0
row_threshold: row_39 >= 2
Bounds
 0 <= sumCol_7 <= 2
 0 <= sumRow_39 <= 2
Generals
sumCol_7
sumRow_39
Binaries
col_7
row_39
End
