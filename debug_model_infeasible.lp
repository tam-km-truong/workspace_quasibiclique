\* max_e_h *\
Maximize
max_sum_rows_columns: sumCol_34 + sumCol_44 + sumRow_24
Subject To
col_slack_34: - 5 col_34 - row_24 + sumCol_34 >= -5
col_slack_44: - 5 col_44 - row_24 + sumCol_44 >= -5
col_sum_34: - row_24 + sumCol_34 = 0
col_sum_44: - row_24 + sumCol_44 = 0
col_threshold: col_34 + col_44 >= 2
row_slack_24: - col_34 - col_44 - 5 row_24 + sumRow_24 >= -5
row_sum_24: - col_34 - col_44 + sumRow_24 = 0
row_threshold: row_24 >= 2
Bounds
 0 <= sumCol_34 <= 2
 0 <= sumCol_44 <= 2
 0 <= sumRow_24 <= 3
Generals
sumCol_34
sumCol_44
sumRow_24
Binaries
col_34
col_44
row_24
End
