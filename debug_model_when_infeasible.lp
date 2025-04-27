\* max_e_h *\
Maximize
max_sum_rows_columns: sumCol_10 + sumCol_18 + sumCol_27 + sumCol_42
 + sumRow_23
Subject To
col_slack_10: - 7 col_10 - row_23 + sumCol_10 >= -7
col_slack_18: - 7 col_18 - row_23 + sumCol_18 >= -7
col_slack_27: - 7 col_27 - row_23 + sumCol_27 >= -7
col_slack_42: - 7 col_42 - row_23 + sumCol_42 >= -7
col_sum_10: - row_23 + sumCol_10 = 0
col_sum_18: - row_23 + sumCol_18 = 0
col_sum_27: - row_23 + sumCol_27 = 0
col_sum_42: - row_23 + sumCol_42 = 0
col_threshold: col_10 + col_18 + col_27 + col_42 >= 2
row_slack_23: - col_10 - col_18 - col_27 - col_42 - 7 row_23 + sumRow_23 >= -7
row_sum_23: - col_10 - col_18 - col_27 - col_42 + sumRow_23 = 0
row_threshold: row_23 >= 2
Bounds
 0 <= sumCol_10 <= 2
 0 <= sumCol_18 <= 2
 0 <= sumCol_27 <= 2
 0 <= sumCol_42 <= 2
 0 <= sumRow_23 <= 5
Generals
sumCol_10
sumCol_18
sumCol_27
sumCol_42
sumRow_23
Binaries
col_10
col_18
col_27
col_42
row_23
End
