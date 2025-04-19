\* max_e_h *\
Maximize
max_sum_rows_columns: sumCol_1 + sumRow_2 + sumRow_3 + sumRow_4
Subject To
col_slack_1: - 6 col_1 - row_2 - row_3 - row_4 + sumCol_1 >= -6
col_sum_1: - row_2 - row_3 - row_4 + sumCol_1 = 0
col_threshold: col_1 >= 2
row_slack_2: - col_1 - 6 row_2 + sumRow_2 >= -6
row_slack_3: - col_1 - 6 row_3 + sumRow_3 >= -6
row_slack_4: - col_1 - 6 row_4 + sumRow_4 >= -6
row_sum_2: - col_1 + sumRow_2 = 0
row_sum_3: - col_1 + sumRow_3 = 0
row_sum_4: - col_1 + sumRow_4 = 0
row_threshold: row_2 + row_3 + row_4 >= 2
Bounds
 0 <= sumCol_1 <= 4
 0 <= sumRow_2 <= 2
 0 <= sumRow_3 <= 2
 0 <= sumRow_4 <= 2
Generals
sumCol_1
sumRow_2
sumRow_3
sumRow_4
Binaries
col_1
row_2
row_3
row_4
End
