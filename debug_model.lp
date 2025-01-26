\* AB_E_r *\
Maximize
max_sum_rows_columns: sumCol_0 + sumCol_1 + sumCol_2 + sumRow_0 + sumRow_1
 + sumRow_2
Subject To
col_err_rate_1_0: - 8 col_0 + slackCol_0 >= -8
col_err_rate_1_1: - 8 col_1 + slackCol_1 >= -8
col_err_rate_1_2: - 8 col_2 + slackCol_2 >= -8
col_slack_0: 0.9 row_0 + 0.9 row_1 + 0.9 row_2 + slackCol_0 - sumCol_0 = 0
col_slack_1: 0.9 row_0 + 0.9 row_1 + 0.9 row_2 + slackCol_1 - sumCol_1 = 0
col_slack_2: 0.9 row_0 + 0.9 row_1 + 0.9 row_2 + slackCol_2 - sumCol_2 = 0
col_sum_0: - row_0 - row_1 + sumCol_0 = 0
col_sum_1: - row_0 - row_1 - row_2 + sumCol_1 = 0
col_sum_2: - row_2 + sumCol_2 = 0
col_threshold: col_0 + col_1 + col_2 >= 2
obj_sum_col: obj_colsum - sumCol_0 - sumCol_1 - sumCol_2 = 0
obj_sum_row: obj_rowsum - sumRow_0 - sumRow_1 - sumRow_2 = 0
row_err_rate_1_0: - 8 row_0 + slackRow_0 >= -8
row_err_rate_1_1: - 8 row_1 + slackRow_1 >= -8
row_err_rate_1_2: - 8 row_2 + slackRow_2 >= -8
row_slack_0: 0.9 col_0 + 0.9 col_1 + 0.9 col_2 + slackRow_0 - sumRow_0 = 0
row_slack_1: 0.9 col_0 + 0.9 col_1 + 0.9 col_2 + slackRow_1 - sumRow_1 = 0
row_slack_2: 0.9 col_0 + 0.9 col_1 + 0.9 col_2 + slackRow_2 - sumRow_2 = 0
row_sum_0: - col_0 - col_1 + sumRow_0 = 0
row_sum_1: - col_0 - col_1 + sumRow_1 = 0
row_sum_2: - col_1 - col_2 + sumRow_2 = 0
row_threshold: row_0 + row_1 + row_2 >= 2
Bounds
 obj_colsum free
 obj_rowsum free
 slackCol_0 free
 slackCol_1 free
 slackCol_2 free
 slackRow_0 free
 slackRow_1 free
 slackRow_2 free
 sumCol_0 free
 sumCol_1 free
 sumCol_2 free
 sumRow_0 free
 sumRow_1 free
 sumRow_2 free
Binaries
col_0
col_1
col_2
row_0
row_1
row_2
End
