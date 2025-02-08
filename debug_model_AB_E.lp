\* AB_E *\
Maximize
maximize_weight: cell_1_2 + cell_1_3 + cell_1_4 + cell_1_6 + cell_1_9
 + cell_4_1 + cell_4_2 + cell_4_3 + cell_4_4 + cell_4_6 + cell_4_9 + cell_5_1
 + cell_5_2 + cell_5_3 + cell_5_5 + cell_5_6 + cell_6_1 + cell_6_2 + cell_6_3
 + cell_6_4 + cell_6_6 + cell_6_9 + cell_7_1 + cell_7_2 + cell_7_3 + cell_7_4
 + cell_7_6 + cell_7_9
Subject To
cell_1_2_3: - cell_1_2 + 0.5 col_2 + 0.5 row_1 >= 0
cell_1_2_4: - cell_1_2 + col_2 + row_1 <= 1
cell_1_3_3: - cell_1_3 + 0.5 col_3 + 0.5 row_1 >= 0
cell_1_3_4: - cell_1_3 + col_3 + row_1 <= 1
cell_1_4_3: - cell_1_4 + 0.5 col_4 + 0.5 row_1 >= 0
cell_1_4_4: - cell_1_4 + col_4 + row_1 <= 1
cell_1_6_3: - cell_1_6 + 0.5 col_6 + 0.5 row_1 >= 0
cell_1_6_4: - cell_1_6 + col_6 + row_1 <= 1
cell_1_9_3: - cell_1_9 + 0.5 col_9 + 0.5 row_1 >= 0
cell_1_9_4: - cell_1_9 + col_9 + row_1 <= 1
cell_4_1_3: - cell_4_1 + 0.5 col_1 + 0.5 row_4 >= 0
cell_4_1_4: - cell_4_1 + col_1 + row_4 <= 1
cell_4_2_3: - cell_4_2 + 0.5 col_2 + 0.5 row_4 >= 0
cell_4_2_4: - cell_4_2 + col_2 + row_4 <= 1
cell_4_3_3: - cell_4_3 + 0.5 col_3 + 0.5 row_4 >= 0
cell_4_3_4: - cell_4_3 + col_3 + row_4 <= 1
cell_4_4_3: - cell_4_4 + 0.5 col_4 + 0.5 row_4 >= 0
cell_4_4_4: - cell_4_4 + col_4 + row_4 <= 1
cell_4_6_3: - cell_4_6 + 0.5 col_6 + 0.5 row_4 >= 0
cell_4_6_4: - cell_4_6 + col_6 + row_4 <= 1
cell_4_9_3: - cell_4_9 + 0.5 col_9 + 0.5 row_4 >= 0
cell_4_9_4: - cell_4_9 + col_9 + row_4 <= 1
cell_5_1_3: - cell_5_1 + 0.5 col_1 + 0.5 row_5 >= 0
cell_5_1_4: - cell_5_1 + col_1 + row_5 <= 1
cell_5_2_3: - cell_5_2 + 0.5 col_2 + 0.5 row_5 >= 0
cell_5_2_4: - cell_5_2 + col_2 + row_5 <= 1
cell_5_3_3: - cell_5_3 + 0.5 col_3 + 0.5 row_5 >= 0
cell_5_3_4: - cell_5_3 + col_3 + row_5 <= 1
cell_5_5_3: - cell_5_5 + 0.5 col_5 + 0.5 row_5 >= 0
cell_5_5_4: - cell_5_5 + col_5 + row_5 <= 1
cell_5_6_3: - cell_5_6 + 0.5 col_6 + 0.5 row_5 >= 0
cell_5_6_4: - cell_5_6 + col_6 + row_5 <= 1
cell_6_1_3: - cell_6_1 + 0.5 col_1 + 0.5 row_6 >= 0
cell_6_1_4: - cell_6_1 + col_1 + row_6 <= 1
cell_6_2_3: - cell_6_2 + 0.5 col_2 + 0.5 row_6 >= 0
cell_6_2_4: - cell_6_2 + col_2 + row_6 <= 1
cell_6_3_3: - cell_6_3 + 0.5 col_3 + 0.5 row_6 >= 0
cell_6_3_4: - cell_6_3 + col_3 + row_6 <= 1
cell_6_4_3: - cell_6_4 + 0.5 col_4 + 0.5 row_6 >= 0
cell_6_4_4: - cell_6_4 + col_4 + row_6 <= 1
cell_6_6_3: - cell_6_6 + 0.5 col_6 + 0.5 row_6 >= 0
cell_6_6_4: - cell_6_6 + col_6 + row_6 <= 1
cell_6_9_3: - cell_6_9 + 0.5 col_9 + 0.5 row_6 >= 0
cell_6_9_4: - cell_6_9 + col_9 + row_6 <= 1
cell_7_1_3: - cell_7_1 + 0.5 col_1 + 0.5 row_7 >= 0
cell_7_1_4: - cell_7_1 + col_1 + row_7 <= 1
cell_7_2_3: - cell_7_2 + 0.5 col_2 + 0.5 row_7 >= 0
cell_7_2_4: - cell_7_2 + col_2 + row_7 <= 1
cell_7_3_3: - cell_7_3 + 0.5 col_3 + 0.5 row_7 >= 0
cell_7_3_4: - cell_7_3 + col_3 + row_7 <= 1
cell_7_4_3: - cell_7_4 + 0.5 col_4 + 0.5 row_7 >= 0
cell_7_4_4: - cell_7_4 + col_4 + row_7 <= 1
cell_7_6_3: - cell_7_6 + 0.5 col_6 + 0.5 row_7 >= 0
cell_7_6_4: - cell_7_6 + col_6 + row_7 <= 1
cell_7_9_3: - cell_7_9 + 0.5 col_9 + 0.5 row_7 >= 0
cell_7_9_4: - cell_7_9 + col_9 + row_7 <= 1
col_err_rate_1_1: - 14 col_1 + slackCol_1 >= -14
col_err_rate_1_2: - 14 col_2 + slackCol_2 >= -14
col_err_rate_1_3: - 14 col_3 + slackCol_3 >= -14
col_err_rate_1_4: - 14 col_4 + slackCol_4 >= -14
col_err_rate_1_5: - 14 col_5 + slackCol_5 >= -14
col_err_rate_1_6: - 14 col_6 + slackCol_6 >= -14
col_err_rate_1_9: - 14 col_9 + slackCol_9 >= -14
col_slack_1: 0.9 row_1 + 0.9 row_4 + 0.9 row_5 + 0.9 row_6 + 0.9 row_7
 + slackCol_1 - sumCol_1 = 0
col_slack_2: 0.9 row_1 + 0.9 row_4 + 0.9 row_5 + 0.9 row_6 + 0.9 row_7
 + slackCol_2 - sumCol_2 = 0
col_slack_3: 0.9 row_1 + 0.9 row_4 + 0.9 row_5 + 0.9 row_6 + 0.9 row_7
 + slackCol_3 - sumCol_3 = 0
col_slack_4: 0.9 row_1 + 0.9 row_4 + 0.9 row_5 + 0.9 row_6 + 0.9 row_7
 + slackCol_4 - sumCol_4 = 0
col_slack_5: 0.9 row_1 + 0.9 row_4 + 0.9 row_5 + 0.9 row_6 + 0.9 row_7
 + slackCol_5 - sumCol_5 = 0
col_slack_6: 0.9 row_1 + 0.9 row_4 + 0.9 row_5 + 0.9 row_6 + 0.9 row_7
 + slackCol_6 - sumCol_6 = 0
col_slack_9: 0.9 row_1 + 0.9 row_4 + 0.9 row_5 + 0.9 row_6 + 0.9 row_7
 + slackCol_9 - sumCol_9 = 0
col_sum_1: - row_4 - row_5 - row_6 - row_7 + sumCol_1 = 0
col_sum_2: - row_1 - row_4 - row_5 - row_6 - row_7 + sumCol_2 = 0
col_sum_3: - row_1 - row_4 - row_5 - row_6 - row_7 + sumCol_3 = 0
col_sum_4: - row_1 - row_4 - row_6 - row_7 + sumCol_4 = 0
col_sum_5: - row_5 + sumCol_5 = 0
col_sum_6: - row_1 - row_4 - row_5 - row_6 - row_7 + sumCol_6 = 0
col_sum_9: - row_1 - row_4 - row_6 - row_7 + sumCol_9 = 0
col_threshold: col_1 + col_2 + col_3 + col_4 + col_5 + col_6 + col_9 >= 1
obj_sum_col: obj_colsum - sumCol_1 - sumCol_2 - sumCol_3 - sumCol_4 - sumCol_5
 - sumCol_6 - sumCol_9 = 0
obj_sum_row: obj_rowsum - sumRow_1 - sumRow_4 - sumRow_5 - sumRow_6 - sumRow_7
 = 0
row_err_rate_1_1: - 14 row_1 + slackRow_1 >= -14
row_err_rate_1_4: - 14 row_4 + slackRow_4 >= -14
row_err_rate_1_5: - 14 row_5 + slackRow_5 >= -14
row_err_rate_1_6: - 14 row_6 + slackRow_6 >= -14
row_err_rate_1_7: - 14 row_7 + slackRow_7 >= -14
row_slack_1: 0.9 col_1 + 0.9 col_2 + 0.9 col_3 + 0.9 col_4 + 0.9 col_5
 + 0.9 col_6 + 0.9 col_9 + slackRow_1 - sumRow_1 = 0
row_slack_4: 0.9 col_1 + 0.9 col_2 + 0.9 col_3 + 0.9 col_4 + 0.9 col_5
 + 0.9 col_6 + 0.9 col_9 + slackRow_4 - sumRow_4 = 0
row_slack_5: 0.9 col_1 + 0.9 col_2 + 0.9 col_3 + 0.9 col_4 + 0.9 col_5
 + 0.9 col_6 + 0.9 col_9 + slackRow_5 - sumRow_5 = 0
row_slack_6: 0.9 col_1 + 0.9 col_2 + 0.9 col_3 + 0.9 col_4 + 0.9 col_5
 + 0.9 col_6 + 0.9 col_9 + slackRow_6 - sumRow_6 = 0
row_slack_7: 0.9 col_1 + 0.9 col_2 + 0.9 col_3 + 0.9 col_4 + 0.9 col_5
 + 0.9 col_6 + 0.9 col_9 + slackRow_7 - sumRow_7 = 0
row_sum_1: - col_2 - col_3 - col_4 - col_6 - col_9 + sumRow_1 = 0
row_sum_4: - col_1 - col_2 - col_3 - col_4 - col_6 - col_9 + sumRow_4 = 0
row_sum_5: - col_1 - col_2 - col_3 - col_5 - col_6 + sumRow_5 = 0
row_sum_6: - col_1 - col_2 - col_3 - col_4 - col_6 - col_9 + sumRow_6 = 0
row_sum_7: - col_1 - col_2 - col_3 - col_4 - col_6 - col_9 + sumRow_7 = 0
row_threshold: row_1 + row_4 + row_5 + row_6 + row_7 >= 1
Bounds
 obj_colsum free
 obj_rowsum free
 slackCol_1 free
 slackCol_2 free
 slackCol_3 free
 slackCol_4 free
 slackCol_5 free
 slackCol_6 free
 slackCol_9 free
 slackRow_1 free
 slackRow_4 free
 slackRow_5 free
 slackRow_6 free
 slackRow_7 free
 sumCol_1 free
 sumCol_2 free
 sumCol_3 free
 sumCol_4 free
 sumCol_5 free
 sumCol_6 free
 sumCol_9 free
 sumRow_1 free
 sumRow_4 free
 sumRow_5 free
 sumRow_6 free
 sumRow_7 free
Binaries
cell_1_2
cell_1_3
cell_1_4
cell_1_6
cell_1_9
cell_4_1
cell_4_2
cell_4_3
cell_4_4
cell_4_6
cell_4_9
cell_5_1
cell_5_2
cell_5_3
cell_5_5
cell_5_6
cell_6_1
cell_6_2
cell_6_3
cell_6_4
cell_6_6
cell_6_9
cell_7_1
cell_7_2
cell_7_3
cell_7_4
cell_7_6
cell_7_9
col_1
col_2
col_3
col_4
col_5
col_6
col_9
row_1
row_4
row_5
row_6
row_7
End
