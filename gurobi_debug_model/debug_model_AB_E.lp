\* AB_E_r *\
Maximize
maximize_weight: cell_0_0 + cell_0_1 + cell_1_0 + cell_1_1
Subject To
ce_ll_0_0_2: - cell_0_0 + col_0 >= 0
ce_ll_0_1_2: - cell_0_1 + col_1 >= 0
ce_ll_1_0_2: - cell_1_0 + col_0 >= 0
ce_ll_1_1_2: - cell_1_1 + col_1 >= 0
cell_0_0_1: - cell_0_0 + row_0 >= 0
cell_0_1_1: - cell_0_1 + row_0 >= 0
cell_1_0_1: - cell_1_0 + row_1 >= 0
cell_1_1_1: - cell_1_1 + row_1 >= 0
col_err_rate_1_0: - 6 col_0 + 0.1 row_0 + 0.1 row_1 >= -6
col_err_rate_1_1: - 6 col_1 + 0.1 row_0 + 0.1 row_1 >= -6
col_threshold: col_0 + col_1 >= 2
row_err_rate_1_0: 0.1 col_0 + 0.1 col_1 - 6 row_0 >= -6
row_err_rate_1_1: 0.1 col_0 + 0.1 col_1 - 6 row_1 >= -6
row_threshold: row_0 + row_1 >= 2
Bounds
 cell_0_0 <= 1
 cell_0_1 <= 1
 cell_1_0 <= 1
 cell_1_1 <= 1
Binaries
col_0
col_1
row_0
row_1
End
