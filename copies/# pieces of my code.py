# pieces of my code

def get_data(path:str, model):

    rows_data = []
    cols_data = []
    edges = []
    model_name = model
    df = pd.read_csv(path, header=0 ,index_col=0 )
    df[df == -1] = 0

    rows = df.sum(axis=1)
    row_names = rows.index
    #without resetting the indexes
    #rows_data = list(zip(rows.index, rows))
    #reset the indexes
    rows_data = list(zip(range(len(row_names)), rows))

    cols = df.sum(axis=0)
    col_names = cols.index
    #cols_data = list(zip(cols.index, cols))
    cols_data = list(zip(range(len(col_names)), cols))

    df = df.reset_index(drop=True)
    df = df.T.reset_index(drop=True).T
    if model_name == 'AB_V'  or model_name == 'AB_V_h'  or model_name == 'max_Surface' or model_name == 'max_Vertices' or model_name == 'max_Ones_comp'  or model_name == 'max_Ones' or model_name == 'AB_E' or model_name == 'AB_E_r' or model_name == 'AB_E_c_r':

        edges = list(df[df == 1].stack().index)
    else:
        edges = list(df[df == 0].stack().index)
     
    # print('-get data-' * 40)
    # print('edges =', edges)
    # print('rows =')
    # print(rows)
    # print('cols =')
    # print( cols)
    # print('rows_data =', rows_data)
    # print('cols_data =', cols_data)
    # print('rows_names =', row_names)
    # print('col_names =', col_names)
    # print()
    # print('-' * 40)
    

    return rows_data, cols_data, edges, row_names, col_names, df


def get_data_txt_file(path, model):
    file = open(path,'r')
    content = file.readlines()
    name = content[0][2:-1]
    num_row = int(content[1][7:-1])
    num_col = int(content[2][7:-1])
    num_edge = int(content[3][7:-1])
    deg_row = [0]*num_row
    deg_col = [0]*num_col
    edges = []
    df = pd.DataFrame([[0]*num_col]*num_row)
    for line in content[4:]:
        #regex split with mult delimiter
        splitted_line = re.split('\t|\n',line)
        u, v = int(splitted_line[0]),int(splitted_line[1])
        edges = edges + [(u,v)]
        deg_row[u] = deg_row[u] + 1
        deg_col[v] = deg_col[v] + 1
        df.iloc[u,v] = 1 
        
    rows_data = list(zip(range(num_row), deg_row))
    cols_data = list(zip(range(num_col), deg_col))

    return rows_data, cols_data, edges, range(num_row), range(num_col), df


#----------------------------------------
Input data :
Rows Data: [(0, 3), (1, 5), (2, 0), (3, 1), (4, 6), (5, 5), (6, 6), (7, 6), (8, 1), (9, 3)]
Columns Data: [(0, 0), (1, 5), (2, 6), (3, 7), (4, 4), (5, 1), (6, 6), (7, 1), (8, 1), (9, 5)]
Original Edges: [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 9), (1, 0), (1, 1), (1, 5), (1, 7), (1, 8), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (4, 0), (4, 5), (4, 7), (4, 8), (5, 0), (5, 4), (5, 7), (5, 8), (5, 9), (6, 0), (6, 5), (6, 7), (6, 8), (7, 0), (7, 5), (7, 7), (7, 8), (8, 0), (8, 1), (8, 2), (8, 4), (8, 5), (8, 6), (8, 7), (8, 8), (8, 9), (9, 0), (9, 4), (9, 5), (9, 6), (9, 7), (9, 8), (9, 9)]
Adjacency Matrix:
    0  1  2  3  4  5  6  7  8  9
0  0  0  0  0  0  0  1  1  1  0
1  0  0  1  1  1  0  1  0  0  1
2  0  0  0  0  0  0  0  0  0  0
3  0  0  0  0  0  0  0  0  0  1
4  0  1  1  1  1  0  1  0  0  1
5  0  1  1  1  0  1  1  0  0  0
6  0  1  1  1  1  0  1  0  0  1
7  0  1  1  1  1  0  1  0  0  1
8  0  0  0  1  0  0  0  0  0  0
9  0  1  1  1  0  0  0  0  0  0
rows_names = Index(['r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10'], dtype='object')
col_names = Index(['c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8', 'c9', 'c10'], dtype='object')


-after rows deletion ****
len_rows_res= 6
row_res= ['1', '4', '5', '6', '7', '9']
len_rows_del= 4
rows_del= ['0', '2', '3', '8']
len_cols_res= 10
len_cols_del= 0
