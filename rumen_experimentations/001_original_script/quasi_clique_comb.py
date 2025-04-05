# -*- coding=utf-8 -*-

"""
             Searching Maximum Quasi-Bicliques

Install Gurobi Solver with free academic licence using your University account at: https://www.gurobi.com/features/academic-named-user-license/

Install Pandas, using PyPI or Conda. A detailed installation tutorial can be find at: https://pandas.pydata.org/docs/getting_started/install.html

an instance 
 python3 quasi_clique_comb.py --filepath data/data3.csv   --model  max_e_c  --rho 0.1  --delta 0.1 --debug 0  --threshold 0.83 --dec_conq 0 --timelimit 100 


"""
from pulp import value, LpStatus, GUROBI_CMD
import sys
import csv 
import re
import argparse 
import itertools 
parser = argparse.ArgumentParser()
from argparse import ArgumentParser

import pandas as pd
import itertools
import gurobipy as gp
import gurobipy as GRB 
import numpy as np 
import heapq 
#from tabulate import tabulate
from pulp import (
    GUROBI_CMD,
    PULP_CBC_CMD,
    LpBinary,
    LpMaximize,
    LpMinimize,
    LpProblem,
    LpStatus,
    LpVariable,
    LpContinuous,
    lpSum,
)
from pulp import value
from gurobipy import Model, GRB


def max_v(rows_data, cols_data, edges, delta):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    rho: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    """
    # FORCING 
    #rho = 0.0
    # FORCING 
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M= Big_R+Big_C 

    model = LpProblem(name="max_v", sense=LpMaximize)

    # Variables for rows and columns

    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}
    # Objective function: Maximize sum of selected row and column variables
    model += lpSum(
        [lpvar for lpvar, _ in lpRows.values()] +
        [lpvar for lpvar, _ in lpCols.values()]), "max_sum_vertices"
    

    # Constraints for row and column thresholds
    row_threshold = 1
    col_threshold = 1
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"

    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)

    return model 

def AB_V_h(rows_data, cols_data, edges, rho):
    """<
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    rho: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.

    This is a variant of max_v !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    """
    debug = 0
    model = LpProblem(name="AB_V_h", sense=LpMaximize)

    # Variables for rows and columns

    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}             
    # Objective function: Maximize sum of selected row and column variables
    model += lpSum(
        [lpvar for lpvar, _ in lpRows.values()] +
        [lpvar for lpvar, _ in lpCols.values()]), "max_sum_vertices"

    # Constraints for row and column thresholds
    row_threshold = 2
    col_threshold = 2
    if debug >= 2:
        print()
        print('-' * 40)
        print('row_threshold=', row_threshold )
        print('col_threshold=', col_threshold )
        print()
        print('-' * 40)
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    # Add row density constraints
    #__row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, rho)
    #__col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, rho)
    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, rho)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, rho)

    return model 


def max_e(rows_data, cols_data, edges, delta):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    """


    model = LpProblem(name="max_e", sense=LpMaximize)
    # Variables for rows and columns

    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpSum_row = {row: (LpVariable(f'sumRow_{row}', cat='Continuous')) for row, _ in rows_data}
    lpSum_col= {col: (LpVariable(f'sumCol_{col}', cat='Continuous')) for col, _ in cols_data}
    lpCells = {(row, col): LpVariable(
                    f'cell_{row}_{col}', cat="Binary")
                    #f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1) 
                    for (row, col) in edges
    }

    model += lpSum(lpCells), 'maximize_sum'
    # Constraints for row and column thresholds
    row_threshold = 2
    col_threshold = 2
    if debug >= 3:
        print()
        print('-' * 40) 
        print(f"******** Solving model  ********", model.name, " with delta = ", delta)
        print(' # rows_data =', len(rows_data),' # cols_data =', len(cols_data), ' # edges =', len(edges)) 
        print('row_threshold=', row_threshold )
        print('col_threshold=', col_threshold )
        print()
        print('-' * 40)
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    for row, col in edges:
    #for row, col in lpCells:
        model += (lpRows[row][0]+lpCols[col][0])/2 >= lpCells[(row, col)], f'cell_{row}_{col}_3'
        model += lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)], f'cell_{row}_{col}_4'
    # Add row/cols density constraints
    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    return model 


def max_e_r(rows_data, cols_data, edges, delta):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    
    This is an improved model from Chnag et al. 
    """
    model = LpProblem(name="max_e_r", sense=LpMaximize)

    # Variables for rows and columns

    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpCells = {(row, col): LpVariable(
                    #f'cell_{row}_{col}', cat="Binary")
                    f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1) 
                    for (row, col) in edges
    }

    #
    model += lpSum(lpCells), 'maximize_sum'
    # Constraints for row and column thresholds
    row_threshold = 1
    col_threshold = 1
    if debug >= 3:
        print()
        print('-' * 40) 
        print(f"******** Solving model  ********", model.name, " with delta = ", delta)
        print(' # rows_data =', len(rows_data),' # cols_data =', len(cols_data), ' # edges =', len(edges)) 
        print(' # rows_data =', rows_data, '\n, # cols_data =', cols_data, '\n, # edges =', edges)
        print('row_threshold=', row_threshold )
        print('col_threshold=', col_threshold )
        print()
        print('-' * 40)
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    # for row, col in edges:
    # #for row, col in lpCells:
    #     model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
    #     model += (lpCols[col][0] >= lpCells[(row, col)]), f'ce ll_{row}_{col}_2'
    # Constraints for matrix structure
    # Check for missing columns in cols_data
    missing_cols = {col for _, col in edges} - {col for col, _ in cols_data}
    if missing_cols:
        print(f"Warning: These columns are missing from cols_data but appear in edges: {missing_cols}")

    for row, col in edges:
        if row in lpRows and col in lpCols and (row, col) in lpCells:
            model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
            model += (lpCols[col][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_2'
        else:
            print(f"Warning: row={row}, col={col} not found in lpRows/lpCols!")

    # for row, col in edges:
    #     if (row, col) in lpCells:
    #         model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
    #         model += (lpCols[col][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_2'
    #     else:
    #         print(f"Warning: ({row}, {col}) not found in lpCells!")

    # Add row density constraints
    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)

    return model 

def max_e_wr(rows_data, cols_data, edges, rows_res, cols_res, prev_obj, delta):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    rows_res : list of indexies of rows (input solution)
    cols_res : list of indexies of cols (input solution)
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    
    """
    model = LpProblem(name="max_e_wr", sense=LpMaximize)

    # Variables for rows and columns

    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpCells = {(row, col): LpVariable(
                    #f'cell_{row}_{col}', cat="Binary")
                    f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1) 
                    for (row, col) in edges
    }

    cols_res_set = set(map(int, cols_res)) 
    rows_res_set = set(map(int, rows_res))
    # print("I am in warm start!!!!!!!!!!!!!!!!")
    #print("rows_res_set: bis", rows_res_set)
    #print("cols_res_set: bis", cols_res_set )

    print(" !!!!!!!!!!!!!!!!!! I got a lower bound ", prev_obj)
    # print(" !!!!!!!!!!!!!!!!!! nb col =", len(cols_data))

# Warm start (PuLP way: use setInitialValue if needed)
    # Set and store initial values
    #Attention: check what happens when the lower bound is given by previous bound !!!!!!!!!!!!!!!!!!!!!!!!
    ####################################################
    #row_indices: list[int] = []
    #column_indices: list[int] = []
    cells_indices = []
    row_initial_values = {}
    for row, (var, _) in lpRows.items():
        if row in rows_res_set:
            val = 1 
            var.setInitialValue(val)
            row_initial_values[row] = val  # Store for later use
            #row_indices.append(row)
        else:
            val = 0
            row_initial_values[row] = val  # Store for later use 

    col_initial_values = {}
    for col, (var, _) in lpCols.items():
        if col in cols_res_set:
            val = 1  
            var.setInitialValue(val)
            col_initial_values[col] = val
            #column_indices.append(col)
        else: 
            val = 0
            col_initial_values[col] = val

    cell_initial_values = {}
    for (row, col), var in lpCells.items():
        if (row, col) in itertools.product(row_initial_values, col_initial_values):
        #if (row, col) in itertools.product(row_indices, column_indices): 
            var.setInitialValue(1)
            val = 1  
            #cells_indices.append((row,col))
            cell_initial_values[(row, col)] = val
        else:
            var.setInitialValue(0)
            val = 0
            col_initial_values[col] = val 
    if debug >=3:
        print("\n Initial point before solving:")
        print("ROWS:", row_initial_values)
        print("COLS:", col_initial_values)
        #print("CELLS:", cell_initial_values)
    if debug >=4:
        print("\n Initial point (cells) before solving:")
        print("CELLS:", cell_initial_values)



    # Objective 
    model += lpSum(lpCells), 'current value'


    # row_initial_values = {}

    # for row, (var, _) in lpRows.items():
    #     initial_val = 1 if row in rows_res_set else 0
    #     var.setInitialValue(initial_val)
    #     row_initial_values[row] = initial_val  # Store separately

    # col_initial_values = {}

    # for col, (var, _) in lpCols.items():
    #     initial_val = 1 if col in cols_res_set else 0
    #     var.setInitialValue(initial_val)
    #     col_initial_values[col] = initial_val  # Store separately

    # Print initial values
    # for row, val in row_initial_values.items():
    #     print(f"Row {row}: {val}")

    # for col, val in col_initial_values.items():
    #     print(f"Col {col}: {val}")

    # for row in rows_res:
    #     lpRows[row][0].setInitialValue(1)
    # for col in cols_res:
    #     lpCols[col][0].setInitialValue(1)

    # Constraints for row and column thresholds
    #sys.exit("Terminating program to check warm start value ")
    # row_threshold = 2
    # col_threshold = 2
    # if debug >= 3:
    #     print()
    #     print('-' * 40) 
    #     print(f"******** Solving model  ********", model.name, " with delta = ", delta)
    #     print(' # rows_data =', len(rows_data),' # cols_data =', len(cols_data), ' # edges =', len(edges)) 
    #     print('row_threshold=', row_threshold )
    #     print('col_threshold=', col_threshold )
    #     print()
    #     print('-' * 40)
    # model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    # model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    # #constraint for obj improvement
    # model += lpSum(lpCells) >= prev_obj + 1, "impovement????"
    # #end constraint for obj improvement
    # for row, col in edges:
    # #for row, col in lpCells:
    #     model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
    #     model += (lpCols[col][0] >= lpCells[(row, col)]), f'ce ll_{row}_{col}_2'

    row_threshold = 2
    col_threshold = 2

    if debug >= 3:
        print()
        print('-' * 40) 
        print(f"******** Solving model ******** {model.name} with delta = {delta}")
        print(f"# rows_data = {len(rows_data)}, # cols_data = {len(cols_data)}, # edges = {len(edges)}") 
        print(f"row_threshold = {row_threshold}")
        print(f"col_threshold = {col_threshold}")
        print()
        print('-' * 40)

    # Add constraints for row and column thresholds
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"

    # Constraint for objective improvement
    model += lpSum(lpCells.values()) >= prev_obj + 1, "improvement"
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Constraints for matrix structure
    for row, col in edges:
        if (row, col) in lpCells:
            model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
            model += (lpCols[col][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_2'
        else:
            print(f"Warning: ({row}, {col}) not found in lpCells!")



    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)

    return model 

def max_e_h(rows_data, cols_data, edges, delta):
    """
    A heuristic proposed to max edge biclique problem.  The particularity is to use an objective that does not require  quadratic variables 
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    
    Goal : compute heuristically a good solution to be used as warm start 
    Particularities : no quadratic variables
    """
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M= Big_R+Big_C 
    model = LpProblem(name="max_e_h", sense=LpMaximize)

    # Variables for rows and columns

    lpRows = {row: (LpVariable(f'row_{row}', cat='Binary'), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Binary'), degree) for col, degree in cols_data}
    
    #lpRows = {row: (LpVariable(f'row_{row}', cat='Integer', lowBound=0, upBound=1), degree) for row, degree in rows_data}
    #lpCols = {col: (LpVariable(f'col_{col}', cat='Integer', lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpSum_row = {row: (LpVariable(f'sumRow_{row}', cat='Integer',lowBound=0, upBound=Big_C)) for row, _ in rows_data}
    lpSum_col= {col: (LpVariable(f'sumCol_{col}', cat='Integer',lowBound=0, upBound=Big_R)) for col, _ in cols_data}

    # Objective function: Maximize sum of selected row and column variables
    #model +=  lpSum(lpSum_row), "max_sum_rows_columns" 
    model +=  lpSum(lpSum_col) + lpSum(lpSum_row), "max_sum_rows_columns"  # lpSum(lpSum_row) +
    # #model += lpSum(
    #     [lpvar*deg for lpvar, deg in lpRows.values()] +
    #     [lpvar*deg for lpvar, deg in lpCols.values()]), "max_sum_vertices*degree"
    
    # Constraints for row and column thresholds
    row_threshold = 2
    col_threshold = 2
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    #model += obj_rowsum == lpSum(lpSum_row), "obj_sum_row" # why ????
    #model += obj_colsum == lpSum(lpSum_col), "obj_sum_col" # why ???? 
    #Add row density constraints
    # __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, rho)
    # __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, rho)

    __row_density_h(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_row, delta)
    __col_density_h(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_col, delta)

    return model


def max_e_h_slack(rows_data, cols_data, edges, delta):
    """
    A heuristic proposed to max edge biclique problem.  The particularity is to use an objective that does not require  quadratic variables 
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    
    Goal : compute heuristically a good solution to be used as warm start 
    Particularities : no quadratic variables
    """
    model = LpProblem(name="max_e_h", sense=LpMaximize)

    # Variables for rows and columns

    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpSum_row = {row: (LpVariable(f'sumRow_{row}', cat='Continuous')) for row, _ in rows_data}
    lpSum_col= {col: (LpVariable(f'sumCol_{col}', cat='Continuous')) for col, _ in cols_data}
    lpSlack_row = {row: (LpVariable(f'slackRow_{row}', cat='Continuous')) for row, _ in rows_data}
    lpSlack_col= {col: (LpVariable(f'slackCol_{col}', cat='Continuous')) for col, _ in cols_data}


    obj_rowsum =  LpVariable(f"obj_rowsum", cat=LpContinuous)
    obj_colsum =  LpVariable(f"obj_colsum", cat=LpContinuous)

    # # print()
    # print('-' * 40)
    # print('lpCells=')
    # print(lpCells)
    # print()
    # print('-' * 40)
    # Objective function: Maximize sum of selected row and column variables
    #
    model += lpSum(lpSum_row) + lpSum(lpSum_col), "max_sum_rows_columns" 
    # Constraints for row and column thresholds
    #row_threshold = 2
    #col_threshold = 2
    #model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    #model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += obj_rowsum == lpSum(lpSum_row), "obj_sum_row" # why ????
    model += obj_colsum == lpSum(lpSum_col), "obj_sum_col" # why ???? 
    # Add row density constraints
    #__row_density(rows_data, cols_data, edges, model, lpRows, lpCols, rho)
    #_col_density(rows_data, cols_data, edges, model, lpRows, lpCols, rho)

    __row_density_slack(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_row, lpSlack_row, delta)
    __col_density_slack(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_col, lpSlack_col, delta)

    return model

# def max_e_c(rows_data, cols_data, edges, delta):
#     """
#     Arguments:
#     ----------
#     rows_data: list of tuples (row, degree) of rows in the matrix.
#     cols_data: list of tuples (col, degree) of columns in the matrix.
#     edges: list of tuples (row, col) corresponding to the ones of the matrix.
#     delta: float, error tolerance for density constraints.

#     Returns:
#     --------
#     LpProblem:
#         The ILP model.
    
#     This is an improved model from Chnag et al. 
#     """
#     model = LpProblem(name="max_e_c", sense=LpMaximize)

#     # Variables for rows and columns

#     lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
#                     lowBound=0, upBound=1), degree) for row, degree in rows_data}
#     lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
#                                lowBound=0, upBound=1), degree) for col, degree in cols_data}
#     lpCells = {(row, col): LpVariable(
#                     #f'cell_{row}_{col}', cat="Binary")
#                     f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1) 
#                     for (row, col) in edges
#     }

#     #
#     model += lpSum(lpCells), 'maximize_sum'
#     # Constraints for row and column thresholds
#     row_threshold = 2
#     col_threshold = 2
#     if debug >= 3:
#         print()
#         print('-' * 40) 
#         print(f"******** Solving model  ********", model.name, " with delta = ", delta)
#         print(' # rows_data =', len(rows_data),' # cols_data =', len(cols_data), ' # edges =', len(edges)) 
#         print('row_threshold=', row_threshold )
#         print('col_threshold=', col_threshold )
#         print()
#         print('-' * 40)
#     model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
#     model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
#     for row, col in edges:
#     #for row, col in lpCells:
#         model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
#         model += (lpCols[col][0] >= lpCells[(row, col)]), f'ce ll_{row}_{col}_2'

#         #########################################
#         #compacting  with degree 
#         #########################################
      
#     # for col, _ in cols_data:
#     #     col_edges = [u for u, v in edges if v == col]           
#     #     model += (
#     #         lpSum(lpCells[(row, col)][0] for row in col_edges) <= lpCols[col][1]*lpCols[col][0]
#     #     ), f"col_degre_{col}"
#     # for row, _ in rows_data:
#     #     row_edges = [v for u, v in edges if u == row]           
#     #     model += (
#     #         lpSum(lpCells[(row, col)][0] for col in row_edges) <= lpRows[row][1]*lpRows[row][0]
#     #     ), f"row_degre_{row}"
        
#          #########################################


#     # Add row density constraints
#     __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)
#     __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta)

#     return model 

def __col_density_h(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_col, delta):
    """
    Adds col density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    model: LpProblem to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    rho: float, error tolerance for density constraints.
    """
    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M =Big_R + Big_C
    #print('Big_R=', Big_R)
    #print('Big_C=', Big_C)

    for col, _ in cols_data:
        #print(f"Adding col density constraints for col {col}:")
        col_edges = [u for u, v in edges if v == col]
        model += (lpSum_col[col] ==
                    lpSum(lpRows[row][0] for row in col_edges),f"col_sum_{col}"
        )        
        #print(f"Col edges: {col_edges}") 
        #print(f"lpCols[col][0]: {lpCols[col][0]}")             
        # Constraint for col density upper bound
        #print(f"Sum of edge variables: {lpSum(lpRows[row][0] for row in col_edges)}")
        #print(f"Sum of row variables: {lpSum(lpRows[row][0] for  row, _ in rows_data)}")
        model += (
            lpSum_col[col] - (1 - delta) * lpSum(lpRows[row][0] for row, _ in rows_data)  >= (lpCols[col][0]-1) * Big_M, f"col_slack_{col}" 
        )
        #    >= (lpCols[col][0]-1) * Big_M, f"col_slack_{col}"
        #)
        # #model += (
        #     lpSum_col[col] - (1 - delta) * lpSum(lpRows[row][0] for row, _ in rows_data  >= (lpCols[col][0]-1) * Big_M), f"col_slack_{col}" 
        # #    >= (lpCols[col][0]-1) * Big_M, f"col_slack_{col}" 
        # )
        #model += (lpSlack_col[col] >= (lpCols[col][0]-1) * Big_M, f"col_err_rate_1_{col}"
            # lpSum_col[col] - (1 - rho) * lpSum(lpRows[row][0] for row, _ in rows_data) 
            # >= (lpCols[col][0]-1) * Big_M, f"col_err_rate_1_{col}"
        #)


def __col_density_slack(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_col, lpSlack_col, rho):
    """
    Adds col density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    model: LpProblem to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    rho: float, error tolerance for density constraints.
    """
    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M =Big_R + Big_C
    #print('Big_R=', Big_R)
    #print('Big_C=', Big_C)

    for col, _ in cols_data:
        #print(f"Adding col density constraints for col {col}:")
        col_edges = [u for u, v in edges if v == col]
        model += (lpSum_col[col] ==
                   lpSum(lpRows[row][0] for row in col_edges),f"col_sum_{col}"
        )        
        #print(f"Col edges: {col_edges}") 
        #print(f"lpCols[col][0]: {lpCols[col][0]}")             
        # Constraint for col density upper bound
        #print(f"Sum of edge variables: {lpSum(lpRows[row][0] for row in col_edges)}")
        #print(f"Sum of row variables: {lpSum(lpRows[row][0] for  row, _ in rows_data)}")
        model += (lpSlack_col[col] == 
            lpSum_col[col] - (1 - rho) * lpSum(lpRows[row][0] for row, _ in rows_data), f"col_slack_{col}" 
        #    >= (lpCols[col][0]-1) * Big_M, f"col_slack_{col}"
        )
        model += (lpSlack_col[col] >= (lpCols[col][0]-1) * Big_M, f"col_err_rate_1_{col}"
            # lpSum_col[col] - (1 - rho) * lpSum(lpRows[row][0] for row, _ in rows_data) 
            # >= (lpCols[col][0]-1) * Big_M, f"col_err_rate_1_{col}"
        )

def __row_density_h(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_row, delta):
    """
    Adds row density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    model: LpProblem to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    delta: float, error tolerance for density constraints.
    """

    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M = Big_R + Big_C 
    if debug >= 3:
        print()
        print('-' * 40)
        print(f"row data in row_density_h = {rows_data}")
        print(f"col data in row_density_h = {cols_data}")

    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row]

        # Debugging output
        if debug >= 4:
            print(f"Row_edges for row in row_density_h {row}: {row_edges}")
            print(f"Keys in lpCols in row_density_h: {list(lpCols.keys())}")

        # Avoid KeyError by filtering missing keys
        model += (lpSum_row[row] == lpSum(lpCols[col][0] for col in row_edges if col in lpCols), f"row_sum_{row}")

        model += (
            lpSum_row[row] - (1 - delta) * lpSum(lpCols[col][0] for col, _ in cols_data) >= (lpRows[row][0] - 1) * Big_M,
            f"row_slack_{row}"
        )


def __row_density_slack(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_row, lpSlack_row, rho):
    """
    Adds row density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    model: LpProblem to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    rho: float, error tolerance for density constraints.
    """
    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M= Big_R+Big_C 
    for row, _ in rows_data:
        #print(f"Adding row density constraints for row {row}:")
        row_edges = [v for u, v in edges if u == row]
        model += (lpSum_row[row] == 
                  lpSum(lpCols[col][0] for col in row_edges), f"row_sum_{row}"
        )
        #print(f"Row edges: {row_edges}") 
        #print(f"lpRows[row][0]: {lpRows[row][0]}")             
        # Constraint for row density upper bound
        #print(f"Sum of edge variables: {lpSum(lpCols[col][0] for col in row_edges)}")
        #print(f"Sum of column variables: {lpSum(lpCols[col][0] for  col, _ in cols_data)}")
        model += (lpSlack_row[row]==
            lpSum_row[row] - (1 - rho) * lpSum(lpCols[col][0] for col, _ in cols_data), f"row_slack_{row}" 
           # >= (lpRows[row][0]-1) * Big_M, f"row_err_rate_1_{row}"
        )

        model += (
            lpSlack_row[row] >= (lpRows[row][0]-1) * Big_M, f"row_err_rate_1_{row}"
        )


def __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, delta):
    """
    Adds col density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    model: LpProblem to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    delta: float, error tolerance for density constraints. 
    Attention : in the text it is written as ALPHA TO DO: (ALPHA : = 1-DELTA) !!!!
    """
    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M= Big_R+Big_C 
    for col, _ in cols_data:
        col_edges = [u for u, v in edges if v == col]
        model += (
            lpSum(lpRows[row][0] for row in col_edges) - (1 -  delta) * lpSum(lpRows[row][0] for row, _ in rows_data) >= 
            (lpCols[col][0]-1) * Big_M
        ), f"col_err_rate_1_{col}"


def __row_density(rows_data, cols_data, edges, model, lpRows, lpCols,  delta):
    """
    Adds row density constraints to the model.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    model: LpProblem to add constraints to.
    lpRows: dict of row variables and their degrees.
    lpCols: dict of column variables and their degrees.
    delta: float, error tolerance for density constraints.
    Attention : in the text it is written as ALPHA TO DO: (ALPHA : = 1-DELTA) !!!!
    """
    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M = Big_R+Big_C

    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row]
        model += (
            lpSum(lpCols[col][0] for col in row_edges) - (1 - delta) * lpSum(lpCols[col][0] for col, _ in cols_data) >= 
            (lpRows[row][0]-1) * Big_M
        ), f"row_err_rate_1_{row}"


#============================================================================ #
#                     LP MODEL - MAXIMIZE DESIRED CELL                         #
# ============================================================================ #


def max_Ones(rows_data, cols_data, edges, rho):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros  of the matrix.
    * rho: percentage of accepted zeros 
    """

    # ------------------------------------------------------------------------ #
    # Model with maximization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='maximize_ones', sense=LpMaximize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            if (row, col) in edges:
                lpCells[(row, col)] = (LpVariable(
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)
            else:
                lpCells[(row, col)] = (LpVariable(  
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([cellValue*lpvar for lpvar,
                   cellValue in lpCells.values()]), 'maximize_ones_in_matrix'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #

    for row, col in lpCells:
        model += (lpRows[row][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_1'
        model += (lpCols[col][0] >= lpCells[(row, col)][0]), f'ce ll_{row}_{col}_2'
        #if (row, col) in edges:
        model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'

    model += (lpSum([(1-cellValue)*lpvar for lpvar, cellValue in lpCells.values()]) <= rho *
              lpSum([lpvar for lpvar, _ in lpCells.values()])), f'err_rate'

    return model



def max_Ones_comp(rows_data, cols_data, edges, rho):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros  of the matrix.
    * rho: percentage of accepted zeros 
    """

    # ------------------------------------------------------------------------ #
    # Model with maximization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='maximize_ones_compacted', sense=LpMaximize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            if (row, col) in edges:
                lpCells[(row, col)] = (LpVariable(
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)
            else:
                lpCells[(row, col)] = (LpVariable(  
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)
    # y = LpVariable('yyy', cat=LpBinary,)
    # y_edge = LpVariable('y_edge', cat=LpBinary,)
    # big_m = 1000

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([cellValue*lpvar for lpvar,
                   cellValue in lpCells.values()]), 'maximize_ones_compacted'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    row_threshold = 2
    col_threshold = 2
    model += (lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold), "row_threshold"
    model += (lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold), "col_threshold"
    #
    for row, col in lpCells:
    #for row, col in edges:
        #model += (lpRows[row][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_1'
        #model += (lpCols[col][0] >= lpCells[(row, col)][0]), f'ce ll_{row}_{col}_2'
        model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'
      #  if (row, col) in edges:
    # for (row, col) in edges:
      #      model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'
              #compacting  with degree 
        #########################################
    for col, _ in cols_data:
        col_edges = [u for u, v in edges if v == col]           
        model += (
            lpSum(lpCells[(row, col)][0] for row in col_edges) <= lpCols[col][1]*lpCols[col][0]
        ), f"col_degre_{col}"
    for row, _ in rows_data:
        row_edges = [v for u, v in edges if u == row]           
        model += (
            lpSum(lpCells[(row, col)][0] for col in row_edges) <= lpRows[row][1]*lpRows[row][0]
        ), f"row_degre_{row}"
    #      #########################################

    model += (lpSum([(1-cellValue)*lpvar for lpvar, cellValue in lpCells.values()]) <= rho *
              lpSum([lpvar for lpvar, _ in lpCells.values()])), f'err_rate'

    return model


def KP_QBr(rows_data, col_length, nb_edges_0, debug, rho=0.1):
    """
    In this model, we use knapsack model
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * col_length: number of columns 
    * nb_edges_0: number of zeros in the matrix. Coresponds to the capacity of the knapsack
    * rho: percentage of non deleted zeros 
    """
    row_length = len(rows_data) #number of rows 
    if debug >=3:
        total_degree_0 = sum((col_length - degree) for _, degree in rows_data)
        print('I am currently solving row KP_QBr with Input data : ***************')
        print()
        print('edges =',nb_edges_0, "rho=", rho, "number of columns =", col_length, "number of rows =", row_length,"total_degree_0=", total_degree_0, "RHS=" , rho* nb_edges_0)
        print(' rows_data =', rows_data)
        print('-' * 70) 
    
    # ------------------------------------------------------------------------ #
    # Model with minimization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='row_knapsack_problem', sense=LpMinimize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #

    lpRows = [(LpVariable(f'row_{row}', cat='Continuous',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data]
    # lpCols = [(LpVariable(f'col_{col}', cat='Integer',
    #                 lowBound=0, upBound=1), degree) for col, degree in cols_data]
    #------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([degree*lpvar for lpvar, degree in lpRows]), 'knapsack' #


    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #

    # # delete  at least (1-eps)*nb_edges  zeros
    model += lpSum([(col_length - degree) * lpvar for lpvar, degree in lpRows]) >=  rho * nb_edges_0, 'sensitivity'

    return model

def KP_QBc(cols_data, row_length, nb_edges_0, debug, rho=0.1):
    """
    In this model, we use knapsack model
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    * rho: percentage of accepted undesired cell over accepted desired cell
    """
    col_length = len(cols_data) # # of columns 
    if debug >=2:
        total_degree_0 = sum((row_length - degree) for _, degree in cols_data)
        print('-' * 40)
        print('I will solve column KP_QBc with Input data : ***************')
        print()
        print(" input cols_data =", cols_data)
        print('nb_edges_0 =',nb_edges_0, "rho=", rho,"# of rows =", row_length, "# of columns=", col_length, "total_degree_0=", total_degree_0, "RHS=" , rho* nb_edges_0)
    #print('cols_data =', cols_data)
        print('-' * 40)
    
    # ------------------------------------------------------------------------ #
    # Model with minimization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='column_knapsack_problem', sense=LpMinimize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #

    # lpRows = [(LpVariable(f'row_{row}', cat='Integer',
    #                 lowBound=0, upBound=1), degree) for row, degree in rows_data]
    lpCols = [(LpVariable(f'col_{col}', cat='Continuous',
                    lowBound=0, upBound=1), degree) for col, degree in cols_data]
    # Assign initial values correctly
    # for var, _ in lpCols:  # Unpack tuple (variable, degree)
    #    var.varValue = 0  # Set initial value
    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([degree*lpvar for lpvar, degree in lpCols]), 'knapsack'


    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #

    model += lpSum([(row_length - degree) * lpvar for lpvar, degree in lpCols]) >= rho * nb_edges_0, 'sensitivity'


    return model

#============================================================================
#  greedy approaches
#============================================================================


def zero_cleaner(rows, cols, row_names, col_names, edges_1, nb_edges_0, iter, rho, debug, KP_time): 
    # here we first call KP by rows, then KP by columns 
    rows_res = []
    cols_res = []
    rows_res_name = []
    cols_res_name = []
    rows_del = []
    rows_del_name = []
    cols_del = []
    cols_del_name = []
    col_length = len(cols)
    model = KP_QBr(rows, col_length, nb_edges_0, debug, rho)
    #model.solve(PULP_CBC_CMD(msg=True, timeLimit= 3600, gapRel = 0.5),)
    model.solve(GUROBI_CMD(msg=False, timeLimit= 600)#, options=[("Heuristics", 0.0), ("NoRelHeurTime", 0)] )#,gapRel=0.3)
    )
    #model.solve(GUROBI_CMD(msg=True, timeLimit= 60, MIPGap = 0.03),)
    # Check status
    if model.status == 9:  # GRB.TIME_LIMIT:
        print("Gurobi stopped due to time limit!")
        print("Model is . Exporting LP file for debugging...")
    model.writeLP("debug_model_row.lp")
    obj_value =  value(model.objective)
    if debug >=2:
        print(f"Model status: {LpStatus[model.status]}")
        print(f"Solver Status after model with status :{model.name}, {model.status}")  
        print(f"Computed Objective Value : {obj_value:.3f}")
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
         sys.exit("Terminating program due to infeasibility. EXIT 1")
    if debug >=2:
        print('I solved model name =', model.name, "for iteration i = ", iter, "with obj value:", obj_value, "KP time =", KP_time, 'debug :',  debug)
    if model.name == "row_knapsack_problem":
            for var in model.variables():
                if debug >= 3:
                    print('row_knapsack var name =', var.name,'var value =', var.varValue)
                if var.varValue == 0:
                    if var.name[:3] == "row":
                        rows_res = rows_res + [var.name[4:]]
                    elif var.name[:3] == "col":
                        print('Something wrong. var name =', var.name,'var value =', var.varValue)
                        sys.exit("Terminating program. EXIT 7 ")
                else: # i.e. var.varValue == 1:
                    if var.name[:3] == "row":
                        rows_del = rows_del + [var.name[4:]]
                    elif var.name[:3] == "col":
                        print('Something wrong. var name =', var.name,'var value =', var.varValue)
                        sys.exit("Terminating program. EXIT 7 ")
    if model.name == "column_knapsack_problem":
            sys.exit("Terminating program. Non expected case column_knapsack_problem instead of row_knapsack_problem. Exit 10")  
    if model.name == "row_knapsack_problem":
         cols_res = [str(c) for c, _ in cols]
         cols_res = [int(c) for c in cols_res]
         cols_res_name = [col_names[c] for c in cols_res]
    if model.name == "column_knapsack_problem":
         sys.exit("Terminating program. Non expected case column_knapsack_problem instead of row_knapsack_problem. Exit 10")  
        #  rows_res = [str(c) for c, _ in rows]
        #  rows_res = [int(r) for r in rows_res]
        #  rows_res_name = [row_names[r] for r in rows_res]
    #rows_res = [str(c) for c, _ in rows_res]
    rows_res = [int(r) for r in rows_res]
    if len(rows_res) == 0: 
         #sys.exit("Terminating program due to matrix degeneration. All rows deleted. Decrease the value of rho. EXIT 101")
         if debug >= 1:
            print(f"Matrix degeneration. All rows deleted. Continue with previous rows")
         rows_res = [str(c) for c, _ in rows]
         rows_res = [int(r) for r in rows_res]            
    rows_res_name = [row_names[r] for r in rows_res]
    #cols_res = [str(c) for c, _ in cols]
    #cols_res = [int(c) for c in cols_res]
    #cols_res_name = [col_names[c] for c in cols_res]
    KP_time =  print_log_output(model, KP_time, obj_value, len(rows_res), len(cols_res))
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density  = update_data(rows, cols, edges_1, rows_res, cols_res, debug)
    nb_edges_1_rem = len(edges_1_rem)
    if debug >=2:
        print('-' * 40)
        print(f"""
        Updated data after iteration: {iter} . We solved KP model :  {model.name} with delta =  {delta} 
        Found matrix of size : ({len(rows_rem)}, {len(cols_rem)}) and density : {density:.3f} and number of remaining ones : {nb_edges_1_rem}
        """)
        if debug >= 3:
            print("Remaining  Rows before colling KP col :", rows_rem)
            print("Remaining Columns  before colling KP col :", cols_rem)
        print('-' * 40)
        print()
    nb_rows_rem = len(rows_rem)
       
    model = KP_QBc(cols_rem, nb_rows_rem, nb_edges_0_rem, debug, rho) 

    model.solve(GUROBI_CMD(msg=False, timeLimit= 600))
    # Check status
    if model.status == 9 : #GRB.TIME_LIMIT:
        print("Gurobi stopped due to time limit!")
    #print("Gurobi Model Status:", model.status, "-", gp.GRB.Status.__dict__.get(model.status, "Unknown"))
    model.writeLP("debug_model_col.lp")
    if debug >= 2:
            print(f"Model status: {LpStatus[model.status]}")
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
         sys.exit("Terminating program due to infeasibility. EXIT 7")
    #read the result from the solver
    if debug >=2:
        print('I solved model name =', model.name, "for iteration i = ", iter, "with obj value:", obj_value, "KP time =", KP_time, 'debug :',  debug)
    cols_res=[]
    cols_del=[]
    if model.name == "row_knapsack_problem":
            sys.exit("Terminating program. Non expected case row_knapsack_problem instead of column_knapsack_problem Exit 11")  
    if model.name == "column_knapsack_problem":
            for var in model.variables():
                if debug >= 3:
                    print('column_knapsack var name =', var.name,'var value =', var.varValue)
                if var.varValue == 0:
                    if var.name[:3] == "col":
                        cols_res = cols_res + [var.name[4:]]
                    elif var.name[:3] == "row":
                #    cols_res = cols_res + [var.name[4:]]
                        print('Something wrong. var name =', var.name,'var value =', var.varValue)
                        return model 
                else: # i.e. var.varValue == 1:
                    if var.name[:3] == "col":
                        cols_del = cols_del + [var.name[4:]]
                    elif var.name[:3] == "row":
                        #cols_del = cols_del + [var.name[4:]]
                        print('Something wrong. var name =', var.name,'var value =', var.varValue)
                        return model 
    # update  edges      
    KP_time =  print_log_output(model, KP_time, obj_value, len(rows_res), len(cols_res)) 
    if len(cols_res) == 0 : 
         sys.exit("Terminating program due to matrix degeneration. All columns deleted.")               
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density = update_data(rows_rem, cols_rem, edges_1_rem, rows_res, cols_res, debug)
    if debug >=2:
        print('-' * 40)
        print(f"""
        Updated data after iteration: {iter} . We solved KP model :  {model.name} with delta =  {delta} 
        Found matrix of size : ({len(rows_rem)}, {len(cols_rem)}) and density : {density:.3f} and number of remaining ones : {nb_edges_1_rem}
        """)
    if debug >=2:
        print()
        print('-' * 40)
        # print(f"Updated data after solving model = {model.name}", " Density :", density , " delta =  ", delta, "iteration:", iter)
        # nb_edges_1_rem = len(edges_1_rem)
        # print("Number of Remaining  Rows  :", len(rows_rem), "Number of Remaining Columns :", len(cols_rem))
        # print("Number of Remaining Edges_0 : ", nb_edges_0_rem, "Number of Remaining Edges_1 : ", nb_edges_1_rem)
        if debug >= 3:
            print("Remaining  Rows after solving KP col :", rows_rem)
            print("Remaining Columns after solving KP col :", cols_rem)
        print('-' * 40)
        print()
    #sys.exit("Terminating program after KP col") 
    # end of zero_cleaner  
    return rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density, KP_time 


# ============================================================================ #
#                                    SOLVING                                   #
# ============================================================================ #


def solve(prev_lower_bound, dec_conq, matrix_name, rows, cols, edges_1, model, KP_time, QBC_time, rho=1.0, delta=0.1, density_threshold=0.87):
    """
    Function to solve the maximum biclique problem, this function reads the data,
    create a LP model using the data and return a list of rows and a list of 
    columns as a result. 
    ARGUMENTS:
    ----------
    * dec_conq: the level of decomposition and conquering.
    * matrix_name: the name of the matrix.
    * rows (list of tuples): List of (row_index, degree) for rows.
    * cols (list of tuples): List of (col_index, degree) for columns.
    * edges_1 (list of tuples): List of existing edges (row_index, col_index).
    * path_to_data: the path to the csv file.
    * model: the model to be used.
    * rho: zero deletion value used in all greedy approaches 
    * delta: error tolerance (percentage of accepted zeros) in the final result submatrix
    """ 
    nbi_0, nbi_1, sparsity, density = density_calcul(rows, cols)
    if debug >= 1: 
        print()
        print('-' * 70)
        print(f"***Stats Current Input for matrix {matrix_name} in  {file_path} at level {dec_conq} and with model:  {model}***")
        print("Size of current matrix : ", len(rows), "*", len(cols), "=", len(rows) * len(cols), "; dec_conq:", dec_conq)
        print("Prev_lower_bound : ",prev_lower_bound,"number input zeros : ",nbi_0, "; number input ones : ",nbi_1)            
        print("rho = ",rho, "; delta : ", delta)
        # Safe printing: Convert None to "N/A" or provide a default value
        print(f"Input density : {density:.3f}" if density is not None else "Input density: N/A",  f"; density_threshold: {density_threshold:.5f}")
        #print(f"Input density : {density:.3f}, ; density_threshold: {density_threshold}")
        if debug >= 3 : 
            print("rows =", rows)
            print("cols =", cols)
            print("row_names =", row_names)
            print("col_names =", col_names)
            print("edges_1 =", edges_1)
        print('-' * 70)
        print('End Current  Stats')
        print('-' * 70)
        print()

    # end fetching input data
    # start solving the problem  
    nb_edges_0 = len(edges_compl)
    model_name = model
    rows_in = rows
    cols_in = cols
    edges_1_in = edges_1
    nb_edges_0_in =  nb_edges_0 
    row_names_in = row_names
    col_names_in = col_names 
    matrix_limit = 2 
    iter = 0
    kp_density = 0
    if density == None:
                sys.exit("Terminating program due to density == None")
    while density < density_threshold and len(rows_in) > matrix_limit and len(cols_in) > matrix_limit:
        iter+=1
        if debug >= 2:
            print()
            print(f"calling greedy approaches for zero deletion, density= {density:.3f} density  and density_threshold= {density_threshold:.3f}")
            print()
        #sys.exit("Terminating program because not done EXIT 0")
        if debug >= 2:
            print()
            print("I am in the while loop with i=", iter , "density=", density, "and fixed_threshold=", density_threshold)
            print()
        rows_res, cols_res, edges_1_res, nb_edges_0_res, density, KP_time  = zero_cleaner(rows_in, cols_in, row_names, col_names, edges_1_in, nb_edges_0_in, iter, rho, debug, KP_time)
        rows_in = rows_res
        cols_in = cols_res
        edges_1_in = edges_1_res
        nb_edges_0_in =  nb_edges_0_res
        rows_res_in = [str(c) for c, _ in rows_in]
        rows_res_in = [int(r) for r in rows_res_in]
        row_names_in = [row_names[r] for r in rows_res_in]
        cols_res_in = [str(c) for c, _ in cols_in]
        cols_res_in = [int(r) for r in cols_res_in]
        col_names_in = [col_names[r] for r in cols_res_in]
    else:
        if iter > 0:
            kp_density = density 
            nb_kp_rows = len(rows_res)
            nb_kp_cols = len(cols_res)
            nb_kp_ones = len(edges_1_res)
        if debug >= 1 and iter > 0 :
            print('-' * 40)
            #print()
            print(f"""
            End of greedy approaches. Exit while loop. I did: {iter} iterations
            Density : {density:.3f} > {density_threshold:.3f}
            It took me : {KP_time:.3f} time
            Found matrix of size : ({len(rows_res)}, {len(cols_res)}) and density : {density:.3f}
            """)
            if debug == 3 and iter > 0: 
                print('-' * 40)
                print("rows_in =", rows_in)
                print("cols_in =", cols_in)
                print("row_names_in, =", row_names_in)
                print("col_names_in, =", col_names_in)
            if debug == 4 and iter > 0: 
                print('-' * 40)
                print("edges_1_in =", edges_1_in)
        #sys.exit("Terminating program in while loop and because of density. we printed the current data")  
    #else:
    if debug >= 1: 
        print('-' * 40)
        print()
        print(' Calling exact approaches for QB clique discovery with delta =', delta)
        print('-' * 40)
        print()
    if debug >= 3: 
            print(' Start exact approaches with :')
            print("rows_in =", rows_in)
            print("cols_in =", cols_in)
            print("row_names_in =", row_names_in)
            print("col_names_in =", col_names_in)
            print("edges_1_in =", edges_1_in)
    if debug == 4: 
            print("edges_1_in =", edges_1_in)
    #sys.exit("Terminating program before calling exact exit 10 !!!!")
    if model == "max_e_c" or dec_conq >= 1: 
        rows_res, cols_res, density, nb_edges_1, QBC_time_h, QBC_time_g = warm_exact(prev_lower_bound,dec_conq, matrix_name,model_name, rows_in, cols_in, row_names, col_names, edges_1_in, delta, debug, QBC_time)
        #sys.exit("Terminating program before calling exact exit 10 !!!!")  
    else:
        rows_res, cols_res, density, nb_edges_1, QBC_time_g = exact(dec_conq, matrix_name, model_name, rows_in, cols_in, row_names, col_names, edges_1_in, delta, debug, QBC_time)
        QBC_time_h = 0.0 
    if debug >=1:
        print('-' * 40)
        print(f"""
        Exit from the exact approach: {model_name} for matrix {matrix_name}
        with delta =  {delta}
        Found matrix of size : ({len(rows_res)}, {len(cols_res)})
        and density : {density:.3f}
        and # of ones : {nb_edges_1}
        and Heuristic QBC time   : {QBC_time_h:.5f}
        and Global QBC  time  : {QBC_time_g:.5f}
        """)
    if debug >=2:
        print(" Remaining Rows  :", rows_res )
        print(" Remaining  Cols  :", cols_res )
        #sys.exit("Terminating program because not done EXIT -1 ")
    if iter == 0:
        kp_density = 0
        nb_kp_rows = 0
        nb_kp_cols = 0
        nb_kp_ones = 0 
    return  rows_res, cols_res, density, nb_edges_1, iter, KP_time, kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g 
#######################################################################
    #END OF SOLVE
##################################################################

def exact(dec_conq, matrix_name, model_name, rows, cols, row_names, col_names, edges_1, delta, debug, QBC_time):
    """
    Arguments:
    ----------
    rows: list of tuples (row, degree) of rows in the matrix.
    cols: list of tuples (col, degree) of columns in the matrix.
    edges_1: list of tuples (row, col) corresponding to the ones of the matrix.
    delta: float, error tolerance for density constraints.

    Returns:
    --------
    rows_res: list of tuples (row, degree) of rows in the submatrix.
    cols_res: list of tuples (col, degree) of columns in the submatrix.
    density: float, density of the submatrix.
    nb_edges_1: number of ones in the submatrix.
    QBC_time: float, time taken to solve the model.
    """  
    obj_value = None  # ✅ Initialize to None
    rows_res = []
    cols_res = []

    # # Convert rows and cols into dictionaries for fast lookup
    # row_degree_map = {row: degree for row, degree in rows}
    # col_degree_map = {col: degree for col, degree in cols}

    # Select the appropriate model
    if model_name == 'max_e':
        model = max_e(rows, cols, edges_1, delta)
    elif model_name == 'max_e_h':
        model = max_e_h(rows, cols, edges_1, delta)
    elif model_name == 'max_e_r':
        model = max_e_r(rows, cols, edges_1, delta)
    elif model_name == 'max_v':
        model = max_v(rows, cols, edges_1, delta)
    
    try:
        # Solve the model with Gurobi
        model.solve(GUROBI_CMD(msg=False, timeLimit=timelimit))

        if debug >= 2:
            print('-' * 70)
            print(f"Model status: {LpStatus[model.status]}")    
        # Check if the model has an optimal/feasible solution
        if model.status in [1, 2, 9]:  # 1 = Optimal, 2 = Feasible, 9 = Time limit reached
            obj_value = value(model.objective)  # ✅ Extract objective value           
            if debug >= 2:
                print('-' * 70)
                print(f"Computed Objective Value: {obj_value:.3f}")          
            # Extract solution values (only nonzero rows and columns)
            solution = {
                var.name: var.varValue
                for var in model.variables()
                if var.varValue >= 0.5
            }

            # Convert extracted indices back to tuples (index, degree)
            row_indices = [int(var_name[4:]) for var_name in solution if var_name.startswith("row")]
            col_indices = [int(var_name[4:]) for var_name in solution if var_name.startswith("col")]

            rows_res = [var_name[4:] for var_name in solution if var_name.startswith("row")]
            cols_res = [var_name[4:] for var_name in solution if var_name.startswith("col")]
            # rows_res = [(index, row_degree_map[index]) for index in row_indices]
            # cols_res = [(index, col_degree_map[index]) for index in col_indices]

    except Exception as e:
        print(f"Error during solving: {e}")
    
    # Print model status
    if debug >= 1:
        print(f"Model status: {LpStatus.get(model.status, 'Unknown Status')}")
    
    # ✅ Ensure obj_value is valid before using it
    if obj_value is None:
        print("Warning: obj_value is None. Assigning default value 0.")
        obj_value = 0  # ✅ Fallback value if model failed

    # ✅ Only proceed with saving the solution if the model found a feasible solution
    if model.status in [1, 2, 9]: 
        if model.status == 9: 
            print()
            print(f" !!!!! When soving matrix {matrix_name} with  model {model_name} the time limit of {timelimit} has been reached ")
            print()  
        # Save the solution to a CSV file 
        if p.lower().endswith('.csv'):
            file_path_no_ext = file_path.replace(".csv", "").replace("data/", "")
        if p.lower().endswith('.txt'):
         file_path_no_ext = file_path.replace(".txt", "").replace("data/", "")
        solution_file = f"Experiments/{file_path_no_ext}/results_{dec_conq}_M_{matrix_name}.csv"

        with open(solution_file, mode="w", newline="") as file:
            writer = csv.writer(file)

            # Write the objective value
            if model_name == 'max_e_h':
                lower_bound = len(rows_res)*len(cols_res)
                writer.writerow(["Objective Value", lower_bound])
            else:
                writer.writerow(["Objective Value", obj_value])

            # Write number of rows
            writer.writerow(["# of rows", len(rows_res)])
            # Write number of columns
            writer.writerow(["# of columns", len(cols_res)])

            # Write variable values
            writer.writerow(["Variable", "Value"])
            for var_name, var_value in solution.items():
                if var_name.startswith("row") or var_name.startswith("col"): 
                    writer.writerow([var_name, var_value])

        print(f"Solution saved to {solution_file}")
    
    #return rows_res, cols_res, obj_value

    else:
        print("No feasible solution found.")
        print("Model is infeasible. Exporting LP file for debugging...")
        model.writeLP("debug_model_when_infeasible.lp")

    #if not rows_res or not cols_res :  # Equivalent to checking len(rows_res) == 0
    if not cols_res :  # Equivalent to checking len(cols_res) == 0
        if debug >= 1:
            print(f"Matrix degeneration. All columns deleted. Continue with previous columns")
        cols_res = [str(c) for c, _ in cols]
        cols_res = [int(r) for r in cols_res]
    rows_res = [str(c) for c in rows_res]
    rows_res = [int(r) for r in rows_res]
    cols_res = [str(c) for c in cols_res]
    cols_res = [int(c) for c in cols_res]

    # Convert from 0-based to the original row/column names
    # row_names_res = [row_names[r] for r in row_indices if 0 <= r < len(row_names)]
    # col_names_res = [col_names[c] for c in row_indices if 0 <= c < len(col_names)]
    row_names_res = [row_names[r] for r in rows_res if 0 <= r < len(row_names)]
    col_names_res = [col_names[c] for c in cols_res if 0 <= c < len(col_names)]
    row_names_res = [int(r) if isinstance(r, (np.int64, np.int32)) else r for r in row_names_res]
    if debug >= 3: 
        print("\n-- Debugging Step: checking extracted solution after solving model **** --", model.name )
        print('len_rows_res=', len(rows_res))
        print('row_res=', rows_res)
        #print('len_rows_del=', len(rows_del))
        #print('rows_del=',len(rows_del))
        print('len_cols_res=', len(cols_res))
        print('cols_res=', cols_res)
        #print('len_cols_del=', len(cols_del))
        #print('cols_del=', len(cols_del))
        print("nb row_names, =", len(row_names))
        print("row_names_res =", row_names_res)
        print("col_names_res =", col_names_res)
        print(" nb col_names, =", len(col_names)) 
        print()
        print('-' * 40)
    # cols_res_set = set(map(int, cols_res))
    # rows_res_set = set(map(int, rows_res))
    QBC_time = print_log_output(model, QBC_time, obj_value, len(rows_res), len(cols_res))
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density  = update_data(rows, cols, edges_1, rows_res, cols_res, debug)
    # rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density  = update_data(rows, cols, edges_1, row_indices, col_indices, debug) 
    nb_edges_1 = len(edges_1_rem)
    if debug >=1:
        print()
        print('-' * 40)
        print(f"Results from update_data after solving model = {model.name}", " with delta =  ", delta, "and dec_conq= ", dec_conq)
        print("Number of Remaining  Rows  :", len(rows_rem))
        print("Number of Remaining Columns :", len(cols_rem))
        print("Remaining  number Edges_0 P:", nb_edges_0_rem, "Remaining  number Edges_1 :", nb_edges_1, "Density :", density, "current obj value", obj_value )
        print('-' * 40)
        print()
    if debug >=2:
          print('-' * 40)
          print()
          print('Exit from the heuristic (or exact  approach) ',  model_name,' with delta=', delta, 'Found matrix with rows_res of lenght =',len(rows_res), ' and cols_res of lenght =',len(cols_res)) #, 'and density = ', density)
          if debug >= 3:
                    #print(" Density of the found matrix =  :", density )
                    print(" Original Rows  :", rows )
                    print(" Original Cold  :", cols )
                    print(" Remaining Rows  :", rows_res )
                    print(" Remaining  Cols  :", cols_res )
                    print(" Remaining Rows with degree :", rows_rem )
                    print(" Remaining  Cols with degree :", cols_rem )
    
    return rows_rem, cols_rem, density,  nb_edges_1, QBC_time 
###############################################################################
# END OF EXACT 
# ########################################################################  
def warm_exact(prev_lower_bound, dec_conq, matrix_name,model_name, rows, cols, row_names, col_names, edges_1, delta, debug, QBC_time):
    """
    Arguments:
    ----------
    rows: list of tuples (row, degree) of rows in the matrix.
    cols: list of tuples (col, degree) of columns in the matrix.
    edges_1: list of tuples (row, col) corresponding to the ones of the matrix.

    Returns:
    --------
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density 
        
    """     
    # elif model_name == 'max_e_h':
    obj_value = None  # ✅ Initialize to None
    rows_res = []
    cols_res = []
    # obj_total =  0.0
    #############################
    if debug >= 1:
        print('-' * 40)
        print()
        print("I am in warm_exact before calling max_e_h $$$$$$$$$$$$$$$$$$")
        print()
    model = max_e_h(rows, cols, edges_1, delta)
    try:
        # Solve the model  max_e_h  with Gurobi
        model.solve(GUROBI_CMD(msg=False, timeLimit= timelimit))#, options=[("MIPGap", 0.03)]))
        # Print model status
        if debug >= 2:
            print('-' * 70)
            print(f"Model status: {LpStatus[model.status]}")
        # Extract solution values (only nonzero rows and columns)
        if model.status in [1, 2, 9]:  # 1 = Optimal, 2 = Feasible, 9 = Time limit reached
            obj_value = value(model.objective) # ✅ Extract objective value
            if debug >= 2:
                print('-' * 70)
                print(f"Computed by max_e_h Objective Value : {obj_value:.3f} ")
                #print(f"Computed by max_e_h Objective Value : {obj_value}")
            # Extract solution values (only nonzero rows and columns)
            solution = {
                var.name: var.varValue
                for var in model.variables()
                if var.varValue >= 0.5#and var_name.startswith("row") or var_name.startswith("col")
            }
            # Read results and classify into rows and columns
            # row_indices = [int(var_name[4:]) for var_name in solution if var_name.startswith("row")]
            # col_indices = [int(var_name[4:]) for var_name in solution if var_name.startswith("col")]

            rows_res = [var_name[4:] for var_name in solution if var_name.startswith("row")]
            cols_res = [var_name[4:] for var_name in solution if var_name.startswith("col")]
    except Exception as e:
        print(f"Error during solving: {e}")
    # Print model status
    if debug >= 1:
         print(f"Model status: {LpStatus.get(model.status, 'Unknown Status')}")
    # ✅ Ensure obj_value is valid before using it
    if obj_value is None:
        print("Warning: obj_value is None. Assigning default value 0.")
        obj_value = 0  # ✅ Fallback value if model failed

    # ✅ Only proceed with saving the solution if the model found a feasible solution
    if model.status in [1, 2, 9]:
        # Save the solution to a CSV file
        file_path_no_ext = file_path.replace(".csv", "").replace("data/", "")
        solution_file = f"Experiments/{file_path_no_ext}/results_h_{dec_conq}_M_{matrix_name}.csv"
        with open(solution_file, mode="w", newline="") as file:
            writer = csv.writer(file)
            # Write the objective value
            #if model_name == 'max_e_h':
            lower_bound = len(rows_res)*len(cols_res)    
            writer.writerow(["Obj value",lower_bound ])
            # else:
            #     writer.writerow(["Objective Value", obj_value])
            # Write number of rows
            writer.writerow(["# of rows", len(rows_res)])
            # Write number of columns
            writer.writerow(["# of columns", len(cols_res)])
            # Write variable values
            writer.writerow(["Variable", "Value"])
            for var_name, var_value in solution.items():
                 if var_name.startswith("row") or var_name.startswith("col"): 
                    writer.writerow([var_name, var_value])

        print(f"Solution saved to {solution_file}")

    else:
        print("No feasible solution found.")
        print("Model is infeasible. Exporting LP file for debugging...")
        model.writeLP("debug_model_infeasible.lp")
        #sys.exit("Terminating program due to infeasibility. EXIT 1")

    if not rows_res or not cols_res :  # Equivalent to checking len(rows_res) == 0
       if debug >= 1:
            print("Attention!  max_e_h returns degenerated matrix. Continue without max_e_h")
       rows_res = [str(c) for c,_ in rows]
       rows_res = [int(r) for r in rows_res]
       cols_res = [str(c) for c,_ in cols]
       cols_res = [int(c) for c in cols_res]
       #sys.exit("Terminating program due to matrix degeneration (all rows or all columns have been deleted). EXIT 6 ")
    else:
        rows_res = [str(c) for c in rows_res]
        rows_res = [int(r) for r in rows_res]
        cols_res = [str(c) for c in cols_res]
        cols_res = [int(c) for c in cols_res]
    # Convert from 0-based to the original row/column names
    row_names_res = [row_names[r] for r in rows_res if 0 <= r < len(row_names)]
    col_names_res = [col_names[c] for c in cols_res if 0 <= c < len(col_names)]
    row_names_res = [int(r) if isinstance(r, (np.int64, np.int32)) else r for r in row_names_res]
    if debug >= 3: 
        print("\n-- Debugging Step: checking extracted solution after solving model**** --", model.name )
        print('len_rows_res=', len(rows_res))
        print('row_res=', rows_res)
        print('len_cols_res=', len(cols_res))
        print('cols_res=', cols_res)
        print("nb row_names_res, =", len(row_names_res))
        print("row_names_res =", row_names_res)
        print("col_names_res =", col_names_res)
        print(" nb col_names_res, =", len(col_names_res)) 
        print()
        print('-' * 40)
    # now we call max_e_r with rows_res and cols_res

    # cols_res_set = set(map(int, cols_res))
    # rows_res_set = set(map(int, rows_res))
 
    QBC_time_h = print_log_output(model, QBC_time, obj_value, len(rows_res), len(cols_res))
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density  = update_data(rows, cols, edges_1, rows_res, cols_res, debug) 
    if delta == 0:
        nb_edges_1 = len(rows_res)*len(cols_res)
    else:
        nb_edges_1 = len(edges_1_rem)
    if debug >=1:
        print()
        print('-' * 40)
        print(f"Results from updating data after solving model = {model.name}", " delta =  ", delta)
        print("Number of Remaining number  Rows  :", len(rows_res))
        print("Number of Remaining number Columns :", len(cols_rem))
        print("Remaining  number Edges_0 P:", nb_edges_0_rem, "Remaining  number Edges_1 :", nb_edges_1, "Density :", density , "current obj value", obj_value)
        print('-' * 40)
        print()
    #nb_rows_rem = len(rows_rem) 
    if debug >=2:
          print('-' * 40)
          print()
          print('Exit from the heuristic (or exact  approach) ',  model_name,' with delta=', delta, 'Found matrix with rows_res of lenght =',len(rows_res), ' and cols_res of lenght =',len(cols_res)) #, 'and density = ', density)
          if debug >= 3:
                    #print(" Density of the found matrix =  :", density )
                    print(" Original Rows  :", rows )
                    print(" Original Cold  :", cols )
                    print(" Remaining Rows  :", rows_res )
                    print(" Remaining  Cols  :", cols_res )
                    print(" Remaining Rows with degree :", rows_rem )
                    print(" Remaining  Cols with degree :", cols_rem )

    #sys.exit("Terminating program before callong  max_e_wr EXIT 11")
    if dec_conq >= 1: # the task is not going to be solved but to saved in a file for further solving
        #density= 1
        #nb_edges_1 = len(rows_res)*len(cols_res)
        #return rows_res, cols_res, density, nb_edges_1, QBC_time_h, QBC_time_h 
        #return rows_res, cols_res, density, nb_edges_1, QBC_time_h, QBC_time_h
        return rows_res, cols_res, density, nb_edges_1, QBC_time_h, QBC_time_h  
    #the task is going to be solved 
    #model.Params.Cutoff = lower_bound
    #if prev_lower_bound != None: 
    if prev_lower_bound >= 1+ nb_edges_1:
        lower_bound = prev_lower_bound
        if debug >= 1:
            print(f" A lower_bound = {lower_bound} has been previously found and given as input. ")
    else:   
        lower_bound = nb_edges_1 
        if debug >= 1:
            print(f" A lower_bound = {lower_bound} has been found by max_e_h ")
    if debug >= 1:
        print(f" Start computing max_er with lower_bound =  {lower_bound}. ")
    model = max_e_wr(rows, cols, edges_1, rows_res, cols_res, lower_bound, delta)
    try:
        # Solve the model with Gurobi
        model.solve(GUROBI_CMD(msg=False, timeLimit= timelimit))#, options=[("MIPGap", 0.03)]))
        #, options=[("MIPStart", 1), ("Heuristics", 0.0), ("NoRelHeurTime", 0)] )#,gapRel=0.3))
        # Print model status
        # if model.status == GRB.CUTOFF:
        #     print('-' * 70)
        #     print("No solution found above the found by max_e_h lower bound!!! This lower_bound is considered as the solution")
        #     print('-' * 70) 
        if debug >= 2:
            print('-' * 70)
            print(f"Model status: {LpStatus[model.status]}")
        # Extract solution values (only nonzero rows and columns)
        if model.status in [1, 2, 9]:
            obj_value = value(model.objective) # ✅ Extract objective value
            if debug >= 2:
                print('-' * 70)
                print(f"Computed Objective Value: {obj_value:.3f}")
            # Extract solution values (only nonzero rows and columns)
            solution = {
                var.name: var.varValue
                for var in model.variables()
                    if var.varValue >= 0.5 # and var_name.startswith("row") or var_name.startswith("col")
            }
            # Read results and classify into rows and columns
            rows_res = [var_name[4:] for var_name in solution if var_name.startswith("row")]
            cols_res = [var_name[4:] for var_name in solution if var_name.startswith("col")]
    except Exception as e:
        print(f"Error during solving: {e}")
    if debug >= 2:
            print('-' * 70)
            print(f"Model status: {LpStatus[model.status]}")
    # Check if the model has an optimal/feasible solution
    QBC_time_g = print_log_output(model, QBC_time_h, obj_value, len(rows_res), len(cols_res))
    if  model.status == 0:
        print("*****Model in warm start is infeasible.!!! The given lower bound cannot be improved. Return the below data*** ")
        #the solution found by the heuristic cannot be improved. Return it as the final solution
        if debug >=3:
            print(f"rows_rem: {rows_rem}, \n  cols_rem: {cols_rem},  \n  rows_res: {rows_res},  \n  cols_res: {cols_res},  \n  lower_bound  {lower_bound}")
        return rows_rem, cols_rem, density, lower_bound, QBC_time_h, QBC_time_g 
    if debug >= 1:
        print("*****Model in warm start is feasible. The solution has been improved !!!*** ")
        print(f"model status = {model.status}, LPstatus {LpStatus[model.status]}")
    if model.status in [1, 2, 9]:  # 1 = Optimal, 2 = Feasible, 9 = Time limit reached
        #if model.status == 9:
        if LpStatus[model.status] == 9:
            print()
            print(f" !!!!! When soving matrix {matrix_name} with  model {model_name} the time limit of {timelimit} has been reached ")
            print() 
    # Save the solution to a file
        file_path_no_ext = file_path.replace(".csv", "").replace("data/", "")
        solution_file = f"Experiments/{file_path_no_ext}/results_wstart_{dec_conq}_M_{matrix_name}.csv"
        with open(solution_file, mode="w", newline="") as file:
            writer = csv.writer(file)
            # Write the objective value
            writer.writerow(["Objective Value", obj_value])
            # Write number of rows
            writer.writerow(["# of rows", len(rows_res)])
            # Write number of columns
            writer.writerow(["# of columns", len(cols_res)])
            # Write variable values
            writer.writerow(["Variable", "Value"])
            for var_name, var_value in solution.items():
                if var_name.startswith("row") or var_name.startswith("col"): 
                    writer.writerow([var_name, var_value])

        print(f"Solution saved to {solution_file}")

    else:
        print("No feasible solution found.")
        print("Model is infeasible. Exporting LP file for debugging...")
        model.writeLP("debug_model_if_infeasible.lp")
        #sys.exit("Terminating program due to infeasibility. EXIT 10")

   #Reading and analyse of the obtained results")
    #read the results from the solver

    if not rows_res or not cols_res :  # Equivalent to checking len(rows_res) == 0
       print("matrix degenerated (all rows or all columns have been deleted)")
       sys.exit("Terminating program due to matrix degeneration (all rows or all columns have been deleted). EXIT 7 ")

    rows_res = [str(c) for c in rows_res]
    rows_res = [int(r) for r in rows_res]
    cols_res = [str(c) for c in cols_res]
    cols_res = [int(c) for c in cols_res]
    # Convert from 0-based to the original row/column names
    row_names_res = [row_names[r] for r in rows_res if 0 <= r < len(row_names)]
    col_names_res = [col_names[c] for c in cols_res if 0 <= c < len(col_names)]
    row_names_res = [int(r) if isinstance(r, (np.int64, np.int32)) else r for r in row_names_res]
    #row_names_res = [int(r) for r in row_names_res]
    # Print results
    #print("row_names_res =", row_names_res)
    #print("col_names_res =", col_names_res)

    #rows_res_name = [row_names[r] for r in rows_res]
    #cols_res_name = [col_names[c] for c in cols_res]
    print()
    if debug >= 2: 
    #print("\n-- Debugging Step: -after row KP solving **** --")
    # print('rows=', rows)
    # print('cols=', cols)
        print("\n-- Debugging Step: checking extracted solution after solving model**** --", model.name )
        print('len_rows_res=', len(rows_res))
        print('row_res=', rows_res)
        #print('len_rows_del=', len(rows_del))
        #print('rows_del=',len(rows_del))
        print('len_cols_res=', len(cols_res))
        print('cols_res=', cols_res)
        #print('len_cols_del=', len(cols_del))
        #print('cols_del=', len(cols_del))
        print("nb row_names_res, =", len(row_names_res))
        print("row_names_res =", row_names_res)
        print("col_names_res =", col_names_res)
        print(" nb col_names_res, =", len(col_names_res)) 
        print()
        print('-' * 40)
    # now we call max_e_r with rows_res and cols_res
    #sys.exit("Terminating program due undone tasks. EXIT 01")
    cols_res_set = set(map(int, cols_res))
    rows_res_set = set(map(int, rows_res))
    #QBC_time_g = print_log_output(model, QBC_time_h, obj_value, len(rows_res), len(cols_res))
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density  = update_data(rows, cols, edges_1, rows_res, cols_res, debug) 
    nb_edges_1 = len(edges_1_rem)
    # recently added 
    nb_edges_0_rem, nb_edges_1_rem, sparsity, density = density_calcul(rows_rem, cols_rem)
    if debug >=1:
        print()
        print('-' * 40)
        print(f" results from updating data after solving model = {model.name}", " delta =  ", delta)
        print("Number of Remaining  Rows  :", len(rows_rem))
        print("Number of Remaining number Columns :", len(cols_rem))
        print("Remaining  number Edges_0 P:", nb_edges_0_rem, "Remaining  number Edges_1 :", nb_edges_1, "Density :", density )
        print('-' * 40)
        print()
    nb_rows_rem = len(rows_rem) # - len(rows_del) 
    if debug >=2:
          print('-' * 40)
          print()
          print('Exit from the exact  approach ',  model_name,' with delta=', delta, 'Found matrix with rows_res of lenght =',len(rows_rem), ' and cols_res of lenght =',len(cols_rem), 'and density =', density)
          if debug >= 3:
                    print(" Density of the found matrix =  :", density )
                    print(" Remaining Rows  :", rows_res )
                    print(" Remaining  Cols  :", cols_res )
                    print(" Remaining Rows with degree :", rows_rem )
                    print(" Remaining  Cols with degree :", cols_rem )

    return rows_rem, cols_rem, density, nb_edges_1, QBC_time_h, QBC_time_g 
###############################################################################
# END OF WARM EXACT 
####################################################################

def update_data(rows_data, cols_data, edges_1, rows_res, cols_res, debug):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges_1: list of tuples (row, col) corresponding to the ones of the matrix.
    Updates  rows_data, cols_data, edges_1 and returns them in rows_rem, cols_rem, edges_1_rem. To compute that we use rows_res, cols_res information. The last one corresponds to the remaining rows and columns in the matrix (i.e., the solution of the ILP model).

    Returns:
    --------
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density 
        
    """
    
    # first step
    # Convert rows_del to a set for faster lookup
    # print("\n-- Update Debugging Step: edges  --")
    # print("input edges_0 =", edges_0)
    # print("input edges_1 =", edges_1)

    #rows_del_set = set(map(int, rows_del))
    #cols_del_set = set(map(int, cols_del))

    cols_res_set = set(map(int, cols_res))
    rows_res_set = set(map(int, rows_res))

    # Step 1: Compute deleted edges (edges connected to rows or cols in rows_del_set and cols_del_set)
    #edges_0_del= [edge for edge in edges_0 if edge[0] in rows_del_set or edge[1] in cols_del_set]
    edges_1_del= [edge for edge in edges_1 if edge[0] not in rows_res_set or edge[1] not in cols_res_set]

    #print("\n-- Debugging Step: Edges to be removed  --")
    # print("edges_0_del =", edges_0_del)
    #print(" # Deleted edges_1_del =", len(edges_1_del))
    

    # Step 2: Compute edges_new (remaining edges) 
    #edges_0_rem = [edge for edge in edges_0 if edge[0] in rows_res_set and edge[1] in cols_res_set]
    edges_1_rem = [edge for edge in edges_1 if edge[0] in rows_res_set and edge[1] in cols_res_set]
    #print("\n-- Debugging Step: Remaining Edges --")
    # print("Remainig Edges_0 =",  edges_0_rem)
    #print("Remainig edges_1_rem =",  len(edges_1_rem))
    if debug >= 4:
        print("\n-- Update Debugging Step 1 : input data  and first updates --")
    # print("rows_del_set =", rows_del_set)
    # print("cols_del_set =", cols_del_set)
    # print("\n-- Update Debugging Step: input rows_col_remaining_sets --")
        print("rows_res =", rows_res)
        print("cols_res =", cols_res) 
       #print("edges 1 =", edges_1)  
        print("rows_res_set =", rows_res_set) 
        print("cols_res_set =", cols_res_set)
        #print("  edges_1_del  =",  edges_1_del) 
        #print("  edges_1_rem  =",  edges_1_rem)
        print("rows_data =", rows_data) 
        print("cols_data =", cols_data) 

 # Step 3: Update column degrees in cols_data
    col_degree_map = {col: degree for col, degree in cols_data}
    row_degree_map = {row: degree for row, degree in rows_data}
    # print("input cols_data_degree_map =", col_degree_map)
    # print("input rows_data_degree_map =", row_degree_map)
    if debug >= 4:
        print("\n-- Update Debugging Step 2 : input data  and later  updates --")
        print()
        print("col_degree_map =", col_degree_map) 
        print("row_degree_map =", row_degree_map) 
    # Reduce the count of rows/cols degrees based on removed edges 1 
    for row, col in edges_1_del :  # Loop only through removed edges
        if row in row_degree_map:  # Ensure the row exists in the dictionary
            row_degree_map[row] = max(0, row_degree_map[row] - 1)
        else:
            print(f"!!!!!!Warning: Row {row} not found in row_degree_map")  # Debugging line
        #row_degree_map[row] = max(0, row_degree_map[row] - 1)  # Reduce degree only if that row was involved
        if col in col_degree_map:  # Ensure the column exists in the dictionary
            col_degree_map[col] = max(0, col_degree_map[col] - 1)
        else:
            print(f"!!!!!!Warning: col {col} not found in col_degree_map")  # 
        #col_degree_map[col] = max(0, col_degree_map[col] - 1)  # Reduce degree only if that column was involved

    #print("\n-- Debugging Step: Updated ROW/Column Degrees --")
    # print("Updated row_degree_map =", row_degree_map)
    # print("Updated col_degree_map =", col_degree_map)

    # Rebuild cols_data with updated degrees
    cols_rem = [(col, degree) for col, degree in col_degree_map.items()] 
    cols_rem = [item for item in cols_rem if item[1] != 0] 
    #cols_last = [tup for tup in cols_rem if tup[0] in cols_res_set] # ??? should be  cols_last = cols_res 
    rows_rem = [(row, degree) for row, degree in row_degree_map.items()]
    rows_rem = [item for item in rows_rem if item[1] != 0] 
    col_rem_length = len(rows_rem)
    row_rem_length = len(cols_rem)
    #print("row_rem_length=",row_rem_length, "col_rem_length=",col_rem_length  )
    # print("rows_rem =",rows_rem, "row_rem_length=",row_rem_length )
    # print("cols_rem =",cols_rem, "col_rem_length=",col_rem_length )
    #nb_edges_1_rem = sum(degree for _, degree in rows_rem)
    size = row_rem_length * col_rem_length
    #nb_edges_0_rem = size - nb_edges_1_rem
    nb_edges_0_rem, nb_edges_1_rem, sparsity, density = density_calcul(rows_rem, cols_rem)
    if size != 0:
        #density = nb_edges_1_rem/size
        if debug >= 3:
            print("\n-- Update Debugging Step 3 : input data  and later  updates --")
            print()
            print("col_degree_map reduced =", col_degree_map) 
            print("row_degree_map reduced =", row_degree_map) 
            print("cols_rem   =", cols_rem ) 
            print("rows_rem   =", rows_rem ) 
            print(" nb_edges_1_rem   =", nb_edges_1_rem ) 
            print(" nb_edges_0_rem   =", nb_edges_0_rem )
            print("row_rem_length=", row_rem_length,"col_rem_length =", col_rem_length)
            print("Stats in updata_data : row_rem_length =",row_rem_length,"col_rem_length =",col_rem_length, "nb_edges_0_rem=", nb_edges_0_rem, "nb_edges_1_rem=", nb_edges_1_rem, " !!!!! density =", density)
            print()
    if size == 0:
        if debug >= 4:
            print("\n-- Update Debugging Step 3 : input data  and later  updates --")
            print()
            print("col_degree_map reduced =", col_degree_map) 
            print("row_degree_map reduced =", row_degree_map) 
            print("cols_rem   =", cols_rem ) 
            print("rows_rem   =", rows_rem ) 
            print(" nb_edges_1_rem   =", nb_edges_1_rem ) 
            print(" nb_edges_0_rem   =", nb_edges_0_rem )
            print("row_rem_length=", row_rem_length,"col_rem_length =", col_rem_length)
            print("Stats in updata_data : row_rem_length =",row_rem_length,"col_rem_length =",col_rem_length, "nb_edges_0_rem=", nb_edges_0_rem, "nb_edges_1_rem=", nb_edges_1_rem, "  density  non determined,  !!! size =", size)
        print()
        sys.exit("Terminating program due to matrix degeneration EXIT 200.")

    return rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density 
 
def update_data_first(rows_data, cols_data, edges_0, edges_1, rows_del, cols_del, rows_res, cols_res):
    # Convert rows_del to a set for faster lookup
    # print("\n-- Update Debugging Step: edges  --")
    # print("input edges_0 =", edges_0)
    # print("input edges_1 =", edges_1)
    rows_del_set = set(map(int, rows_del))
    cols_del_set = set(map(int, cols_del))
    # print("\n-- Update Debugging Step: input rows_del_set --")
    # print("rows_del =", rows_del)
    # print("rows_del_set =", rows_del_set)
    # print("rows_res =", rows_res)
    # print("\n-- Update Debugging Step: input cols_del_set --")
    # print("cols_del =", cols_del)
    # print("cols_del_set =", cols_del_set)
    # print("cols_res =", cols_res)

    # Step 1: Compute edges_row (edges connected to rows in rows_del)
    edges_0_row_del= [edge for edge in edges_0 if edge[0] in rows_del_set]
    edges_1_row_del= [edge for edge in edges_1 if edge[0] in rows_del_set]
    edges_0_col_del= [edge for edge in edges_0 if edge[1] in cols_del_set]
    edges_1_col_del= [edge for edge in edges_1 if edge[1] in cols_del_set]

    # print("\n-- Debugging Step: Edges to be removed NEW NEW  --")
    # print("edges_0_row_del =", edges_0_row_del)
    # print("edges_1_row_del =", edges_1_row_del)
    # print("edges_0_col_del =", edges_0_col_del)
    # print("edges_1_col_del =", edges_1_col_del)
    

    # Step 2: Compute edges_new (remaining edges)
    edges_0_row_rem = [edge for edge in edges_0 if edge[0] not in rows_del_set]
    edges_1_row_rem = [edge for edge in edges_1 if edge[0] not in rows_del_set]
    edges_0_col_rem = [edge for edge in edges_0 if edge[0] not in cols_del_set]
    edges_1_col_rem = [edge for edge in edges_1 if edge[0] not in cols_del_set]

    # print("\n-- Debugging Step: Remaining Edges --")
    # print("Remainig Edges_0_row =",  edges_0_row_rem)
    # print("Remainig Edges_1_row =",  edges_1_row_rem)
    # print("Remainig Edges_0_col =",  edges_0_col_rem)
    # print("Remainig Edges_1_col =",  edges_1_col_rem) 

 # Step 3: Update column degrees in cols_data
    col_degree_map = {col: degree for col, degree in cols_data}
    row_degree_map = {row: degree for row, degree in rows_data}
    # print("input cols_data =", cols_data)
    # print("input col_degree_map =", col_degree_map)
    # print("input rows_data =", rows_data)
    # print("input row_degree_map =", row_degree_map)


    # Reduce the count of columns based on removed edges (only affected columns)
    for row, _ in edges_1_row_del :  # Loop only through removed edges
        row_degree_map[row] = max(0, row_degree_map[row] - 1)  # Reduce degree only if that row was involved
    for _, col in edges_1_row_del :  # Loop only through removed edges
        col_degree_map[col] = max(0, col_degree_map[col] - 1)  # Reduce degree only if that column was involved

    # print("\n-- Debugging Step: Updated ROW/Column Degrees --")
    # print("Updated row_degree_map =", row_degree_map)
    # print("Updated col_degree_map =", col_degree_map)

    # Rebuild cols_data with updated degrees
    cols_rem = [(col, degree) for col, degree in col_degree_map.items()]
    cols_res_int = set(map(int, cols_res))
    cols_last = [tup for tup in cols_rem if tup[0] in cols_res_int]
    rows_rem = [(row, degree) for row, degree in row_degree_map.items()]
    rows_res_int = set(map(int, rows_res))
    rows_last = [tup for tup in rows_rem if tup[0] in rows_res_int]
    edges_0_row_last = [edge for edge in edges_0_row_rem if edge[0] in rows_rem]
    edges_1_row_last = [edge for edge in edges_1_row_rem if edge[0] in rows_rem]
    edges_0_col_last = [edge for edge in edges_0_col_rem  if edge[0] in cols_rem]
    edges_1_col_last = [edge for edge in edges_1_col_rem  if edge[0]  in cols_rem]
    # print("cols_data_new =", cols_rem)
    # print("rows_data_new =", rows_rem)

    return rows_last, cols_last, edges_0_row_rem, edges_1_row_rem  

# # Given Input Data
# rows_data = [(0, 2), (1, 2), (2, 2)]
# cols_data = [(0, 2), (1, 3), (2, 1)]
# edges = [(0, 2), (1, 2), (2, 0)]
# rows_del = ['1', '2']

# # Run Function
# cols_data_new, edges_new = update_data(rows_data, cols_data, edges, rows_del)

# # Print Final Output
# print("Updated Columns Data:", cols_data_new)
# print("Updated Edges:", edges_new)


###############################################################################

def get_data(path: str):
    rows_data = []
    cols_data = []
    edges_0 = []
    edges_1 = []

    df = pd.read_csv(path, header=0, index_col=0)
    df[df == -1] = 0  # Ensure all values are 0 or 1

    rows = df.sum(axis=1)
    row_names = rows.index
    rows_data = list(zip(range(len(row_names)), rows))

    cols = df.sum(axis=0)
    col_names = cols.index
    cols_data = list(zip(range(len(col_names)), cols))

    df = df.reset_index(drop=True)
    df = df.T.reset_index(drop=True).T
    edges_1 = list(df[df == 1].stack().index)
    edges_0 = list(df[df == 0].stack().index)

    # Compute the complement matrix
    comp_df = 1 - df  # Flipping 0s to 1s and 1s to 0s

    return rows_data, cols_data, edges_0, edges_1, row_names, col_names, df, comp_df

# ###################################################################################

def get_complement_edges(num_row, num_col, edges):  # Working version
    # Create a set of all possible edges using Cartesian product
    all_edges = set(itertools.product(range(num_row), range(num_col)))

    # Convert the original edges list to a set for efficient subtraction
    original_edges = set(edges)

    # Complement edges = all_edges - original_edges
    complement_edges = list(all_edges - original_edges)

    return complement_edges

####################################################################################
def get_complement_rowcols(rows, cols, edges):  # Working version
    cols_compl_map = {col: len(rows)-degree for col, degree in cols}
    rows_compl_map = {row: len(cols) -degree for row, degree in rows}
    cols_compl = [(col,degree) for col, degree in cols_compl_map.items()]
    rows_compl = [(row, degree) for row, degree in rows_compl_map.items()]
    rows_compl = [(row, degree) for row, degree in rows_compl_map.items()]
    edges_compl = [(r, c) for r,_ in rows_compl for c,_ in cols_compl if (r, c) not in edges]
    return rows_compl, cols_compl,  edges_compl 

################################################################

def get_data_txt_file(path):  # Working version
    with open(path, 'r') as file:
        content = file.readlines()
    
    name = content[0][2:-1].strip()
    num_row = int(content[1].split(":")[1].strip())
    num_col = int(content[2].split(":")[1].strip())
    num_edge = int(content[3].split(":")[1].strip())

    deg_row = [0] * num_row
    deg_col = [0] * num_col
    edges = []
    df = pd.DataFrame(0, index=range(num_row), columns=range(num_col))

    for line in content[4:]:
        splitted_line = line.strip().split()
        u, v = int(splitted_line[0]), int(splitted_line[1])
        
        edges.append((u, v))
        deg_row[u] += 1
        deg_col[v] += 1
        #df.iloc[u, v] = 1 

    rows_data = list(zip(range(num_row), deg_row))
    cols_data = list(zip(range(num_col), deg_col))

    # Compute the complement matrix
    #comp_df = 1 - df  # Flip 0s to 1s and 1s to 0s

    return rows_data, cols_data, edges, range(num_row), range(num_col) #, df, comp_df


def print_log_output(prob, in_run_time, obj, nb_rows, nb_cols):
    """Print the log output and problem  solutions and incerease  the run time by the time taken to solve the problem.
    ARGUMENTS:
    ----------
    * prob: an solved LP model (pulp.LpProblem)
    """
    out_run_time = in_run_time + prob.solutionCpuTime
    if debug >= 2: 
        print()
        print('-' * 70)
        print('Stats')
    #print('-' * 70) 
    #print(f'Number variables: {prob.numVariables()}')
    #print(f'Number constraints: {prob.numConstraints()}')
        print(f' Model solved :', prob.name, f' Local Time:  - (real) {prob.solutionTime:.5f}', f'- (Local CPU) {prob.solutionCpuTime:.5f}', f'- (Global CPU) {out_run_time:.5f}')
    #print(f'Number variables: {prob.numVariables()}')) 
    #print(f'- (real) {prob.solutionTime}')
    #print(f'- (CPU) {prob.solutionCpuTime}')
        print()
        print(f' Solve status: {LpStatus[prob.status]}', f'Objective value: {prob.objective.value():.5f}', "nb_rows", nb_rows, "nb_cols", nb_cols )
        print('-' * 40)

    # print()
    # print('-' * 40)
    # print("Variables' values")
    # print('-' * 40)
    # print()
    # for v in prob.variables():
    #       if v.varValue == 1 : 
    #         print(v.name, v.varValue)
    # #breakpoint()
    # print('-' * 40)
    return out_run_time 

def parse_arguments():
    """Parse the input arguments and retrieve the chosen resolution method and
    the instance that must be solved."""
    global debug, timelimit  # Declare as global variables

    argparser = ArgumentParser()

    argparser.add_argument(
        '--filepath', dest='filepath', required=True, default='',
        help='Select the data',
    )

    argparser.add_argument(
        '--model', dest='model', required=False, default='max_e_r',
        help='Select the model to use',
    )

    argparser.add_argument(
        '--rho', dest='rho', required=False, default=1.0, type=float,
        help='Select the zero deletion rho value',
    )

    argparser.add_argument(
        '--delta', dest='delta', required=False, default=0.1, type=float,
        help='Select the error rho value (tolerance)',
    )

    argparser.add_argument(
        '--threshold', dest='threshold', required=False, default=0.87, type=float,
        help='Above this value greedy approaches are used; below, exact methods',
    )

    argparser.add_argument(
        '--debug', nargs='?', required=False, const=1, type=int, default=0, 
        help='Select the debug value',
    )

    argparser.add_argument(
        '--dec_conq', nargs='?', required=False, const=1, type=int, default=0, 
        help='Select the dec_conq value',
    )

    argparser.add_argument(
        '--timelimit', nargs='?', required=False, const=1, type=int, default=1800, 
        help='Select the time limit for all models except knapsacks',
    )

    # Parse arguments
    arg = argparser.parse_args()

    # Assign to global variables
    debug = arg.debug
    timelimit = arg.timelimit

    if arg.model not in ['max_e', 'max_e_h', 'max_e_r','max_e_c','max_v', 'KP_QBr', 'KP_QBc']:
        argparser.print_help()
        sys.exit(1)

    return (arg.filepath, arg.model, arg.rho, arg.delta,  arg.threshold, arg.dec_conq)
    
def get_submatrices(rows_data, cols_data, edges, rows_res, cols_res):       
    """
    Computes submatrices M1 and M2 from the original matrix M.

    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    rows_res: list of row indices or list of tuples (row, degree).
    cols_res: list of column indices or list of tuples (col, degree).

    Returns:
    --------
    (rows_M1, cols_M1, edges_M1), (rows_M2, cols_M2, edges_M2)
    """

    # Extract only indices if tuples are given
    rows_res_set = {row if isinstance(row, int) else row[0] for row in rows_res}
    cols_res_set = {col if isinstance(col, int) else col[0] for col in cols_res}

    # Compute M1 (Excludes selected rows)
    rows_M1 = [(row, deg) for row, deg in rows_data if row not in rows_res_set]
    edges_M1 = [(r, c) for r, c in edges if r not in rows_res_set]
    cols_M1 = [(col, sum(1 for r, c in edges_M1 if c == col)) for col, _ in cols_data]

    # Compute M2 (Excludes selected columns)
    cols_M2 = [(col, deg) for col, deg in cols_data if col not in cols_res_set]
    edges_M2 = [(r, c) for r, c in edges if c not in cols_res_set]
    rows_M2 = [(row, sum(1 for r, c in edges_M2 if r == row)) for row, _ in rows_data]

    return (rows_M1, cols_M1, edges_M1), (rows_M2, cols_M2, edges_M2)


def get_submatrices_old(rows_data, cols_data, edges, rows_res, cols_res):
    # suppose that rows_data and cols_data are the rows and columns of the original matrix M. Let rows_res and cols_res are the rows and columns of a found submatrix M0.  We want to compute the submatrices M1 and M2 such that M1=[rows_M - rows_res, cols_M] and M2=[rows_M, cols_M - cols_res]
    # Compute M1 and M2s
    # Convert row and column data to sets for easy removal
    # rows_res_set = set(row for row, _ in rows_res)
    # cols_res_set = set(col for col, _ in cols_res)
    rows_res_set = set(row for row in rows_res)
    cols_res_set = set(col for col in cols_res)

    # Compute M1
    rows_M1 = [(row, deg) for row, deg in rows_data if row not in rows_res_set]
    edges_M1 = [(r, c) for r, c in edges if r not in rows_res_set]
    cols_M1 = [(col, sum(1 for r, c in edges_M1 if c == col)) for col, _ in cols_data]

    # Compute M2
    cols_M2 = [(col, deg) for col, deg in cols_data if col not in cols_res_set]
    edges_M2 = [(r, c) for r, c in edges if c not in cols_res_set]
    rows_M2 = [(row, sum(1 for r, c in edges_M2 if r == row)) for row, _ in rows_data]

    return (rows_M1, cols_M1, edges_M1), (rows_M2, cols_M2, edges_M2)

def density_calcul(rows, cols): 
    """
    Computes the density and sparsity of a matrix given its row and column degrees.
    
    Args:
        rows (list of tuples): Each tuple (row_index, degree) for row degrees.
        cols (list of tuples): Each tuple (col_index, degree) for column degrees.

    Returns:
        tuple: (nbi_0, nbi_1, sparsity, density)
    """
    if not rows or not cols:  # ✅ Avoid division by zero
        return 0, 0, None, None  # None indicates sparsity and density are undefined
    
    # Count ones
    nbi_1 = sum(degree for _, degree in rows)

    # Total number of elements
    total_elements = len(rows) * len(cols)

    # Count zeros
    nbi_0 = total_elements - nbi_1

    # Compute sparsity and density safely
    sparsity = nbi_0 / total_elements if total_elements > 0 else None
    density = nbi_1 / total_elements if total_elements > 0 else None

    return nbi_0, nbi_1, sparsity, density


def affichage(dec_conq,matrix_name, rows_res, cols_res, density, nb_ones, iter, KP_time,  kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g):
    # pretty print the results
    global_time = KP_time + QBC_time_g  #+ QBC_time_h 
    global_time_c = QBC_time_g #+ QBC_time_h 
    if selected_model == "max_e_c":
        time_e_wr = QBC_time_g - QBC_time_h
    else:
         time_e_wr = 0
    if global_time != 0:
        percentage_greedy = (KP_time / global_time) * 100
    else:
        print("Percentage: Undefined (division by zero). Something wrong. Global time is zero.")
        sys.exit("Terminating program due to division by zero")
    if QBC_time_g != 0:
        percentage_c = (QBC_time_h / QBC_time_g) * 100
    else:
        percentage_c = 0
        print(" Global QBC time is zero. Only greedy approach has been used. ")
        #sys.exit("Terminating program due to division by zero")
    # with Size of initial matrix : {len(rows)} * {len(cols)}
    print()
    print('-' * 70)
    print(f""" 
    End of computations for matrix {matrix_name} in  {file_path} at level {dec_conq} and debug  {debug}
    With  model: {selected_model} and quasi-biclique error: {delta} 
    Size of Remaining matrix : ({ len(rows_res)},{len(cols_res)}), with  density : {density} and number of ones: {nb_ones}
    Global Time (in sec): {global_time:.3f}
    Total Time in QBC approaches: {global_time_c:.3f}
    Time in greedy approach: {KP_time:.3f},  size of matrix found by greedy : ({nb_kp_rows},{nb_kp_cols}) 
    With density : {kp_density:.3f} and # ones : {nb_kp_ones} 
    Time in heuristic only : {QBC_time_h:.3f} 
    Time in warm start only : {time_e_wr :.3f}
    Number of iterations in greedy approach: {iter}
    Percentage of greedy approach run time in global run time : {percentage_greedy:.2f}% 
    Percentage of heuristic run time in QBC run time : {percentage_c:.2f}%
    With zero deletion rate (rho): {rho} and threshold: {threshold}
    """)
    print('-' * 70)
    print()
    if debug >=1:
        # rows_last = [str(c) for c, _ in rows_res]
        # rows_last = [int(r) for r in rows_last]
        # cols_last = [str(c) for c, _ in cols_res]
        # cols_last = [int(c) for c in cols_last]             
        # print(" ***Solution Rows indices :", rows_last )
        # print(" ***Solution Cols indices :", cols_last )
        if debug >= 3:
                print(" Remaining Rows with degree :", rows_res )
                print(" Remaining Cols with degree :", cols_res )
    return matrix_name, rows_res, cols_res, density, nb_ones, global_time

def final_print(dec_conq, rows, cols, edges, model,solved_count, fathomed_count,nb_skipped, rho, delta):
    nbi_0, nbi_1, sparsity, density = density_calcul(rows, cols)
    print('-' * 70)
    print('-' * 70)
    nb_ext = 2**dec_conq
    nb_int = nb_ext - 1 
    print(f""" 
    End of computations for matrix  {file_path} with rows: {len(rows)} and columns {len(cols)}
    with  input density : {density:.3f} and number of ones: {len(edges)}
    time limit has been set to {timelimit}
    using  model: {selected_model}  with quasi-biclique error: {delta} 
    nd density_threshold: {threshold:.3f} and zero deletion rate (rho): {rho} and debug: {debug}
    Decrease and conquer levels:  {dec_conq}, # ext task: {nb_ext}, int task : {nb_int} 
    The solution has been found in matrix : {best_task}  with 
    size max clique  {best_obj}, # rows: {len(best_rows)} # columns: {len(best_cols)},
    # solved  tasks : {solved_count},  # fathomed tasks : {fathomed_count} # skipped tasks : {nb_skipped}
    """)
    print('-' * 70)
    print()
    # print(f"***We solved instance {file_path} at level {dec_conq} and with model:  {model}***")
    # print("Size of current matrix : ", len(rows), "*", len(cols), "=", len(rows) * len(cols), "; dec_conq:", dec_conq)
    # print("number input zeros : ",nbi_0, "; number input ones : ",nbi_1, )            
    # print("rho = ",rho, "; delta : ", delta)
    #     # Safe printing: Convert None to "N/A" or provide a default value
    # print(f"Input density : {density:.3f}" if density is not None else "Input density: N/A",  f"; density_threshold: {threshold:.5f}")

def last_affichage(winning_matrix, matrix_name, rows_res, cols_res, density, nb_ones, global_time):
    # pretty print the results
    print()
    print('-' * 70)
    print(f""" 
    End of computations for matrix {matrix_name} in  {file_path} at level {dec_conq} and debug  {debug}
    The solution has been found in matrix : {winning_matrix} !!!
    Global Time (in sec): {global_time:.3f}
    With  model: {selected_model} and quasi-biclique error: {delta} 
    Size of Remaining matrix : ({ len(rows_res)},{len(cols_res)}), with  density : {density} and number of ones: {nb_ones}
    With zero deletion rate (rho): {rho} and threshold: {threshold}
    """)
    print('-' * 70)
    print()
    if debug >=1:
        rows_last = [str(c) for c, _ in rows_res]
        rows_last = [int(r) for r in rows_last]
        cols_last = [str(c) for c, _ in cols_res]
        cols_last = [int(c) for c in cols_last]             
        print(" ***Solution Rows indices :", rows_last )
        print(" ***Solution Cols indices :", cols_last )
        if debug >= 3:
                print(" Remaining Rows with degree :", rows_res )
                print(" Remaining Cols with degree :", cols_res )


def get_complement_edges(rows, cols, edges):
    # Ensure rows are in (index, degree) format
    #if isinstance(rows[0], tuple):
    row_indices = {row for row, _ in rows}
    #row_indices = set(map(int, rows))
    # #else:
    #     row_indices = set(rows)  # Handle case where rows are just indices

    #if isinstance(cols[0], tuple):
    col_indices = {col for col, _ in cols}
    # else:
    #     col_indices = set(cols)

    # Compute complement edges
    edges_compl = [(r, c) for r in row_indices for c in col_indices if (r, c) not in edges]
    
    return edges_compl 


# Function to add tasks to the priority queue
def add_task(matrix_name, rows, cols, edges, obj):
    edge_count = len(edges)  # Use number of edges instead of size
    heapq.heappush(QUEUE, (-edge_count, edge_count, (matrix_name, rows, cols, edges, obj)))  # Negative for max-heap
# Function to add tasks to the priority queue

def process_tasks(selected_model, global_time):
    global QUEUE, EVALUATED_QUEUE

    best_matrix = None  # Track the best task number
    best_obj = float('-inf')  # Track the highest objective value
    #best_obj = given_lower_bound  # Track the highest objective value
    #prev_lower_bound  = float('-inf') 
    #prev_lower_bound  = given_lower_bound 
    solved_count = 0  # Count of tasks that were solved
    skipped_count = 0  # Count of tasks that were skipped
    fathomed_count = 0 # Count for fathomed tasks
    dec_conq = 0
    KP_time = 0
    all_checked = False
    while QUEUE: 
        _, edge_count, (matrix_name, rows, cols, edges, obj_val) = heapq.heappop(QUEUE)  # Extract highest priority task
        if best_obj >= edge_count :
            print(f"Task {matrix_name} with edges count {edge_count}) has been skipped by the best task  {best_matrix} with obj  : {best_obj}.")
            print(f"All other tasks are also skipped because the queue is sorted")
            return best_matrix, best_obj, best_rows, best_cols, best_density, EVALUATED_QUEUE, fathomed_count, solved_count
        # Solve the problem
        start_solving_task_count = time.time()
        #prev_lower_bound = obj_val 
        print() 
        print(f"***QUEUE We currently process task number {matrix_name} with (edges {len(edges)}) selected_model {selected_model} dec_conq {dec_conq} delta {delta} threshold {threshold} rho {rho} QBC_time {QBC_time} ***")
        print() 
        results = solve(best_obj,dec_conq, matrix_name, rows, cols, edges, selected_model, KP_time, QBC_time, rho, delta, threshold)   
        # Unpack results
        (rows_res, cols_res, density, nb_ones, iter, KP_time, 
        kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, 
        QBC_time_h, QBC_time_g) = results
        # Display results
        view = affichage(dec_conq, matrix_name, rows_res, cols_res, density, nb_ones, iter, KP_time, kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g)     
        (matrix_name, rows_res, cols_res, density, nb_ones, QBC_time_g) = view 
        # Compute objective function
        obj = len(rows_res) * len(cols_res)  # Replace with actual computation
        end_solving_task_count = time.time()
        solved_count += 1  # Increment solved task count
        print(f"I solved TASK NUMBER {matrix_name} with (edges: {edge_count}) and obj: {obj}  with solving TIME : {end_solving_task_count - start_solving_task_count:.4f} sec" )
        # Count fathomed tasks
        if best_obj >= obj:
            print(f"Task {matrix_name} with obj  {obj} has been fathomed by the best task {best_matrix} with obj  : {best_obj}.")
            fathomed_count  += 1 
            continue 
        if obj > best_obj: 
            best_obj = obj
            best_matrix = matrix_name
            best_rows = rows_res
            best_cols = cols_res
            best_density = density
            print(f"Task {best_matrix} with obj {best_obj}) is the current record.")
        # Store the evaluated task
        EVALUATED_QUEUE.append((matrix_name, rows, cols, edge_count, obj, len(rows_res), len(cols_res)))
    # Return the best task, best objective, evaluated queue, fathomed count, solved count, and skipped count
    return best_matrix, best_obj, best_rows, best_cols, best_density, EVALUATED_QUEUE, fathomed_count, solved_count 

def process_tasks_PREV(selected_model, QBC_time):
    global QUEUE, PROCESSED_OBJS, EVALUATED_QUEUE

    best_task = None  # Track the best task number
    best_obj = float('-inf')  # Track the highest objective value
    prev_lower_bound  = None 
    solved_count = 0  # Count of tasks that were solved
    skipped_count = 0  # Count of tasks that were skipped

    while QUEUE:
        _, edge_count, (matrix_name, rows, cols, edges, obj_val) = heapq.heappop(QUEUE)  # Extract highest priority task

        # Check fathoming condition: If a processed obj >= edge_count of this task, skip it
        if any(obj >= edge_count for obj in PROCESSED_OBJS):
            print(f"Skipping {matrix_name} (edges {edge_count}) - Fathomed by an earlier task.")
            skipped_count += 1  # Increment skipped task count
            continue

        # Solve the problem
        start_solving_task_count = time.time()
        prev_lower_bound = obj_val 
        KP_time = 0
        dec_conq = 0
        print() 
        print(f"***QUEUE Processing of task number {matrix_name} with (edges {edge_count}) selected_model {selected_model} dec_conq {dec_conq} delta {delta} threshold {threshold} rho {rho} QBC_time {QBC_time} ***")
        print() 

        results = solve(prev_lower_bound,dec_conq, matrix_name, rows, cols, edges, selected_model, KP_time, QBC_time, rho, delta, threshold)
        
        # Unpack results
        (rows_res, cols_res, density, nb_ones, iter, KP_time, 
        kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, 
        QBC_time_h, QBC_time_g) = results
        
        # Display results
        view = affichage(dec_conq, matrix_name, rows_res, cols_res, density, nb_ones, iter, KP_time,  
                         kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g)
        
        (matrix_name, rows_res, cols_res, density, nb_ones, QBC_time_g) = view 

        # Compute objective function
        obj = len(rows_res) * len(cols_res)  # Replace with actual computation
        PROCESSED_OBJS.append(obj)  # Store obj value
        end_solving_task_count = time.time()
        solved_count += 1  # Increment solved task count

        print(f"PROCESSED TASK NUMBER {matrix_name} with (edges {edge_count}) -> obj = {obj}  with solving TIME : {end_solving_task_count - start_solving_task_count:.4f} sec" )

        # Store the evaluated task
        EVALUATED_QUEUE.append((matrix_name, rows, cols, edge_count, obj, len(rows_res), len(cols_res)))

        # Update best task if this one is better
        if obj > best_obj:
            best_obj = obj
            best_task = matrix_name
            best_rows = rows_res
            best_cols = cols_res
            best_density = density

    # Count fathomed tasks
    fathomed_count = sum(1 for matrix_name, rows, cols, edge_count, obj, nb_rows, nb_cols in EVALUATED_QUEUE if obj < best_obj)

    # Return the best task, best objective, evaluated queue, fathomed count, solved count, and skipped count
    return best_task, best_obj, best_rows, best_cols, best_density, EVALUATED_QUEUE, fathomed_count, solved_count, skipped_count

def process_tasks_OLD(selected_model, QBC_time):
    global QUEUE, PROCESSED_OBJS, EVALUATED_QUEUE

    best_task = None  # Track the best task number
    best_obj = float('-inf')  # Track the highest objective value

    while QUEUE:
        _, edge_count, (matrix_name, rows, cols, edges, obj_val) = heapq.heappop(QUEUE)  # Extract highest priority task

        # Check fathoming condition: If a processed obj >= edge_count of this task, skip it
        if any(obj >= edge_count for obj in PROCESSED_OBJS):
            print(f"SKIPPING task number {matrix_name} (edges {edge_count}) - Fathomed by an earlier task.")
            continue

        # Solve the problem
        KP_time = 0
        dec_conq = 0
        print() 
        print(f"***QUEUE Processing of task number {matrix_name} with (edges {edge_count}) selected_model {selected_model} dec_conq {dec_conq} delta {delta} threshold {threshold} rho {rho} QBC_time {QBC_time} ***")
        print() 

        results = solve(dec_conq, matrix_name, rows, cols, edges, selected_model, KP_time, QBC_time, rho, delta, threshold)
        
        # Unpack results
        (rows_res, cols_res, density, nb_ones, iter, KP_time, 
        kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, 
        QBC_time_h, QBC_time_g) = results
        
        # Display results
        view = affichage(dec_conq, matrix_name, rows_res, cols_res, density, nb_ones, iter, KP_time,  
                         kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g)
        
        (matrix_name, rows_res, cols_res, density, nb_ones, QBC_time_g) = view 

        # Compute objective function
        obj = len(rows_res) * len(cols_res)  # Replace with actual computation
        PROCESSED_OBJS.append(obj)  # Store obj value

        print(f"PROCESSED TASK NUMBER {matrix_name} with (edges {edge_count}) -> obj = {obj}")

        # Store the evaluated task
        EVALUATED_QUEUE.append((matrix_name, rows, cols, edge_count, obj, len(rows_res), len(cols_res)))

        # Update best task if this one is better
        if obj > best_obj:
            best_obj = obj
            best_task = matrix_name
            best_rows = rows_res
            best_cols = cols_res
            best_density = density

    # Count fathomed tasks
    fathomed_count = sum(1 for matrix_name, rows, cols, edge_count, obj, nb_rows, nb_cols in EVALUATED_QUEUE if obj < best_obj)

    # Return the best task, best objective, evaluated queue, and fathomed count
    return best_task, best_obj, best_rows, best_cols, best_density, EVALUATED_QUEUE, fathomed_count

# Function to add tasks to the priority queue
def add_task_OLD(matrix_name, rows, cols, edges):
    size = len(rows) * len(cols)  # Compute the size
    heapq.heappush(QUEUE, (-size, size, (matrix_name, rows, cols, edges)))  # Negative for max-heap


def process_tasks_OLD(selected_model, QBC_time):
    global QUEUE, PROCESSED_OBJS

    best_task = None  # Track the best task number
    best_obj = float('-inf')  # Track the highest objective value

    while QUEUE:
        _, size, (matrix_name, rows, cols, edges) = heapq.heappop(QUEUE)  # Extract highest priority task

        # Check fathoming condition: If a processed obj >= size of this task, skip it
        if any(obj >= size for obj in PROCESSED_OBJS):
            print(f"Skipping {matrix_name} (size {size}) - Fathomed by an earlier task.")
            continue

        # Solve the problem
        KP_time = 0
        dec_conq = 0
        print() 
        print(f"***QUEUE Processing of task number {matrix_name} with (size {size}) selected_model {selected_model} dec_conq {dec_conq} delta {delta} threshold {threshold} rho {rho} QBC_time {QBC_time} ***")
        print() 

        results = solve(dec_conq, matrix_name, rows, cols, edges, selected_model, KP_time, QBC_time, rho, delta, threshold)
        
        # Unpack results
        (rows_res, cols_res, density, nb_ones, iter, KP_time, 
        kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, 
        QBC_time_h, QBC_time_g) = results
        
        # Display results
        view = affichage(dec_conq, matrix_name, rows_res, cols_res, density, nb_ones, iter, KP_time,  
                         kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g)
        
        (matrix_name, rows_res, cols_res, density, nb_ones, QBC_time_g) = view 

        # Compute objective function
        obj = len(rows_res) * len(cols_res)  # Replace with actual computation
        PROCESSED_OBJS.append(obj)  # Store obj value

        print(f"PROCESSED TASK NUMBER {matrix_name} with (size {size}) -> obj = {obj}")

        # Update best task if this one is better
        if obj > best_obj:
            best_obj = obj
            best_task = matrix_name

    # Return the best task and its objective value
    return best_task, best_obj if best_task is not None else (None, None)


def process_tasks_OLD(selected_model, QBC_time):
    global QUEUE
    
    while QUEUE:
        _, size, (matrix_name, rows, cols, edges) = heapq.heappop(QUEUE)  # Extract highest priority task
        
        # Check fathoming condition: If a processed obj >= size of this task, skip it
        if any(obj >= size for obj in PROCESSED_OBJS):
            print(f"Skipping {matrix_name} (size {size}) - Fathomed by an earlier task.")
            continue
        # Solve the problem
        KP_time = 0
        dec_conq = 0
        print() 
        print(f"***QUEUE Processing of task number {matrix_name} with (size {size}) selected_model {selected_model} dec_conq {dec_conq} delta {delta} threshold {threshold} rho {rho} QBC_time {QBC_time} ***")
        print() 
        results = solve(dec_conq, matrix_name, rows, cols, edges, selected_model, KP_time, QBC_time, rho, delta, threshold)
        # Unpack results
        (rows_res, cols_res, density, nb_ones, iter, KP_time, 
        kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, 
        QBC_time_h, QBC_time_g) = results
        # Display results
        view = affichage(dec_conq, matrix_name, rows_res, cols_res, density, nb_ones, iter, KP_time,  kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g)
        (matrix_name, rows_res, cols_res, density, nb_ones, QBC_time_g) = view 
        # Process the task (dummy computation)
        obj = len(rows_res) * len(cols_res)  # Replace with your actual computation
        #obj = compute_obj(matrix_name, rows, cols, edges)  # Replace with your actual computation
        PROCESSED_OBJS.append(obj)  # Store obj value
        
        print(f"PROCESSED TASK NUMBER  {matrix_name} with  (size {size}) -> obj = {obj}")


def decrease_and_conquer(dec_conq, matrix_name, rows, cols, edges_1, KP_time, QBC_time):
    """
    Implements a decrease-and-conquer approach to solve the problem.

    Args:
        matrix_name (int): number of the corresponding matrix. 
        dec_conq : level of the decrease and conquer approach. If dec_conq = 0, the decrease and conquer approach is not applied and the original matrix is used directly for computatiion. When dec_conq >= 1, the decrease and conquer approach is utilized. As a result, the matrix is reduced and divided into two smaller submatrices that certainly contain the maximun size clique. The process is repeated recursively until dec_conq = 0.
        rows (list of tuples): List of (row_index, degree) for rows.
        cols (list of tuples): List of (col_index, degree) for columns.
        edges_1 (list of tuples): List of existing edges (row_index, col_index).
        KP_time (float): Time taken for the greedy approach.
        QBC_time (float): Time taken for the quasi-biclique approach.

    Returns:
        M1, M2: Two submatrices.
        KP_time: Updated KP_time.
        QBC_time_g: Updated QBC_time_g.
    """
    # Compute complementary row and column indices
    if dec_conq == 0: # No decrease-and-conquer
        temp_obj =  float('-inf') 
        add_task(matrix_name, rows, cols, edges_1, temp_obj)  # Add the task to the priority queue
        # Compute the density and  number of ones in the matrix
        nb_zeros, nb_ones, sparsity, density = density_calcul(rows, cols)
        if debug >= 2:
            print() 
            print(f"Task with matrix {matrix_name} with size ({len(rows)},{len(cols)}) and density {density} and number of ones {nb_ones}  and number of zeros {nb_zeros} has been added to the queue.")
        return  matrix_name, rows, cols, density, nb_ones, QBC_time
    if dec_conq >= 1:
        # Compute complementary row and column indices
        rows_compl, cols_compl, edges_compl = get_complement_rowcols(rows, cols, edges_1)
        # Solve the problem
        results = solve(None,dec_conq, matrix_name, rows_compl, cols_compl, edges_compl, selected_model, KP_time, QBC_time, rho, delta, threshold)
        # Unpack results
        (rows_res, cols_res, density, nb_ones, iter, KP_time, 
        kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, 
        QBC_time_h, QBC_time_g) = results
        # Display results
        view = affichage(dec_conq, matrix_name, rows_res, cols_res, density, nb_ones, iter, KP_time,  kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g)
        ( matrix_name, rows_res, cols_res, density, nb_ones, QBC_time_g) = view 
        # print("rows =", rows)
        # print("cols =", cols)
        # print("edges_1 =", edges_1)
        # print("rows_res =", rows_res)
        # print("cols_res =", cols_res)

        # Compute the density and  number of ones in the matrix
        # Compute submatrices M1 and M2
        print("dec_conq=", dec_conq)
        #sys.exit("Terminating program for cheking . EXIT 108.")
        ML, MR = get_submatrices(rows, cols, edges_1, rows_res, cols_res)
        #ML, MR = get_submatrices(rows_compl, cols_compl, edges_compl, rows_res, cols_res)
    else:
        print("dec_conq=", dec_conq)
        sys.exit("Terminating program due to invalid value for dec_conq. EXIT 1.")
    # Debugging output
    node = 2*matrix_name
    node1= node + 1
    if debug >= 1:
        print(f"\n Level {dec_conq-1}, Matrix {node}:")
        print("Size Rows:", len(ML[0]))
        print("Size Cols:", len(ML[1]))
        print(f"\n Level {dec_conq-1} Matrix  {node1}:")
        print("Size Rows:", len(MR[0]))
        print("Size Cols:", len(MR[1]))
    if debug >= 3:
        print(f"\n Level {dec_conq-1}  Matrix {node} Rows:", ML[0])
        print(f"Level  {dec_conq-1}  Matrix {node} Cols:", ML[1])
        print(f"Level  {dec_conq-1}  Matrix {node} Edges :", ML[2])
        print(f"\n Level  {dec_conq-1}  Matrix {node1} Rows:", MR[0])
        print(f"Level {dec_conq-1}  Matrix {node1} Cols:", MR[1])
        print(f"Level {dec_conq-1}  Matrix {node1} Edges :", MR[2])
    #sys.exit("Terminating program for cheking . EXIT 108.")
    if len(ML[0]) <=  min_number_rows or len(ML[1]) <=  min_number_cols:
        nb_ones_left = 0
        global_time_left = 0
        node_left = None
        rows_ind_left = None 
        cols_ind_left = None
        density_left = None
        if debug >= 1:
            print() 
            print(f" Node {node} has been fathomed because of min number rows = {min_number_rows} or min number columns = {min_number_cols} "  )
            #sys.exit("Terminating program for cheking . EXIT 111.")
            print() 
        #sys.exit("Terminating program for cheking . EXIT 111.")
    else:
        result_left = decrease_and_conquer(dec_conq-1, node, ML[0], ML[1], ML[2], KP_time, QBC_time)
        node_left, rows_ind_left, cols_ind_left, density_left, nb_ones_left, global_time_left  = result_left
        if debug >= 2:
            print()  
            print(f" Return from {node} with winning node = {node_left}" )
    if len(MR[0]) <=  min_number_rows or len(MR[1]) <=  min_number_cols:
        nb_ones_right = 0
        global_time_right = 0
        node_right = None
        rows_ind_right = None
        cols_ind_right = None
        density_right = None
        if debug >= 1:
            print() 
            print(f" Node {node1} has been fathomed because of min number rows = {min_number_rows} or min number columns = {min_number_cols} "  )
            #sys.exit("Terminating program for cheking . EXIT 111.")
            print() 
    else:
        if debug >= 2:
            print() 
            print(f"calling decrease_and_conquer for node {node+1}" )
        result_right = decrease_and_conquer(dec_conq-1, node+1, MR[0], MR[1], MR[2], KP_time, QBC_time)
        node_right, rows_ind_right, cols_ind_right, density_right, nb_ones_right, global_time_right = result_right
        if debug >= 2:
            print() 
            print(f"return in decrease_and_conquer from {node+1} with winning node = {node_right}" )
    if nb_ones_left > nb_ones_right:
        rows_ind = rows_ind_left
        cols_ind = cols_ind_left
        density = density_left
        nb_ones = nb_ones_left
        winning_node = node_left
    if nb_ones_left <= nb_ones_right:
        rows_ind = rows_ind_right
        cols_ind = cols_ind_right
        density = density_right
        nb_ones = nb_ones_right
        winning_node = node_right
    global_time = global_time_left + global_time_right
    if debug >= 2:
            print() 
            print(f"return in decrease_and_conquer from {matrix_name} with winning node = {winning_node}" )

    return winning_node, rows_ind, cols_ind, density, nb_ones, global_time


def decrease_and_conquer_BIS(dec_conq, matrix_name, rows, cols, edges_1, KP_time, QBC_time):
    """
    Implements a decrease-and-conquer approach to solve the problem.

    Args:
        matrix_name (int): Name of the matrix, e.i. nulber of the corresponding node. 
        dec_conq : level of the decrease and conquer approach. If dec_conq = 0, the decrease and conquer approach is not applied and the original matrix is used directly for computatiion. When dec_conq >= 1, the decrease and conquer approach is utilized and the matrix is reduced and divided into two smaller submatrices. One of them contains the maximun size clique. The process is repeated recursively until dec_conq = 0.
        rows (list of tuples): List of (row_index, degree) for rows.
        cols (list of tuples): List of (col_index, degree) for columns.
        edges_1 (list of tuples): List of existing edges (row_index, col_index).
        KP_time (float): Time taken for the greedy approach.
        QBC_time (float): Time taken for the quasi-biclique approach.

    Returns:
        M1, M2: Two submatrices.
        KP_time: Updated KP_time.
        QBC_time_g: Updated QBC_time_g.
    """
    # Compute complementary row and column indices
    if dec_conq == 0: # No decrease-and-conquer
        # add_task(matrix_name, rows, cols, edges_1)  # Add the task to the priority queue
        # nb_ones = len(rows) * len(cols)  # Compute the number of ones in the matrix
        # density = nb_ones / (len(rows) * len(cols))  # Compute the density
        # return  matrix_name, rows, cols, density, nb_ones, QBC_time
        # Solve the problem
        results = solve(dec_conq, matrix_name, rows, cols, edges_1, selected_model, KP_time, QBC_time, rho, delta, threshold)
        # Unpack results
        (rows_res, cols_res, density, nb_ones, iter, KP_time, 
        kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, 
        QBC_time_h, QBC_time_g) = results
        # Display results
        view = affichage(dec_conq, matrix_name, rows_res, cols_res, density, nb_ones, iter, KP_time,  kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g)
        (matrix_name, rows_res, cols_res, density, nb_ones, QBC_time_g) = view 
        return matrix_name, rows_res, cols_res, density, nb_ones, QBC_time_g
        #sys.exit(" Terminating program because dec_conq == 0. End of computations!  EXIT 0")
    if dec_conq >= 1:
        # Compute complementary row and column indices
        rows_compl, cols_compl, edges_compl = get_complement_rowcols(rows, cols, edges_1)
        # Solve the problem
        results = solve(dec_conq, matrix_name, rows_compl, cols_compl, edges_compl, selected_model, KP_time, QBC_time, rho, delta, threshold)
        # Unpack results
        (rows_res, cols_res, density, nb_ones, iter, KP_time, 
        kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, 
        QBC_time_h, QBC_time_g) = results
        # Display results
        view = affichage(dec_conq, matrix_name, rows_res, cols_res, density, nb_ones, iter, KP_time,  kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g)
        ( matrix_name, rows_res, cols_res, density, nb_ones, QBC_time_g) = view 
        # Compute submatrices M1 and M2
        ML, MR = get_submatrices(rows, cols, edges_1, rows_res, cols_res)
    else:
        print("dec_conq=", dec_conq)
        sys.exit("Terminating program due to invalid value for dec_conq. EXIT 1.")
    # Debugging output
    node = 2*matrix_name
    node1= node + 1
    if debug >= 1:
        print(f"\n Level {dec_conq-1}, Matrix {node}:")
        print("Size Rows:", len(ML[0]))
        print("Size Cols:", len(ML[1]))
        print(f"\n Level {dec_conq-1} Matrix  {node1}:")
        print("Size Rows:", len(MR[0]))
        print("Size Cols:", len(MR[1]))
    if debug >= 2:
        print(f"\n Level {dec_conq-1}  Matrix {node} Rows:", ML[0])
        print(f"Level  {dec_conq-1}  Matrix {node} Cols:", ML[1])
        print(f"\n Level  {dec_conq-1}  Matrix {node1} Rows:", MR[0])
        print(f"Level {dec_conq-1}  Matrix {node1} Cols:", MR[1])
    if len(ML[0]) <=  min_number_rows or len(ML[1]) <=  min_number_cols:
        nb_ones_left = 0
        global_time_left = 0
        node_left = None
        rows_ind_left = None 
        cols_ind_left = None
        density_left = None
        print(f" Node {node} has been fathomed because of min number rows = {min_number_rows} or min number columns = {min_number_cols} "  )
    else:
        result_left = decrease_and_conquer(dec_conq-1, node, ML[0], ML[1], ML[2], KP_time, QBC_time)
        node_left, rows_ind_left, cols_ind_left, density_left, nb_ones_left, global_time_left  = result_left 
        print(f" Return from {node} with winning node = {node_left}" )
    if len(MR[0]) <=  min_number_rows or len(MR[1]) <=  min_number_cols:
        nb_ones_right = 0
        global_time_right = 0
        node_right = None
        rows_ind_right = None
        cols_ind_right = None
        density_right = None
        print(f" Node {node1} has been fathomed because of min number rows = { min_number_rows} or min number columns = {min_number_cols} "  )
    else:
        result_right = decrease_and_conquer(dec_conq-1, node+1, MR[0], MR[1], MR[2], KP_time, QBC_time)
        node_right, rows_ind_right, cols_ind_right, density_right, nb_ones_right, global_time_right = result_right
        print(f"return from {node+1} with winning node = {node_right}" )
    if nb_ones_left > nb_ones_right:
        rows_ind = rows_ind_left
        cols_ind = cols_ind_left
        density = density_left
        nb_ones = nb_ones_left
        winning_node = node_left
    if nb_ones_left <= nb_ones_right:
        rows_ind = rows_ind_right
        cols_ind = cols_ind_right
        density = density_right
        nb_ones = nb_ones_right
        winning_node = node_right
    global_time = global_time_left + global_time_right
    print(f"return from {matrix_name} with winning node = {winning_node}" )
    return winning_node, rows_ind, cols_ind, density, nb_ones, global_time

def write_matrix(rows_data, cols_data, edges):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.

    Returns:
    --------
    a csv file with the matrix data.
    """
        
    # Extract row and column indices
    row_indices = [row for row, _ in rows_data]
    col_indices = [col for col, _ in cols_data]

    # Sort indices to ensure correct ordering
    row_indices.sort()
    col_indices.sort()

    # Create a matrix initialized with zeros
    matrix = {row: {col: 0 for col in col_indices} for row in row_indices}

    # Fill the matrix with ones based on edges
    for row, col in edges:
        if row in row_indices and col in col_indices:
            matrix[row][col] = 1

    # Write to CSV
    csv_filename = "matrix.csv"
    with open(csv_filename, mode="w", newline="") as file:
        writer = csv.writer(file)

        # Write header (column labels)
        writer.writerow([""] + [f"c{col}" for col in col_indices])

        # Write rows with data
        for row in row_indices:
            writer.writerow([f"r{row}"] + [matrix[row][col] for col in col_indices])

    print(f"Matrix saved to {csv_filename}")


def count_ones_in_submatrix(edges, row_indices, col_indices):
    """
    Counts the number of ones in the submatrix B of the given matrix M.

    :param matrix: List of tuples representing the edges (row, col) in M
    :param row_indices: Set of row indices defining submatrix B
    :param col_indices: Set of column indices defining submatrix B
    :return: Number of ones in the submatrix B
    """
    count = sum(1 for row, col in edges if row in row_indices and col in col_indices)
    return count

# # Define the original matrix M as a list of (row, col) pairs
# matrix_M = [
#     (0, 0),
#     (0, 1),
#     (1, 0),
#     (1, 1),
#     (2, 1),
#     (2, 2)
# ]

    # # Define the row and column indices for submatrix B
    # row_indices_B = {0, 1}
    # col_indices_B = {0, 1}

    # # Count ones in submatrix B
    # num_ones = count_ones_in_submatrix(matrix_M, row_indices_B, col_indices_B)

    # # Print the result
    # print(f"Number of ones in submatrix B: {num_ones}")



if __name__ == '__main__':
    import time
    min_number_rows = 3
    min_number_cols = 3
    # Define a priority queue (max-heap using negative size)
    QUEUE = []
    COPY_QUEUE = []
    #PROCESSED_OBJS = []  # Store processed obj values
    # List to store evaluated tasks
    EVALUATED_QUEUE = []
    import heapq

    start_model_building_and_solving = time.time()
    # Read the arguments
    file_path, selected_model, rho, delta,  threshold, dec_conq= parse_arguments()
    p = file_path 
    if p.lower().endswith('.txt'):
        rows, cols, edges_1, row_names, col_names  = get_data_txt_file(file_path)
        # for test only to be canceled later
        rows_compl,  cols_compl, edges_compl =  get_complement_rowcols(rows, cols, edges_1) # for test only to be canceled later
        if debug >=4:
            print("rows:\n", rows)
            print("cols:\n", cols)
            print(" row_names:\n",  row_names)
            print(" col_names:\n",  col_names)
            print(" edges_1 :\n", edges_1)
            print(" edges_0 :\n", edges_compl)
            #print("Original Matrix M:\n", M)
    elif p.lower().endswith('.csv'):
        rows, cols, edges_compl, edges_1, row_names, col_names, df, comp_df = get_data(file_path) 
        # for test only to be canceled late
        rows_compl, cols_compl, edges_compl = get_complement_rowcols(rows, cols, edges_1) # for test only to be canceled later
    else:
        raise ValueError('Input need to be a matrix csv file, or a text file with a specific layout')
    ###################################################################
    nb_eges_0 = len(edges_compl)
    nb_eges_1 = len(edges_1)
    if debug >= 2:
            print('-' * 40)
            if p.lower().endswith('.txt'):
                print(f" Input Data in txt files : {file_path}")
            if p.lower().endswith('.csv'):
                print(f" Input Data in csv files : {file_path}")
            print("Number Rows Data :", len(rows) )
            print("Number Cols Data :",  len(cols) )
            print("Number Edges_1 :", len(edges_1))
            print("Number Edges_0 :", len(edges_compl))
            if debug >= 3:
                print(" Rows Data :", rows )
                print(" Cols Data :", cols )
                print("  row_names :",  row_names)
                print("  col_names :",  col_names)
                if debug >= 4:
                    print(" Edges_1 :", edges_1)
                    print(" Edges_0 :", edges_compl)
                    #print("Adjacency Matrix:\n", df)
                    #print("Adjacency Matrix:\n", comp_df)
            print('-' * 40)
    # # end fetching input data
    # start tasks_generation  and computations
    in_rows = rows
    in_cols = cols 
    in_edges_1 = edges_1
    start_tasks_generation = time.time()
    KP_time = 0.0
    QBC_time = 0.0 
    #dec_conq = 2
    node = 1
    winning_node, rows, cols, density, nb_ones, global_time = decrease_and_conquer(dec_conq, node , rows, cols, edges_1, KP_time, QBC_time)
    end_tasks_generation = time.time()
    print()
    print('-' * 70)
    print() 
    print(f"End of tasks generation stage. Last generated task from matrix {node} with winning node {winning_node}")
    print(f"Tasks_generation time: {end_tasks_generation  - start_tasks_generation :.4f} sec")     
    print('-' * 70)
    print(f" Size of the queue: {len(QUEUE)}")
    for _, size, (matrix_name, rows, cols, edges, obj) in QUEUE:
        print(f" Matrix: {matrix_name}, Size: {size}, Number Rows: {len(rows)},  Number Cols: {len(cols)},  Number Edges: {len(edges)}, obj {obj}")
    print()
    print('-' * 70)
    print()
    import heapq
    COPY_QUEUE = list(QUEUE)  # Copy the heap
    heapq.heapify(COPY_QUEUE)  # Ensure it maintains heap properties
    start_tasks_solving = time.time()
    best_task, best_obj, best_rows, best_cols, best_density, EVALUATED_QUEUE, fathomed_count, solved_count = process_tasks(selected_model, global_time)
    end_tasks_solving = time.time()
    nb_skipped = len(COPY_QUEUE)- solved_count
        # Print evaluated queue
    print("\nEvaluated Queue:")
    # for task in EVALUATED_QUEUE:
    #     print(task)
    # print()
    print('-' * 70)
    print(f" Size of the evaluated  queue: {len(EVALUATED_QUEUE)}")
    for matrix_name, rows, cols, edge_count, obj, nb_rows_res, nb_cols_res in EVALUATED_QUEUE:
        #size = len(rows)*len(cols)
        print(f" Matrix: {matrix_name}, # Rows: {len(rows)},  # Cols: {len(cols)},  # Edges: {edge_count}, size max clique  {obj}, # rows: {nb_rows_res} # columns: {nb_cols_res}")
    print()
    print('-' * 70)
    print()
    print('-' * 70)
    sorted_copy_queue = sorted(COPY_QUEUE, key=lambda x: x[1], reverse=True)
    print(f"Size of the COPY_QUEUE: {len(sorted_copy_queue)}")
    for _, size, (matrix_name, task_rows, task_cols, edges, obj) in sorted_copy_queue:
        print(f"Matrix: {matrix_name}, Size: {size}, Number Rows: {len(task_rows)}, Number Cols: {len(task_cols)}, Number Edges: {len(edges)}")
    print()
    print('-' * 70)
    print()
    print(f"Best task: {best_task}, Best objective: {best_obj} with # rows {len(best_rows)} and # cols {len(best_cols)}")
    print(f"Number of solved tasks: {solved_count}, Number of fathomed tasks: {fathomed_count}, Number of skipped tasks: {nb_skipped}")
    print()
    print('-' * 70)
    print("***End of computations !!!")
    print()
    final_print(dec_conq, in_rows, in_cols, in_edges_1,selected_model, solved_count, fathomed_count,nb_skipped, rho, delta)
    print(f"Best task: {best_task}, Best objective: {best_obj} with # rows {len(best_rows)} and # cols {len(best_cols)}")
    print()
    print('-' * 70)
    end_model_building_and_solving = time.time()
    print(f" Model building and solving global time: {end_model_building_and_solving - start_model_building_and_solving:.4f} sec")
    print(f" Tasks_generation time: {end_tasks_generation  - start_tasks_generation :.4f} sec")
    print(f" Tasks_pure_solving time: {end_tasks_solving  - start_tasks_solving:.4f} sec")
    # Your matrix data
    # rows_data = [(0, 7), (1, 5), (2, 10), (3, 9), (4, 4), (5, 5), (6, 4), (7, 4), (8, 9), (9, 7)] 
    # cols_data = [(0, 10), (1, 5), (2, 4), (3, 3), (4, 6), (5, 9), (6, 4), (7, 9), (8, 9), (9, 5)] 
    # edges = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 9), (1, 0), (1, 1), (1, 5), (1, 7), (1, 8), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (4, 0), (4, 5), (4, 7), (4, 8), (5, 0), (5, 4), (5, 7), (5, 8), (5, 9), (6, 0), (6, 5), (6, 7), (6, 8), (7, 0), (7, 5), (7, 7), (7, 8), (8, 0), (8, 1), (8, 2), (8, 4), (8, 5), (8, 6), (8, 7), (8, 8), (8, 9), (9, 0), (9, 4), (9, 5), (9, 6), (9, 7), (9, 8), (9, 9)]
    # rows_data = [(0, 5), (1, 1), (2, 6), (3, 5), (4, 0), (5, 2), (6, 0), (7, 0), (8, 5), (9, 3)] 
    # cols_data = [(1, 5), (2, 4), (3, 3), (4, 6), (6, 4), (9, 5)] 
    # edges = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 9), (1, 1), (2, 1), (2, 2), (2, 3), (2, 4), 
    #         (2, 6), (2, 9), (3, 1), (3, 2), (3, 3), (3, 4), (3, 6), (5, 4), (5, 9), (8, 1), 
    #         (8, 2), (8, 4), (8, 6), (8, 9), (9, 4), (9, 6), (9, 9)]
    #write_matrix(rows_data, cols_data, edges)

    print()
    print(f"size edges_1: {len(in_edges_1)}")
    print()
    #print(f"egdes_1: {in_edges_1}")

    row_indices = [row for row, _ in best_rows]
    col_indices = [col for col, _ in best_cols]

    row_set = set(row_indices)
    col_set = set(col_indices)
    if debug >= 2:
        print(f"row_set: {row_set}")
        print(f"col_set: {col_set}")
    num_ones = count_ones_in_submatrix(in_edges_1,  row_set,  col_set)

    # Print the result
    print(f"Number of ones in the best submatrix : {num_ones}")
    sys.exit("***Computation done EXIT 108 !!!")
    # ################################################
 
