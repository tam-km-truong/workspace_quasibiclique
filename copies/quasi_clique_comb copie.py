# -*- coding=utf-8 -*-

"""
Skeleton for the project: Searching Maximum Quasi-Bicliques

Install Gurobi Solver with free academic licence using your University account at: https://www.gurobi.com/features/academic-named-user-license/

Install Pandas, using PyPI or Conda. A detailed installation tutorial can be find at: https://pandas.pydata.org/docs/getting_started/install.html

an instance 
 python3 quasi_clique_comb.py --filepath data/data4.csv   --model  max_e_h  --rho 0.1  --delta 0.1 --debug 0  --threshold 0.83

"""

import sys
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
        model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
        model += (lpCols[col][0] >= lpCells[(row, col)]), f'ce ll_{row}_{col}_2'

        #########################################
        #compacting  with degree 
        #########################################
      
    # for col, _ in cols_data:
    #     col_edges = [u for u, v in edges if v == col]           
    #     model += (
    #         lpSum(lpCells[(row, col)][0] for row in col_edges) <= lpCols[col][1]*lpCols[col][0]
    #     ), f"col_degre_{col}"
    # for row, _ in rows_data:
    #     row_edges = [v for u, v in edges if u == row]           
    #     model += (
    #         lpSum(lpCells[(row, col)][0] for col in row_edges) <= lpRows[row][1]*lpRows[row][0]
    #     ), f"row_degre_{row}"
        
         #########################################


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
    #constraint for obj improvement
    model += lpSum(lpCells) >= prev_obj + 1, "impovement????"
    #end constraint for obj improvement
    for row, col in edges:
    #for row, col in lpCells:
        model += (lpRows[row][0] >= lpCells[(row, col)]), f'cell_{row}_{col}_1'
        model += (lpCols[col][0] >= lpCells[(row, col)]), f'ce ll_{row}_{col}_2'


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

    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}
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
    # Add row density constraints
    #__row_density(rows_data, cols_data, edges, model, lpRows, lpCols, rho)
    #_col_density(rows_data, cols_data, edges, model, lpRows, lpCols, rho)

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
    Big_M= Big_R+Big_C 
    for row, _ in rows_data:
        #print(f"Adding row density constraints for row {row}:")
        row_edges = [v for u, v in edges if u == row]
        model += (lpSum_row[row] == 
                  lpSum(lpCols[col][0] for col in row_edges), f"row_sum_{row}"
        )
        model += (
            lpSum_row[row] - (1 - delta) * lpSum(lpCols[col][0] for col, _ in cols_data) >= (lpRows[row][0]-1) * Big_M, f"row_slack_{row}" 
           # >= (lpRows[row][0]-1) * Big_M, f"row_err_rate_1_{row}"
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
    #print('rows_data =', rows_data)
        print('-' * 40) 
    
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
    model.solve(GUROBI_CMD(msg=False, timeLimit= 100)#, options=[("Heuristics", 0.0), ("NoRelHeurTime", 0)] )#,gapRel=0.3)
    )
    #model.solve(GUROBI_CMD(msg=True, timeLimit= 60, MIPGap = 0.03),)
    # Check status
    if model.status == 9:  # GRB.TIME_LIMIT:
        print("Gurobi stopped due to time limit!")
        print("Model is . Exporting LP file for debugging...")
    model.writeLP("debug_model_row.lp")
    print(f"Model status: {LpStatus[model.status]}")
    print(f"Solver Status after model with status :{model.name}, {model.status}")  
    obj_value =  value(model.objective)
    print(f"Computed Objective Value : ", obj_value)
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
         sys.exit("Terminating program due to infeasibility. EXIT 1")
    if debug >=1:
        print('I solved model name =', model.name, "for iteration i = ", iter, "KP time =", KP_time, 'debug :',  debug)
    if model.name == "row_knapsack_problem":
            for var in model.variables():
                #print('row_knapsack var name =', var.name,'var value =', var.varValue)
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
         sys.exit("Terminating program due to matrix degeneration. All rows deleted. Increase the value of rho.")    
    rows_res_name = [row_names[r] for r in rows_res]
    #cols_res = [str(c) for c, _ in cols]
    #cols_res = [int(c) for c in cols_res]
    #cols_res_name = [col_names[c] for c in cols_res]
    KP_time =  print_log_output(model, KP_time, obj_value, len(rows_res), len(cols_res))
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density  = update_data(rows, cols, edges_1, rows_res, cols_res, debug)
    nb_edges_1_rem = len(edges_1_rem)
    if debug >=1:
        print('-' * 40)
        print(f"""
        Updated data after iteration: {iter} . We solved KP model :  {model.name} with delta =  {delta} 
        Found matrix of size : ({len(rows_rem)}, {len(cols_rem)}) and density : {density} and number of remaining ones : {nb_edges_1_rem}
        """)
        if debug >= 2:
            print("Remaining  Rows before colling KP col :", rows_rem)
            print("Remaining Columns  before colling KP col :", cols_rem)
        print('-' * 40)
        print()
    nb_rows_rem = len(rows_rem)
       
    model = KP_QBc(cols_rem, nb_rows_rem, nb_edges_0_rem, debug, rho) 

    model.solve(GUROBI_CMD(msg=False, timeLimit= 100))
    # Check status
    if model.status == 9 : #GRB.TIME_LIMIT:
        print("Gurobi stopped due to time limit!")
    #print("Gurobi Model Status:", model.status, "-", gp.GRB.Status.__dict__.get(model.status, "Unknown"))
    model.writeLP("debug_model_col.lp")
    print(f"Model status: {LpStatus[model.status]}")
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
         sys.exit("Terminating program due to infeasibility. EXIT 7")
    #read the result from the solver
    if debug >=1:
        print('I solved model name =', model.name, "for iteration i = ", iter)
    cols_res=[]
    cols_del=[]
    if model.name == "row_knapsack_problem":
            sys.exit("Terminating program. Non expected case row_knapsack_problem instead of column_knapsack_problem Exit 11")  
    if model.name == "column_knapsack_problem":
            for var in model.variables():
                if var.varValue == 0:
                    #print('column_knapsack var name =', var.name,'var value =', var.varValue)
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
    if debug >=0:
        print('-' * 40)
        print(f"""
        Updated data after iteration: {iter} . We solved KP model :  {model.name} with delta =  {delta} 
        Found matrix of size : ({len(rows_rem)}, {len(cols_rem)}) and density : {density} and number of remaining ones : {nb_edges_1_rem}
        """)
    if debug >=2:
        print()
        print('-' * 40)
        # print(f"Updated data after solving model = {model.name}", " Density :", density , " delta =  ", delta, "iteration:", iter)
        # nb_edges_1_rem = len(edges_1_rem)
        # print("Number of Remaining  Rows  :", len(rows_rem), "Number of Remaining Columns :", len(cols_rem))
        # print("Number of Remaining Edges_0 : ", nb_edges_0_rem, "Number of Remaining Edges_1 : ", nb_edges_1_rem)
        if debug >= 2:
            print("Remaining  Rows before colling KP col :", rows_rem)
            print("Remaining Columns  before colling KP col :", cols_rem)
        print('-' * 40)
        print()
    #sys.exit("Terminating program after KP col") 
    # end of zero_cleaner  
    return rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density, KP_time 


# ============================================================================ #
#                                    SOLVING                                   #
# ============================================================================ #


def solve(path_to_data, model, rho=1.0, delta=0.1, density_threshold=0.87, debug=0):
    """
    Function to solve the maximum biclique problem, this function reads the data,
    create a LP model using the data and return a list of rows and a list of 
    columns as a result. 
    ARGUMENTS:
    ----------
    * path_to_data: the path to the csv file.
    * model: the model to be use.
    * rho: rho of zero deletion used in all greedy approaches 
    * delta: error tolerance (percentage of accepted zeros) in the final result submatrix
    """
    #density_threshold =  0.8
    print()
    print(" ****** Density_threashold = ", density_threshold )
    print(" DEBUG  = !!!!!!!!!!!!!!!!", debug)
    KP_time = 0.0
    QBC_time = 0.0  
    #obj_total =  0.0 
    #slack_total =  0.0 
    #total_error = 0.0 
    # start reading input data
    complementary = 1
    p = path_to_data
    if p.lower().endswith('.txt'):
        #rows, cols, edges_1, row_names, col_names, M, comp_M = get_data_txt_file(path_to_data)
        rows, cols, edges_1, row_names, col_names  = get_data_txt_file(path_to_data)
        edges_0 = get_complement_edges(len(rows), len(cols), edges_1)
        if complementary == 1: 
            cols_compl_map = {col: len(cols)-degree for col, degree in cols}
            rows_compl_map = {row: len(rows)-degree for row, degree in rows}
            cols_compl = [(col,degree) for col, degree in cols_compl_map.items()]
            rows_compl = [(row, degree) for row, degree in rows_compl_map.items()]
            if debug >=4:
        # Print original and complement matrix
                print("rows_compl:\n", rows_compl)
                print("cols_compl:\n", cols_compl)
                print("Complement Matrix comp(M):\n", comp_M)
        if debug >=4:
            print("rows:\n", rows)
            print("cols:\n", cols)
            print(" row_names:\n",  row_names)
            print(" col_names:\n",  col_names)
            print(" edges_1 :\n", edges_1)
            print(" edges_0 :\n", edges_0)
            print("Original Matrix M:\n", M)

        #sys.exit("Terminating program because not done EXIT 0")
    elif p.lower().endswith('.csv'):
        rows, cols, edges_0, edges_1, row_names, col_names, df = get_data(path_to_data) 
    else:
        raise ValueError('Input need to be a matrix csv file, or a text file with a specific layout')
    ########################################################
    if complementary == 1: 
        rows = rows_compl
        cols = cols_compl 
        edge_temp = edges_1
        edges_1 = edges_0
        edges_0 = edge_temp
    ###############################################################
    nb_eges_0 = len(edges_0)
    nb_eges_1 = len(edges_1)
    col_length = len(cols)
    row_length = len(rows)
    
    if debug >= 2:
            print('-' * 40)
            if p.lower().endswith('.txt'):
                print(" Input Data in txt files :")
            if p.lower().endswith('.csv'):
                print(" Input Data in csv files :")
            print("Number Rows Data :", row_length )
            print("Number Cols Data :", col_length )
            print("Number Edges_1 :", len(edges_1))
            print("Number Edges_0 :", len(edges_0))
            if debug >= 3:
                print(" Rows Data :", rows )
                print(" Cols Data :", cols )
                print("  row_names :",  row_names)
                print("  col_names :",  col_names)
                if debug >= 4:
                    print(" Edges_1 :", edges_1)
                    print(" Edges_0 :", edges_0)
                    print("Adjacency Matrix:\n", M)
            print('-' * 40)
    # end  reading input data

    nbi_0 = 0
    nbi_1 = 0
    for _, degree in rows:
        nbi_1 = nbi_1 + degree
    nbi_0 = len(rows)*len(cols) - nbi_1
    sparsity = nbi_0/(len(rows) * len(cols))
    density = nbi_1/(len(rows) * len(cols))
    print('Initial Stats')
    print(f"Input data = {path_to_data}")
    print(f"model.name ={model}")
    print("Size of initial matrix : ", len(rows), "*", len(cols), "=", len(rows) * len(cols))
    print("number input zeros : ",nbi_0)     
    print("number input ones : ",nbi_1)            
    print("rho = ",rho, " delta : ", delta)
    print(f"Input sparsity : {sparsity:.3f}")
    print(f"Input density : {density:.3f}, density_threshold: {density_threshold}")
    print('End Initial Stats')
    print('-' * 40)
    print()
    # end fetching input data
    # start solving the problem  
    nb_edges_0 = len(edges_0)
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
    while density < density_threshold and len(rows_in) > matrix_limit and len(cols_in) > matrix_limit:
        iter+=1
        if debug >= 1:
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
            Density : {density} > {density_threshold}
            It took me : {KP_time:.3f} time
            Found matrix of size : ({len(rows_res)}, {len(cols_res)}) and density : {density}
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
    if debug == 4: 
            print("edges_1_in =", edges_1_in)
    #sys.exit("Terminating program before calling exact exit 10 !!!!")
    if model == "max_e_c": 
        rows_res, cols_res, density, nb_edges_1, QBC_time_h, QBC_time_g = warm_exact(model_name, rows_in, cols_in, row_names, col_names, edges_1_in, delta, debug, QBC_time)
        #sys.exit("Terminating program before calling exact exit 10 !!!!")  
    else:
        rows_res, cols_res, density, nb_edges_1, QBC_time_g = exact(model_name, rows_in, cols_in, row_names, col_names, edges_1_in, delta, debug, QBC_time)
        QBC_time_h = 0.0 
    if debug >=0:
        print('-' * 40)
        print(f"""
        Exit from the exact approach: {model_name}
        with delta =  {delta}
        Found matrix of size : ({len(rows_res)}, {len(cols_res)})
        and density : {density}
        and # of ones : {nb_edges_1}
        and QBC time heuristic  : {QBC_time_h:.5f}
        and QBC time  : {QBC_time_g:.5f}
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
def exact(model_name, rows, cols, row_names, col_names, edges_1, delta, debug, QBC_time):     
    if model_name == 'max_e':
        model = max_e(rows, cols, edges_1, delta)
    elif model_name == 'max_e_h':
        model = max_e_h(rows, cols, edges_1, delta)
    elif model_name == 'max_e_r':
        model = max_e_r(rows, cols, edges_1, delta)
    #elif model_name == 'max_e_c':
        #model = max_e_c(rows, cols, edges_1, delta)
    elif model_name == 'max_v':
        model = max_v(rows, cols, edges_1, delta)
    #model.solve(PULP_CBC_CMD(msg=True, timeLimit= 3600, gapRel = 0.5),)
    model.solve(GUROBI_CMD(msg=True, timeLimit= 600))#, options=[("Heuristics", 0.0), ("NoRelHeurTime", 0)] )#,gapRel=0.3)
    #model.solve(GUROBI_CMD(msg=False, timeLimit= 60, MIPGap = 0.03),)
    # Check status
    if model.status == 9:  # GRB.TIME_LIMIT:
        print("Gurobi stopped due to time limit!")
    print("Model is . Exporting LP file for debugging...")
    model.writeLP("debug_model_exact.lp")
    print(f"Model status: {LpStatus[model.status]}")
    print(f"Solver Status after model with status :{model.name}, {model.status}")  
    obj_value =  value(model.objective)
    print(f"Computed Objective Value :", obj_value)
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
         sys.exit("Terminating program due to infeasibility. EXIT 1")
    
   #Reading and analyse of the obtained results")
    #read the results from the solver
    rows_res = []
    cols_res = []
    rows_res_name = []
    cols_res_name = []
    rows_del = []
    rows_del_name = []
    cols_del = []
    cols_del_name = []
    obj_total =  0.0

    for var in model.variables():
            #if var.varValue == 1:
            if var.varValue !=  0:
                if var.name[:3] == "row":
                    rows_res = rows_res + [var.name[4:]]
                elif var.name[:3] == "col":
                    cols_res = cols_res + [var.name[4:]]

    if not rows_res or not cols_res :  # Equivalent to checking len(rows_res) == 0
       print("matrix degenerated (all rows or all columns have been deleted)")
       sys.exit("Terminating program due to matrix degeneration (all rows or all columns have been deleted). EXIT 6 ")
       #return rows_res, cols_res
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
    # print("row_names_res =", row_names_res)
    # print("col_names_res =", col_names_res)

    #rows_res_name = [row_names[r] for r in rows_res]
    #cols_res_name = [col_names[c] for c in cols_res]
    print()
    if debug >= 2: 
    #print("\n-- Debugging Step: -after row KP solving **** --")
    # print('rows=', rows)
    # print('cols=', cols)
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
    cols_res_set = set(map(int, cols_res))
    rows_res_set = set(map(int, rows_res))
    QBC_time = print_log_output(model, QBC_time, obj_value, len(rows_res), len(cols_res))
    # if model_name == 'max_e_h':
    #     QBC_time_h = print_log_output(model, QBC_time, obj_value, len(rows_res), len(cols_res))
    # else:
    #     QBC_time_g = print_log_output(model, QBC_time, obj_value, len(rows_res), len(cols_res))
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density  = update_data(rows, cols, edges_1, rows_res, cols_res, debug) 
    nb_edges_1 = len(edges_1_rem)
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
    # if model_name == 'max_v'  or model_name == 'AB_V_h'  or model_name == 'max_Surface' or model_name == 'max_Vertices' or model_name == 'max_Ones_comp'  or model_name == 'max_Ones' or model_name == 'max_e' or model_name == 'max_e_h' or model_name == 'max_e_r': 
    if debug >=2:
          print('-' * 40)
          print()
          print('Exit from the exact  approach ',  model_name,' with delta=', delta, 'Found matrix with rows_res of lenght =',len(rows_rem), ' and cols_res of lenght =',len(cols_rem), 'and density =', density)
          if debug >= 2:
                    print(" Density of the found matrix =  :", density )
                    print(" Remaining Rows  :", rows_res )
                    print(" Remaining  Cols  :", cols_res )
                    print(" Remaining Rows with degree :", rows_rem )
                    print(" Remaining  Cols with degree :", cols_rem )
    
    return rows_rem, cols_rem, density,  nb_edges_1, QBC_time 
###############################################################################
# END OF EXACT 
# ########################################################################  
def warm_exact(model_name, rows, cols, row_names, col_names, edges_1, delta, debug, QBC_time):     
    # elif model_name == 'max_e_h':
    print('-' * 40)
    print()
    print("I am in warm_exact before calling max_e_h $$$$$$$$$$$$$$$$$$")
    print()
    model = max_e_h(rows, cols, edges_1, delta)
    #model.solve(PULP_CBC_CMD(msg=True, timeLimit= 3600, gapRel = 0.5),)
    model.solve(GUROBI_CMD(msg=True, timeLimit= 600)#, options=[("Heuristics", 0.0), ("NoRelHeurTime", 0)] )#,gapRel=0.3)
    )
    #model.solve(GUROBI_CMD(msg=False, timeLimit= 60, MIPGap = 0.03),)
    # Check status
    if model.status == 9:  # GRB.TIME_LIMIT:
        print("Gurobi stopped due to time limit!")
    print("Model is . Exporting LP file for debugging...")
    model.writeLP("debug_model_heuristic.lp")
    print(f"Solver Status after model with status :{model.name}, {model.status}")  
    obj_value_h =  value(model.objective)
    print(f"Computed Objective Value in heuristic :", obj_value_h)
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
         sys.exit("Terminating program due to infeasibility. EXIT 1")
   #Reading and analyse of the obtained results")
    #read the results from the solver
    rows_res = []
    cols_res = []
    rows_res_name = []
    cols_res_name = []
    rows_del = []
    rows_del_name = []
    cols_del = []
    cols_del_name = []
    obj_total =  0.0

    for var in model.variables():
            #if var.varValue == 1:
            if var.varValue !=  0:
                if var.name[:3] == "row":
                    rows_res = rows_res + [var.name[4:]]
                elif var.name[:3] == "col":
                    cols_res = cols_res + [var.name[4:]]

    if not rows_res or not cols_res :  # Equivalent to checking len(rows_res) == 0
       print("matrix degenerated (all rows or all columns have been deleted)")
       sys.exit("Terminating program due to matrix degeneration (all rows or all columns have been deleted). EXIT 6 ")
    rows_res = [str(c) for c in rows_res]
    rows_res = [int(r) for r in rows_res]
    cols_res = [str(c) for c in cols_res]
    cols_res = [int(c) for c in cols_res]
    # Convert from 0-based to the original row/column names
    row_names_res = [row_names[r] for r in rows_res if 0 <= r < len(row_names)]
    col_names_res = [col_names[c] for c in cols_res if 0 <= c < len(col_names)]
    row_names_res = [int(r) if isinstance(r, (np.int64, np.int32)) else r for r in row_names_res]

    print()
    if debug >= 3: 
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

    cols_res_set = set(map(int, cols_res))
    rows_res_set = set(map(int, rows_res))
 
    QBC_time_h = print_log_output(model, QBC_time, obj_value_h, len(rows_res), len(cols_res))
    rows_rem, cols_rem, edges_1_rem, nb_edges_0_rem, density  = update_data(rows, cols, edges_1, rows_res, cols_res, debug) 
    nb_edges_1 = len(edges_1_rem)
    if debug >=1:
        print()
        print('-' * 40)
        print(f" results from updating data after solving model = {model.name}", " delta =  ", delta)
        print("Number of Remaining  Rows  :", len(rows_rem))
        print("Number of Remaining number Columns :", len(cols_rem))
        print("Remaining  number Edges_0 P:", nb_edges_0_rem, "Remaining  number Edges_1 :", nb_edges_1, "Density :", density , "current obj value", obj_value_h)
        print('-' * 40)
        print()
    nb_rows_rem = len(rows_rem) 
    if debug >=2:
          print('-' * 40)
          print()
          print('Exit from the exact  approach ',  model_name,' with delta=', delta, 'Found matrix with rows_res of lenght =',len(rows_rem), ' and cols_res of lenght =',len(cols_rem), 'and density =', density)
          if debug >= 1:
                    print(" Density of the found matrix =  :", density )
                    print(" Original Rows  :", rows )
                    print(" Original Cold  :", cols )
                    print(" Remaining Rows  :", rows_res )
                    print(" Remaining  Cols  :", cols_res )
                    print(" Remaining Rows with degree :", rows_rem )
                    print(" Remaining  Cols with degree :", cols_rem )

    #sys.exit("Terminating program before collong  max_e_wr EXIT 11")
    cols_res_set = set(map(int, cols_res))
    rows_res_set = set(map(int, rows_res))
    # print("I am in warm start!!!!!!!!!!!!!!!!")
    #from pulp import LpProblem, LpMaximize, LpVariable, lpSum, value

    # Initialize the model
    #model = LpProblem(name="max_e_wr", sense=LpMaximize)

    #sys.exit("Terminating program for checking")
    model = max_e_wr(rows, cols, edges_1, rows_res, cols_res, nb_edges_1, delta)

    #model.solve(PULP_CBC_CMD(msg=True, timeLimit= 3600, gapRel = 0.5),)
    #model.update()
    # for row in rows_res:
    #     lpRows[row][0].start = 1
    # for col in cols_res:
    #     lpCols[col][0].start = 1
    model.solve(GUROBI_CMD(msg=True, timeLimit= 600)) 
    #, options=[("MIPStart", 1), ("Heuristics", 0.0), ("NoRelHeurTime", 0)] )#,gapRel=0.3))
    #model.solve(GUROBI_CMD(msg=False, timeLimit= 60, MIPGap = 0.03),)
    # Check status
    #status = model.solve()
    obj_value =  value(model.objective)
    print()
    print('-' * 40)
    #print(f" !!!!!!!!!! Value of 'current': {value(lpCurrent)} and obj value {value(model.objective)} ")
    print(f"Solver Status after solving model with status and computed objective : {model.name}, {model.status}", obj_value)
    QBC_time_g = print_log_output(model, QBC_time_h, obj_value, len(rows_res), len(cols_res))
    if   model.status == 0:
         print("Model in warm start is infeasible. ")
         #QBC_time_g = print_log_output(model, QBC_time_h, obj_value, len(rows_res), len(cols_res))
         return rows_rem, cols_rem, density, nb_edges_1, QBC_time_h, QBC_time_g 
    #print(f"Solver Status after max_e_wr: {LpStatus[status]}")  
    #print(f"Objective Value in max_e_wr: {value(model.objective)}")  
    if model.status == 9:  # GRB.TIME_LIMIT:
        print("Gurobi stopped due to time limit!")
    print("Model is . Exporting LP file for debugging...")
    model.writeLP("debug_model_exact_warm_start.lp")
    #print(f"Model status: {LpStatus[model.status]}")
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
         sys.exit("Terminating program due to infeasibility. EXIT 1")
   #Reading and analyse of the obtained results")
    #read the results from the solver
    rows_res = []
    cols_res = []
    rows_res_name = []
    cols_res_name = []
    rows_del = []
    rows_del_name = []
    cols_del = []
    cols_del_name = []
    obj_total =  0.0

    for var in model.variables():
            #if var.varValue == 1:
            if var.varValue !=  0:
                if var.name[:3] == "row":
                    rows_res = rows_res + [var.name[4:]]
                elif var.name[:3] == "col":
                    cols_res = cols_res + [var.name[4:]]

    if not rows_res or not cols_res :  # Equivalent to checking len(rows_res) == 0
       print("matrix degenerated (all rows or all columns have been deleted)")
       sys.exit("Terminating program due to matrix degeneration (all rows or all columns have been deleted). EXIT 6 ")

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
    # if model_name == 'max_v'  or model_name == 'AB_V_h'  or model_name == 'max_Surface' or model_name == 'max_Vertices' or model_name == 'max_Ones_comp'  or model_name == 'max_Ones' or model_name == 'max_e' or model_name == 'max_e_h' or model_name == 'max_e_r': 
    if debug >=2:
          print('-' * 40)
          print()
          print('Exit from the exact  approach ',  model_name,' with delta=', delta, 'Found matrix with rows_res of lenght =',len(rows_rem), ' and cols_res of lenght =',len(cols_rem), 'and density =', density)
          if debug >= 2:
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
    #updates  rows_data, cols_data, edges_1 and returs them in rows_rem, cols_rem, edges_1_rem. To compute that we use rows_del, cols_del, rows_res, cols_res information
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
        print("rows_data =", rows_data) 
        print("cols_data =", cols_data) 
        #print("edges 1 =", edges_1)  
        print("rows_res_set =", rows_res_set) 
        print("cols_res_set =", cols_res_set)
        #print("  edges_1_del  =",  edges_1_del) 
        #print("  edges_1_rem  =",  edges_1_rem)

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
        row_degree_map[row] = max(0, row_degree_map[row] - 1)  # Reduce degree only if that row was involved
        col_degree_map[col] = max(0, col_degree_map[col] - 1)  # Reduce degree only if that column was involved

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
    nb_edges_1_rem = sum(degree for _, degree in rows_rem)
    size = row_rem_length * col_rem_length
    nb_edges_0_rem = size - nb_edges_1_rem
    density = nb_edges_1_rem/size
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
        print("Stats in updata_data : row_rem_length =",row_rem_length,"col_rem_length =",col_rem_length, "nb_edges_0_rem=", nb_edges_0_rem, "nb_edges_1_rem=", nb_edges_1_rem, " !!!!! density =", density)
        print()
    if size == 0:
        sys.exit("Terminating program due to matrix degeneration EXIT 2.")


    #rows_last = [tup for tup in rows_rem if tup[0] in rows_res_set] # ??? should be  rows_last = rows_res 
    #edges_0_row_last = [edge for edge in edges_0_row_rem if edge[0] in rows_rem]
    #edges_1_last = [edge for edge in edges_1_row_rem if edge[0] in rows_rem]
    #edges_0_col_last = [edge for edge in edges_0_col_rem  if edge[0] in cols_rem]
    #edges_1_col_last = [edge for edge in edges_1_col_rem  if edge[0]  in cols_rem]
    # print("cols_data_new =", cols_rem)
    # print("rows_data_new =", rows_rem)

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

# def get_data_txt_file(path):   
#     with open(path, 'r') as file:
#         content = file.readlines()

#     # Extract the number of rows, columns, and edges
#     num_row = int(content[1].split(":")[1].strip())  # |U|
#     num_col = int(content[2].split(":")[1].strip())  # |V|

#     # Initialize structures
#     deg_row = [0] * num_row
#     deg_col = [0] * num_col
#     edges_1 = []

#     # Read edges (starting from line 5, skipping headers)
#     for line in content[4:]:
#         parts = line.strip().split()
#         if len(parts) == 2:  # Ensure it's a valid edge
#             u, v = map(int, parts)
#             edges_1.append((u, v))
#             deg_row[u] += 1
#             deg_col[v] += 1

#     # Create original binary matrix M
#     df_M = pd.DataFrame(0, index=range(num_row), columns=range(num_col))
#     for u, v in edges_1:
#         df_M.iloc[u, v] = 1  

#     # Compute complement edges
#     all_edges = set(itertools.product(range(num_row), range(num_col)))
#     original_edges = set(edges_1)
#     edges_0 = list(all_edges - original_edges)

#     # Create complement matrix
#     df_comp_M = pd.DataFrame(1, index=range(num_row), columns=range(num_col))
#     for u, v in edges_1:
#         df_comp_M.iloc[u, v] = 0  

#     # Prepare row and column data
#     rows = list(zip(range(num_row), deg_row))
#     cols = list(zip(range(num_col), deg_col))
#     row_names = list(range(num_row))
#     col_names = list(range(num_col))

#     return rows, cols, edges_0, edges_1, row_names, col_names, df_M, df_comp_M


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


# def get_data_txt_file(path):  # working version 
#     with open(path, 'r') as file:  # Use 'with' to properly handle file closing
#         content = file.readlines()
    
#     name = content[0][2:-1].strip()  # Strip removes unwanted whitespace
#     num_row = int(content[1].split(":")[1].strip())  # Extract numbers correctly
#     num_col = int(content[2].split(":")[1].strip())
#     num_edge = int(content[3].split(":")[1].strip())

#     deg_row = [0] * num_row
#     deg_col = [0] * num_col
#     edges = []
#     df = pd.DataFrame(0, index=range(num_row), columns=range(num_col))  # Corrected dataframe initialization

#     for line in content[4:]:
#         splitted_line = line.strip().split()  # Automatically splits on any whitespace
#         u, v = int(splitted_line[0]), int(splitted_line[1])
        
#         edges.append((u, v))
#         deg_row[u] += 1
#         deg_col[v] += 1
#         df.iloc[u, v] = 1 

#     rows_data = list(zip(range(num_row), deg_row))
#     cols_data = list(zip(range(num_col), deg_col))

#     return rows_data, cols_data, edges, range(num_row), range(num_col), df



def print_log_output(prob, run_time, obj, nb_rows, nb_cols):
    """Print the log output and problem solutions.
    ARGUMENTS:
    ----------
    * prob: an solved LP model (pulp.LpProblem)
    """
    print()
    print('-' * 70)
    print('Stats')
    #print('-' * 70) 
    #print(f'Number variables: {prob.numVariables()}')
    #print(f'Number constraints: {prob.numConstraints()}')
    run_time = run_time + prob.solutionCpuTime
    print(f' Model solved :', prob.name, f' Local Time:  - (real) {prob.solutionTime:.5f}', f'- (Local CPU) {prob.solutionCpuTime:.5f}', f'- (Global CPU) {run_time:.5f}')
    #print(f'Number variables: {prob.numVariables()}')) 
    #print(f'- (real) {prob.solutionTime}')
    #print(f'- (CPU) {prob.solutionCpuTime}')
    print()
    print(f' Solve status: {LpStatus[prob.status]}', f'Objective value: {prob.objective.value()}', "nb_rows", nb_rows, "nb_cols", nb_cols )
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
    return run_time 

def parse_arguments():
    """Parse the input arguments and retrieve the choosen resolution method and
    the instance that must be solve."""
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
        help='Select the error rho value (tolerance) ',
    )

    argparser.add_argument(
        '--threshold', dest='threshold', required=False, default=0.87, type=float,
        help='Above this value greedy  approaches (KP zero clreaners)  are used? Below, exact ',
    )


    argparser.add_argument(
        '--debug', nargs='?', required=False, const=1, type=int, default=0, 
        help='Select the debug value  ',
    )

    #parser.add_argument('--debug', nargs='?', const=1, type=int, default=0)  


    arg = argparser.parse_args()

    if arg.model not in ['max_e', 'max_e_h', 'max_e_r','max_e_c','max_v', 'KP_QBr', 'KP_QBc']:
        argparser.print_help()
        sys.exit(1)

    return (arg.filepath, arg.model, arg.rho, arg.delta,  arg.threshold, arg.debug)


if __name__ == '__main__':

    # Read the arguments
    file_path, selected_model, rho, delta,  threshold, debug, = parse_arguments()


    rows_res, cols_res, density, nb_ones, iter, KP_time,  kp_density, nb_kp_rows, nb_kp_cols, nb_kp_ones, QBC_time_h, QBC_time_g  =  solve(file_path,selected_model,rho,delta, threshold, debug) 
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
    End of computations for = {file_path} 
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

    # print(" Computed density : ", density, " Number computed ones:",nb_ones )
    # print(" Size of Remaining matrix : (", len(rows_res), ",", len(cols_res), ")") 
    #print(" Number Remaining  Cols  :", len(cols_res) )
    if debug >=1:
        rows_last = [str(c) for c, _ in rows_res]
        rows_last = [int(r) for r in rows_last]
        cols_last = [str(c) for c, _ in cols_res]
        cols_last = [int(c) for c in cols_last]   
        print(" Remaining Rows last :", rows_last )
        print(" Remaining  Cols  last :", cols_last )




