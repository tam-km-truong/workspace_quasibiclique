# -*- coding=utf-8 -*-

"""
Skeleton for TP project: Searching Maximum Quasi-Bicliques
Student:

Install Gurobi Solver with free academic licence using your University account at: https://www.gurobi.com/features/academic-named-user-license/

Install Pandas, using PyPI or Conda. A detailed installation tutorial can be find at: https://pandas.pydata.org/docs/getting_started/install.html
"""

import sys
import re

from argparse import ArgumentParser

import pandas as pd
import itertools
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


# ============================================================================ #
#                     LP MODEL - König theorem classical                        #
# ============================================================================ #

def AB_V(rows_data, cols_data, edges, epsilon):
    """<
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    epsilon: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    """
    model = LpProblem(name="AB_V", sense=LpMaximize)

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

    lpCells = {} #can be improved because we dont need quadratic variable. This state is temperary 
    for row, _ in rows_data:
        for col, _ in cols_data:
            if (row, col) in edges:
                lpCells[(row, col)] = (LpVariable(
                    #f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1), 1)
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)
            else:
                lpCells[(row, col)] = (LpVariable(  
                    #f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1), 0)
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)
    # Objective function: Maximize sum of selected row and column variables
    model += lpSum(
        [lpvar for lpvar, _ in lpRows.values()] +
        [lpvar for lpvar, _ in lpCols.values()]), "max_sum_vertices"
    

    # Constraints for row and column thresholds
    row_threshold = 2
    col_threshold = 2
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    model += obj_rowsum == lpSum(lpSum_row), "obj_sum_row"
    model += obj_colsum == lpSum(lpSum_col), "obj_sum_col"
    # print()
    # print('-' * 40)
    # print('row_threshold=', row_threshold )
    # print('col_threshold=', col_threshold )
    # print()
    # print('-' * 40)
    #just for testing 
    # for row, col in lpCells:
    #     #model += (lpRows[row][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_1'
    #     #model += (lpCols[col][0] >= lpCells[(row, col)][0]), f'ce ll_{row}_{col}_2'
    #     model += ((lpRows[row][0]+lpCols[col][0])/2 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'
    #     model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_4'
    #just for testing     
    # Add row density constraints
    #__row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    #__col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_row, lpSlack_row, epsilon)
    __col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_col, lpSlack_col, epsilon)

    return model 

def AB_V_h(rows_data, cols_data, edges, epsilon):
    """<
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
    epsilon: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    """
    model = LpProblem(name="AB_V_h", sense=LpMaximize)

    # Variables for rows and columns

    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}             
    # Objective function: Maximize sum of selected row and column variables
    model += lpSum(
        [lpvar for lpvar, _ in lpRows.values()] +
        [lpvar for lpvar, _ in lpCols.values()], "max_sum_vertices"
    )

    # Constraints for row and column thresholds
    row_threshold = 2
    col_threshold = 2
    print()
    print('-' * 40)
    print('row_threshold=', row_threshold )
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    print('col_threshold=', col_threshold )
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    print()
    print('-' * 40)
    # Add row density constraints
    #__row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    #__col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)

    return model 


def AB_E(rows_data, cols_data, edges, epsilon):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    epsilon: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    """
    model = LpProblem(name="AB_E", sense=LpMaximize)

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

    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            if (row, col) in edges:
                lpCells[(row, col)] = (LpVariable(
                    #f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1), 1)
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)
            else:
                lpCells[(row, col)] = (LpVariable(  
                    #f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1), 0)
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)
  
    # Objective function: Maximize sum of selected row and column variables
    #
    model += lpSum([cellValue*lpvar for lpvar,
                   cellValue in lpCells.values()]), 'maximize_weight'
    # 
    #model += lpSum(
    #    [lpvar for lpvar, _ in lpRows.values()] +
    #    [lpvar for lpvar, _ in lpCols.values()]
    #), "max_sum_vertices"

    # Constraints for row and column thresholds
    row_threshold = 2
    col_threshold = 2
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    # print()
    # print('-' * 40)
    # print('row_threshold=', row_threshold )
    #print('col_threshold=', col_threshold )
    # print()
    # print('-' * 40)
    for row, col in edges:
    #for row, col in lpCells:
        #model += (lpRows[row][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_1'
        #model += (lpCols[col][0] >= lpCells[(row, col)][0]), f'ce ll_{row}_{col}_2'
        model += (lpRows[row][0]+lpCols[col][0])/2 >= lpCells[(row, col)][0], f'cell_{row}_{col}_3'
        model += lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0], f'cell_{row}_{col}_4'
    # Add row/cols density constraints
    # __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    # __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    model += obj_rowsum == lpSum(lpSum_row), "obj_sum_row"
    model += obj_colsum == lpSum(lpSum_col), "obj_sum_col"
    # Add row density constraints
    __row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_row, lpSlack_row, epsilon)
    __col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_col, lpSlack_col, epsilon)


    return model 

def AB_E_r(rows_data, cols_data, edges, epsilon):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    epsilon: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    """
    model = LpProblem(name="AB_E_r", sense=LpMaximize)

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

    lpCells = {} #can be improved because we dont need quadratic variable. This state is temporary 
    for row, _ in rows_data:
        for col, _ in cols_data:
            if (row, col) in edges:
                lpCells[(row, col)] = (LpVariable(
                    #f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1), 1)
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)
            else:
                lpCells[(row, col)] = (LpVariable(  
                    #f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1), 0)
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)
    # print()
    # print('-' * 40)
    # print('lpCells=')
    # print(lpCells)
    # print()
    # print('-' * 40)
    # Objective function: Maximize sum of selected row and column variables
    #
    model += lpSum(lpSum_row) + lpSum(lpSum_col), "max_sum_rows_columns"
    # +lpSum(
    #     [lpvar for lpvar, _ in lpRows.values()] +
    #     [lpvar for lpvar, _ in lpCols.values()]), "max_sum_vertices_rows_columns"
    #model += lpSum([cellValue*lpvar for lpvar, cellValue in lpCells.values()]), 'maximize_weight'
    # 
    #model += lpSum(
    #    [lpvar for lpvar, _ in lpRows.values()] +
    #    [lpvar for lpvar, _ in lpCols.values()]
    #), "max_sum_vertices"

    # Constraints for row and column thresholds
    row_threshold = 2
    col_threshold = 2
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    model += obj_rowsum == lpSum(lpSum_row), "obj_sum_row"
    model += obj_colsum == lpSum(lpSum_col), "obj_sum_col"
    # Add row density constraints
    __row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_row, lpSlack_row, epsilon)
    __col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_col, lpSlack_col, epsilon)

    return model

def AB_E_c_r(rows_data, cols_data, edges, epsilon):
    """
    Arguments:
    ----------
    rows_data: list of tuples (row, degree) of rows in the matrix.
    cols_data: list of tuples (col, degree) of columns in the matrix.
    edges: list of tuples (row, col) corresponding to the ones of the matrix.
    epsilon: float, error tolerance for density constraints.

    Returns:
    --------
    LpProblem:
        The ILP model.
    """
    model = LpProblem(name="AB_E_c_r", sense=LpMaximize)

    # Variables for rows and columns

    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                               lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            if (row, col) in edges:
                lpCells[(row, col)] = (LpVariable(
                    #f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1), 1)
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)
            else:
                lpCells[(row, col)] = (LpVariable(  
                    #f'cell_{row}_{col}', cat='Continuous', lowBound=0, upBound=1), 0)
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)
    # print()
    # print('-' * 40)
    # print('lpCells=')
    # print(lpCells)
    # print()
    # print('-' * 40)
    # Objective function: Maximize sum of selected row and column variables
    #
    model += lpSum([cellValue*lpvar for lpvar,
                   cellValue in lpCells.values()]), 'maximize_weight'
    # 
    #model += lpSum(
    #    [lpvar for lpvar, _ in lpRows.values()] +
    #    [lpvar for lpvar, _ in lpCols.values()]
    #), "max_sum_vertices"

    # Constraints for row and column thresholds
    row_threshold = 2
    col_threshold = 2
    print()
    print('-' * 40)
    print('row_threshold=', row_threshold )
    model += lpSum(lpvar for lpvar, _ in lpRows.values()) >= row_threshold, "row_threshold"
    print('col_threshold=', col_threshold )
    model += lpSum(lpvar for lpvar, _ in lpCols.values()) >= col_threshold, "col_threshold"
    print()
    print('-' * 40)
    for row, col in edges:
    #for row, col in lpCells:
        #model += (lpRows[row][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_1'
        #model += (lpCols[col][0] >= lpCells[(row, col)][0]), f'ce ll_{row}_{col}_2'
        #model += ((lpRows[row][0]+lpCols[col][0])/2 >= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'
        model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_4'
        #########################################
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
        
         #########################################


    # Add row density constraints
    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)

    return model 


def __col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_col, lpSlack_col, epsilon):
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
    epsilon: float, error tolerance for density constraints.
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
            lpSum_col[col] - (1 - epsilon) * lpSum(lpRows[row][0] for row, _ in rows_data), f"col_slack_{col}" 
        #    >= (lpCols[col][0]-1) * Big_M, f"col_slack_{col}"
        )
        model += (lpSlack_col[col] >= (lpCols[col][0]-1) * Big_M, f"col_err_rate_1_{col}"
            # lpSum_col[col] - (1 - epsilon) * lpSum(lpRows[row][0] for row, _ in rows_data) 
            # >= (lpCols[col][0]-1) * Big_M, f"col_err_rate_1_{col}"
        )
        # Constraint for col density lower bound
        
        # model += (
        #     lpSum_col[col] - (1 - epsilon) * lpSum(lpRows[row][0] for row, _ in rows_data) <=
        #     -mu + lpCols[col][0] * Big_M
        # ), f"col_err_rate_0_{col}"

def __row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, lpSum_row, lpSlack_row, epsilon):
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
    epsilon: float, error tolerance for density constraints.
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
            lpSum_row[row] - (1 - epsilon) * lpSum(lpCols[col][0] for col, _ in cols_data), f"row_slack_{row}" 
           # >= (lpRows[row][0]-1) * Big_M, f"row_err_rate_1_{row}"
        )

        model += (
            lpSlack_row[row] >= (lpRows[row][0]-1) * Big_M, f"row_err_rate_1_{row}"
        )

        # Constraint for row density lower bound
        
        # model += (
        #     lpSum(lpCols[col][0] for col in row_edges) - (1 - epsilon) * lpSum(lpCols[col][0] for col, _ in cols_data) <=
        #     -mu + lpRows[row][0] * Big_M
        # ), f"row_err_rate_0_{row}"

def __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon):
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
    epsilon: float, error tolerance for density constraints.
    """
    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M= Big_R+Big_C 
    #print('Big_R=', Big_R)
    #print('Big_C=', Big_C)

    for col, _ in cols_data:
        #print(f"Adding col density constraints for col {col}:")
        col_edges = [u for u, v in edges if v == col]
        #print(f"Col edges: {col_edges}") 
        #print(f"lpCols[col][0]: {lpCols[col][0]}")             
        # Constraint for col density upper bound
        #print(f"Sum of edge variables: {lpSum(lpRows[row][0] for row in col_edges)}")
        #print(f"Sum of row variables: {lpSum(lpRows[row][0] for  row, _ in rows_data)}")
        model += (
            lpSum(lpRows[row][0] for row in col_edges) - (1 - epsilon) * lpSum(lpRows[row][0] for row, _ in rows_data) >= 
            (lpCols[col][0]-1) * Big_M
        ), f"col_err_rate_1_{col}"

        # Constraint for col density lower bound
        
        # model += (
        #     lpSum(lpRows[row][0] for row in col_edges) - (1 - epsilon) * lpSum(lpRows[row][0] for row, _ in rows_data) <=
        #     -mu + lpCols[col][0] * Big_M
        # ), f"col_err_rate_0_{col}"

def __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon):
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
    epsilon: float, error tolerance for density constraints.
    """
    mu = 0.0001
    Big_R = len(rows_data) + 1
    Big_C = len(cols_data) + 1
    Big_M = Big_R+Big_C
    #print('Big_R=', Big_R)
    #print('Big_C=', Big_C)

    for row, _ in rows_data:
        #print(f"Adding row density constraints for row {row}:")
        row_edges = [v for u, v in edges if u == row]
        #print(f"Row edges: {row_edges}") 
        #print(f"lpRows[row][0]: {lpRows[row][0]}")             
        # Constraint for row density upper bound
        #print(f"Sum of edge variables: {lpSum(lpCols[col][0] for col in row_edges)}")
        #print(f"Sum of column variables: {lpSum(lpCols[col][0] for  col, _ in cols_data)}")
        model += (
            lpSum(lpCols[col][0] for col in row_edges) - (1 - epsilon) * lpSum(lpCols[col][0] for col, _ in cols_data) >= 
            (lpRows[row][0]-1) * Big_M
        ), f"row_err_rate_1_{row}"

        # Constraint for row density lower bound
        
        # model += (
        #     lpSum(lpCols[col][0] for col in row_edges) - (1 - epsilon) * lpSum(lpCols[col][0] for col, _ in cols_data) <=
        #     -mu + lpRows[row][0] * Big_M
        # ), f"row_err_rate_0_{row}"

def König_V(rows_data, cols_data, edges):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuples (row,col) corresponding to the zeros of the matrix.
    """

    # ------------------------------------------------------------------------ #
    # Model with minimization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='König_V', sense=LpMinimize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {row: (LpVariable(f'row_{row}', 
    cat='Continous',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data} 
    #cat='Binary', lowBound=0), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Continous',
                    lowBound=0, upBound=1), degree) for col, degree in cols_data}
    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    # model += lpSum([degree*lpvar for lpvar, degree in lpRows.values()] +
    #                [degree*lpvar for lpvar, degree in lpCols.values()]), 'min_weighted_cover'

    model += lpSum([lpvar for lpvar, _ in lpRows.values()] +
                   [lpvar for lpvar, _ in lpCols.values()]), 'min_vertices'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    for row, col in edges:
        model += (lpRows[row][0]+lpCols[col][0] >= 1), f'edge_{row}_{col}'
    #model += (lpRows['r1'] + lpRows['r2'] <= 1,'constraint uncomp')

    return model
# ============================================================================ #
#                     LP MODEL - König with degrees                       #
# ============================================================================ #


def König_E(rows_data, cols_data, edges):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuples (row,col) corresponding to the zeros  of the matrix.
    """

    # ------------------------------------------------------------------------ #
    # Model with minimization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='König_E', sense=LpMinimize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {row: (LpVariable(f'row_{row}', 
    cat='Continous',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', 
    cat='Continous', lowBound=0, upBound=1), degree) for col, degree in cols_data}
    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([degree*lpvar for lpvar, degree in lpRows.values()] +
                   [degree*lpvar for lpvar, degree in lpCols.values()]), 'min_degrees'

    # model += lpSum([lpvar for lpvar, _ in lpRows.values()] +
    #                [lpvar for lpvar, _ in lpCols.values()]), 'min_vertices'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    for row, col in edges:
        model += (lpRows[row][0]+lpCols[col][0] >= 1), f'edge_{row}_{col}'

    #model += (lpSum(lpRows[row][0] for row, _ in rows_data)  <= 3, 'row_limits')
              
    #model += (lpRows['row_0'] + lpRows['row_2'] <= 1,'constraint uncomp')
    #model += (lpRows[0][0] + lpRows[2][0] <= 1,'constraint uncomp')
    #model += (lpSum([lpvar for lpvar, _ in lpRows.values()])  <= 2,'row_limits_bis',    )
    #model += (lpSum([lpvar for lpvar, _ in lpCols.values()])  <= 2,'col_limits_bis',    )
    
    
    return model


# ============================================================================ #
#                LP MODEL - Deleting rows/columns for zeros eliminating                    #
# ============================================================================ #

def  minDel_RC(rows_data, cols_data, edges, epsilon=0.3): 
#     bigraph: ZerosBiGraph,
#     epsilon: float,
# ) -> tuple[LpProblem, VertexChoicesLP, EdgeChoicesLP]:
    """Minimizing the number of deleted rows.
    Implement the LP model for deleting minimum rows and columns. 
    In this model, we introduce epsilon with the goal of alowing errors in the result. 
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros  of the matrix.
    * epsilon: a percentage of the zeros that will be accepted in the final matrix 
    """


    model = LpProblem(name='minDel_rows', sense=LpMinimize)

    # 
    #
    # Variables
    #
    # u_vertex_choices = __u_vertex_choices_variables(bigraph)
    # edge_choices = __edge_choices_variables(bigraph)
    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                    lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpEdges = {(row, col): LpVariable(f'edge_{row}_{col}', cat='Integer',
                                      lowBound=0, upBound=1) for row, col in edges}

    #
    # Objective function
    #
    model += lpSum([lpvar for lpvar, _ in lpRows.values()] ), 'min_sum_vertices'

    # 
    #
    # Constraints
    #
    for row, col in edges:
        model += lpRows[row][0] >= lpEdges[(row, col)], f'edge_{row}_{col}'
        
    #__select_edge_select_u_endpoint(prob, bigraph, u_vertex_choices, edge_choices)

    #__at_most_epsilon_remaining_zeros(prob, bigraph, epsilon, edge_choices)
    model += (
            lpSum(
                [degree*lpvar for lpvar, degree in lpRows.values()]
                    # var * deg
                    # for var, deg in zip(u_vertex_choices, bigraph.u_bipart_degrees())
                
                )
            # + lpSum(
            #     (
            #         var * deg
            #         for var, deg in zip(v_vertex_choices, bigraph.v_bipart_degrees())
            #     ),
            # )
            >= (1 - epsilon) * len(edges)
        ), "sensitivity"

    return model 

def minDel_RC_original(rows_data, cols_data, edges, epsilon=0.3):
    """
    Implement the LP model for deleting minimum rows and columns. 
    In this model, we introduce epsilon with the goal of alowing errors in the result. 
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros  of the matrix.
    * epsilon: a percentage of the zeros that will be accepted in the final matrix 
    """

    # ------------------------------------------------------------------------ #
    # Model with minimization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='minDel_rows_cols', sense=LpMinimize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                    lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpEdges = {(row, col): LpVariable(f'edge_{row}_{col}', cat='Integer',
                                      lowBound=0, upBound=1) for row, col in edges}

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    # model += lpSum([degree*lpvar for lpvar, degree in lpRows.values()] +
    #                [degree*lpvar for lpvar, degree in lpCols.values()]), 'precise_min_weighted_cover'
    model += lpSum([lpvar for lpvar, _ in lpRows.values()] +
                   [lpvar for lpvar, _ in lpCols.values()]), 'min_sum_vertices'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    for row, col in edges:
        model += (lpRows[row][0]+lpCols[col][0] >=
                  lpEdges[(row, col)]), f'edge_{row}_{col}'
        
    adjustment_ = 1

    model += (lpSum(lpEdges) >= (1-epsilon*adjustment_) * len(edges)), f'sensitivity'

    # nbi0 = len(edges)
    # nbi1 = sum([degree for _,degree in rows_data])
    # nbi2 = len(rows_data)*len(cols_data)-nbi0
    # print()
    # print('-' * 40)
    # print("nbi0 = ",nbi0)
    # print("nbi1 = ",nbi1)
    # print("nbi2 = ",nbi2)
    # print(nbi0, nbi1, nbi2)
    # print('-' * 40)
    # if nbi0>nbi1:
    #     diff = nbi0-nbi1
    #     model += (lpSum(lpEdges) >= diff+(1-epsilon)*nbi1), f'sensitivity'
    # else:
    #     model += (lpSum(lpEdges) >= (1-epsilon) * len(edges)), f'sensitivity'

    return model
#==========================================================================#
#                LP MODEL - Deleting rows/columns for zeros eliminating                    #
# ============================================================================ #


def minDel_Ones(rows_data, cols_data, edges, epsilon):
    """
    Implement the LP model for deleting minimum rows and columns. 
    In this model, we introduce epsilon with the goal of alowing errors in the result. 
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros  of the matrix. 
    * epsilon: a percentage of zeros that will be tolerated in the final matrix 
    """

    # ------------------------------------------------------------------------ #
    # Model with minimization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='minDel_Ones', sense=LpMinimize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {row: (LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Integer',
                    lowBound=0, upBound=1), degree) for col, degree in cols_data}
    lpEdges = {(row, col): LpVariable(f'edge_{row}_{col}', cat='Integer',
                                      lowBound=0, upBound=1) for row, col in edges}

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([degree*lpvar for lpvar, degree in lpRows.values()]
               +   [degree*lpvar for lpvar, degree in lpCols.values()]), 'min_weighted_rows/cols'
    #model += lpSum([lpvar for lpvar, _ in lpRows.values()]), 'min_weighted_rows/cols'
                   
    #               + [lpvar for lpvar, _ in lpCols.values()]), 'min_weighted_rows/cols'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ ##
    # for row, col in edges:
    #     model += (lpRows[row][0]+lpCols[col][0] >=
    #               lpEdges[(row, col)]), f'edge_{row}_{col}'
        
    # adjustment_ = 1

    # model += (lpSum(lpEdges) >= (1-epsilon*adjustment_) * len(edges)), f'sensitivity'

    # return model

    for row, col in edges:
        model += (lpRows[row][0]+lpCols[col][0] >=
                  lpEdges[(row, col)]), f'edge_{row}_{col}'
        

    nbi0 = len(edges)
    nbi1 = sum([degree for _,degree in rows_data])
    nbi2 = len(rows_data)*len(cols_data)-nbi0
    # print()
    # print('-' * 40)
    # print("nbi0 = ",nbi0)
    # print("nbi1 = ",nbi1)
    # print("nbi2 = ",nbi2)
    # print(nbi0, nbi1, nbi2)
    # print('-' * 40)
    # if nbi0>nbi1:
    #     diff = nbi0-nbi1
    #     model += (lpSum(lpEdges) >= diff+(1-epsilon)*nbi1), f'sensitivity'
    # else:
    model += (lpSum(lpEdges) >= (1-epsilon) * len(edges)), f'sensitivity'

    return model
# # ============================================================================ #
#                     LP MODEL - MAXIMIZE DESIRED CELL                         #
# ============================================================================ #


def max_Ones(rows_data, cols_data, edges, epsilon):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros  of the matrix.
    * epsilon: percentage of accepted zeros 
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

    model += (lpSum([(1-cellValue)*lpvar for lpvar, cellValue in lpCells.values()]) <= epsilon *
              lpSum([lpvar for lpvar, _ in lpCells.values()])), f'err_rate'

    return model



def max_Ones_comp(rows_data, cols_data, edges, epsilon):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros  of the matrix.
    * epsilon: percentage of accepted zeros 
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

    model += (lpSum([(1-cellValue)*lpvar for lpvar, cellValue in lpCells.values()]) <= epsilon *
              lpSum([lpvar for lpvar, _ in lpCells.values()])), f'err_rate'

    return model


def max_Surface(rows_data, cols_data, edges, epsilon=0.1):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple (row,col) corresponding to the zeros  of the matrix. 
    * epsilon: percentage of accepted zeros 
    """

    # ------------------------------------------------------------------------ #
    # Model with maximization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='maximize_surface', sense=LpMaximize)

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
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)
            else:
                lpCells[(row, col)] = (LpVariable(  
                    f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([lpvar for lpvar,_ in lpCells.values()]), 'maximize_surface_matrix'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
 
    for row, col in lpCells:
        model += (lpRows[row][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_1'
        model += (lpCols[col][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_2'
        model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'

    model += (lpSum([(1-cellValue)*lpvar for lpvar, cellValue in lpCells.values()]) <= epsilon *
              lpSum([lpvar for lpvar, _ in lpCells.values()])), f'err_rate'

    #model += (lpRows['r1'] + lpRows['r2'] <= 1,'constraint uncomp')

    return model

def max_Vertices(rows_data, cols_data, edges, epsilon=0.1):
    """
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * edges: list of tuple of the coordination (row,col) of zeros of the matrix.
    * epsilon: percentage of accepted zeros 
    """

    # ------------------------------------------------------------------------ #
    # Model with maximization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='maximize_sum_of_rows_cols', sense=LpMaximize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #
    lpRows = {row: (LpVariable(f'row_{row}', cat='Binary'), degree) for row, degree in rows_data}
    lpCols = {col: (LpVariable(f'col_{col}', cat='Binary'), degree) for col, degree in cols_data}
    lpCells = {}
    for row, _ in rows_data:
        for col, _ in cols_data:
            if (row, col) in edges:
                lpCells[(row, col)] = (LpVariable(f'cell_{row}_{col}', cat='Binary'), 0)
            else:
                lpCells[(row, col)] = (LpVariable(f'cell_{row}_{col}', cat='Binary'), 1)

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([lpvar for lpvar, _ in lpRows.values()] +
                   [lpvar for lpvar, _ in lpCols.values()]), 'max_sum_vertices'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
 
    for row, col in lpCells:
        model += (lpRows[row][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_1'
        model += (lpCols[col][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_2'
        model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'

    model += (lpSum([(1-cellValue)*lpvar for lpvar, cellValue in lpCells.values()]) <= epsilon *
              lpSum([lpvar for lpvar, _ in lpCells.values()])), f'err_rate'

    #model += (lpRows['r1'] + lpRows['r2'] <= 1,'constraint uncomp')

    return model

# ============================================================================ #
#                         LP MODEL - KNAPSACK MODEL                            #
# ============================================================================ #
# def KP_QBr(rows_data,  edges, epsilon=0.1):
#     """
#     Knapsack problem model
#     """
#     print('-' * 40)
#     print('INput data in KP ***************')
#     print()
#     print('edges =', edges)
#     print('rows_data =', rows_data)
#     print('cols_data =', cols_data)
#     print('-' * 40)
    
#     # ------------------------------------------------------------------------ #
#     # Model with minimization
#     # ------------------------------------------------------------------------ #
#     model = LpProblem(name='knapsack_problem', sense=LpMinimize)

#     # ------------------------------------------------------------------------ #
#     # Variables
#     # ------------------------------------------------------------------------ #
#     lpRows = [(LpVariable(f'row_{row}', cat='Integer',
#                           lowBound=0, upBound=1), degree) for row, degree in rows_data]
#     lpCols = [(LpVariable(f'col_{col}', cat='Integer',
#                           lowBound=0, upBound=1), degree) for col, degree in cols_data]

#     # ------------------------------------------------------------------------ #
#     # Objective function (ONLY includes rows, NOT columns)
#     # ------------------------------------------------------------------------ #
#     model += lpSum([degree * lpvar for lpvar, degree in lpRows]), 'knapsack'

#     # ------------------------------------------------------------------------ #
#     # Constraints (NO `lpCols` in "sensitivity" constraint)
#     # ------------------------------------------------------------------------ #
#     nb_edges = len(edges)
#     row_length = len(rows_data)

#     # ✅ Only `lpRows` variables appear in this constraint
#     model += lpSum([(row_length - degree) * lpvar for lpvar, degree in lpRows]) >= (1 - epsilon) * nb_edges, 'sensitivity'

#     return model

def KP_QBr(rows_data, edges, epsilon=0.1):
    """
    In this model, we use knapsack model
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * epsilon: percentage of accepted undesired cell over accepted desired cell
    """
    # print('-' * 40)
    # print('INput data in KP ***************')
    # print()
    # print('edges =', edges)
    # print('rows_data =', rows_data)
    # #print('cols_data =', cols_data)
    # print('-' * 40)
    
    # ------------------------------------------------------------------------ #
    # Model with minimization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='row_knapsack_problem', sense=LpMinimize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #

    lpRows = [(LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data]
    # lpCols = [(LpVariable(f'col_{col}', cat='Integer',
    #                 lowBound=0, upBound=1), degree) for col, degree in cols_data]
    # Assign initial values correctly
    # for var, _ in lpCols:  # Unpack tuple (variable, degree)
    #    var.varValue = 0  # Set initial value
    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([degree*lpvar for lpvar, degree in lpRows]), 'knapsack' # + lpSum([0*lpvar for lpvar, degree in lpCols]), 'knapsack'


    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    nb_edges = len(edges)
    print("# of edges in KP=", len(edges))
    print("edges in KP=", edges)
    row_length = len(rows_data)
    #col_length = len(cols_data)
    # Ensure col_* variables are included in the model
    # for var, _ in lpCols:
    #     model += var == 0, f"force_col_{var.name}"  # Dummy constraint

    model += lpSum([(row_length - degree) * lpvar for lpvar, degree in lpRows]) >= (1 - epsilon) * nb_edges, 'sensitivity'

    #model += lpSum([(row_length-degree)*lpvar for lpvar, degree in lpRows]) >= (1-epsilon) * nb_edges, f'sensitivity'
    # +       lpSum([(col_length-degree)*lpvar for lpvar, degree in lpCols])  >= (1-epsilon) * nb_edges), f'sensitivity'
    #


    return model

def KP_QBc(cols_data, edges, epsilon=0.1):
    """
    In this model, we use knapsack model
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * epsilon: percentage of accepted undesired cell over accepted desired cell
    """
    print('-' * 40)
    print('INput data in KPc ***************')
    print()
    print('edges =', edges)
    print('cols_data =', cols_data)
    print('-' * 40)
    
    # ------------------------------------------------------------------------ #
    # Model with minimization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='column knapsack_problem', sense=LpMinimize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #

    # lpRows = [(LpVariable(f'row_{row}', cat='Integer',
    #                 lowBound=0, upBound=1), degree) for row, degree in rows_data]
    lpCols = [(LpVariable(f'col_{col}', cat='Integer',
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
    nb_edges = len(edges)
    col_length = len(cols_data)
    # Ensure col_* variables are included in the model
    # for var, _ in lpCols:
    #     model += var == 0, f"force_col_{var.name}"  # Dummy constraint

    model += lpSum([(col_length - degree) * lpvar for lpvar, degree in lpCols]) >= (1 - epsilon) * nb_edges, 'sensitivity'

    #model += lpSum([(row_length-degree)*lpvar for lpvar, degree in lpRows]) >= (1-epsilon) * nb_edges, f'sensitivity'
    # +       lpSum([(col_length-degree)*lpvar for lpvar, degree in lpCols])  >= (1-epsilon) * nb_edges), f'sensitivity'
    #


    return model
# ============================================================================ #
#                                    SOLVING                                   #
# ============================================================================ #


def solve(path_to_data, model, epsilon=0.1):
    """
    Function to solve the maximum biclique problem, this function reads the data,
    create a LP model using the data and return a list of rows and a list of 
    columns as a result. 
    ARGUMENTS:
    ----------
    * path_to_data: the path to the csv file.
    * model: the model to be use.
    * epsilon: percentage of errors (accepted edges) in the final result submatrix
    """
    obj_total =  0.0 
    slack_total =  0.0 
    total_error = 0.0 
    p = path_to_data
    if p.lower().endswith('.txt'):
        rows, cols, edges, row_names, col_names, df = get_data_txt_file(path_to_data, model)
        # Get complement edges
        num_row = len(row_names)
        num_col = len(col_names)
        CompEdges = get_complement_edges(num_row, num_col, edges)
        print("Complement Edges:", CompEdges)      
    elif p.lower().endswith('.csv'):
        rows, cols, edges, row_names, col_names, df = get_data(path_to_data, model)
    else:
        raise ValueError('Input need to be a matrix csv file, or a text file with a specific layout')
    
    print('-' * 40)
    print("Input data :")
    print("Rows Data:", rows)
    print("Columns Data:", cols)
    print("Original Edges:", edges) 
    print("Adjacency Matrix:\n", df)
    print('rows_names =', row_names)
    print('col_names =', col_names) 
    print('-' * 40)
          
    nbi_0 = 0
    nbi_1 = 0
    for _, degree in rows:
        nbi_1 = nbi_1 + degree
    nbi_0 = len(rows)*len(cols) - nbi_1

    model_name = model
    if model == 'König_V':
        model = König_V(rows, cols, edges)
    elif model == 'König_E':
        model = König_E(rows, cols, edges)
    elif model == 'AB_E':
        model = AB_E(rows, cols, edges, epsilon)
    elif model == 'AB_E_r':
        model = AB_E_r(rows, cols, edges, epsilon)
    elif model == 'AB_E_c_r':
        model = AB_E_c_r(rows, cols, edges, epsilon)
    elif model == 'AB_V':
        model = AB_V(rows, cols, edges, epsilon)
    elif model == 'AB_V_h':
        model = AB_V_h(rows, cols, edges, epsilon)
    elif model == 'max_Ones':
        model = max_Ones(rows, cols, edges, epsilon)
    elif model == 'max_Ones_comp':
        model = max_Ones_comp(rows, cols, edges, epsilon)
    elif model == 'max_Surface':
        model = max_Surface(rows, cols, edges, epsilon)
    elif model == 'max_Vertices':
        model = max_Vertices(rows, cols, edges, epsilon)
    elif p.lower().endswith('.txt'):
        if   model == 'minDel_RC':
             model = minDel_RC(rows, cols, CompEdges, epsilon)
        elif model == 'minDel_Ones':
             model = minDel_Ones(rows, cols, CompEdges, epsilon)
        elif model == 'KP_QBr':
             model = KP_QBr(rows, CompEdges, epsilon)
        elif model == 'KP_QBc':
             model = KP_QBc(cols, CompEdges, epsilon)  
    elif p.lower().endswith('.csv'):
        if   model == 'minDel_RC':
             model = minDel_RC(rows, cols, edges, epsilon)
        elif model == 'minDel_Ones':
             model = minDel_Ones(rows, cols, edges, epsilon)
        elif model == 'KP_QBr':
             model = KP_QBr(rows, edges, epsilon)
        elif model == 'KP_QBc':
             model = KP_QBc(cols, edges, epsilon)
                    
    #create the model using one of the previously implemented models

    #solve the model using GUROBI_CMD. it is possible for the solver to take a long time
    #the time limit is set to 1 hour. The solver will be automatically stop after 1h.
    #model.solve(PULP_CBC_CMD(msg=True, timeLimit= 3600, gapRel = 0.5),)
    model.solve(GUROBI_CMD(msg=True, timeLimit= 600)#, options=[("Heuristics", 0.0), ("NoRelHeurTime", 0)] )#,gapRel=0.3)
    )
    #model.solve(GUROBI_CMD(msg=True, timeLimit= 60, MIPGap = 0.03),)
    # Check status
    print("Model is . Exporting LP file for debugging...")
    model.writeLP("debug_model_row.lp")
    print(f"Model status: {LpStatus[model.status]}")
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
    #read the result from the solver
    rows_res = []
    cols_res = []
    rows_res_name = []
    cols_res_name = []
    rows_del = []
    rows_del_name = []
    cols_del = []
    cols_del_name = []
    obj_total =  0.0

    if model_name == 'AB_V'  or model_name == 'AB_V_h'  or model_name == 'max_Surface' or model_name == 'max_Vertices' or model_name == 'max_Ones_comp'  or model_name == 'max_Ones' or model_name == 'AB_E' or model_name == 'AB_E_r' or model_name == 'AB_E_c_r':
        for var in model.variables():
            #if var.varValue == 1:
            if var.varValue !=  0:
                # if var.name[:10] == "obj_colsum":
                #     print(f"Value: {var.varValue}, Name: {var.name}")
                #     obj_total =  obj_total + var.varValue
                # if var.name[:10] == "obj_rowsum":
                #     print(f"Value: {var.varValue}, Name: {var.name}")
                #     obj_total =  obj_total + var.varValue
                # if var.name[:3] == "row":
                #     print(f"Value: {var.varValue}, Name: {var.name}")
                # if var.name[:3] == "col":
                #     print(f"Value: {var.varValue}, Name: {var.name}")
                # if var.name[:3] == "sum":
                #     print(f"Value: {var.varValue}, Name: {var.name}")  
                # if var.name[:5] == "slack":
                #     print(f"Value: {var.varValue}, Name: {var.name}")    
                #     if   var.varValue <= 0:
                #          slack_total = slack_total + var.varValue
                if var.name[:3] == "row":
                    rows_res = rows_res + [var.name[4:]]
                elif var.name[:3] == "col":
                    cols_res = cols_res + [var.name[4:]]
    else:
        print('I solved model name =', model.name)
        if model.name == "row_knapsack_problem":
            for var in model.variables():
                print('row_knapsack var name =', var.name,'var value =', var.varValue)
                if var.varValue == 0:
                    if var.name[:3] == "row":
                        rows_res = rows_res + [var.name[4:]]
                    elif var.name[:3] == "col":
                #    cols_res = cols_res + [var.name[4:]]
                        print('Something wrong. var name =', var.name,'var value =', var.varValue)
                        return model 
                else: # i.e. var.varValue == 1:
                    if var.name[:3] == "row":
                        rows_del = rows_del + [var.name[4:]]
                    elif var.name[:3] == "col":
                        #cols_del = cols_del + [var.name[4:]]
                        print('Something wrong. var name =', var.name,'var value =', var.varValue)
                        return model 
        if model.name == "column_knapsack_problem":
            for var in model.variables():
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
                    

    print_log_output(model)

    if obj_total != 0.0: 
        total_error = -slack_total/obj_total
    
    print()
    print('-' * 40)
    print('obj_total=', obj_total)
    print('slack_total=', slack_total)
    print('total_error=', total_error)
    print()
    print('-' * 40)
    #
    if model.name == "row_knapsack_problem":
         cols_res = [str(c) for c, _ in cols]
         cols_res = [int(c) for c in cols_res]
         cols_res_name = [col_names[c] for c in cols_res]
    if model.name == "column_knapsack_problem":
         rows_res = [str(c) for c, _ in rows]
         rows_res = [int(r) for r in rows_res]
         rows_res_name = [row_names[r] for r in rows_res]
    

    print()
    print('-after KP ****')
   #print('rows=', rows)
    #print('cols=', cols)
    print('len_rows_res=', len(rows_res))
    print('row_res=', rows_res)
    print('len_rows_del=', len(rows_del))
    print('rows_del=',rows_del)
    print('len_cols_res=', len(cols_res))
    #print('cols_res=', cols_res)
    print('len_cols_del=', len(cols_del))
    #print('cols_del=',cols_del)
    print()
    print('-' * 40)
    

    if not rows_res:  # Equivalent to checking len(rows_res) == 0
       print("matrix degenerated (all rows or all columns have been deleted)")
       return rows_res, cols_res
    #rows_res = [str(r) for r, _ in rows]
    # if len(rows_res) == 0:
    #     row_res = [str(r) for r, _ in rows]

    # print()
    # print('-after KP ****')
    # print('cols_res')
    # print(cols_res)
    # print()
    # print('-' * 40)
    
    # #
    #rows_res = [int(r) for r in rows_res]

    # rows_res = [r.varValue for r in rows_res if r.varValue is not None]
    # cols_res = [c.varValue for c in cols_res if c.varValue is not None]
    # rows_res = [ r.varValue for r in rows_res if isinstance(r, LpVariable) and r.varValue is not None]
    # cols_res = [ c.varValue for c in cols_res if isinstance(c, LpVariable) and c.varValue is not None]

    print(f"input data = {path_to_data}")
    print('Initial Stats')
    print("Size of initial matrix : ", len(rows), "*", len(cols), "=", len(rows) * len(cols))
    
    nbi_0 = 0
    nbi_1 = 0
    for _, degree in rows:
        nbi_1 = nbi_1 + degree
    
    nbi_0 = len(rows)*len(cols) - nbi_1

    print("number zero initial = ",nbi_0)     
    print("number one initial = ",nbi_1)            
    print("epsilon = ",epsilon)
    print("Initial sparsity = ", nbi_0/(len(rows) * len(cols)))
    print("Initial density = ", nbi_1/(len(rows) * len(cols)))
    print('-' * 40)
    print()

    rows_res = [int(r) for r in rows_res]
    cols_res = [int(c) for c in cols_res]
    rows_res_name = [row_names[r] for r in rows_res]
    cols_res_name = [col_names[c] for c in cols_res]

    print()
    print('-' * 40)
    print("solution as remaining row and column names =  ")

    print('row_res_name')
    print(rows_res_name)
    print('cols_res_name')
    print(cols_res_name)
    print("solution as row and column indices = ")
    print('row_res_indices')
    print(rows_res)
    print('cols_res_indices')
    print(cols_res)

 # Example usage suggested by AI :
    print("******Example usage suggested by AI*********")
    #
    print('cols=')
    print(cols)
    print('rows=')
    print(rows)
    print('rows_del=')
    print(rows_del)
    print('edges for KP=')
    print(edges)


    cols_data_new, edges_new = update_data(rows, cols, edges, rows_del)

    print("Updated Columns Data:", cols_data_new)
    print("Updated  Edges:", edges_new)

    model = KP_QBc(cols_data_new, edges_new, epsilon)
    model.solve(GUROBI_CMD(msg=True, timeLimit= 600))
    # Check status
    print("Model is . Exporting LP file for debugging...")
    model.writeLP("debug_model_col.lp")
    print(f"Model status: {LpStatus[model.status]}")
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
    #read the result from the solver
    if model.name == "column_knapsack_problem":
            for var in model.variables():
                print('column_knapsack var name =', var.name,'var value =', var.varValue)
                # if var.varValue == 0:
                #     if var.name[:3] == "col":
                #         cols_res = cols_res + [var.name[4:]]
                #     elif var.name[:3] == "row":
                # #    cols_res = cols_res + [var.name[4:]]
                #         print('Something wrong. var name =', var.name,'var value =', var.varValue)
                #         return model 
                # else: # i.e. var.varValue == 1:
                #     if var.name[:3] == "col":
                #         cols_del = cols_del + [var.name[4:]]
                #     elif var.name[:3] == "row":
                #         #cols_del = cols_del + [var.name[4:]]
                #         print('Something wrong. var name =', var.name,'var value =', var.varValue)
                #         return model 
                    


    nb_0 = len(rows_res)*len(cols_res)-df.iloc[rows_res,cols_res].sum().sum()
    nb_1 =  df.iloc[rows_res,cols_res].sum().sum()
    
    # if len(rows_res) == 0:
    #     nb_1 = df.iloc[cols_res].values.sum()
    #     nb_0 = len(cols_res)* df.shape[1] - nb_1
    # else:
    #     if len(cols_res) == 0: 
    #         nb_1 = df.iloc[rows_res].values.sum()
    #         nb_0 = len(rows_res) * df.shape[1] - nb_1
    #     else: 
    #         nb_0 = len(rows_res)*len(cols_res)-df.iloc[rows_res,cols_res].sum().sum()
    #         nb_1 =  nbi_1 - df.iloc[rows_res,cols_res].sum().sum()
    #len(rows_res)*len(cols_res)-nb_0 
    if model_name == 'AB_V'  or model_name == 'AB_V_h'  or model_name == 'max_Surface' or model_name == 'max_Vertices' or model_name == 'max_Ones_comp'  or model_name == 'max_Ones' or model_name == 'AB_E' or model_name == 'AB_E_r' or model_name == 'AB_E_c_r': 

        print("Final Perimeter : ", len(rows_res), "+", len(cols_res), "=", len(rows_res) + len(cols_res))
        print("number zero after solving = ",nb_0) 
        print("number ones after solving = ",nb_1) 
        print("epsilon = ",epsilon)
        print(f"model = {model_name}")
        if len(rows_res)== 0 or len(cols_res)==0 : 
            print("final matrix degenerated (all rows or all columns have been deleted)")
        else : 
            print("Final sparsity", nb_0/(len(rows_res) * len(cols_res)))
            print("Final density", 1-(nb_0/(len(rows_res) * len(cols_res))))
    else: 
        print("Size of found matrix  = ", len(rows_res), "*",  len(cols_res),  " = ", len(rows_res)* len(cols_res))      
        print("nbi_0, nb_0 = ", nbi_0, nb_0) 
        print("nbi_1, nb_1 = ", nbi_1, nb_1) 
       # print("nbi_1, nb_1, df.iloc  = ",nbi_1,  nb_1, df.iloc[rows_res,cols_res].sum().sum() )
        print("number zero after solving = ", nb_0) 
        print("number ones after solving = ", nb_1)
        if nb_0 + nb_1 !=0:
            print("Final sparsity",  nb_0/(nb_0 + nb_1))
            print("Final density",  nb_1/(nb_0 + nb_1))
        print("epsilon = ",epsilon)
        print(f"model = {model_name}")
        if len(rows_res)== 0 or len(cols_res)==0 : 
            print("final matrix degenerated (all rows or all columns have been deleted)")


    return rows_res, cols_res

def update_data(rows_data, cols_data, edges, rows_del):
    # Convert rows_del to a set for faster lookup
    rows_del_set = set(map(int, rows_del))
    print("\n-- Debugging Step: rows in rows_del_set --")
    print("rows_del_set =", rows_del_set)

    # Step 1: Compute edges_row (edges connected to rows in rows_del)
    edges_row_del= [edge for edge in edges if edge[0] in rows_del_set]

    print("\n-- Debugging Step: Edges to be removed --")
    print("Edges_row_deleted =", edges_row_del)

    # Step 2: Compute edges_new (remaining edges)
    edges_new = [edge for edge in edges if edge[0] not in rows_del_set]

    print("\n-- Debugging Step: Remaining Edges --")
    print("Remainig Edges_new =", edges_new)

 # Step 3: Update column degrees in cols_data
    print("cols_data =", cols_data)
    col_degree_map = {col: degree for col, degree in cols_data}
    print("OLD col_degree_map =", col_degree_map)

    # Reduce the count of columns based on removed edges (only affected columns)
    for _, col in edges_row_del :  # Loop only through removed edges
        col_degree_map[col] = max(0, col_degree_map[col] - 1)  # Reduce degree only if that column was involved

    print("\n-- Debugging Step: Updated Column Degrees --")
    print("Updated col_degree_map =", col_degree_map)

    # Rebuild cols_data with updated degrees
    cols_data_new = [(col, degree) for col, degree in col_degree_map.items()]
    print("cols_data_new =", cols_data_new)


    # # Step 3: Update column degrees in cols_data
    # print("cols_data =", cols_data)
    # col_degree_map = {col: degree for col, degree in cols_data}
    # print("cols_degremap =", col_degree_map)

    # # Reduce the count of columns based on removed edges
    # for r  in rows_del_set:
    #     for col in col_degree_map:
    #         col_degree_map[col] = max(0, col_degree_map[col] - 1)  # Ensure non-negative

    return cols_data_new, edges_new

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

# def get_data_txt_file(path, model=None):
#     with open(path, 'r') as file:
#         content = file.readlines()
    
#     # Extract metadata
#     num_row = int(content[0].split(":")[1].strip())  # Number of rows (|U|)
#     num_col = int(content[1].split(":")[1].strip())  # Number of columns (|V|)
#     num_edge = int(content[2].split(":")[1].strip())  # Number of edges (|E|)

#     # Initialize row and column degrees
#     deg_row = [0] * num_row
#     deg_col = [0] * num_col
#     edges = []

#     # Initialize a DataFrame for the adjacency matrix
#     df = pd.DataFrame(0, index=range(num_row), columns=range(num_col))

#     # Parse the edges
#     for line in content[3:]:
#         # Split the line by whitespace or tabs
#         splitted_line = re.split(r'\s+', line.strip())
#         if len(splitted_line) < 2:
#             continue  # Skip empty or invalid lines
#         u, v = int(splitted_line[0]), int(splitted_line[1])
#         edges.append((u, v))
#         deg_row[u] += 1
#         deg_col[v] += 1
#         df.iloc[u, v] = 1

#     # Prepare row and column data
#     rows_data = list(zip(range(num_row), deg_row))
#     cols_data = list(zip(range(num_col), deg_col))

#     return rows_data, cols_data, edges, range(num_row), range(num_col), df

def get_complement_edges(num_row, num_col, edges):
    # Create a set of all possible edges using Cartesian product
    all_edges = set(itertools.product(range(num_row), range(num_col)))

    # Convert the original edges list to a set for efficient subtraction
    original_edges = set(edges)

    # Complement edges = all_edges - original_edges
    complement_edges = list(all_edges - original_edges)

    return complement_edges


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

def print_log_output(prob):
    """Print the log output and problem solutions.
    ARGUMENTS:
    ----------
    * prob: an solved LP model (pulp.LpProblem)
    """
    print()
    print('-' * 40)
    print('Stats')
    print('-' * 40)
    print()
    print(f'Number variables: {prob.numVariables()}')
    print(f'Number constraints: {prob.numConstraints()}')
    print()
    print('Time:')
    print(f'- (real) {prob.solutionTime}')
    print(f'- (CPU) {prob.solutionCpuTime}')
    print()

    print(f'Solve status: {LpStatus[prob.status]}')
    print(f'Objective value: {prob.objective.value()}')

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


def parse_arguments():
    """Parse the input arguments and retrieve the choosen resolution method and
    the instance that must be solve."""
    argparser = ArgumentParser()

    argparser.add_argument(
        '--filepath', dest='filepath', required=True, default='',
        help='Select the data',
    )
    # argparser.add_argument(
    #     '--model', dest='model', required=False, default='AB_V',
    #     help='Select the model to use',
    # )
    argparser.add_argument(
        '--model', dest='model', required=False, default='AB_V_h',
        help='Select the model to use',
    )

    argparser.add_argument(
        '--epsilon', dest='epsilon', required=False, default=0.1, type=float,
        help='Select the error rate value',
    )

    arg = argparser.parse_args()

    if arg.model not in ['König_V', 'König_E', 'AB_E', 'AB_E_r', 'AB_E_c_r','AB_V','AB_V_h','max_Ones','max_Ones_comp','max_Surface','max_Vertices','minDel_RC', 'minDel_Ones', 'KP_QBr', 'KP_QBc']:
        argparser.print_help()
        sys.exit(1)

    return (arg.filepath, arg.model, arg.epsilon)


if __name__ == '__main__':

    # Read the arguments
    file_path, selected_model, epsilon = parse_arguments()

    solve(file_path,selected_model,epsilon)



