# -*- coding=utf-8 -*-

"""
Skeleton for TP project: Searching Maximum Quasi-Bicliques
Student:

Install Gurobi Solver with free academic licence using your University account at: https://www.gurobi.com/features/academic-named-user-license/

Install Pandas, using PyPI or Conda. A detailed installation tutorial can be find at: https://pandas.pydata.org/docs/getting_started/install.html
"""

import sys
from argparse import ArgumentParser

import pandas as pd
from pulp import (
    GUROBI_CMD,
    PULP_CBC_CMD,
    LpBinary,
    LpMaximize,
    LpMinimize,
    LpProblem,
    LpStatus,
    LpVariable,
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
    # lpCells = {}
    # for row, _ in rows_data:
    #     for col, _ in cols_data:
    #         if (row, col) in edges:
    #             lpCells[(row, col)] = (LpVariable(
    #                 f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 1)
    #         else:
    #             lpCells[(row, col)] = (LpVariable(  
    #                 f'cell_{row}_{col}', cat='Integer', lowBound=0, upBound=1), 0)
                
    # Objective function: Maximize sum of selected row and column variables
    model += lpSum(
        [lpvar for lpvar, _ in lpRows.values()] +
        [lpvar for lpvar, _ in lpCols.values()]
    ), "max_sum_vertices"

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
    #just for testing 
    # for row, col in lpCells:
    #     #model += (lpRows[row][0] >= lpCells[(row, col)][0]), f'cell_{row}_{col}_1'
    #     #model += (lpCols[col][0] >= lpCells[(row, col)][0]), f'ce ll_{row}_{col}_2'
    #     model += ((lpRows[row][0]+lpCols[col][0])/2 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'
    #     model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_4'
    #just for testing     
    # Add row density constraints
    __row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    #__row_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    #__col_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)

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
        [lpvar for lpvar, _ in lpCols.values()]
    ), "max_sum_vertices"

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
    edges: list of tuples (row, col) corresponding to the zeros of the matrix.
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
        model += ((lpRows[row][0]+lpCols[col][0])/2 >= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'
        model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_4'
        #########################################
        #compacting  with degree 
        #########################################
        '''        
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
          '''
    # Add row density constraints
    __row_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __col_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    #__row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    #__col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)

    return model 

def AB_E_r(rows_data, cols_data, edges, epsilon):
    """
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
    model = LpProblem(name="AB_E_r", sense=LpMaximize)

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
        model += ((lpRows[row][0]+lpCols[col][0])/2 >= lpCells[(row, col)][0]), f'cell_{row}_{col}_3'
        model += (lpRows[row][0]+lpCols[col][0] -1 <= lpCells[(row, col)][0]), f'cell_{row}_{col}_4'
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
        
    #      #########################################


    # Add row density constraints
    #__row_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    #__col_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)

    return model 

def AB_E_c_r(rows_data, cols_data, edges, epsilon):
    """
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
    #__row_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    #__col_density(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)
    __col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon)

    return model 


def __col_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon):
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
        
        model += (
            lpSum(lpRows[row][0] for row in col_edges) - (1 - epsilon) * lpSum(lpRows[row][0] for row, _ in rows_data) <=
            -mu + lpCols[col][0] * Big_M
        ), f"col_err_rate_0_{col}"

def __row_density_iff(rows_data, cols_data, edges, model, lpRows, lpCols, epsilon):
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
        
        model += (
            lpSum(lpCols[col][0] for col in row_edges) - (1 - epsilon) * lpSum(lpCols[col][0] for col, _ in cols_data) <=
            -mu + lpRows[row][0] * Big_M
        ), f"row_err_rate_0_{row}"

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
    model += (lpSum([lpvar for lpvar, _ in lpRows.values()])  <= 2,'row_limits_bis',    )
    #model += (lpSum([lpvar for lpvar, _ in lpCols.values()])  <= 2,'col_limits_bis',    )
    
    
    return model


# ============================================================================ #
#                LP MODEL - Deleting rows/columns for zeros eliminating                    #
# ============================================================================ #


def minDel_RC(rows_data, cols_data, edges, epsilon=0.3):
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


def minDel_Ones(rows_data, cols_data, edges, epsilon=0.3):
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
    model += lpSum([degree*lpvar for lpvar, degree in lpRows.values()] +
                   [degree*lpvar for lpvar, degree in lpCols.values()]), 'min_weighted_rows/cols'
    # model += lpSum([lpvar for lpvar, _ in lpRows.values()] +
    #                [lpvar for lpvar, _ in lpCols.values()]), 'min_weighted_rows/cols'

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
    if nbi0>nbi1:
        diff = nbi0-nbi1
        model += (lpSum(lpEdges) >= diff+(1-epsilon)*nbi1), f'sensitivity'
    else:
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

def KP_QB(rows_data, cols_data, edges, epsilon=0.1):
    """
    In this model, we use knapsack model
    ARGUMENTS:
    ----------
    * rows_data: list of the tuples (row, degree) of rows the matrix.
    * cols_data: list of the tuples (col, degree) of columns the matrix.
    * epsilon: percentage of accepted undesired cell over accepted desired cell
    """
    # ------------------------------------------------------------------------ #
    # Model with minimization
    # ------------------------------------------------------------------------ #
    model = LpProblem(name='knapsack_problem', sense=LpMinimize)

    # ------------------------------------------------------------------------ #
    # Variables
    # ------------------------------------------------------------------------ #

    lpRows = [(LpVariable(f'row_{row}', cat='Integer',
                    lowBound=0, upBound=1), degree) for row, degree in rows_data]
    lpCols = [(LpVariable(f'col_{col}', cat='Integer',
                    lowBound=0, upBound=1), degree) for col, degree in cols_data]

    # ------------------------------------------------------------------------ #
    # Objective function
    # ------------------------------------------------------------------------ #
    model += lpSum([degree*lpvar for lpvar, degree in lpRows] +
                   [degree*lpvar for lpvar, degree in lpCols]), 'knapsack'

    # ------------------------------------------------------------------------ #
    # Constraints
    # ------------------------------------------------------------------------ #
    nb_edges = len(edges)
    row_length = len(rows_data)
    col_length = len(cols_data)
    model += (lpSum([(row_length-degree)*lpvar for lpvar, degree in lpRows]) + 
              lpSum([(col_length-degree)*lpvar for lpvar, degree in lpCols])  >= (1-epsilon) * nb_edges), f'sensitivity'
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
    

    rows, cols, edges, row_names, col_names, df = get_data(path_to_data)


    # print()
    # print('-' * 40)
    # print('Initial Stats')
    # print("Perimeter of initial matrix : ", len(rows), "+", len(cols), "=", len(rows) + len(cols))
    
    nbi_0 = 0
    nbi_1 = 0
    for _, degree in rows:
        nbi_1 = nbi_1 + degree
    
    nbi_0 = len(rows)*len(cols) - nbi_1

    # print("number zero initial = ",nbi_0)     
    # print("number one initial = ",nbi_1)            
    # print("epsilon = ",epsilon)
    # print("Initial sparsity = ", nbi_0/(len(rows) * len(cols)))
    # print("Initial density = ", 1-(nbi_0/(len(rows) * len(cols))))
    # print('-' * 40)
    # print()


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
    elif model == 'minDel_RC':
        model = minDel_RC(rows, cols, edges, epsilon)
    elif model == 'minDel_Ones':
        model = minDel_Ones(rows, cols, edges, epsilon)
    elif model == 'KP_QB':
        model = KP_QB(rows, cols, edges, epsilon)
    #create the model using one of the previously implemented models

    #solve the model using GUROBI_CMD. it is possible for the solver to take a long time
    #the time limit is set to 1 hour. The solver will be automatically stop after 1h.
    #model.solve(PULP_CBC_CMD(msg=True, timeLimit= 3600, gapRel = 0.5),)
    model.solve(GUROBI_CMD(msg=True, timeLimit= 600)#, options=[("Heuristics", 0.0), ("NoRelHeurTime", 0)] )#,gapRel=0.3)
    )
    #model.solve(GUROBI_CMD(msg=True, timeLimit= 60, MIPGap = 0.03),)
# Check status
    #print("Model is . Exporting LP file for debugging...")
    #model.writeLP("debug_model.lp")
    print(f"Model status: {LpStatus[model.status]}")
    if model.status == -1:
         print("Model is infeasible. Exporting LP file for debugging...")
         model.writeLP("debug_model.lp")
    #read the result from the solver
    rows_res = []
    cols_res = []
    if model_name == 'AB_V'  or model_name == 'AB_V_h'  or model_name == 'max_Surface' or model_name == 'max_Vertices' or model_name == 'max_Ones_comp'  or model_name == 'max_Ones' or model_name == 'AB_E' or model_name == 'AB_E_r' or model_name == 'AB_E_c_r':
        for var in model.variables():
            if var.varValue == 1:
                if var.name[:3] == "row":
                    rows_res = rows_res + [var.name[4:]]
                elif var.name[:3] == "col":
                    cols_res = cols_res + [var.name[4:]]
    else:
        for var in model.variables():
            if var.varValue == 0:
                if var.name[:3] == "row":
                    rows_res = rows_res + [var.name[4:]]
                elif var.name[:3] == "col":
                    cols_res = cols_res + [var.name[4:]]

    print_log_output(model)

    rows_res = [int(r) for r in rows_res]
    cols_res = [int(c) for c in cols_res]
    #print('row_res=', rows_res)
    #print('cols_res=', cols_res)
    print()
    print('-' * 40)
    print(f"input data = {path_to_data}")
    print('Initial Stats')
    print("Perimeter of initial matrix : ", len(rows), "+", len(cols), "=", len(rows) + len(cols))
    
    nbi_0 = 0
    nbi_1 = 0
    for _, degree in rows:
        nbi_1 = nbi_1 + degree
    
    nbi_0 = len(rows)*len(cols) - nbi_1

    print("number zero initial = ",nbi_0)     
    print("number one initial = ",nbi_1)            
    print("epsilon = ",epsilon)
    print("Initial sparsity = ", nbi_0/(len(rows) * len(cols)))
    print("Initial density = ", 1-(nbi_0/(len(rows) * len(cols))))
    print('-' * 40)
    print()
    
    print("Final Perimeter : ", len(rows_res), "+", len(cols_res), "=", len(rows_res) + len(cols_res))

    nb_0 = len(rows_res)*len(cols_res)-df.iloc[rows_res,cols_res].sum().sum()
    nb_1 =  len(rows_res)*len(cols_res)-nb_0 

    print("number zero after solving = ",nb_0) 
    print("number ones after solving = ",nb_1) 
    print("epsilon = ",epsilon)
    print(f"model = {model_name}")
    if len(rows_res)== 0 or len(cols_res)==0 : 
        print("final matrix degenerated (all rows or all columns have been deleted)")
    else : 
        print("Final sparsity", nb_0/(len(rows_res) * len(cols_res)))
        print("Final density", 1-(nb_0/(len(rows_res) * len(cols_res))))

    rows_res = [row_names[r] for r in rows_res]
    cols_res = [col_names[c] for c in cols_res]
    print()
    print('-' * 40)
    print("solution as rows and columns = ")

    print('row_res')
    print(rows_res)
    print('cols_res')
    print(cols_res)
    #print(rows_res,cols_res)
    return rows_res, cols_res

  
def get_data(path:str):

    rows_data = []
    cols_data = []
    edges = []

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
    edges = list(df[df == 1].stack().index)
    #edges = list(df[df == 0].stack().index)
    '''
    print('-' * 40)
    print('edges =', edges)
    print('rows_data =', rows_data)
    print('cols_data =', cols_data)
    print()
    print('-' * 40)
    '''
    return rows_data, cols_data, edges, row_names, col_names, df

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
    #breakpoint()
    #print('-' * 40)


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

    if arg.model not in ['König_V', 'König_E', 'AB_E', 'AB_E_r', 'AB_E_c_r','AB_V','AB_V_h','max_Ones','max_Ones_comp','max_Surface','max_Vertices','minDel_RC', 'minDel_Ones', 'KP_QB']:
        argparser.print_help()
        sys.exit(1)

    return (arg.filepath, arg.model, arg.epsilon)


if __name__ == '__main__':

    # Read the arguments
    file_path, selected_model, epsilon = parse_arguments()

    solve(file_path,selected_model,epsilon)



