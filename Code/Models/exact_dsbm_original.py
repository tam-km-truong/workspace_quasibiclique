"""ILP for the exact densest sub-binary matrix (DSBM) problem."""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

from pulp import LpInteger, LpMaximize, LpProblem, LpVariable, lpSum

from qbc_dsbm import QBCExactModels

if TYPE_CHECKING:
    from qbc_dsbm.binary_matrix import BinMatrix

# ============================================================================ #
#                                     TYPES                                    #
# ============================================================================ #
RowColumnChoicesLP = list[LpVariable]
CellChoicesLP = list[list[LpVariable]]


# ============================================================================ #
#                     LP MODEL - MAXIMIZE DESIRED CELL                         #
# ============================================================================ #
def max_ones(
    bin_matrix: BinMatrix,
    epsilon: float,
) -> tuple[LpProblem, RowColumnChoicesLP, RowColumnChoicesLP, CellChoicesLP]:
    """Maximize the number of ones in the matrix.

    Parameters
    ----------
    bin_matrix : BinMatrix
        Binary matrix
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        ILP model
    RowColumnChoicesLP
        Row choices
    RowColumnChoicesLP
        Column choices
    CellChoicesLP
        Cell choices

    """
    prob = LpProblem(name=QBCExactModels.MAX_ONES.value, sense=LpMaximize)

    #
    # Variables
    #
    (
        row_choices,
        column_choices,
        cell_choices,
    ) = __row_column_cell_choices_variables(bin_matrix)

    #
    # Objective function
    #
    prob += (
        lpSum(
            bin_matrix[u, v] * cell_choices[u][v]
            for u in range(bin_matrix.number_of_rows())
            for v in range(bin_matrix.number_of_columns())
        ),
        "maximize_ones_in_matrix",
    )

    #
    # Constraints
    #
    __select_cell_select_its_row_column(
        bin_matrix,
        prob,
        row_choices,
        column_choices,
        cell_choices,
    )
    __select_no_more_epsilon_zeros(prob, bin_matrix, epsilon, cell_choices)

    return prob, row_choices, column_choices, cell_choices


# TODO add these options to the other models
def max_ones_compact(
    bin_matrix: BinMatrix,
    epsilon: float,
    min_number_of_rows: int = 0,
    min_number_of_columns: int = 0,
    warm_solution: BinMatrix | None = None,
    check_satisfiability: bool = True,  # noqa: FBT001, FBT002 # TODO use this option to the parent functions
) -> tuple[LpProblem, RowColumnChoicesLP, RowColumnChoicesLP, CellChoicesLP]:
    """Maximize the number of ones in the matrix (compact version).

    Parameters
    ----------
    bin_matrix : BinMatrix
        Binary matrix
    epsilon : float
        Error rate
    min_number_of_rows : int
        Minimum number of rows
    min_number_of_columns : int
        Minimum number of columns
    warm_solution : BinMatrix
        Warm solution
    check_satisfiability : bool
        Check satisfiability, do not solve an optimization problem

    Returns
    -------
    LpProblem
        ILP model
    RowColumnChoicesLP
        Row choices
    RowColumnChoicesLP
        Column choices
    CellChoicesLP
        Cell choices

    """
    prob = LpProblem(name=QBCExactModels.MAX_ONES_COMPACT.value, sense=LpMaximize)

    #
    # Variables
    #
    (
        row_choices,
        column_choices,
        cell_choices,
    ) = __row_column_cell_choices_variables(bin_matrix)

    #
    # Objective function
    #
    if check_satisfiability:
        prob += 0, "dummy_objective"
    else:
        prob += (
            lpSum(
                bin_matrix[u, v] * cell_choices[u][v]
                for u in range(bin_matrix.number_of_rows())
                for v in range(bin_matrix.number_of_columns())
            ),
            "maximize_ones_in_matrix",
        )

    #
    # Constraints
    #
    __select_cell_select_its_row_column_compact(
        bin_matrix,
        prob,
        row_choices,
        column_choices,
        cell_choices,
    )
    __select_no_more_epsilon_zeros(prob, bin_matrix, epsilon, cell_choices)
    __respect_the_minimum_dimensions(
        prob,
        min_number_of_rows,
        min_number_of_columns,
        row_choices,
        column_choices,
    )

    if warm_solution is not None:
        __set_initial_values_from_warm_solution(
            bin_matrix,
            warm_solution,
            row_choices,
            column_choices,
            cell_choices,
        )

    return prob, row_choices, column_choices, cell_choices


def max_surface(
    bin_matrix: BinMatrix,
    epsilon: float,
) -> tuple[LpProblem, RowColumnChoicesLP, RowColumnChoicesLP, CellChoicesLP]:
    """Maximize the surface of the matrix.

    Parameters
    ----------
    bin_matrix : BinMatrix
        Binary matrix
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        ILP model
    RowColumnChoicesLP
        Row choices
    RowColumnChoicesLP
        Column choices
    CellChoicesLP
        Cell choices

    """
    prob = LpProblem(name=QBCExactModels.MAX_SURFACE.value, sense=LpMaximize)
    #
    # Variables
    #
    (
        row_choices,
        column_choices,
        cell_choices,
    ) = __row_column_cell_choices_variables(bin_matrix)
    #
    # Objective function
    #
    prob += (
        lpSum(
            cell_choices[u][v]
            for u in range(bin_matrix.number_of_rows())
            for v in range(bin_matrix.number_of_columns())
        ),
        "maximize_surface_matrix",
    )
    #
    # Constraints
    #
    __select_cell_select_its_row_column(
        bin_matrix,
        prob,
        row_choices,
        column_choices,
        cell_choices,
    )
    __select_no_more_epsilon_zeros(prob, bin_matrix, epsilon, cell_choices)

    return prob, row_choices, column_choices, cell_choices


def max_perimeter(
    bin_matrix: BinMatrix,
    epsilon: float,
) -> tuple[LpProblem, RowColumnChoicesLP, RowColumnChoicesLP, CellChoicesLP]:
    """Maximize the perimeter of the matrix.

    Parameters
    ----------
    bin_matrix : BinMatrix
        Binary matrix
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        ILP model
    RowColumnChoicesLP
        Row choices
    RowColumnChoicesLP
        Column choices
    CellChoicesLP
        Cell choices

    """
    prob = LpProblem(name=QBCExactModels.MAX_PERIMETER.value, sense=LpMaximize)
    #
    # Variables
    #
    (
        row_choices,
        column_choices,
        cell_choices,
    ) = __row_column_cell_choices_variables(bin_matrix)
    #
    # Objective function
    #
    prob += lpSum(row_choices) + lpSum(column_choices), "max_sum_vertices"
    #
    # Constraints
    #
    __select_cell_select_its_row_column(
        bin_matrix,
        prob,
        row_choices,
        column_choices,
        cell_choices,
    )
    __select_no_more_epsilon_zeros(prob, bin_matrix, epsilon, cell_choices)

    return prob, row_choices, column_choices, cell_choices


# ============================================================================ #
#                                   VARIABLES                                  #
# ============================================================================ #
def __row_column_cell_choices_variables(
    bin_matrix: BinMatrix,
) -> tuple[RowColumnChoicesLP, RowColumnChoicesLP, CellChoicesLP]:
    row_choices = [
        LpVariable(f"xu_{u}", cat=LpInteger, lowBound=0, upBound=1)
        for u in range(bin_matrix.number_of_rows())
    ]
    column_choices = [
        LpVariable(f"xv_{v}", cat=LpInteger, lowBound=0, upBound=1)
        for v in range(bin_matrix.number_of_columns())
    ]
    cell_choices = [
        [
            LpVariable(f"cell_{u}_{v}", cat=LpInteger, lowBound=0, upBound=1)
            for v in range(bin_matrix.number_of_columns())
        ]
        for u in range(bin_matrix.number_of_rows())
    ]

    return row_choices, column_choices, cell_choices


# ============================================================================ #
#                                  CONSTRAINTS                                 #
# ============================================================================ #
def __select_cell_select_its_row_column(
    bin_matrix: BinMatrix,
    prob: LpProblem,
    row_choices: RowColumnChoicesLP,
    column_choices: RowColumnChoicesLP,
    cell_choices: CellChoicesLP,
) -> None:
    # OPTIMIZE forall or sum? Perhaps have a look at robust ILP
    for u in range(bin_matrix.number_of_rows()):
        for v in range(bin_matrix.number_of_columns()):
            prob += (row_choices[u] >= cell_choices[u][v]), f"cell_{u}_{v}_1"

            prob += (column_choices[v] >= cell_choices[u][v]), f"cell_{u}_{v}_2"

            if not bin_matrix[u, v]:
                prob += (
                    (row_choices[u] + column_choices[v] - 1 <= cell_choices[u][v]),
                    f"cell_{u}_{v}_3",
                )


def __select_cell_select_its_row_column_compact(
    bin_matrix: BinMatrix,
    prob: LpProblem,
    row_choices: RowColumnChoicesLP,
    column_choices: RowColumnChoicesLP,
    cell_choices: CellChoicesLP,
) -> None:
    for u in range(bin_matrix.number_of_rows()):
        prob += (
            (
                sum(
                    1 for v in range(bin_matrix.number_of_columns()) if bin_matrix[u, v]
                )
                * row_choices[u]
                >= lpSum(
                    cell_choices[u][v]
                    for v in range(bin_matrix.number_of_columns())
                    if bin_matrix[u, v]
                )
            ),
            f"row_{u}_1",
        )
    for v in range(bin_matrix.number_of_columns()):
        prob += (
            (
                sum(1 for u in range(bin_matrix.number_of_rows()) if bin_matrix[u, v])
                * column_choices[v]
                >= lpSum(
                    cell_choices[u][v]
                    for u in range(bin_matrix.number_of_rows())
                    if bin_matrix[u, v]
                )
            ),
            f"col_{v}_2",
        )
    # OPTIMIZE forall or sum? Perhaps have a look at robust ILP
    for u in range(bin_matrix.number_of_rows()):
        for v in range(bin_matrix.number_of_columns()):
            if not bin_matrix[u, v]:
                prob += (
                    (row_choices[u] + column_choices[v] - 1 <= cell_choices[u][v]),
                    f"cell_{u}_{v}_3",
                )


def __select_no_more_epsilon_zeros(
    prob: LpProblem,
    bin_matrix: BinMatrix,
    epsilon: float,
    cell_choices: CellChoicesLP,
) -> None:
    prob += (
        (
            lpSum(
                (1 - bin_matrix[u, v]) * cell_choices[u][v]
                for u in range(bin_matrix.number_of_rows())
                for v in range(bin_matrix.number_of_columns())
            )
            <= epsilon
            * lpSum(
                cell_choices[u][v]
                for u in range(bin_matrix.number_of_rows())
                for v in range(bin_matrix.number_of_columns())
            )
        ),
        "err_rate",
    )


def __respect_the_minimum_dimensions(
    prob: LpProblem,
    min_number_of_rows: int,
    min_number_of_columns: int,
    row_choices: RowColumnChoicesLP,
    column_choices: RowColumnChoicesLP,
) -> None:
    prob += (
        lpSum(row_choices) >= min_number_of_rows,
        f"The result must have at least {min_number_of_rows} rows",
    )
    prob += (
        lpSum(column_choices) >= min_number_of_columns,
        f"The result must have at least {min_number_of_columns} columns",
    )


def __set_initial_values_from_warm_solution(
    bin_matrix: BinMatrix,
    warm_solution: BinMatrix,
    row_choices: RowColumnChoicesLP,
    column_choices: RowColumnChoicesLP,
    cell_choices: CellChoicesLP,
) -> None:
    row_indices: list[int] = []
    column_indices: list[int] = []
    rows_in_warm_solution = set(warm_solution.row_header())
    for row_idx, row_id in enumerate(bin_matrix.row_header()):
        if row_id in rows_in_warm_solution:
            row_choices[row_idx].setInitialValue(1)
            row_indices.append(row_idx)
        else:
            row_choices[row_idx].setInitialValue(0)

    columns_in_warm_solution = set(warm_solution.column_header())
    for column_idx, column_id in enumerate(bin_matrix.column_header()):
        if column_id in columns_in_warm_solution:
            column_choices[column_idx].setInitialValue(1)
            column_indices.append(column_idx)
        else:
            column_choices[column_idx].setInitialValue(0)

    for i in range(bin_matrix.number_of_rows()):
        for j in range(bin_matrix.number_of_columns()):
            cell_choices[i][j].setInitialValue(0)
    for row_idx, column_idx in itertools.product(row_indices, column_indices):
        cell_choices[row_idx][column_idx].setInitialValue(1)
