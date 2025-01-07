"""ILP models module."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pulp import LpContinuous, LpInteger, LpMinimize, LpProblem, LpVariable, lpSum

from qbc_dsbm import QBCHeuristicModels

if TYPE_CHECKING:
    from qbc_dsbm.bigraph import ZerosBiGraph

# ============================================================================ #
#                                     TYPES                                    #
# ============================================================================ #
VertexChoicesLP = list[LpVariable]
EdgeChoicesLP = dict[tuple[int, int], LpVariable]


# ============================================================================ #
#            LP MODEL - DELETING ROWS/COLUMNS FOR ZEROS ELIMINATING            #
# ============================================================================ #
def min_del_rc(
    bigraph: ZerosBiGraph,
    epsilon: float,
) -> tuple[LpProblem, VertexChoicesLP, VertexChoicesLP, EdgeChoicesLP]:
    """Minimizing the number of deleted rows/columns.

    Parameters
    ----------
    bigraph : BiGraph
        Bigraph
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        PuLP problem
    VertexChoicesLP
        U vertex choices
    VertexChoicesLP
        V vertex choices
    EdgeChoicesLP
        Edge choices

    """
    prob = LpProblem(name=QBCHeuristicModels.MIN_DEL_RC.value, sense=LpMinimize)
    #
    # Variables
    #
    u_vertex_choices = __u_vertex_choices_variables(bigraph)
    v_vertex_choices = __v_vertex_choices_variables(bigraph)
    edge_choices = __edge_choices_variables(bigraph)
    #
    # Objective function
    #
    prob += lpSum(u_vertex_choices) + lpSum(v_vertex_choices), "min_vertices"
    #
    # Constraints
    #
    __select_edge_select_at_least_one_endpoint(
        prob,
        bigraph,
        u_vertex_choices,
        v_vertex_choices,
        edge_choices,
    )
    __at_most_epsilon_remaining_zeros(prob, bigraph, epsilon, edge_choices)

    return prob, u_vertex_choices, v_vertex_choices, edge_choices


def min_del_rows(
    bigraph: ZerosBiGraph,
    epsilon: float,
) -> tuple[LpProblem, VertexChoicesLP, EdgeChoicesLP]:
    """Minimizing the number of deleted rows.

    Parameters
    ----------
    bigraph : BiGraph
        Bigraph
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        PuLP problem
    VertexChoicesLP
        U vertex choices
    EdgeChoicesLP
        Edge choices

    """
    prob = LpProblem(name=QBCHeuristicModels.MIN_DEL_ROWS.value, sense=LpMinimize)
    #
    # Variables
    #
    u_vertex_choices = __u_vertex_choices_variables(bigraph)
    edge_choices = __edge_choices_variables(bigraph)
    #
    # Objective function
    #
    prob += lpSum(u_vertex_choices), "min_u_vertices"
    #
    # Constraints
    #
    __select_edge_select_u_endpoint(prob, bigraph, u_vertex_choices, edge_choices)
    __at_most_epsilon_remaining_zeros(prob, bigraph, epsilon, edge_choices)

    return prob, u_vertex_choices, edge_choices


def min_del_cols(
    bigraph: ZerosBiGraph,
    epsilon: float,
) -> tuple[LpProblem, VertexChoicesLP, EdgeChoicesLP]:
    """Minimizing the number of deleted rows.

    Parameters
    ----------
    bigraph : BiGraph
        Bigraph
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        PuLP problem
    VertexChoicesLP
        V vertex choices
    EdgeChoicesLP
        Edge choices

    """
    prob = LpProblem(name=QBCHeuristicModels.MIN_DEL_COLS.value, sense=LpMinimize)
    #
    # Variables
    #
    v_vertex_choices = __v_vertex_choices_variables(bigraph)
    edge_choices = __edge_choices_variables(bigraph)
    #
    # Objective function
    #
    prob += lpSum(v_vertex_choices), "min_v_vertices"
    #
    # Constraints
    #
    __select_edge_select_v_endpoint(prob, bigraph, v_vertex_choices, edge_choices)
    __at_most_epsilon_remaining_zeros(prob, bigraph, epsilon, edge_choices)

    return prob, v_vertex_choices, edge_choices


def min_del_rows_relax(
    bigraph: ZerosBiGraph,
    epsilon: float,
) -> tuple[LpProblem, VertexChoicesLP, EdgeChoicesLP]:
    """Minimizing the number of deleted rows.

    Parameters
    ----------
    bigraph : BiGraph
        Bigraph
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        PuLP problem
    VertexChoicesLP
        U vertex choices
    EdgeChoicesLP
        Edge choices

    """
    prob = LpProblem(name=QBCHeuristicModels.MIN_DEL_ROWS_RELAX.value, sense=LpMinimize)
    #
    # Variables
    #
    u_vertex_choices = __u_vertex_choices_variables(bigraph, var_type=LpContinuous)
    edge_choices = __edge_choices_variables(bigraph)
    #
    # Objective function
    #
    prob += lpSum(u_vertex_choices), "min_u_vertices"
    #
    # Constraints
    #
    __select_edge_select_u_endpoint(prob, bigraph, u_vertex_choices, edge_choices)
    __at_most_epsilon_remaining_zeros(prob, bigraph, epsilon, edge_choices)

    return prob, u_vertex_choices, edge_choices


def min_del_cols_relax(
    bigraph: ZerosBiGraph,
    epsilon: float,
) -> tuple[LpProblem, VertexChoicesLP, EdgeChoicesLP]:
    """Minimizing the number of deleted rows.

    Parameters
    ----------
    bigraph : BiGraph
        Bigraph
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        PuLP problem
    VertexChoicesLP
        V vertex choices
    EdgeChoicesLP
        Edge choices

    """
    prob = LpProblem(name=QBCHeuristicModels.MIN_DEL_COLS_RELAX.value, sense=LpMinimize)
    #
    # Variables
    #
    v_vertex_choices = __v_vertex_choices_variables(bigraph, var_type=LpContinuous)
    edge_choices = __edge_choices_variables(bigraph)
    #
    # Objective function
    #
    prob += lpSum(v_vertex_choices), "min_v_vertices"
    #
    # Constraints
    #
    __select_edge_select_v_endpoint(prob, bigraph, v_vertex_choices, edge_choices)
    __at_most_epsilon_remaining_zeros(prob, bigraph, epsilon, edge_choices)

    return prob, v_vertex_choices, edge_choices


def min_del_ones(
    bigraph: ZerosBiGraph,
    epsilon: float,
) -> tuple[LpProblem, VertexChoicesLP, VertexChoicesLP, EdgeChoicesLP]:
    """Minimizing the number of deleted ones.

    Parameters
    ----------
    bigraph : BiGraph
        Bigraph
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        PuLP problem
    VertexChoicesLP
        U vertex choices
    VertexChoicesLP
        V vertex choices
    EdgeChoicesLP
        Edge choices

    """
    prob = LpProblem(name=QBCHeuristicModels.MIN_DEL_ONES.value, sense=LpMinimize)
    #
    # Variables
    #
    u_vertex_choices = __u_vertex_choices_variables(bigraph)
    v_vertex_choices = __v_vertex_choices_variables(bigraph)
    edge_choices = __edge_choices_variables(bigraph)
    #
    # Objective function
    #
    prob += (
        (
            lpSum(
                (
                    var * deg / bigraph.card_v_bipart()
                    for var, deg in zip(
                        u_vertex_choices,
                        bigraph.u_bipart_anti_degrees(),
                    )
                ),
            )
            + lpSum(
                (
                    var * deg / bigraph.card_u_bipart()
                    for var, deg in zip(
                        v_vertex_choices,
                        bigraph.v_bipart_anti_degrees(),
                    )
                ),
            )
        ),
        "min_weighted_rows/cols",
    )
    #
    # Constraints
    #
    __select_edge_select_at_least_one_endpoint(
        prob,
        bigraph,
        u_vertex_choices,
        v_vertex_choices,
        edge_choices,
    )
    __at_most_epsilon_remaining_zeros(prob, bigraph, epsilon, edge_choices)

    return prob, u_vertex_choices, v_vertex_choices, edge_choices


# ============================================================================ #
#                           LP MODEL - KNAPSACK MODEL                          #
# ============================================================================ #
def kp_qb(
    bigraph: ZerosBiGraph,
    epsilon: float,
) -> tuple[LpProblem, VertexChoicesLP, VertexChoicesLP]:
    """Knapsack model with zero degrees uncertainty.

    Parameters
    ----------
    bigraph : BiGraph
        Bigraph
    epsilon : float
        Error rate

    Returns
    -------
    LpProblem
        PuLP problem
    VertexChoicesLP
        U vertex choices
    VertexChoicesLP
        V vertex choices

    """
    prob = LpProblem(name=QBCHeuristicModels.KP_QB.value, sense=LpMinimize)
    #
    # Variables
    #
    u_vertex_choices = __u_vertex_choices_variables(bigraph)
    v_vertex_choices = __v_vertex_choices_variables(bigraph)
    #
    # Objective function
    #
    prob += (
        (
            lpSum(
                (
                    var * deg / bigraph.card_v_bipart()
                    for var, deg in zip(
                        u_vertex_choices,
                        bigraph.u_bipart_anti_degrees(),
                    )
                ),
            )
            + lpSum(
                (
                    var * deg / bigraph.card_u_bipart()
                    for var, deg in zip(
                        v_vertex_choices,
                        bigraph.v_bipart_anti_degrees(),
                    )
                ),
            )
        ),
        "knapsack",
    )
    #
    # Constraints
    #
    __at_most_epsilon_zero_degrees(
        prob,
        bigraph,
        epsilon,
        u_vertex_choices,
        v_vertex_choices,
    )

    return prob, u_vertex_choices, v_vertex_choices


# ============================================================================ #
#                                   VARIABLES                                  #
# ============================================================================ #
def __u_vertex_choices_variables(
    bigraph: ZerosBiGraph,
    var_type: str = LpInteger,
) -> VertexChoicesLP:
    return [
        LpVariable(f"xu_{u}", cat=var_type, lowBound=0, upBound=1)
        for u in range(bigraph.card_u_bipart())
    ]


def __v_vertex_choices_variables(
    bigraph: ZerosBiGraph,
    var_type: str = LpInteger,
) -> VertexChoicesLP:
    return [
        LpVariable(f"xv_{v}", cat=var_type, lowBound=0, upBound=1)
        for v in range(bigraph.card_v_bipart())
    ]


def __edge_choices_variables(
    bigraph: ZerosBiGraph,
    var_type: str = LpInteger,
) -> EdgeChoicesLP:
    return {
        (u, v): LpVariable(f"edge_{u}_{v}", cat=var_type, lowBound=0, upBound=1)
        for u, v in bigraph.edges()
    }


# ============================================================================ #
#                                  CONSTRAINTS                                 #
# ============================================================================ #
def __select_edge_select_at_least_one_endpoint(
    prob: LpProblem,
    bigraph: ZerosBiGraph,
    u_vertex_choices: VertexChoicesLP,
    v_vertex_choices: VertexChoicesLP,
    edge_choices: EdgeChoicesLP,
) -> None:
    for u, v in bigraph.edges():
        prob += (
            (u_vertex_choices[u] + v_vertex_choices[v] >= edge_choices[u, v]),
            f"edge_{u}_{v}",
        )


def __select_edge_select_u_endpoint(
    prob: LpProblem,
    bigraph: ZerosBiGraph,
    u_vertex_choices: VertexChoicesLP,
    edge_choices: EdgeChoicesLP,
) -> None:
    for u, v in bigraph.edges():
        prob += (u_vertex_choices[u] >= edge_choices[u, v]), f"edge_{u}_{v}"


def __select_edge_select_v_endpoint(
    prob: LpProblem,
    bigraph: ZerosBiGraph,
    v_vertex_choices: VertexChoicesLP,
    edge_choices: EdgeChoicesLP,
) -> None:
    for u, v in bigraph.edges():
        prob += (v_vertex_choices[v] >= edge_choices[u, v]), f"edge_{u}_{v}"


def __at_most_epsilon_remaining_zeros(
    prob: LpProblem,
    bigraph: ZerosBiGraph,
    epsilon: float,
    edge_choices: EdgeChoicesLP,
) -> None:
    prob += (lpSum(edge_choices) >= (1 - epsilon) * len(bigraph.edges())), "sensitivity"


def __at_most_epsilon_zero_degrees(
    prob: LpProblem,
    bigraph: ZerosBiGraph,
    epsilon: float,
    u_vertex_choices: VertexChoicesLP,
    v_vertex_choices: VertexChoicesLP,
) -> None:
    prob += (
        (
            lpSum(
                (
                    var * deg
                    for var, deg in zip(u_vertex_choices, bigraph.u_bipart_degrees())
                ),
            )
            + lpSum(
                (
                    var * deg
                    for var, deg in zip(v_vertex_choices, bigraph.v_bipart_degrees())
                ),
            )
            >= (1 - epsilon) * len(bigraph.edges())
        ),
        "sensitivity",
    )
