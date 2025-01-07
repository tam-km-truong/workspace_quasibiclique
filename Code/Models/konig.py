"""ILP König models."""

from __future__ import annotations

from typing import TYPE_CHECKING

from pulp import LpContinuous, LpMinimize, LpProblem, LpVariable, lpSum

from qbc_dsbm import BCModels

if TYPE_CHECKING:
    from qbc_dsbm.bigraph import ZerosBiGraph

# ============================================================================ #
#                                     TYPES                                    #
# ============================================================================ #
VertexChoicesLP = list[LpVariable]


# ============================================================================ #
#                      LP MODEL - KÖNIG THEOREM CLASSICAL                      #
# ============================================================================ #
def konig_v(
    bigraph: ZerosBiGraph,
) -> tuple[LpProblem, VertexChoicesLP, VertexChoicesLP]:
    """König_V model.

    Parameters
    ----------
    bigraph : BiGraph
        Bigraph

    Returns
    -------
    LpProblem
        PuLP problem
    VertexChoicesLP
        U vertex choices
    VertexChoicesLP
        V vertex choices

    """
    prob = LpProblem(name=BCModels.KONIG_V.value, sense=LpMinimize)
    #
    # Variables
    #
    u_vertex_choices, v_vertex_choices = __vertex_choices_variables(bigraph)
    #
    # Objective function
    #
    prob += lpSum(u_vertex_choices) + lpSum(v_vertex_choices), "min_vertices"
    #
    # Constraints
    #
    __delete_row_or_col_contains_zero(prob, bigraph, u_vertex_choices, v_vertex_choices)

    return prob, u_vertex_choices, v_vertex_choices


# ============================================================================ #
#                         LP MODEL - KÖNIG WITH DEGREES                        #
# ============================================================================ #
def konig_e(
    bigraph: ZerosBiGraph,
) -> tuple[LpProblem, VertexChoicesLP, VertexChoicesLP]:
    """König_E model.

    Parameters
    ----------
    bigraph : BiGraph
        Bigraph

    Returns
    -------
    LpProblem
        PuLP problem
    VertexChoicesLP
        U vertex choices
    VertexChoicesLP
        V vertex choices

    """
    prob = LpProblem(name=BCModels.KONIG_E.value, sense=LpMinimize)
    #
    # Variables
    #
    u_vertex_choices, v_vertex_choices = __vertex_choices_variables(bigraph)
    #
    # Objective function
    #
    prob += (
        (
            lpSum(
                (
                    var * deg
                    for var, deg in zip(
                        u_vertex_choices,
                        bigraph.u_bipart_anti_degrees(),
                    )
                ),
            )
            + lpSum(
                (
                    var * deg
                    for var, deg in zip(
                        v_vertex_choices,
                        bigraph.v_bipart_anti_degrees(),
                    )
                ),
            )
        ),
        "min_lost_of_ones",
    )
    #
    # Constraints
    #
    __delete_row_or_col_contains_zero(prob, bigraph, u_vertex_choices, v_vertex_choices)

    return prob, u_vertex_choices, v_vertex_choices


# ============================================================================ #
#                                   VARIABLES                                  #
# ============================================================================ #
def __vertex_choices_variables(
    bigraph: ZerosBiGraph,
) -> tuple[VertexChoicesLP, VertexChoicesLP]:
    return (
        [
            LpVariable(f"xu_{u}", cat=LpContinuous, lowBound=0, upBound=1)
            for u in range(bigraph.card_u_bipart())
        ],
        [
            LpVariable(f"xv_{v}", cat=LpContinuous, lowBound=0, upBound=1)
            for v in range(bigraph.card_v_bipart())
        ],
    )


# ============================================================================ #
#                                  CONSTRAINTS                                 #
# ============================================================================ #
def __delete_row_or_col_contains_zero(
    prob: LpProblem,
    bigraph: ZerosBiGraph,
    u_vertex_choices: VertexChoicesLP,
    v_vertex_choices: VertexChoicesLP,
) -> None:
    for u, v in bigraph.edges():
        prob += (u_vertex_choices[u] + v_vertex_choices[v] >= 1), f"edge_{u}_{v}"
