"""Solvers module."""

from __future__ import annotations

import logging
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, cast

from pulp import (
    GUROBI,
    PULP_CBC_CMD,
    LpProblem,
    LpSolutionInfeasible,
    LpStatus,
    LpStatusNotSolved,
)

from qbc_dsbm import UniqueModels
from qbc_dsbm.models.views import (
    LPStats,
    ModelLPMILPStats,
    ModelLPStats,
    ModelSolveStats,
    WarmMILPStats,
)

_LOGGER = logging.getLogger(__name__)


# ============================================================================ #
#                                    CLASSES                                   #
# ============================================================================ #
class Solver(ABC):
    """Solver base class."""

    _NAME: str = "base solver class"
    KEY_NAME = "name"
    KEY_SOLVE_TIME_LIM = "time_limit"
    KEY_PRINT_LOG = "print_log"

    @classmethod
    def name(cls) -> str:
        """Return the solver name.

        Returns
        -------
        str
            Solver name

        """
        return cls._NAME

    @classmethod
    def from_dict(cls, solver_config_dict: dict[str, Any]) -> Solver:
        """Convert dict to Solver.

        Parameters
        ----------
        solver_config_dict : dict
            Dictionnary of solver config

        Returns
        -------
        Solver
            Solver

        """
        return cls(
            time_limit=solver_config_dict[cls.KEY_SOLVE_TIME_LIM],
            print_log=solver_config_dict[cls.KEY_PRINT_LOG],
        )

    def __init__(
        self,
        time_limit: float | None = None,
        *,
        print_log: bool = False,
    ) -> None:
        """Initialize solver."""
        self._time_limit = time_limit
        self._print_log = print_log

    def time_limit(self) -> float | None:
        """Return the time limit.

        Returns
        -------
        float | None
            Time limit

        """
        return self._time_limit

    def print_log(self) -> bool:
        """Return the print log flag.

        Returns
        -------
        bool
            Print log flag

        """
        return self._print_log

    @abstractmethod
    def solve_lp(
        self,
        prob: LpProblem,
        warm_start: bool = False,  # noqa: FBT001, FBT002
    ) -> ModelLPStats:
        """Solve LP with LP solver.

        Parameters
        ----------
        prob : LpProblem
            Problem to solve
        warm_start : bool, default=False
            Warm start

        Returns
        -------
        ModelLPStats
            LP statistics

        Raises
        ------
        NoSolutionError
            If no solution is found
        SolverTimeLimitReached
            If solver time limit is reached

        """
        raise NotImplementedError

    @abstractmethod
    def solve_milp(
        self,
        prob: LpProblem,
        warm_start: bool = False,  # noqa: FBT001, FBT002
    ) -> ModelLPMILPStats:
        """Solve problem with CBC.

        Parameters
        ----------
        prob : LpProblem
            Problem to solve
        warm_start : bool, default=False
            Warm start

        Returns
        -------
        MIPStats
            MIP statistics

        Raises
        ------
        NoSolutionError
            If no solution is found
        SolverTimeLimitReached
            If solver time limit is reached

        """
        raise NotImplementedError

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict.

        Returns
        -------
        dict[str, Any]
            Dictionary

        """
        return {
            self.KEY_NAME: self._NAME,
            self.KEY_SOLVE_TIME_LIM: self._time_limit,
            self.KEY_PRINT_LOG: self._print_log,
        }


class CBCSolver(Solver):
    """CBC solver class."""

    _NAME = "CBC"
    __LOG_FILE_PATH = Path("tmp_cbc.log")

    def solve_lp(
        self,
        prob: LpProblem,
        warm_start: bool = False,  # noqa: FBT001, FBT002
    ) -> ModelLPStats:
        """Solve LP with CBC LP solver.

        Parameters
        ----------
        prob : LpProblem
            Problem to solve
        warm_start : bool, default=False
            Warm start

        Returns
        -------
        ModelLPStats
            LP statistics

        Raises
        ------
        NoSolutionError
            If no solution is found
        SolverTimeLimitReached
            If solver time limit is reached

        """
        if self._print_log:
            _LOGGER.info("Solve stats:")
        prob.solve(
            PULP_CBC_CMD(
                mip=False,
                timeLimit=self._time_limit,
                msg=self._print_log,
                warmStart=warm_start,
            ),
        )
        lp_stats = LPStats(
            LpStatus[prob.status],
            prob.solutionTime,
            prob.solutionCpuTime,
        )
        solve_stats = ModelLPStats(
            cast(UniqueModels, prob.name),
            prob.numVariables(),
            prob.numConstraints(),
            lp_stats,
        )
        if prob.status == LpSolutionInfeasible:
            raise NoSolutionError(solve_stats)
        if prob.status == LpStatusNotSolved:
            raise SolverTimeLimitReachedError(solve_stats)
        lp_stats.set_objective_value(
            float(prob.objective.value()),
        )
        return solve_stats

    def solve_milp(
        self,
        prob: LpProblem,
        warm_start: bool = False,  # noqa: FBT001, FBT002
    ) -> ModelLPMILPStats:
        """Solve problem with CBC.

        Parameters
        ----------
        prob : LpProblem
            Problem to solve
        warm_start : bool, default=False
            Warm start

        Returns
        -------
        MIPStats
            MIP statistics

        Raises
        ------
        NoSolutionError
            If no solution is found
        SolverTimeLimitReached
            If solver time limit is reached

        """
        solve_stats = ModelLPMILPStats(
            cast(UniqueModels, prob.name),
            prob.numVariables(),
            prob.numConstraints(),
        )
        if self._print_log:
            _LOGGER.info("Solve stats:")

        if not warm_start:
            prob.solve(
                PULP_CBC_CMD(
                    mip=False,
                    timeLimit=self._time_limit,
                    msg=self._print_log,
                    warmStart=warm_start,
                ),
            )
            relax_stats = LPStats(
                LpStatus[prob.status],
                prob.solutionTime,
                prob.solutionCpuTime,
            )
            solve_stats.set_relax_stats(relax_stats)
            if prob.status == LpSolutionInfeasible:
                raise NoSolutionError(solve_stats)
            if prob.status == LpStatusNotSolved:
                raise SolverTimeLimitReachedError(solve_stats)
            relax_stats.set_objective_value(
                float(prob.objective.value()),
            )

        prob.solve(
            PULP_CBC_CMD(
                warmStart=True,
                logPath=self.__LOG_FILE_PATH,
                msg=self._print_log,
                timeLimit=self._time_limit,
            ),
        )
        iterations, node_count = self.__get_iterations_node_count()

        milp_stats = WarmMILPStats(
            LpStatus[prob.status],
            prob.solutionTime,
            prob.solutionCpuTime,
            iterations=iterations,
            node_count=node_count,
        )
        solve_stats.set_milp_stats(milp_stats)
        if prob.status == LpSolutionInfeasible:
            raise NoSolutionError(solve_stats)
        if prob.status == LpStatusNotSolved:
            raise SolverTimeLimitReachedError(solve_stats)
        milp_stats.set_objective_value(
            float(prob.objective.value()),
        )
        return solve_stats

    def __get_iterations_node_count(self) -> tuple[int | None, int | None]:
        """Use regex to retrieve these informations."""
        iterations, node_count = None, None
        with self.__LOG_FILE_PATH.open(encoding="utf-8") as file:
            for line in file:
                match = re.search(
                    r"Cbc0001I .+? (\d+) iterations and (\d+) nodes .+?",
                    line,
                )
                if match is not None:
                    iterations, node_count = (
                        int(match.group(1)),
                        int(match.group(2)),
                    )
        self.__LOG_FILE_PATH.unlink()
        return iterations, node_count


class GurobiSolver(Solver):
    """Gurobi solver class."""

    _NAME = "Gurobi"

    def solve_lp(
        self,
        prob: LpProblem,
        warm_start: bool = False,  # noqa: FBT001, FBT002
    ) -> ModelLPStats:
        """Solve LP with CBC LP solver.

        Parameters
        ----------
        prob : LpProblem
            Problem to solve
        warm_start : bool, default=False
            Warm start

        Returns
        -------
        ModelLPStats
            LP statistics

        Raises
        ------
        NoSolutionError
            If no solution is found
        SolverTimeLimitReached
            If solver time limit is reached

        """
        if self._print_log:
            _LOGGER.info("Solve stats:")
        prob.solve(
            GUROBI(
                mip=False,
                timeLimit=self._time_limit,
                msg=self._print_log,
                warmStart=warm_start,
            ),
        )
        lp_stats = LPStats(
            LpStatus[prob.status],
            prob.solutionTime,
            prob.solutionCpuTime,
        )
        solve_stats = ModelLPStats(
            cast(UniqueModels, prob.name),
            prob.numVariables(),
            prob.numConstraints(),
            lp_stats,
        )
        if prob.status == LpSolutionInfeasible:
            raise NoSolutionError(solve_stats)
        if prob.status == LpStatusNotSolved:
            raise SolverTimeLimitReachedError(solve_stats)
        lp_stats.set_objective_value(
            float(prob.objective.value()),
        )
        return solve_stats

    def solve_milp(
        self,
        prob: LpProblem,
        warm_start: bool = False,  # noqa: FBT001, FBT002
    ) -> ModelLPMILPStats:
        """Solve problem with Gurobi.

        Parameters
        ----------
        prob : LpProblem
            Problem to solve
        warm_start : bool, default=False
            Warm start

        Returns
        -------
        MIPStats
            MIP statistics

        Raises
        ------
        NoSolutionError
            If no solution is found
        SolverTimeLimitReached
            If solver time limit is reached

        """
        if self._print_log:
            _LOGGER.info("Solve stats:")
        solve_stats = ModelLPMILPStats(
            cast(UniqueModels, prob.name),
            prob.numVariables(),
            prob.numConstraints(),
        )
        if not warm_start:
            prob.solve(
                GUROBI(
                    mip=False,
                    timeLimit=self._time_limit,
                    msg=self._print_log,
                    warmStart=warm_start,
                ),
            )
            relax_stats = LPStats(
                LpStatus[prob.status],
                prob.solutionTime,
                prob.solutionCpuTime,
            )
            solve_stats.set_relax_stats(relax_stats)
            if prob.status == LpSolutionInfeasible:
                raise NoSolutionError(solve_stats)
            if prob.status == LpStatusNotSolved:
                raise SolverTimeLimitReachedError(solve_stats)
            relax_stats.set_objective_value(
                float(prob.objective.value()),
            )

        prob.solve(
            GUROBI(
                warmStart=True,
                timeLimit=self._time_limit,
                msg=self._print_log,
            ),
        )
        milp_stats = WarmMILPStats(
            LpStatus[prob.status],
            prob.solutionTime,
            prob.solutionCpuTime,
            iterations=prob.solverModel.IterCount,
            node_count=prob.solverModel.NodeCount,
        )
        solve_stats.set_milp_stats(milp_stats)
        if prob.status == LpSolutionInfeasible:
            raise NoSolutionError(solve_stats)
        if prob.status == LpStatusNotSolved:
            raise SolverTimeLimitReachedError(solve_stats)
        milp_stats.set_objective_value(
            float(prob.objective.value()),
        )
        return solve_stats


# ============================================================================ #
#                                   FUNCTIONS                                  #
# ============================================================================ #
def solver_from_dict(solver_config_dict: dict[str, Any]) -> Solver:
    """Convert dict to solver.

    Parameters
    ----------
    solver_config_dict : dict
        Dictionnary of solver config

    Returns
    -------
    Solver
        Solver

    Raises
    ------
    ValueError
        If the solver is unknown

    """
    if solver_config_dict[Solver.KEY_NAME] == GurobiSolver.name():
        return GurobiSolver.from_dict(solver_config_dict)
    if solver_config_dict[Solver.KEY_NAME] == CBCSolver.name():
        return CBCSolver.from_dict(solver_config_dict)
    msg = f"Unknown solver: {solver_config_dict[Solver.KEY_NAME]}"
    raise ValueError(msg)


# ============================================================================ #
#                                  EXCEPTIONS                                  #
# ============================================================================ #
class NoSolutionError(Exception):
    """No solution exception."""

    def __init__(self, solve_stats: ModelSolveStats) -> None:
        """Initialize the exception."""
        super().__init__()
        self.__solve_stats: ModelSolveStats = solve_stats

    def solve_stats(self) -> ModelSolveStats:
        """Solve statistics.

        Returns
        -------
        SolveStats
            Solve statistics

        """
        return self.__solve_stats

    def __str__(self) -> str:
        """Print the exception message.

        Returns
        -------
        str
            Exception message

        """
        return (
            f"The problem is unfeasible, see the solve statistics:\n\n"
            f"{self.__solve_stats}"
        )


class SolverTimeLimitReachedError(Exception):
    """Solver time limit reached exception."""

    def __init__(self, solve_stats: ModelSolveStats) -> None:
        """Initialize the exception."""
        super().__init__()
        self.__solve_stats: ModelSolveStats = solve_stats

    def solve_stats(self) -> ModelSolveStats:
        """Solve statistics.

        Returns
        -------
        SolveStats
            Solve statistics

        """
        return self.__solve_stats

    def __str__(self) -> str:
        """Print the exception message.

        Returns
        -------
        str
            Exception message

        """
        return (
            f"The solver time limit is reached,"
            " see the solve statistics:\n\n"
            f"{self.__solve_stats}"
        )
