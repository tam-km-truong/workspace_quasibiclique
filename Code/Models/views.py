"""MIP stats module."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from qbc_dsbm import UniqueModels

import enum


class ModelType(enum.Enum):
    """Model type enum."""

    LP = "LP"
    """LP model."""

    MILP = "MILP"
    """MILP model."""


# ============================================================================ #
#                                    CLASSES                                   #
# ============================================================================ #
class LPStats:
    """LP statistics."""

    __STATUS = "status"
    __USER_TIME = "user_time"
    __CPU_TIME = "cpu_time"
    __OBJ_VALUE = "objective_value"

    @classmethod
    def from_dict(cls, relax_stats_dict: dict[str, Any]) -> LPStats:
        """Convert dict to RelaxStats.

        Parameters
        ----------
        relax_stats_dict : dict
            Dictionnary of relax stats

        Returns
        -------
        RelaxStats
            RelaxStats

        """
        return cls(
            status=relax_stats_dict[cls.__STATUS],
            relax_user_time=relax_stats_dict[cls.__USER_TIME],
            relax_cpu_time=relax_stats_dict[cls.__CPU_TIME],
            objective_value=relax_stats_dict[cls.__OBJ_VALUE],
        )

    def __init__(
        self,
        status: str,
        relax_user_time: float,
        relax_cpu_time: float,
        objective_value: float | None = None,
    ) -> None:
        """Initialize the stats."""
        self.__status = status
        self.__user_time = relax_user_time
        self.__cpu_time = relax_cpu_time
        self.__objective_value = objective_value

    def status(self) -> str:
        """Status.

        Returns
        -------
        str
            Status

        """
        return self.__status

    def user_time(self) -> float:
        """User time (in seconds).

        Returns
        -------
        float
            User time

        """
        return self.__user_time

    def cpu_time(self) -> float:
        """CPU time (in seconds).

        Returns
        -------
        float
            CPU time

        """
        return self.__cpu_time

    def objective_value(self) -> float | None:
        """Objective value.

        Returns
        -------
        float or None
            Objective value

        """
        return self.__objective_value

    def set_objective_value(self, value: float) -> None:
        """Set objective value.

        Parameters
        ----------
        value : float
            Objective value

        """
        self.__objective_value = value

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict.

        Returns
        -------
        dict
            Dict

        """
        return {
            self.__STATUS: self.__status,
            self.__USER_TIME: self.__user_time,
            self.__CPU_TIME: self.__cpu_time,
            self.__OBJ_VALUE: self.__objective_value,
        }

    def __str__(self) -> str:
        """Print the statistics.

        Returns
        -------
        str
            Statistics

        """
        return "\n".join(
            [
                f"Status: {self.__status}",
                f"User time: {self.__user_time}",
                f"CPU time: {self.__cpu_time}",
                f"Objective value: {self.__objective_value}",
            ],
        )


class WarmMILPStats:
    """Statistics for MIP problems."""

    __STATUS = "status"
    __USER_TIME = "user_time"
    __CPU_TIME = "cpu_time"
    __OBJ_VALUE = "objective_value"
    __ITERATIONS = "iterations"
    __NODE_COUNT = "node_count"

    @classmethod
    def from_dict(cls, mip_stats_dict: dict[str, Any]) -> WarmMILPStats:
        """Convert dict to MIPStats.

        Parameters
        ----------
        mip_stats_dict : dict
            Dictionnary of MIP stats

        Returns
        -------
        MIPStats
            MIPStats

        """
        return WarmMILPStats(
            status=mip_stats_dict[cls.__STATUS],
            user_time=mip_stats_dict[cls.__USER_TIME],
            cpu_time=mip_stats_dict[cls.__CPU_TIME],
            obj_value=mip_stats_dict[cls.__OBJ_VALUE],
            iterations=mip_stats_dict[cls.__ITERATIONS],
            node_count=mip_stats_dict[cls.__NODE_COUNT],
        )

    def __init__(  # noqa: PLR0913
        self,
        status: str,
        user_time: float,
        cpu_time: float,
        obj_value: float | None = None,
        iterations: int | None = None,
        node_count: int | None = None,
    ) -> None:
        """Initialize the stats."""
        self.__status: str = status
        self.__user_time: float = user_time
        self.__cpu_time: float = cpu_time
        self.__obj_value: float | None = obj_value
        self.__iterations: int | None = iterations
        self.__node_count: int | None = node_count

    def status(self) -> str:
        """Solver status.

        Returns
        -------
        str
            Solver status

        """
        return self.__status

    def user_time(self) -> float:
        """User time (in seconds).

        Returns
        -------
        float
            User time

        """
        return self.__user_time

    def cpu_time(self) -> float:
        """CPU time (in seconds).

        Returns
        -------
        float
            CPU time

        """
        return self.__cpu_time

    def objective_value(self) -> float | None:
        """Objective value.

        Returns
        -------
        float or None
            Objective value

        """
        return self.__obj_value

    def iterations(self) -> int | None:
        """Give the number of iterations used to solve the MIP.

        Returns
        -------
        int
            Number of iterations

        """
        return self.__iterations

    def node_count(self) -> int | None:
        """Give the number of nodes explored by the solver.

        Returns
        -------
        int or None
            Number of nodes

        """
        return self.__node_count

    def set_objective_value(self, value: float) -> None:
        """Set objective value.

        Parameters
        ----------
        value : float
            Objective value

        """
        self.__obj_value = value

    def set_iterations(self, value: int) -> None:
        """Set number of iterations.

        Parameters
        ----------
        value : int
            Number of iterations

        """
        self.__iterations = value

    def set_node_count(self, value: int) -> None:
        """Set number of nodes.

        Parameters
        ----------
        value : int
            Number of nodes

        """
        self.__node_count = value

    def to_dict(self) -> dict:
        """Convert MIPStats to dict.

        Returns
        -------
        dict
            Dictionnary

        """
        return {
            self.__STATUS: self.__status,
            self.__USER_TIME: self.__user_time,
            self.__CPU_TIME: self.__cpu_time,
            self.__OBJ_VALUE: self.__obj_value,
            self.__ITERATIONS: self.__iterations,
            self.__NODE_COUNT: self.__node_count,
        }

    def __str__(self) -> str:
        """Give the string representation.

        Returns
        -------
        str
            String representation

        """
        return "\n".join(
            [
                f"Status: {self.__status}",
                f"User time: {self.__user_time}",
                f"CPU time: {self.__cpu_time}",
                f"Objective value: {self.__obj_value}",
                f"Iterations: {self.__iterations}",
                f"Node count: {self.__node_count}",
            ],
        )


class ModelSolveStats:
    """Solve stats parent class."""

    _MODEL_TYPE = ModelType.MILP

    _MODEL_TYPE_KEY = "model_type"
    _MODEL_NAME_KEY = "model"

    _NUMBER_VARIABLES_KEY = "number_variables"
    _NUMBER_CONSTRAINTS_KEY = "number_constraints"

    @classmethod
    def model_type(cls) -> ModelType:
        """Model type.

        Returns
        -------
        ModelType
            Model type

        """
        return cls._MODEL_TYPE

    @classmethod
    def model_type_from_dict(cls, stats_dict: dict) -> ModelType:
        """Retrun model type from a dictionnary.

        Parameters
        ----------
        stats_dict : dict
            Dictionary

        Returns
        -------
        ModelType
            Model type

        Raises
        ------
        ValueError
            Model type not found

        """
        if cls._MODEL_TYPE_KEY in stats_dict:
            return ModelType(stats_dict[cls._MODEL_TYPE_KEY])
        msg = f"Model type not found in {stats_dict}."
        raise ValueError(msg)

    @classmethod
    def from_dict(cls, stats_dict: dict) -> ModelSolveStats:
        """Convert dict to ModelSolveStats.

        Parameters
        ----------
        stats_dict : dict
            Dictionary

        Returns
        -------
        ModelSolveStats
            ModelSolveStats

        """
        return cls(
            stats_dict[cls._MODEL_NAME_KEY],
            stats_dict[cls._NUMBER_VARIABLES_KEY],
            stats_dict[cls._NUMBER_CONSTRAINTS_KEY],
        )

    def __init__(
        self,
        model: UniqueModels,
        number_variables: int,
        number_constraints: int,
    ) -> None:
        """Initialize the stats."""
        self._model = model
        self._number_variables = number_variables
        self._number_constraints = number_constraints

    def model(self) -> UniqueModels:
        """Model.

        Returns
        -------
        UniqueModels
            Model

        """
        return self._model

    def number_variables(self) -> int:
        """Give the number of variables in the MIP.

        Returns
        -------
        int
            Number of variables

        """
        return self._number_variables

    def number_constraints(self) -> int:
        """Give the number of constraints in the MIP.

        Returns
        -------
        int
            Number of constraints

        """
        return self._number_constraints

    def to_dict(self) -> dict:
        """Convert SolveStats to dict.

        Returns
        -------
        dict
            Dictionnary

        """
        return {
            self._MODEL_TYPE_KEY: self._MODEL_TYPE.value,
            self._MODEL_NAME_KEY: self._model,
            self._NUMBER_VARIABLES_KEY: self._number_variables,
            self._NUMBER_CONSTRAINTS_KEY: self._number_constraints,
        }

    def __str__(self) -> str:
        """Give the string representation.

        Returns
        -------
        str
            String representation

        """
        return "\n".join(
            [
                f"Model type: {self._MODEL_TYPE.value}",
                f"Model: {self._model}",
                f"Number of variables: {self._number_variables}",
                f"Number of constraints: {self._number_constraints}",
            ],
        )


class ModelLPStats(ModelSolveStats):
    """LP solve stats."""

    _MODEL_TYPE = ModelType.LP

    __SOLVE_STATS_KEY = "solve_stats"

    @classmethod
    def from_dict(cls, stats_dict: dict) -> ModelLPStats:
        """Convert dict to ModelLPStats.

        Parameters
        ----------
        stats_dict : dict
            Dictionary

        Returns
        -------
        ModelLPStats
            ModelLPStats

        """
        return cls(
            stats_dict[cls._MODEL_NAME_KEY],
            stats_dict[cls._NUMBER_VARIABLES_KEY],
            stats_dict[cls._NUMBER_CONSTRAINTS_KEY],
            LPStats.from_dict(stats_dict[cls.__SOLVE_STATS_KEY]),
        )

    def __init__(
        self,
        model: UniqueModels,
        number_variables: int,
        number_constraints: int,
        relax_stats: LPStats,
    ) -> None:
        """Initialize the stats."""
        super().__init__(model, number_variables, number_constraints)
        self._relax_stats = relax_stats

    def relax_stats(self) -> LPStats:
        """Relaxation stats.

        Returns
        -------
        LPStats
            Relaxation stats

        """
        return self._relax_stats

    def to_dict(self) -> dict:
        """Convert SolveStats to dict.

        Returns
        -------
        dict
            Dictionnary

        """
        return {
            **super().to_dict(),
            self.__SOLVE_STATS_KEY: self._relax_stats.to_dict(),
        }

    def __str__(self) -> str:
        """Give the string representation.

        Returns
        -------
        str
            String representation

        """
        return "\n".join(
            [
                super().__str__(),
                self._relax_stats.__str__(),
            ],
        )


class ModelLPMILPStats(ModelSolveStats):
    """LP + warm MILP solve stats."""

    _MODEL_TYPE = ModelType.MILP

    __RELAX_STATS_KEY = "relaxation"
    __MILP_STATS_KEY = "milp"

    @classmethod
    def from_dict(cls, stats_dict: dict) -> ModelLPMILPStats:
        """Convert dict to ModelLPMILPStats.

        Parameters
        ----------
        stats_dict : dict
            Dictionary

        Returns
        -------
        ModelLPMILPStats
            ModelLPMILPStats

        """
        return cls(
            stats_dict[cls._MODEL_NAME_KEY],
            stats_dict[cls._NUMBER_VARIABLES_KEY],
            stats_dict[cls._NUMBER_CONSTRAINTS_KEY],
            (
                LPStats.from_dict(stats_dict[cls.__RELAX_STATS_KEY])
                if stats_dict[cls.__RELAX_STATS_KEY] is not None
                else None
            ),
            (
                WarmMILPStats.from_dict(stats_dict[cls.__MILP_STATS_KEY])
                if stats_dict[cls.__MILP_STATS_KEY] is not None
                else None
            ),
        )

    def __init__(
        self,
        model: UniqueModels,
        number_variables: int,
        number_constraints: int,
        relax_stats: LPStats | None = None,
        milp_stats: WarmMILPStats | None = None,
    ) -> None:
        """Initialize the stats."""
        super().__init__(model, number_variables, number_constraints)
        self.__relax_stats: LPStats | None = relax_stats
        self.__milp_stats: WarmMILPStats | None = milp_stats

    def relax_stats(self) -> LPStats | None:
        """Relaxation statistics.

        Returns
        -------
        RelaxStats
            Relaxation statistics

        """
        return self.__relax_stats

    def milp_stats(self) -> WarmMILPStats | None:
        """MILP statistics.

        Returns
        -------
        MILPStats
            MILP statistics

        """
        return self.__milp_stats

    def set_relax_stats(self, relax_stats: LPStats) -> None:
        """Set relaxation statistics.

        Parameters
        ----------
        relax_stats : RelaxStats
            Relaxation statistics

        """
        self.__relax_stats = relax_stats

    def set_milp_stats(self, milp_stats: WarmMILPStats) -> None:
        """Set MILP statistics.

        Parameters
        ----------
        milp_stats : MILPStats
            MILP statistics

        """
        self.__milp_stats = milp_stats

    def to_dict(self) -> dict:
        """Convert ModelLPMILPStats to dict.

        Returns
        -------
        dict
            Dictionnary

        """
        return {
            **super().to_dict(),
            self.__RELAX_STATS_KEY: (
                self.__relax_stats.to_dict() if self.__relax_stats is not None else None
            ),
            self.__MILP_STATS_KEY: (
                self.__milp_stats.to_dict() if self.__milp_stats is not None else None
            ),
        }

    def __str__(self) -> str:
        """Give the string representation.

        Returns
        -------
        str
            String representation

        """
        return "\n".join(
            [
                f"{super().__str__()}",
                "",
                "Relaxation statistics",
                "---------------------",
                f"{self.__relax_stats}",
                "",
                "MILP statistics",
                "---------------",
                f"{self.__milp_stats}",
            ],
        )
