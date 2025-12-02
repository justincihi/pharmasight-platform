from __future__ import annotations
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, Callable, Tuple, List
from enum import Enum
import numpy as np
from scipy.integrate import solve_ivp
import logging

logger = logging.getLogger(__name__)


class AdminRoute(Enum):
    """Supported administration routes."""
    IV = "iv"
    ORAL = "oral"
    IM = "im"
    SC = "sc"
    INHALATION = "inhalation"


class ModelType(Enum):
    """Supported model types."""
    ONE_COMPARTMENT = "one_compartment"
    TWO_COMPARTMENT = "two_compartment"
    THREE_COMPARTMENT = "three_compartment"
    PBPK = "pbpk"
    CUSTOM = "custom"


@dataclass
class PKParameters:
    """Generic container for PK parameters with validation."""
    CL: float  # clearance (L/h)
    Vc: float  # central volume (L)
    Ka: Optional[float] = None  # absorption rate (1/h), if applicable
    Qp: Optional[float] = None  # intercompartmental clearance (L/h), optional
    Vp: Optional[float] = None  # peripheral volume (L), optional
    F: float = 1.0  # bioavailability (0-1)
    extra: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate parameter values."""
        if self.CL <= 0:
            raise ValueError(f"CL must be > 0, got {self.CL}")
        if self.Vc <= 0:
            raise ValueError(f"Vc must be > 0, got {self.Vc}")
        if self.Ka is not None and self.Ka <= 0:
            raise ValueError(f"Ka must be > 0, got {self.Ka}")
        if self.Qp is not None and self.Qp < 0:
            raise ValueError(f"Qp must be >= 0, got {self.Qp}")
        if self.Vp is not None and self.Vp <= 0:
            raise ValueError(f"Vp must be > 0, got {self.Vp}")
        if not (0 <= self.F <= 1):
            raise ValueError(f"F (bioavailability) must be between 0 and 1, got {self.F}")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "CL": self.CL,
            "Vc": self.Vc,
            "Ka": self.Ka,
            "Qp": self.Qp,
            "Vp": self.Vp,
            "F": self.F,
            "extra": self.extra,
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> PKParameters:
        """Create from dictionary representation."""
        return cls(**data)


@dataclass
class SimulationResult:
    """Results container for PK simulations."""
    t: np.ndarray
    y: np.ndarray
    meta: Dict[str, Any] = field(default_factory=dict)
    
    def get_central_compartment(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return time and central compartment concentration."""
        # Assumes first row of y is central compartment
        return self.t, self.y[0]
    
    def get_all_compartments(self) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """Return all compartment profiles."""
        compartment_names = self.meta.get("compartment_names", [f"comp_{i}" for i in range(self.y.shape[0])])
        return {
            name: (self.t, self.y[i])
            for i, name in enumerate(compartment_names)
        }
    
    def get_auc(self, t_start: float = 0, t_end: Optional[float] = None) -> float:
        """Calculate AUC using trapezoidal rule."""
        if t_end is None:
            t_end = self.t[-1]
        
        mask = (self.t >= t_start) & (self.t <= t_end)
        t_subset = self.t[mask]
        y_subset = self.y[0][mask]  # central compartment
        
        if len(t_subset) < 2:
            return 0.0
        
        auc = np.trapz(y_subset, t_subset)
        return float(auc)
    
    def get_cmax(self) -> Tuple[float, float]:
        """Return (Cmax, Tmax) for central compartment."""
        central_conc = self.y[0]
        cmax_idx = np.argmax(central_conc)
        return float(central_conc[cmax_idx]), float(self.t[cmax_idx])


class PKModel(ABC):
    """
    Base class for all PK models in PharmaSight.
    Concrete subclasses must implement `rhs` and `initial_state`.
    
    Features:
    - Parameter validation
    - DDI support with audit trail
    - Covariate adjustments
    - Multiple administration routes
    - Numerical stability controls
    """

    def __init__(self, params: PKParameters, model_type: ModelType = ModelType.CUSTOM):
        self.params = params
        self.model_type = model_type
        self.ddi_history: List[Dict[str, Any]] = []
        self._validate_params()

    def _validate_params(self):
        """Validate that model-specific parameters are present."""
        pass  # Override in subclasses for specific requirements

    @abstractmethod
    def rhs(self, t: float, y: np.ndarray,
            covariates: Optional[Dict[str, Any]] = None,
            ddi_factors: Optional[Dict[str, float]] = None) -> np.ndarray:
        """
        Right-hand side of the ODE system: dy/dt = f(t, y, params, covariates, ddi_factors).
        
        Args:
            t: Current time
            y: State vector
            covariates: Optional covariate dictionary (e.g., weight, age, CYP activity)
            ddi_factors: Optional DDI factors (multiplicative adjustments to parameters)
            
        Returns:
            dy/dt: Time derivatives
        """
        raise NotImplementedError

    @abstractmethod
    def initial_state(self, dose_mg: float,
                      route: str = "iv",
                      covariates: Optional[Dict[str, Any]] = None) -> np.ndarray:
        """
        Return initial state vector for given dose and route.
        
        Args:
            dose_mg: Dose in mg
            route: Administration route (from AdminRoute enum)
            covariates: Optional covariates affecting initial state
            
        Returns:
            Initial state vector y0
        """
        raise NotImplementedError
    
    @abstractmethod
    def get_compartment_names(self) -> List[str]:
        """Return names of compartments in order (e.g., ['central', 'peripheral'])."""
        raise NotImplementedError

    def apply_ddi(self, ddi_factors: Optional[Dict[str, float]]) -> PKParameters:
        """
        Adjust parameters for DDI. Example: {'CL': 0.5} means CL is halved.
        
        Args:
            ddi_factors: Dictionary of parameter adjustments (multiplicative factors)
            
        Returns:
            Modified PKParameters object
        """
        if not ddi_factors:
            return self.params

        # Create copy of original parameters
        original_dict = self.params.to_dict()
        adjusted_dict = dict(original_dict)

        # Track DDI changes
        ddi_record = {"timestamp": len(self.ddi_history), "changes": {}}

        for key, factor in ddi_factors.items():
            if key in original_dict and original_dict[key] is not None:
                original_value = original_dict[key]
                adjusted_value = original_value * factor
                adjusted_dict[key] = adjusted_value
                ddi_record["changes"][key] = {
                    "original": original_value,
                    "factor": factor,
                    "adjusted": adjusted_value,
                }

        self.ddi_history.append(ddi_record)
        
        try:
            return PKParameters(**adjusted_dict)
        except ValueError as e:
            logger.error(f"Invalid DDI adjustment: {e}")
            raise

    def simulate(
        self,
        dose_mg: float,
        route: str,
        t_end: float,
        n_points: int = 200,
        covariates: Optional[Dict[str, Any]] = None,
        ddi_factors: Optional[Dict[str, float]] = None,
        rtol: float = 1e-6,
        atol: float = 1e-8,
        method: str = "RK45",
    ) -> SimulationResult:
        """
        Run a simulation and return concentrations / amounts over time.
        
        Args:
            dose_mg: Dose in mg
            route: Administration route
            t_end: End time in hours
            n_points: Number of output time points
            covariates: Optional covariate adjustments
            ddi_factors: Optional DDI adjustments
            rtol: Relative tolerance for ODE solver
            atol: Absolute tolerance for ODE solver
            method: Integration method (RK45, BDF, etc.)
            
        Returns:
            SimulationResult with time, states, and metadata
            
        Raises:
            RuntimeError: If simulation fails
            ValueError: If parameters are invalid
        """
        if t_end <= 0:
            raise ValueError(f"t_end must be > 0, got {t_end}")
        if n_points < 2:
            raise ValueError(f"n_points must be >= 2, got {n_points}")

        times = np.linspace(0.0, t_end, n_points)
        
        try:
            y0 = self.initial_state(dose_mg, route=route, covariates=covariates)
        except Exception as e:
            logger.error(f"Failed to compute initial state: {e}")
            raise
        
        effective_params = self.apply_ddi(ddi_factors)

        # Create a closure that captures the effective parameters
        def wrapped_rhs(t, y):
            return self.rhs(t, y, covariates=covariates, ddi_factors=ddi_factors)

        try:
            sol = solve_ivp(
                wrapped_rhs,
                t_span=(times[0], times[-1]),
                y0=y0,
                t_eval=times,
                rtol=rtol,
                atol=atol,
                method=method,
            )
        except Exception as e:
            logger.error(f"ODE solver failed: {e}")
            raise RuntimeError(f"PK simulation ODE integration failed: {e}")

        if not sol.success:
            raise RuntimeError(f"PK simulation failed: {sol.message}")

        return SimulationResult(
            t=sol.t,
            y=sol.y,
            meta={
                "dose_mg": dose_mg,
                "route": route,
                "params": effective_params.to_dict(),
                "covariates": covariates or {},
                "ddi_factors": ddi_factors or {},
                "model_type": self.model_type.value,
                "compartment_names": self.get_compartment_names(),
            },
        )
    
    def get_ddi_history(self) -> List[Dict[str, Any]]:
        """Return DDI adjustment history."""
        return self.ddi_history.copy()
    
    def reset_ddi_history(self):
        """Clear DDI adjustment history."""
        self.ddi_history.clear()
