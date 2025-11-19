"""
Three-compartment pharmacokinetic models.

Supports IV bolus and first-order absorption with two peripheral compartments.
Useful for drugs with complex distribution patterns (e.g., lipophilic compounds).
"""

from __future__ import annotations
from typing import Dict, Any, Optional, List
import numpy as np
import logging

from .base import PKModel, PKParameters, ModelType, AdminRoute

logger = logging.getLogger(__name__)


class ThreeCompartmentIV(PKModel):
    """
    3-compartment IV bolus model with two peripheral compartments.
    
    Model equations:
        dA_c/dt  = - (CL/Vc + Qp1/Vc + Qp2/Vc) * A_c + (Qp1/Vp1) * A_p1 + (Qp2/Vp2) * A_p2
        dA_p1/dt =   (Qp1/Vc) * A_c - (Qp1/Vp1) * A_p1
        dA_p2/dt =   (Qp2/Vc) * A_c - (Qp2/Vp2) * A_p2
        C_c = A_c / Vc
    
    where:
        A_c: amount in central compartment (mg)
        A_p1: amount in peripheral compartment 1 (mg)
        A_p2: amount in peripheral compartment 2 (mg)
        Qp1, Qp2: intercompartmental clearances (L/h)
        Vp1, Vp2: peripheral volumes (L)
        CL: systemic clearance (L/h)
        Vc: central volume (L)
    
    Note: This model uses the 'extra' dict to store Qp2 and Vp2:
        params.extra = {"Qp2": <value>, "Vp2": <value>}
    
    Typical uses: Highly lipophilic drugs, oncology agents, complex PK
    State vector: [A_c, A_p1, A_p2]
    """

    def __init__(self, params: PKParameters):
        """Initialize three-compartment IV model."""
        super().__init__(params, model_type=ModelType.THREE_COMPARTMENT)
        self._validate_params()

    def _validate_params(self):
        """Ensure required parameters are present."""
        if self.params.Qp is None or self.params.Vp is None:
            raise ValueError(
                "ThreeCompartmentIV requires Qp (Qp1) and Vp (Vp1). "
                f"Got Qp={self.params.Qp}, Vp={self.params.Vp}"
            )
        
        if "Qp2" not in self.params.extra or "Vp2" not in self.params.extra:
            raise ValueError(
                "ThreeCompartmentIV requires 'Qp2' and 'Vp2' in params.extra dict. "
                f"Got extra={self.params.extra}"
            )
        
        Qp2 = self.params.extra.get("Qp2")
        Vp2 = self.params.extra.get("Vp2")
        
        if Qp2 is None or Qp2 < 0:
            raise ValueError(f"Qp2 must be >= 0, got {Qp2}")
        if Vp2 is None or Vp2 <= 0:
            raise ValueError(f"Vp2 must be > 0, got {Vp2}")

    def rhs(
        self,
        t: float,
        y: np.ndarray,
        covariates: Optional[Dict[str, Any]] = None,
        ddi_factors: Optional[Dict[str, float]] = None,
    ) -> np.ndarray:
        """
        Compute right-hand side of ODE system.
        
        Args:
            t: Current time (h)
            y: State vector [A_c, A_p1, A_p2] (amounts in mg)
            covariates: Optional covariate adjustments
            ddi_factors: Optional DDI factors (applied upstream)
            
        Returns:
            dy/dt: [dA_c/dt, dA_p1/dt, dA_p2/dt] (mg/h)
        """
        A_c, A_p1, A_p2 = y
        
        # Get parameters with covariate adjustments
        CL = self.params.CL
        Vc = self.params.Vc
        Qp1 = self.params.Qp  # First peripheral clearance
        Vp1 = self.params.Vp  # First peripheral volume
        Qp2 = self.params.extra.get("Qp2", 0.0)
        Vp2 = self.params.extra.get("Vp2", 1.0)
        
        if covariates:
            CL = self._adjust_cl_by_covariates(CL, covariates)
            Vc = self._adjust_vc_by_covariates(Vc, covariates)
            Qp1 = self._adjust_qp_by_covariates(Qp1, covariates, "Qp1")
            Vp1 = self._adjust_vp_by_covariates(Vp1, covariates, "Vp1")
            Qp2 = self._adjust_qp_by_covariates(Qp2, covariates, "Qp2")
            Vp2 = self._adjust_vp_by_covariates(Vp2, covariates, "Vp2")
        
        # Rate constants
        k_el = CL / Vc
        k_cp1 = Qp1 / Vc  # central to peripheral1
        k_p1c = Qp1 / Vp1  # peripheral1 to central
        k_cp2 = Qp2 / Vc  # central to peripheral2
        k_p2c = Qp2 / Vp2  # peripheral2 to central
        
        # ODE system
        dA_c_dt = -(k_el + k_cp1 + k_cp2) * A_c + k_p1c * A_p1 + k_p2c * A_p2
        dA_p1_dt = k_cp1 * A_c - k_p1c * A_p1
        dA_p2_dt = k_cp2 * A_c - k_p2c * A_p2
        
        return np.array([dA_c_dt, dA_p1_dt, dA_p2_dt])

    def initial_state(
        self,
        dose_mg: float,
        route: str = "iv",
        covariates: Optional[Dict[str, Any]] = None,
    ) -> np.ndarray:
        """
        Compute initial state for IV bolus.
        
        Args:
            dose_mg: Dose in mg
            route: Administration route (must be 'iv')
            covariates: Optional covariate adjustments
            
        Returns:
            Initial state vector y0 = [A_c, A_p1, A_p2] (mg)
            
        Raises:
            ValueError: If route is not IV
        """
        route_lower = route.lower().strip()
        
        if route_lower not in ["iv", AdminRoute.IV.value]:
            raise ValueError(
                f"ThreeCompartmentIV only supports IV route, got '{route}'. "
                f"For oral, use ThreeCompartmentOral."
            )
        
        initial_amount = dose_mg * self.params.F
        
        return np.array([initial_amount, 0.0, 0.0])

    def get_compartment_names(self) -> List[str]:
        """Return compartment names."""
        return ["central", "peripheral_1", "peripheral_2"]

    def _adjust_cl_by_covariates(self, cl: float, covariates: Dict[str, Any]) -> float:
        """Allometric scaling: CL_adj = CL * (weight/70)^0.75"""
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            cl = cl * (weight / 70.0) ** 0.75
        return cl

    def _adjust_vc_by_covariates(self, vc: float, covariates: Dict[str, Any]) -> float:
        """Linear scaling: Vc_adj = Vc * (weight/70)"""
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            vc = vc * (weight / 70.0)
        return vc

    def _adjust_qp_by_covariates(
        self, qp: float, covariates: Dict[str, Any], param_name: str = "Qp"
    ) -> float:
        """Allometric scaling for intercompartmental clearance."""
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            qp = qp * (weight / 70.0) ** 0.75
        return qp

    def _adjust_vp_by_covariates(
        self, vp: float, covariates: Dict[str, Any], param_name: str = "Vp"
    ) -> float:
        """Linear scaling for peripheral volume."""
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            vp = vp * (weight / 70.0)
        return vp

    def get_concentration(self, amount: float) -> float:
        """Convert central compartment amount to concentration."""
        return amount / self.params.Vc

    def get_steady_state_volume(self) -> float:
        """Calculate Vss (volume of distribution at steady state)."""
        Vp2 = self.params.extra.get("Vp2", 1.0)
        return self.params.Vc + self.params.Vp + Vp2

    def get_lambda_values(self) -> tuple:
        """
        Calculate the three exponential coefficients (lambda1, lambda2, lambda3).
        
        Returns:
            Tuple of three lambda values (eigenvalues of the system)
        """
        k_el = self.params.CL / self.params.Vc
        k_cp1 = self.params.Qp / self.params.Vc
        k_p1c = self.params.Qp / self.params.Vp
        k_cp2 = self.params.extra.get("Qp2", 0.0) / self.params.Vc
        Vp2 = self.params.extra.get("Vp2", 1.0)
        k_p2c = self.params.extra.get("Qp2", 0.0) / Vp2
        
        # Characteristic equation coefficients (simplified for 3-compartment)
        # For full analytical solution, eigenvalues of the system matrix
        a0 = k_el * k_p1c * k_p2c
        a1 = k_el * (k_p1c + k_p2c) + k_p1c * k_p2c + k_cp1 * k_p2c + k_cp2 * k_p1c
        a2 = k_el + k_p1c + k_p2c + k_cp1 + k_cp2
        
        # Cubic formula (approximate eigenvalue calculation)
        # For precise values, use numpy.linalg.eigvals
        import numpy as np
        
        # System matrix
        A = np.array([
            [-(k_el + k_cp1 + k_cp2), k_p1c, k_p2c],
            [k_cp1, -k_p1c, 0],
            [k_cp2, 0, -k_p2c]
        ])
        
        eigenvalues = np.linalg.eigvals(A)
        lambdas = np.sort(np.abs(eigenvalues))[::-1]  # sort descending
        
        return tuple(lambdas)

    def get_phase_half_lives(self) -> Dict[str, float]:
        """
        Calculate half-lives for each exponential phase.
        
        Returns:
            Dictionary with 'alpha', 'beta', 'gamma' half-lives
        """
        lambdas = self.get_lambda_values()
        
        return {
            "alpha": 0.693 / lambdas[0],  # fastest
            "beta": 0.693 / lambdas[1],   # intermediate
            "gamma": 0.693 / lambdas[2],  # slowest
        }


class ThreeCompartmentOral(PKModel):
    """
    3-compartment model with first-order absorption.
    
    Model equations:
        dA_gut/dt = - Ka * A_gut
        dA_c/dt   =  Ka * A_gut - (CL/Vc + Qp1/Vc + Qp2/Vc) * A_c + (Qp1/Vp1) * A_p1 + (Qp2/Vp2) * A_p2
        dA_p1/dt  =  (Qp1/Vc) * A_c - (Qp1/Vp1) * A_p1
        dA_p2/dt  =  (Qp2/Vc) * A_c - (Qp2/Vp2) * A_p2
    
    State vector: [A_gut, A_c, A_p1, A_p2]
    """

    def __init__(self, params: PKParameters):
        """Initialize three-compartment oral model."""
        super().__init__(params, model_type=ModelType.THREE_COMPARTMENT)
        self._validate_params()

    def _validate_params(self):
        """Validate required parameters."""
        if self.params.Ka is None:
            raise ValueError("ThreeCompartmentOral requires Ka parameter")
        if self.params.Qp is None or self.params.Vp is None:
            raise ValueError("ThreeCompartmentOral requires Qp1 and Vp1")
        if "Qp2" not in self.params.extra or "Vp2" not in self.params.extra:
            raise ValueError("ThreeCompartmentOral requires Qp2 and Vp2 in extra dict")

    def rhs(
        self,
        t: float,
        y: np.ndarray,
        covariates: Optional[Dict[str, Any]] = None,
        ddi_factors: Optional[Dict[str, float]] = None,
    ) -> np.ndarray:
        """Compute ODE right-hand side."""
        A_gut, A_c, A_p1, A_p2 = y
        
        Ka = self.params.Ka
        CL = self.params.CL
        Vc = self.params.Vc
        Qp1 = self.params.Qp
        Vp1 = self.params.Vp
        Qp2 = self.params.extra.get("Qp2", 0.0)
        Vp2 = self.params.extra.get("Vp2", 1.0)
        
        if covariates:
            Ka = self._adjust_ka_by_covariates(Ka, covariates)
            CL = self._adjust_cl_by_covariates(CL, covariates)
            Vc = self._adjust_vc_by_covariates(Vc, covariates)
            Qp1 = self._adjust_qp_by_covariates(Qp1, covariates, "Qp1")
            Vp1 = self._adjust_vp_by_covariates(Vp1, covariates, "Vp1")
            Qp2 = self._adjust_qp_by_covariates(Qp2, covariates, "Qp2")
            Vp2 = self._adjust_vp_by_covariates(Vp2, covariates, "Vp2")
        
        k_el = CL / Vc
        k_cp1 = Qp1 / Vc
        k_p1c = Qp1 / Vp1
        k_cp2 = Qp2 / Vc
        k_p2c = Qp2 / Vp2
        
        dA_gut_dt = -Ka * A_gut
        dA_c_dt = Ka * A_gut - (k_el + k_cp1 + k_cp2) * A_c + k_p1c * A_p1 + k_p2c * A_p2
        dA_p1_dt = k_cp1 * A_c - k_p1c * A_p1
        dA_p2_dt = k_cp2 * A_c - k_p2c * A_p2
        
        return np.array([dA_gut_dt, dA_c_dt, dA_p1_dt, dA_p2_dt])

    def initial_state(
        self,
        dose_mg: float,
        route: str = "oral",
        covariates: Optional[Dict[str, Any]] = None,
    ) -> np.ndarray:
        """Initial state for oral dosing."""
        route_lower = route.lower().strip()
        valid_routes = {"po", "oral", "p.o.", AdminRoute.ORAL.value}
        
        if route_lower not in valid_routes:
            raise ValueError(f"ThreeCompartmentOral expects oral route, got '{route}'")
        
        available_dose = dose_mg * self.params.F
        return np.array([available_dose, 0.0, 0.0, 0.0])

    def get_compartment_names(self) -> List[str]:
        """Return compartment names."""
        return ["gut", "central", "peripheral_1", "peripheral_2"]

    def _adjust_ka_by_covariates(self, ka: float, covariates: Dict[str, Any]) -> float:
        """Adjust Ka for fed state."""
        if "fed_state" in covariates and covariates["fed_state"]:
            ka = ka * 0.7
        return ka

    def _adjust_cl_by_covariates(self, cl: float, covariates: Dict[str, Any]) -> float:
        """Allometric scaling for CL."""
        if "weight_kg" in covariates:
            cl = cl * (covariates["weight_kg"] / 70.0) ** 0.75
        return cl

    def _adjust_vc_by_covariates(self, vc: float, covariates: Dict[str, Any]) -> float:
        """Linear scaling for Vc."""
        if "weight_kg" in covariates:
            vc = vc * (covariates["weight_kg"] / 70.0)
        return vc

    def _adjust_qp_by_covariates(
        self, qp: float, covariates: Dict[str, Any], param_name: str = "Qp"
    ) -> float:
        """Allometric scaling for Qp."""
        if "weight_kg" in covariates:
            qp = qp * (covariates["weight_kg"] / 70.0) ** 0.75
        return qp

    def _adjust_vp_by_covariates(
        self, vp: float, covariates: Dict[str, Any], param_name: str = "Vp"
    ) -> float:
        """Linear scaling for Vp."""
        if "weight_kg" in covariates:
            vp = vp * (covariates["weight_kg"] / 70.0)
        return vp

    def get_concentration(self, amount: float) -> float:
        """Convert central compartment amount to concentration."""
        return amount / self.params.Vc

    def get_steady_state_volume(self) -> float:
        """Calculate Vss."""
        Vp2 = self.params.extra.get("Vp2", 1.0)
        return self.params.Vc + self.params.Vp + Vp2

    def get_absorption_half_life(self) -> float:
        """Calculate absorption half-life."""
        return 0.693 / self.params.Ka
