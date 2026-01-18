"""
Two-compartment pharmacokinetic models.

Supports IV bolus and first-order absorption with intercompartmental distribution.
"""

from __future__ import annotations
from typing import Dict, Any, Optional, List
import numpy as np
import logging

from .base import PKModel, PKParameters, ModelType, AdminRoute

logger = logging.getLogger(__name__)


class TwoCompartmentIV(PKModel):
    """
    2-compartment IV bolus model with rapid central and slower peripheral kinetics.
    
    Model equations:
        dA_c/dt = - (CL/Vc + Qp/Vc) * A_c + (Qp/Vp) * A_p
        dA_p/dt =   (Qp/Vc) * A_c - (Qp/Vp) * A_p
        C_c = A_c / Vc
    
    where:
        A_c: amount in central compartment (mg)
        A_p: amount in peripheral compartment (mg)
        C_c: concentration in central (mg/L)
        CL: clearance (L/h, occurs from central)
        Vc: central volume (L)
        Qp: intercompartmental clearance (L/h)
        Vp: peripheral volume (L)
    
    Typical uses: IV bolus with extended elimination, lipophilic drugs
    State vector: [A_c, A_p]
    """

    def __init__(self, params: PKParameters):
        """Initialize two-compartment IV model."""
        super().__init__(params, model_type=ModelType.TWO_COMPARTMENT)
        self._validate_params()

    def _validate_params(self):
        """Ensure required parameters are present for two-compartment model."""
        if self.params.Qp is None or self.params.Vp is None:
            raise ValueError(
                "TwoCompartmentIV requires both Qp and Vp parameters. "
                f"Got Qp={self.params.Qp}, Vp={self.params.Vp}"
            )
        if self.params.Qp < 0:
            raise ValueError(f"Qp must be >= 0, got {self.params.Qp}")
        if self.params.Vp <= 0:
            raise ValueError(f"Vp must be > 0, got {self.params.Vp}")

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
            y: State vector [A_c, A_p] (amounts in mg)
            covariates: Optional covariate adjustments
            ddi_factors: Optional DDI factors (applied upstream)
            
        Returns:
            dy/dt: [dA_c/dt, dA_p/dt] (mg/h)
        """
        A_c, A_p = y
        
        # Get parameters with covariate adjustments
        CL = self.params.CL
        Vc = self.params.Vc
        Qp = self.params.Qp
        Vp = self.params.Vp
        
        if covariates:
            CL = self._adjust_cl_by_covariates(CL, covariates)
            Vc = self._adjust_vc_by_covariates(Vc, covariates)
            Qp = self._adjust_qp_by_covariates(Qp, covariates)
            Vp = self._adjust_vp_by_covariates(Vp, covariates)
        
        # Rate constants
        k_el = CL / Vc
        k_cp = Qp / Vc  # central to peripheral
        k_pc = Qp / Vp  # peripheral to central
        
        # ODE system
        dA_c_dt = -(k_el + k_cp) * A_c + k_pc * A_p
        dA_p_dt = k_cp * A_c - k_pc * A_p
        
        return np.array([dA_c_dt, dA_p_dt])

    def initial_state(
        self,
        dose_mg: float,
        route: str = "iv",
        covariates: Optional[Dict[str, Any]] = None,
    ) -> np.ndarray:
        """
        Compute initial state for IV bolus.
        
        For IV: entire dose enters central compartment.
        
        Args:
            dose_mg: Dose in mg
            route: Administration route (must be 'iv')
            covariates: Optional covariate adjustments
            
        Returns:
            Initial state vector y0 = [A_c, A_p] (mg)
            
        Raises:
            ValueError: If route is not IV
        """
        route_lower = route.lower().strip()
        
        if route_lower not in ["iv", AdminRoute.IV.value]:
            raise ValueError(
                f"TwoCompartmentIV only supports IV route, got '{route}'. "
                f"For oral, use TwoCompartmentOral."
            )
        
        initial_amount = dose_mg * self.params.F
        
        return np.array([initial_amount, 0.0])

    def get_compartment_names(self) -> List[str]:
        """Return compartment names."""
        return ["central", "peripheral"]

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

    def _adjust_qp_by_covariates(self, qp: float, covariates: Dict[str, Any]) -> float:
        """Allometric scaling for intercompartmental clearance: Qp_adj = Qp * (weight/70)^0.75"""
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            qp = qp * (weight / 70.0) ** 0.75
        return qp

    def _adjust_vp_by_covariates(self, vp: float, covariates: Dict[str, Any]) -> float:
        """Linear scaling: Vp_adj = Vp * (weight/70)"""
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            vp = vp * (weight / 70.0)
        return vp

    def get_concentration(self, amount: float) -> float:
        """Convert central compartment amount to concentration."""
        return amount / self.params.Vc

    def get_distribution_half_life(self) -> float:
        """
        Calculate half-life of distribution phase (alpha phase).
        
        For 2-compartment: t_1/2,alpha is rapid
        """
        # Alpha half-life (rapid phase)
        k_el = self.params.CL / self.params.Vc
        k_cp = self.params.Qp / self.params.Vc
        k_pc = self.params.Qp / self.params.Vp
        
        lambda_alpha = (k_el + k_cp + k_pc + np.sqrt((k_el + k_cp + k_pc)**2 - 4*k_el*k_pc)) / 2
        return 0.693 / lambda_alpha

    def get_elimination_half_life(self) -> float:
        """
        Calculate half-life of elimination phase (beta phase).
        
        For 2-compartment: t_1/2,beta is slower
        """
        # Beta half-life (slow phase)
        k_el = self.params.CL / self.params.Vc
        k_cp = self.params.Qp / self.params.Vc
        k_pc = self.params.Qp / self.params.Vp
        
        lambda_beta = (k_el + k_cp + k_pc - np.sqrt((k_el + k_cp + k_pc)**2 - 4*k_el*k_pc)) / 2
        return 0.693 / lambda_beta

    def get_steady_state_volume(self) -> float:
        """
        Calculate steady-state volume of distribution (Vss).
        
        Vss = Vc + Vp (at equilibrium)
        """
        return self.params.Vc + self.params.Vp

    def get_alpha_coefficient(self) -> float:
        """Get alpha coefficient for two-compartment equation."""
        k_el = self.params.CL / self.params.Vc
        k_cp = self.params.Qp / self.params.Vc
        k_pc = self.params.Qp / self.params.Vp
        return (k_el + k_cp + k_pc + np.sqrt((k_el + k_cp + k_pc)**2 - 4*k_el*k_pc)) / 2

    def get_beta_coefficient(self) -> float:
        """Get beta coefficient for two-compartment equation."""
        k_el = self.params.CL / self.params.Vc
        k_cp = self.params.Qp / self.params.Vc
        k_pc = self.params.Qp / self.params.Vp
        return (k_el + k_cp + k_pc - np.sqrt((k_el + k_cp + k_pc)**2 - 4*k_el*k_pc)) / 2


class TwoCompartmentOral(PKModel):
    """
    2-compartment model with first-order absorption.
    
    Model equations:
        dA_gut/dt = - Ka * A_gut
        dA_c/dt   =  Ka * A_gut - (CL/Vc + Qp/Vc) * A_c + (Qp/Vp) * A_p
        dA_p/dt   =  (Qp/Vc) * A_c - (Qp/Vp) * A_p
        C_c = A_c / Vc
    
    State vector: [A_gut, A_c, A_p]
    """

    def __init__(self, params: PKParameters):
        """Initialize two-compartment oral model."""
        super().__init__(params, model_type=ModelType.TWO_COMPARTMENT)
        self._validate_params()

    def _validate_params(self):
        """Validate required parameters."""
        if self.params.Ka is None:
            raise ValueError("TwoCompartmentOral requires Ka parameter")
        if self.params.Qp is None or self.params.Vp is None:
            raise ValueError("TwoCompartmentOral requires Qp and Vp parameters")

    def rhs(
        self,
        t: float,
        y: np.ndarray,
        covariates: Optional[Dict[str, Any]] = None,
        ddi_factors: Optional[Dict[str, float]] = None,
    ) -> np.ndarray:
        """Compute ODE right-hand side."""
        A_gut, A_c, A_p = y
        
        Ka = self.params.Ka
        CL = self.params.CL
        Vc = self.params.Vc
        Qp = self.params.Qp
        Vp = self.params.Vp
        
        if covariates:
            Ka = self._adjust_ka_by_covariates(Ka, covariates)
            CL = self._adjust_cl_by_covariates(CL, covariates)
            Vc = self._adjust_vc_by_covariates(Vc, covariates)
            Qp = self._adjust_qp_by_covariates(Qp, covariates)
            Vp = self._adjust_vp_by_covariates(Vp, covariates)
        
        k_el = CL / Vc
        k_cp = Qp / Vc
        k_pc = Qp / Vp
        
        dA_gut_dt = -Ka * A_gut
        dA_c_dt = Ka * A_gut - (k_el + k_cp) * A_c + k_pc * A_p
        dA_p_dt = k_cp * A_c - k_pc * A_p
        
        return np.array([dA_gut_dt, dA_c_dt, dA_p_dt])

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
            raise ValueError(f"TwoCompartmentOral expects oral route, got '{route}'")
        
        available_dose = dose_mg * self.params.F
        return np.array([available_dose, 0.0, 0.0])

    def get_compartment_names(self) -> List[str]:
        """Return compartment names."""
        return ["gut", "central", "peripheral"]

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

    def _adjust_qp_by_covariates(self, qp: float, covariates: Dict[str, Any]) -> float:
        """Allometric scaling for Qp."""
        if "weight_kg" in covariates:
            qp = qp * (covariates["weight_kg"] / 70.0) ** 0.75
        return qp

    def _adjust_vp_by_covariates(self, vp: float, covariates: Dict[str, Any]) -> float:
        """Linear scaling for Vp."""
        if "weight_kg" in covariates:
            vp = vp * (covariates["weight_kg"] / 70.0)
        return vp

    def get_concentration(self, amount: float) -> float:
        """Convert central compartment amount to concentration."""
        return amount / self.params.Vc

    def get_absorption_half_life(self) -> float:
        """Calculate absorption half-life."""
        return 0.693 / self.params.Ka

    def get_elimination_half_life(self) -> float:
        """Calculate beta (elimination) half-life."""
        k_el = self.params.CL / self.params.Vc
        k_cp = self.params.Qp / self.params.Vc
        k_pc = self.params.Qp / self.params.Vp
        
        lambda_beta = (k_el + k_cp + k_pc - np.sqrt((k_el + k_cp + k_pc)**2 - 4*k_el*k_pc)) / 2
        return 0.693 / lambda_beta

    def get_steady_state_volume(self) -> float:
        """Calculate Vss."""
        return self.params.Vc + self.params.Vp
