"""
One-compartment pharmacokinetic models with first-order absorption.

Supports oral and other extravascular routes with Ka-driven absorption.
"""

from __future__ import annotations
from typing import Dict, Any, Optional, List
import numpy as np
import logging

from .base import PKModel, PKParameters, ModelType, AdminRoute

logger = logging.getLogger(__name__)


class OralOneCompartment(PKModel):
    """
    1-compartment model with first-order absorption from GI tract.
    
    Model equations:
        dA_gut/dt = - Ka * A_gut
        dA_c/dt   =  Ka * A_gut - (CL / Vc) * A_c
        C_c = A_c / Vc
    
    where:
        A_gut: amount in GI tract (absorption compartment) (mg)
        A_c: amount in central compartment (mg)
        C_c: concentration in central compartment (mg/L)
        Ka: first-order absorption rate constant (1/h)
        CL: clearance (L/h)
        Vc: central volume of distribution (L)
        F: bioavailability (0-1)
    
    Typical uses: Oral tablets, capsules, solutions with first-order absorption
    State vector: [A_gut, A_c]
    """

    def __init__(self, params: PKParameters):
        """Initialize oral one-compartment model."""
        super().__init__(params, model_type=ModelType.ONE_COMPARTMENT)
        self._validate_params()

    def _validate_params(self):
        """Ensure required parameters are present for oral model."""
        if self.params.Ka is None:
            raise ValueError(
                "OralOneCompartment requires Ka parameter to be set. "
                "Ka (absorption rate constant) is essential for oral models."
            )
        if self.params.Ka <= 0:
            raise ValueError(f"Ka must be > 0, got {self.params.Ka}")

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
            y: State vector [A_gut, A_c] (amounts in mg)
            covariates: Optional covariate adjustments
            ddi_factors: Optional DDI factors (applied upstream)
            
        Returns:
            dy/dt: [dA_gut/dt, dA_c/dt] (mg/h)
        """
        A_gut, A_c = y
        
        # Get parameters with covariate adjustments
        Ka = self.params.Ka
        CL = self.params.CL
        Vc = self.params.Vc
        
        if covariates:
            Ka = self._adjust_ka_by_covariates(Ka, covariates)
            CL = self._adjust_cl_by_covariates(CL, covariates)
            Vc = self._adjust_vc_by_covariates(Vc, covariates)
        
        # ODE system
        k_el = CL / Vc  # elimination rate constant (1/h)
        
        dA_gut_dt = -Ka * A_gut
        dA_c_dt = Ka * A_gut - k_el * A_c
        
        return np.array([dA_gut_dt, dA_c_dt])

    def initial_state(
        self,
        dose_mg: float,
        route: str = "oral",
        covariates: Optional[Dict[str, Any]] = None,
    ) -> np.ndarray:
        """
        Compute initial state for given dose and route.
        
        For oral dosing: entire dose (adjusted for F) enters GI tract compartment.
        
        Args:
            dose_mg: Dose in mg (nominal dose)
            route: Administration route (oral, po, etc.)
            covariates: Optional covariate adjustments
            
        Returns:
            Initial state vector y0 = [A_gut, A_c] (mg)
            
        Raises:
            ValueError: If route is not oral/po
        """
        route_lower = route.lower().strip()
        
        valid_routes = {"po", "oral", "p.o.", AdminRoute.ORAL.value}
        if route_lower not in valid_routes:
            raise ValueError(
                f"OralOneCompartment expects oral route, got '{route}'. "
                f"Valid routes: {valid_routes}. "
                f"For IV administration, use OneCompartmentIV."
            )
        
        # Apply bioavailability to dose
        available_dose = dose_mg * self.params.F
        
        # All available dose enters GI tract
        A_gut0 = available_dose
        A_c0 = 0.0  # Central compartment starts empty
        
        return np.array([A_gut0, A_c0])

    def get_compartment_names(self) -> List[str]:
        """Return compartment names in state vector order."""
        return ["gut", "central"]

    def _adjust_ka_by_covariates(self, ka: float, covariates: Dict[str, Any]) -> float:
        """
        Adjust Ka by covariates (e.g., GI motility, pH, food).
        
        Example: fed state can reduce Ka
        """
        if "fed_state" in covariates:
            # Approximate: fed reduces Ka by ~30%
            fed = covariates["fed_state"]
            if fed:
                ka = ka * 0.7
        
        return ka

    def _adjust_cl_by_covariates(self, cl: float, covariates: Dict[str, Any]) -> float:
        """
        Adjust clearance by covariates (allometric scaling).
        
        CL_adj = CL * (weight/70)^0.75
        """
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            cl = cl * (weight / 70.0) ** 0.75
        
        return cl

    def _adjust_vc_by_covariates(self, vc: float, covariates: Dict[str, Any]) -> float:
        """
        Adjust volume of distribution by covariates.
        
        Vc_adj = Vc * (weight/70)
        """
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            vc = vc * (weight / 70.0)
        
        return vc

    def get_concentration(self, amount: float) -> float:
        """Convert central compartment amount to concentration (mg/L)."""
        return amount / self.params.Vc

    def get_half_life(self) -> float:
        """Calculate elimination half-life (hours)."""
        k_el = self.params.CL / self.params.Vc
        return 0.693 / k_el

    def get_absorption_half_life(self) -> float:
        """Calculate absorption half-life (hours)."""
        if self.params.Ka is None:
            raise ValueError("Ka not set")
        return 0.693 / self.params.Ka

    def get_tmax_theoretical(self) -> float:
        """
        Calculate theoretical time to maximum concentration (Tmax) in hours.
        
        For 1-compartment first-order absorption:
        Tmax = ln(Ka / k_el) / (Ka - k_el)
        where k_el = CL / Vc
        
        Only valid when Ka != k_el.
        """
        if self.params.Ka is None:
            raise ValueError("Ka not set")
        
        k_el = self.params.CL / self.params.Vc
        
        if abs(self.params.Ka - k_el) < 1e-6:
            logger.warning("Ka ≈ k_el; Tmax calculation may be unstable")
            return 0.0
        
        tmax = np.log(self.params.Ka / k_el) / (self.params.Ka - k_el)
        return max(0.0, float(tmax))  # ensure non-negative

    def get_cmax_theoretical(self, dose_mg: float, covariates: Optional[Dict[str, Any]] = None) -> float:
        """
        Calculate theoretical maximum concentration (Cmax) in mg/L.
        
        For 1-compartment first-order absorption:
        Cmax = (F * dose * Ka) / (Vc * (Ka - k_el))
        where k_el = CL / Vc
        """
        if self.params.Ka is None:
            raise ValueError("Ka not set")
        
        Ka = self.params.Ka
        CL = self.params.CL
        Vc = self.params.Vc
        
        if covariates:
            Ka = self._adjust_ka_by_covariates(Ka, covariates)
            CL = self._adjust_cl_by_covariates(CL, covariates)
            Vc = self._adjust_vc_by_covariates(Vc, covariates)
        
        k_el = CL / Vc
        
        if abs(Ka - k_el) < 1e-6:
            logger.warning("Ka ≈ k_el; Cmax calculation may be unstable")
            return 0.0
        
        cmax = (self.params.F * dose_mg * Ka) / (Vc * (Ka - k_el))
        return max(0.0, float(cmax))

    def get_bioavailability(self) -> float:
        """Get bioavailability factor (0-1)."""
        return self.params.F
