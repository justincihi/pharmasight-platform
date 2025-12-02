"""
One-compartment pharmacokinetic models.

Supports IV bolus administration with optional covariate adjustments.
"""

from __future__ import annotations
from typing import Dict, Any, Optional, List
import numpy as np
import logging

from .base import PKModel, PKParameters, ModelType, AdminRoute

logger = logging.getLogger(__name__)


class OneCompartmentIV(PKModel):
    """
    Classic 1-compartment IV bolus model.
    
    Model equation:
        dA_c/dt = - (CL / Vc) * A_c
        C_c = A_c / Vc
    
    where:
        A_c: amount in central compartment (mg)
        C_c: concentration in central compartment (mg/L)
        CL: clearance (L/h)
        Vc: central volume of distribution (L)
    
    Typical uses: IV bolus drugs with rapid distribution
    """

    def __init__(self, params: PKParameters):
        """Initialize IV one-compartment model."""
        super().__init__(params, model_type=ModelType.ONE_COMPARTMENT)
        self._validate_params()

    def _validate_params(self):
        """Ensure required parameters are present for IV model."""
        if self.params.Ka is not None:
            logger.warning("OneCompartmentIV ignores Ka parameter (set for IV route)")

    def rhs(
        self,
        t: float,
        y: np.ndarray,
        covariates: Optional[Dict[str, Any]] = None,
        ddi_factors: Optional[Dict[str, float]] = None,
    ) -> np.ndarray:
        """
        Compute right-hand side of ODE.
        
        Args:
            t: Current time (h)
            y: State vector [A_c] (amount in central compartment, mg)
            covariates: Optional covariate adjustments
            ddi_factors: Optional DDI factors (applied upstream via apply_ddi)
            
        Returns:
            dy/dt: [dA_c/dt] (mg/h)
        """
        A_c = y[0]
        
        # Apply covariate adjustments if provided
        CL = self.params.CL
        Vc = self.params.Vc
        
        if covariates:
            CL = self._adjust_cl_by_covariates(CL, covariates)
            Vc = self._adjust_vc_by_covariates(Vc, covariates)
        
        k_el = CL / Vc  # elimination rate constant (1/h)
        dA_c_dt = -k_el * A_c
        
        return np.array([dA_c_dt])

    def initial_state(
        self,
        dose_mg: float,
        route: str = "iv",
        covariates: Optional[Dict[str, Any]] = None,
    ) -> np.ndarray:
        """
        Compute initial state for given dose and route.
        
        For IV bolus: entire dose enters central compartment instantaneously.
        
        Args:
            dose_mg: Dose in mg
            route: Administration route (must be 'iv')
            covariates: Optional covariate adjustments
            
        Returns:
            Initial state vector y0 = [A_c] (mg)
            
        Raises:
            ValueError: If route is not IV
        """
        route_lower = route.lower().strip()
        
        if route_lower not in ["iv", AdminRoute.IV.value]:
            raise ValueError(
                f"OneCompartmentIV only supports IV route, got '{route}'. "
                f"For oral administration, use OneCompartmentOral."
            )
        
        # Apply bioavailability (typically 1.0 for IV)
        initial_amount = dose_mg * self.params.F
        
        return np.array([initial_amount])

    def get_compartment_names(self) -> List[str]:
        """Return compartment names."""
        return ["central"]

    def _adjust_cl_by_covariates(self, cl: float, covariates: Dict[str, Any]) -> float:
        """
        Adjust clearance by covariates (e.g., weight, age, renal function).
        
        Simple allometric scaling: CL_adj = CL * (weight/70)^0.75
        """
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            cl = cl * (weight / 70.0) ** 0.75
        
        return cl

    def _adjust_vc_by_covariates(self, vc: float, covariates: Dict[str, Any]) -> float:
        """
        Adjust volume of distribution by covariates (e.g., weight, body composition).
        
        Simple linear scaling: Vc_adj = Vc * (weight/70)
        """
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            vc = vc * (weight / 70.0)
        
        return vc

    def get_concentration(self, amount: float) -> float:
        """Convert amount to concentration."""
        return amount / self.params.Vc

    def get_half_life(self) -> float:
        """Calculate elimination half-life (hours)."""
        k_el = self.params.CL / self.params.Vc
        return 0.693 / k_el

    def get_clearance_unit_adjusted(self, covariates: Optional[Dict[str, Any]] = None) -> float:
        """Get clearance with optional covariate adjustment."""
        cl = self.params.CL
        if covariates:
            cl = self._adjust_cl_by_covariates(cl, covariates)
        return cl
