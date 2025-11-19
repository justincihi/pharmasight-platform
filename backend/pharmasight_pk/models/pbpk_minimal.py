"""
Physiologically-Based Pharmacokinetic (PBPK) models.

Supports organ-based models with tissue compartments, blood flow-based distribution,
and CYP450-mediated metabolism.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Tuple
import numpy as np
import logging

from .base import PKModel, PKParameters, ModelType, AdminRoute

logger = logging.getLogger(__name__)


class MetabolismType:
    """Enumeration of supported metabolism types."""
    LINEAR = "linear"
    MICHAELIS_MENTEN = "michaelis_menten"
    CYP450 = "cyp450"


@dataclass
class CYP450Info:
    """Information about CYP450 enzyme involvement."""
    cyp_name: str  # e.g., "CYP3A4", "CYP2D6"
    km: float  # Michaelis constant (mg/L)
    vmax: float  # Maximum velocity (mg/h)
    relative_activity: float = 1.0  # Activity factor (0-1 for inhibition, >1 for induction)
    
    def __post_init__(self):
        if self.km < 0:
            raise ValueError(f"Km must be >= 0, got {self.km}")
        if self.vmax < 0:
            raise ValueError(f"Vmax must be >= 0, got {self.vmax}")
        if self.relative_activity < 0:
            raise ValueError(f"relative_activity must be >= 0, got {self.relative_activity}")


@dataclass
class PBPKCompartment:
    """
    Represents a tissue/organ compartment in PBPK model.
    
    Attributes:
        name: Compartment identifier (e.g., "liver", "kidney", "heart")
        volume: Tissue volume in liters
        blood_flow: Tissue blood flow rate in L/h
        partition_coefficient: Tissue-to-plasma concentration ratio
        is_elimination_site: Whether drug is metabolized/eliminated here
        metabolism_type: Type of metabolism (linear, Michaelis-Menten, CYP450)
        cyp450_info: CYP450 enzyme information (if applicable)
        relative_activity: Baseline enzymatic activity (0-1 scale)
    """
    name: str
    volume: float  # L
    blood_flow: float  # L/h
    partition_coefficient: float  # tissue/plasma ratio
    is_elimination_site: bool = False
    metabolism_type: str = MetabolismType.LINEAR
    cyp450_info: Optional[CYP450Info] = None
    relative_activity: float = 1.0  # For activity scaling (e.g., age, genetic polymorphisms)
    
    def __post_init__(self):
        if self.volume <= 0:
            raise ValueError(f"Compartment {self.name}: volume must be > 0, got {self.volume}")
        if self.blood_flow < 0:
            raise ValueError(f"Compartment {self.name}: blood_flow must be >= 0, got {self.blood_flow}")
        if self.partition_coefficient <= 0:
            raise ValueError(f"Compartment {self.name}: partition_coefficient must be > 0, got {self.partition_coefficient}")
        if not (0 <= self.relative_activity <= 2.0):
            logger.warning(f"Compartment {self.name}: relative_activity {self.relative_activity} outside typical range [0, 2]")


@dataclass
class PBPKStructure:
    """
    Physiologically-based compartment structure.
    
    Implements well-stirred tissue model with blood flow-based distribution
    and organ-specific elimination pathways.
    
    Attributes:
        compartments: List of tissue/organ compartments
        cardiac_output: Total body blood flow (L/h)
        plasma_volume: Volume of central plasma compartment (L)
        metabolism_scaling: Factor to scale clearance for population differences
    """
    compartments: List[PBPKCompartment] = field(default_factory=list)
    cardiac_output: float = 75.0  # L/h (typical: ~72 L/h at rest)
    plasma_volume: float = 3.0  # L (typical adult)
    metabolism_scaling: float = 1.0  # Population scaling
    
    def __post_init__(self):
        if self.cardiac_output <= 0:
            raise ValueError(f"cardiac_output must be > 0, got {self.cardiac_output}")
        if self.plasma_volume <= 0:
            raise ValueError(f"plasma_volume must be > 0, got {self.plasma_volume}")
        if self.metabolism_scaling <= 0:
            raise ValueError(f"metabolism_scaling must be > 0, got {self.metabolism_scaling}")
        
        # Validate total blood flow doesn't exceed cardiac output
        total_flow = sum(c.blood_flow for c in self.compartments)
        if total_flow > self.cardiac_output * 1.01:  # Allow 1% tolerance
            logger.warning(
                f"Total tissue blood flow ({total_flow} L/h) exceeds cardiac output "
                f"({self.cardiac_output} L/h). Consider adjusting."
            )
    
    def get_elimination_compartments(self) -> List[PBPKCompartment]:
        """Return all compartments where elimination occurs."""
        return [c for c in self.compartments if c.is_elimination_site]
    
    def get_compartment_by_name(self, name: str) -> Optional[PBPKCompartment]:
        """Retrieve compartment by name."""
        for c in self.compartments:
            if c.name.lower() == name.lower():
                return c
        return None


class MinimalPBPK(PKModel):
    """
    Simplified Physiologically-Based Pharmacokinetic (PBPK) model.
    
    Features:
    - Well-stirred tissue model with blood flow-driven distribution
    - Multi-organ compartments (liver, kidney, adipose, muscle, etc.)
    - CYP450-mediated and linear elimination pathways
    - Michaelis-Menten metabolism support
    - Covariate adjustments for organ function
    - Plasma and tissue distribution
    
    State vector: [A_plasma, A_comp1, A_comp2, ..., A_compN]
    where A_i is the amount (mg) in each compartment.
    
    Model assumptions:
    - Rapid venous equilibration (well-mixed central plasma)
    - Tissue equilibration follows first-order rate = Q*(Cp - Ct/Pc)
    - Linear or saturable metabolism in elimination sites
    - No active transport (passive distribution)
    """

    def __init__(self, params: PKParameters, structure: PBPKStructure):
        """
        Initialize PBPK model.
        
        Args:
            params: Base PK parameters (CL, Vc, F, etc.)
            structure: PBPK compartment structure with organ definitions
        """
        super().__init__(params, model_type=ModelType.PBPK)
        self.structure = structure
        
        # Map compartment names to indices (0 = plasma)
        self._name_to_idx = {
            "plasma": 0,
            **{c.name: i + 1 for i, c in enumerate(self.structure.compartments)}
        }
        
        self._validate_params()
        logger.info(f"Initialized PBPK model with {len(self.structure.compartments)} tissue compartments")

    def _validate_params(self):
        """Validate PBPK-specific parameters."""
        if len(self.structure.compartments) == 0:
            raise ValueError("PBPK structure must contain at least one tissue compartment")

    def rhs(
        self,
        t: float,
        y: np.ndarray,
        covariates: Optional[Dict[str, Any]] = None,
        ddi_factors: Optional[Dict[str, float]] = None,
    ) -> np.ndarray:
        """
        Compute right-hand side of the PBPK ODE system.
        
        Implements blood flow-limited tissue distribution and organ-specific elimination.
        
        Args:
            t: Current time (h)
            y: State vector [A_plasma, A_tissue1, A_tissue2, ..., A_tissueN] (mg)
            covariates: Optional covariate adjustments (weight, age, organ function, etc.)
            ddi_factors: Optional DDI factors (applied upstream)
            
        Returns:
            dy/dt: Rate of change for all compartments (mg/h)
        """
        p = self.params
        s = self.structure
        
        n_comps = len(s.compartments) + 1  # +1 for plasma
        if len(y) != n_comps:
            raise ValueError(f"State vector size mismatch: expected {n_comps}, got {len(y)}")
        
        # Calculate concentrations in each compartment
        C = np.zeros(n_comps)
        
        # Plasma concentration
        C[0] = y[0] / s.plasma_volume
        
        # Tissue concentrations (accounting for partition coefficients)
        for i, comp in enumerate(s.compartments):
            C[i + 1] = y[i + 1] / comp.volume
        
        # Initialize rate vector
        dA_dt = np.zeros(n_comps)
        
        # Apply covariate adjustments
        CL = self._adjust_cl_by_covariates(p.CL, covariates)
        cardiac_output = self.structure.cardiac_output
        
        if covariates:
            cardiac_output = self._adjust_cardiac_output_by_covariates(cardiac_output, covariates)
        
        # Plasma dynamics: input from GI, distribution to tissues, elimination
        C_plasma = C[0]
        
        for i, comp in enumerate(s.compartments):
            tissue_idx = i + 1
            C_tissue = C[tissue_idx] / max(comp.partition_coefficient, 1e-6)
            
            # Blood flow-based exchange (well-stirred model)
            Q = comp.blood_flow
            
            # Influx from plasma to tissue
            # dA_tissue/dt += Q * (C_plasma - C_tissue)
            flux_plasma_to_tissue = Q * (C_plasma - C_tissue)
            
            dA_dt[tissue_idx] += flux_plasma_to_tissue
            dA_dt[0] -= flux_plasma_to_tissue
            
            # Elimination/metabolism in organ
            if comp.is_elimination_site:
                dA_dt[tissue_idx] -= self._calculate_elimination_rate(
                    A=y[tissue_idx],
                    C_tissue=C_tissue,
                    compartment=comp,
                    covariates=covariates
                )
        
        # Renal clearance from plasma (if present)
        if CL > 0:
            renal_elimination = (CL / s.plasma_volume) * y[0]
            dA_dt[0] -= renal_elimination
        
        return dA_dt

    def _calculate_elimination_rate(
        self,
        A: float,
        C_tissue: float,
        compartment: PBPKCompartment,
        covariates: Optional[Dict[str, Any]] = None
    ) -> float:
        """
        Calculate elimination/metabolism rate for a tissue compartment.
        
        Supports linear, Michaelis-Menten, and CYP450-based mechanisms.
        
        Args:
            A: Amount in tissue (mg)
            C_tissue: Concentration in tissue (mg/L)
            compartment: Tissue compartment definition
            covariates: Optional covariates for activity adjustments
            
        Returns:
            Elimination rate (mg/h)
        """
        if compartment.metabolism_type == MetabolismType.LINEAR:
            # Linear clearance: dA/dt = -CLint * C
            # Intrinsic clearance distributed across elimination sites
            n_elim_sites = len(self.structure.get_elimination_compartments())
            if n_elim_sites > 0:
                clint = (self.params.CL / n_elim_sites) / compartment.volume
                activity = compartment.relative_activity
                
                if covariates:
                    activity = self._adjust_enzymatic_activity(activity, compartment, covariates)
                
                return activity * clint * A
        
        elif compartment.metabolism_type == MetabolismType.MICHAELIS_MENTEN:
            # Michaelis-Menten: dA/dt = -(Vmax * C) / (Km + C)
            if compartment.cyp450_info:
                cyp_info = compartment.cyp450_info
                vmax = cyp_info.vmax * compartment.relative_activity
                km = cyp_info.km
                
                if covariates:
                    vmax = self._adjust_enzymatic_activity(vmax, compartment, covariates)
                
                return (vmax * C_tissue) / (km + C_tissue + 1e-12)
        
        elif compartment.metabolism_type == MetabolismType.CYP450:
            # CYP450 metabolism with potential for inhibition/induction
            if compartment.cyp450_info:
                cyp_info = compartment.cyp450_info
                vmax = cyp_info.vmax
                km = cyp_info.km
                activity = cyp_info.relative_activity * compartment.relative_activity
                
                if covariates:
                    activity = self._adjust_enzymatic_activity(activity, compartment, covariates)
                
                # Apply activity modulation
                vmax_effective = vmax * activity
                
                return (vmax_effective * C_tissue) / (km + C_tissue + 1e-12)
        
        return 0.0

    def initial_state(
        self,
        dose_mg: float,
        route: str = "iv",
        covariates: Optional[Dict[str, Any]] = None,
    ) -> np.ndarray:
        """
        Compute initial state based on dose and route.
        
        Args:
            dose_mg: Dose in mg
            route: Administration route (iv, oral, etc.)
            covariates: Optional covariates
            
        Returns:
            Initial state vector [A_plasma, A_tissue1, ..., A_tissueN]
        """
        n_comps = len(self.structure.compartments) + 1
        y0 = np.zeros(n_comps)
        
        route_lower = route.lower().strip()
        
        if route_lower in ["iv", AdminRoute.IV.value]:
            # IV bolus: dose enters plasma immediately
            y0[0] = dose_mg * self.params.F
        
        elif route_lower in ["oral", "po", "p.o.", AdminRoute.ORAL.value]:
            # Oral: dose enters GI tract (would need absorption compartment)
            # For now, enter as available in plasma
            y0[0] = dose_mg * self.params.F
            logger.info("Oral route: dose available in plasma (GI absorption not explicitly modeled)")
        
        else:
            raise ValueError(f"Unsupported route: {route}")
        
        return y0

    def get_compartment_names(self) -> List[str]:
        """Return names of all compartments in order."""
        return ["plasma"] + [c.name for c in self.structure.compartments]

    def get_tissue_concentration(self, tissue_name: str, amount: float) -> float:
        """Convert tissue amount to concentration."""
        compartment = self.structure.get_compartment_by_name(tissue_name)
        if compartment is None:
            raise ValueError(f"Tissue '{tissue_name}' not found in structure")
        return amount / compartment.volume

    def get_plasma_concentration(self, amount: float) -> float:
        """Convert plasma amount to concentration."""
        return amount / self.structure.plasma_volume

    def _adjust_cl_by_covariates(self, cl: float, covariates: Optional[Dict[str, Any]]) -> float:
        """
        Adjust clearance by covariates.
        
        Supports weight-based allometric scaling and liver/kidney function.
        """
        if not covariates:
            return cl
        
        # Weight-based scaling
        if "weight_kg" in covariates:
            weight = covariates["weight_kg"]
            cl = cl * (weight / 70.0) ** 0.75
        
        # Liver function adjustment
        if "liver_function" in covariates:
            liver_func = covariates["liver_function"]  # 0-1, where 1 = normal
            cl = cl * liver_func
        
        # Kidney function adjustment (GFR)
        if "kidney_function" in covariates:
            kidney_func = covariates["kidney_function"]  # 0-1, where 1 = normal
            cl = cl * kidney_func
        
        return cl

    def _adjust_cardiac_output_by_covariates(
        self, cardiac_output: float, covariates: Optional[Dict[str, Any]]
    ) -> float:
        """Adjust cardiac output by covariates (e.g., heart failure)."""
        if not covariates:
            return cardiac_output
        
        if "cardiac_output_factor" in covariates:
            factor = covariates["cardiac_output_factor"]  # 0-1, where 1 = normal
            cardiac_output = cardiac_output * factor
        
        return cardiac_output

    def _adjust_enzymatic_activity(
        self,
        base_activity: float,
        compartment: PBPKCompartment,
        covariates: Optional[Dict[str, Any]]
    ) -> float:
        """
        Adjust enzymatic activity by covariates.
        
        Supports CYP450 inhibition/induction and genetic polymorphisms.
        """
        if not covariates:
            return base_activity
        
        activity = base_activity
        
        # CYP450 inhibition (e.g., from DDI)
        if "cyp_inhibition" in covariates and compartment.cyp450_info:
            cyp_name = compartment.cyp450_info.cyp_name
            inhibitors = covariates["cyp_inhibition"]  # Dict[cyp_name, factor]
            if cyp_name in inhibitors:
                activity = activity * inhibitors[cyp_name]
        
        # CYP450 induction
        if "cyp_induction" in covariates and compartment.cyp450_info:
            cyp_name = compartment.cyp450_info.cyp_name
            inducers = covariates["cyp_induction"]  # Dict[cyp_name, factor]
            if cyp_name in inducers:
                activity = activity * inducers[cyp_name]
        
        # Genetic polymorphisms
        if "cyp_phenotype" in covariates and compartment.cyp450_info:
            cyp_name = compartment.cyp450_info.cyp_name
            phenotypes = covariates["cyp_phenotype"]  # Dict[cyp_name, activity_factor]
            if cyp_name in phenotypes:
                activity = activity * phenotypes[cyp_name]
        
        return activity

    def get_organ_clearance(self, organ_name: str) -> float:
        """
        Get intrinsic clearance for a specific organ.
        
        Returns:
            Intrinsic clearance in L/h
        """
        organ = self.structure.get_compartment_by_name(organ_name)
        if organ is None:
            raise ValueError(f"Organ '{organ_name}' not found")
        
        if not organ.is_elimination_site:
            return 0.0
        
        n_elim = len(self.structure.get_elimination_compartments())
        return (self.params.CL / n_elim) if n_elim > 0 else 0.0

    def get_organ_blood_flow_fraction(self, organ_name: str) -> float:
        """
        Get fraction of cardiac output going to organ.
        
        Returns:
            Fraction of cardiac output (0-1)
        """
        organ = self.structure.get_compartment_by_name(organ_name)
        if organ is None:
            raise ValueError(f"Organ '{organ_name}' not found")
        
        return organ.blood_flow / max(self.structure.cardiac_output, 1e-6)

    def get_hepatic_extraction_ratio(self) -> float:
        """
        Calculate hepatic extraction ratio.
        
        ER_H = CLint_H / (Qh + CLint_H)
        where Qh = hepatic blood flow
        
        Returns:
            Extraction ratio (0-1)
        """
        liver = self.structure.get_compartment_by_name("liver")
        if liver is None:
            logger.warning("Liver not found in PBPK structure")
            return 0.0
        
        if not liver.is_elimination_site:
            return 0.0
        
        clint = self.get_organ_clearance("liver")
        qh = liver.blood_flow
        
        return clint / (qh + clint + 1e-12)
