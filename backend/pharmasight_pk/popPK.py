"""
Population Pharmacokinetics (PopPK) models.

Implements covariate-based PK parameter prediction for diverse patient populations.
Links patient characteristics to PK parameters using established models.

Supports:
- Allometric scaling (weight, body surface area)
- Age-based adjustments
- Renal function adjustments (GFR, creatinine clearance)
- Liver function assessment
- CYP450 genetic polymorphisms
- Disease state effects
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Any, Optional, Callable
import numpy as np
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


class Sex(Enum):
    """Biological sex for dose calculations."""
    MALE = "male"
    FEMALE = "female"
    OTHER = "other"


class RenalFunction(Enum):
    """Renal function categories based on eGFR (mL/min/1.73m²)."""
    NORMAL = (90, 120)  # eGFR ≥ 90
    MILD = (60, 89)  # eGFR 60-89
    MODERATE = (30, 59)  # eGFR 30-59
    SEVERE = (15, 29)  # eGFR 15-29
    DIALYSIS = (0, 15)  # eGFR < 15 or on dialysis
    
    @classmethod
    def from_egfr(cls, egfr: float) -> RenalFunction:
        """Determine renal function category from eGFR."""
        if egfr >= 90:
            return cls.NORMAL
        elif egfr >= 60:
            return cls.MILD
        elif egfr >= 30:
            return cls.MODERATE
        elif egfr >= 15:
            return cls.SEVERE
        else:
            return cls.DIALYSIS


class HepaticFunction(Enum):
    """Hepatic function categories (Child-Pugh)."""
    NORMAL = "normal"  # Normal liver function
    MILD = "mild"  # Mild hepatic impairment (Child-Pugh A)
    MODERATE = "moderate"  # Moderate hepatic impairment (Child-Pugh B)
    SEVERE = "severe"  # Severe hepatic impairment (Child-Pugh C)


class CYPPhenotype(Enum):
    """CYP450 metabolizer phenotype."""
    POOR = "poor"  # PM - ~0.5x activity
    INTERMEDIATE = "intermediate"  # IM - ~0.8x activity
    EXTENSIVE = "extensive"  # EM - ~1.0x activity (normal)
    ULTRA_RAPID = "ultra_rapid"  # UM - ~2.0x activity


@dataclass
class CovariateProfile:
    """
    Patient covariate profile for PopPK predictions.
    
    Attributes:
        age_years: Age in years
        weight_kg: Body weight in kg
        height_cm: Height in cm (optional, for BSA calculation)
        sex: Biological sex
        egfr: Estimated glomerular filtration rate (mL/min/1.73m²)
        albumin_g_dl: Serum albumin (g/dL) - affects protein binding
        liver_function: Hepatic function status
        cyp3a4_phenotype: CYP3A4 metabolizer phenotype
        cyp2d6_phenotype: CYP2D6 metabolizer phenotype
        cyp2c9_phenotype: CYP2C9 metabolizer phenotype
        smoking: Cigarette smoking status (affects CYP1A2, CYP2B6)
        pregnancy_trimester: Pregnancy status (1, 2, 3, or None)
        extra: Additional patient-specific factors
    """
    age_years: float
    weight_kg: float
    sex: Sex
    egfr: float = 90.0
    height_cm: Optional[float] = None
    albumin_g_dl: float = 4.0
    liver_function: HepaticFunction = HepaticFunction.NORMAL
    cyp3a4_phenotype: CYPPhenotype = CYPPhenotype.EXTENSIVE
    cyp2d6_phenotype: CYPPhenotype = CYPPhenotype.EXTENSIVE
    cyp2c9_phenotype: CYPPhenotype = CYPPhenotype.EXTENSIVE
    smoking: bool = False
    pregnancy_trimester: Optional[int] = None
    extra: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate covariate values."""
        if self.age_years < 0:
            raise ValueError(f"Age must be >= 0, got {self.age_years}")
        if self.weight_kg <= 0:
            raise ValueError(f"Weight must be > 0, got {self.weight_kg}")
        if self.height_cm is not None and self.height_cm <= 0:
            raise ValueError(f"Height must be > 0, got {self.height_cm}")
        if self.egfr < 0:
            raise ValueError(f"eGFR must be >= 0, got {self.egfr}")
        if self.albumin_g_dl < 0:
            raise ValueError(f"Albumin must be >= 0, got {self.albumin_g_dl}")
        if self.pregnancy_trimester is not None and self.pregnancy_trimester not in [1, 2, 3]:
            raise ValueError(f"Pregnancy trimester must be 1, 2, 3, or None")
    
    def get_bsa(self) -> float:
        """
        Calculate Body Surface Area (BSA) in m².
        
        Uses Mosteller formula if height available, otherwise uses DuBois approximation.
        
        BSA (Mosteller) = sqrt(weight_kg * height_cm / 3600)
        BSA (DuBois) = 0.007184 * weight^0.425 * height^0.725
        """
        if self.height_cm is not None:
            # Mosteller formula
            return np.sqrt((self.weight_kg * self.height_cm) / 3600.0)
        else:
            # DuBois approximation (uses weight only, assumes average height)
            # Estimate height from weight (rough correlation)
            estimated_height = 1.3 + 0.6 * (self.weight_kg / 70.0)  # in decimeters for DuBois
            return 0.007184 * (self.weight_kg ** 0.425) * (estimated_height * 10 ** 0.725)
    
    def get_renal_function_category(self) -> RenalFunction:
        """Get renal function category from eGFR."""
        return RenalFunction.from_egfr(self.egfr)
    
    def get_cyp_activity_factor(self, cyp_name: str) -> float:
        """
        Get enzymatic activity factor for CYP (0-2.0, where 1.0 = normal).
        
        Args:
            cyp_name: CYP enzyme name (e.g., "CYP3A4", "CYP2D6")
            
        Returns:
            Activity factor relative to extensive metabolizer
        """
        phenotype_factors = {
            CYPPhenotype.POOR: 0.5,
            CYPPhenotype.INTERMEDIATE: 0.8,
            CYPPhenotype.EXTENSIVE: 1.0,
            CYPPhenotype.ULTRA_RAPID: 2.0,
        }
        
        if cyp_name == "CYP3A4":
            return phenotype_factors.get(self.cyp3a4_phenotype, 1.0)
        elif cyp_name == "CYP2D6":
            return phenotype_factors.get(self.cyp2d6_phenotype, 1.0)
        elif cyp_name == "CYP2C9":
            return phenotype_factors.get(self.cyp2c9_phenotype, 1.0)
        elif cyp_name in ["CYP1A2", "CYP2B6"] and self.smoking:
            # Smoking induces CYP1A2 and CYP2B6 by ~1.5-2x
            return 1.5
        else:
            return 1.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize covariate profile to dictionary."""
        return {
            "age_years": self.age_years,
            "weight_kg": self.weight_kg,
            "height_cm": self.height_cm,
            "sex": self.sex.value,
            "egfr": self.egfr,
            "albumin_g_dl": self.albumin_g_dl,
            "liver_function": self.liver_function.value,
            "cyp3a4_phenotype": self.cyp3a4_phenotype.value,
            "cyp2d6_phenotype": self.cyp2d6_phenotype.value,
            "cyp2c9_phenotype": self.cyp2c9_phenotype.value,
            "smoking": self.smoking,
            "pregnancy_trimester": self.pregnancy_trimester,
            "bsa": self.get_bsa(),
            "renal_function": self.get_renal_function_category().name,
        }


class PopPKModel:
    """
    Population PK model relating covariates to PK parameters.
    
    Implements allometric scaling, renal/hepatic adjustments, and CYP phenotype effects.
    """

    @staticmethod
    def predict_clearance(
        base_cl: float,
        covariates: CovariateProfile,
        renal_excretion_fraction: float = 0.0,
        hepatic_extraction_ratio: float = 0.7
    ) -> float:
        """
        Predict total clearance based on covariates.
        
        CL_pred = CL_base * (weight/70)^0.75 * (age_factor) * (renal_factor) * (hepatic_factor) * (cyp_factor)
        
        Args:
            base_cl: Baseline clearance at 70 kg, normal function (L/h)
            covariates: Patient covariate profile
            renal_excretion_fraction: Fraction of drug cleared renally (0-1)
            hepatic_extraction_ratio: Hepatic extraction ratio (0-1)
            
        Returns:
            Predicted clearance (L/h)
        """
        cl = base_cl
        
        # Allometric scaling with weight
        cl = cl * (covariates.weight_kg / 70.0) ** 0.75
        
        # Age adjustment (slight decline with age)
        # Typical: ~1% decline per year after 30 years
        if covariates.age_years > 30:
            age_factor = 1.0 - 0.01 * (covariates.age_years - 30) / 70.0
            cl = cl * max(age_factor, 0.5)  # Floor at 50%
        
        # Renal function adjustment
        if renal_excretion_fraction > 0:
            egfr_factor = covariates.egfr / 90.0  # 90 = normal eGFR
            renal_cl = base_cl * renal_excretion_fraction * egfr_factor
            cl = cl * (1 - renal_excretion_fraction) + renal_cl
        
        # Hepatic function adjustment
        hepatic_factor = PopPKModel._get_hepatic_function_factor(covariates.liver_function)
        hepatic_cl = base_cl * hepatic_extraction_ratio * hepatic_factor
        cl = cl * (1 - hepatic_extraction_ratio) + hepatic_cl
        
        # CYP3A4 activity (if major metabolic pathway)
        if hepatic_extraction_ratio > 0.5:
            cyp_factor = covariates.get_cyp_activity_factor("CYP3A4")
            cl = cl * cyp_factor
        
        # Pregnancy adjustments (increased clearance, especially 3rd trimester)
        if covariates.pregnancy_trimester is not None:
            pregnancy_factors = {1: 1.1, 2: 1.3, 3: 1.5}
            cl = cl * pregnancy_factors.get(covariates.pregnancy_trimester, 1.0)
        
        return cl

    @staticmethod
    def predict_volume(
        base_vc: float,
        covariates: CovariateProfile,
        volume_scaling: str = "weight"  # "weight" or "bsa"
    ) -> float:
        """
        Predict central volume of distribution.
        
        Vc_pred = Vc_base * (weight/70) or (BSA/1.73)
        
        Args:
            base_vc: Baseline Vc at 70 kg (L)
            covariates: Patient covariate profile
            volume_scaling: Scaling method ("weight" or "bsa")
            
        Returns:
            Predicted Vc (L)
        """
        if volume_scaling == "bsa":
            bsa = covariates.get_bsa()
            return base_vc * (bsa / 1.73)
        else:  # weight
            return base_vc * (covariates.weight_kg / 70.0)

    @staticmethod
    def predict_absorption_rate(
        base_ka: float,
        covariates: CovariateProfile,
        fed_state: bool = False
    ) -> float:
        """
        Predict first-order absorption rate constant.
        
        Ka_pred = Ka_base * (age_factor) * (fed_factor) * (gi_motility_factor)
        
        Args:
            base_ka: Baseline Ka at standard conditions (1/h)
            covariates: Patient covariate profile
            fed_state: Whether patient is in fed state
            
        Returns:
            Predicted Ka (1/h)
        """
        ka = base_ka
        
        # Age effect (slower GI motility in elderly)
        if covariates.age_years > 65:
            age_factor = 1.0 - 0.005 * (covariates.age_years - 65)
            ka = ka * max(age_factor, 0.7)
        
        # Fed state (decreases Ka)
        if fed_state:
            ka = ka * 0.7
        
        # Pregnancy effects (faster gastric emptying in some trimesters)
        if covariates.pregnancy_trimester == 3:
            ka = ka * 1.2
        
        return ka

    @staticmethod
    def _get_hepatic_function_factor(liver_function: HepaticFunction) -> float:
        """Get clearance factor for hepatic dysfunction."""
        factors = {
            HepaticFunction.NORMAL: 1.0,
            HepaticFunction.MILD: 0.8,  # 20% reduction
            HepaticFunction.MODERATE: 0.6,  # 40% reduction
            HepaticFunction.SEVERE: 0.3,  # 70% reduction
        }
        return factors.get(liver_function, 1.0)

    @staticmethod
    def predict_bioavailability(
        base_f: float,
        covariates: CovariateProfile,
        substrate_of_pgp: bool = False
    ) -> float:
        """
        Predict oral bioavailability.
        
        Factors: hepatic metabolism, gut wall metabolism, P-gp efflux
        
        Args:
            base_f: Baseline bioavailability (0-1)
            covariates: Patient covariate profile
            substrate_of_pgp: Whether drug is P-gp substrate
            
        Returns:
            Predicted F (0-1)
        """
        f = base_f
        
        # Hepatic dysfunction increases F for high-extraction drugs
        if covariates.liver_function != HepaticFunction.NORMAL:
            factor = PopPKModel._get_hepatic_function_factor(covariates.liver_function)
            # Higher factor → lower hepatic CL → higher F
            f = f * (1.0 + (1.0 - factor) * 0.5)
        
        # P-gp activity (reduced in elderly, increased with inducers)
        if substrate_of_pgp:
            if covariates.age_years > 65:
                f = f * 1.15  # Reduced P-gp activity
        
        # Ensure F stays in valid range
        f = max(0.0, min(1.0, f))
        
        return f


@dataclass
class PopPKPrediction:
    """Results of PopPK prediction for a patient."""
    covariates: CovariateProfile
    predicted_cl: float  # L/h
    predicted_vc: float  # L
    predicted_ka: Optional[float]  # 1/h (for oral models)
    predicted_f: float  # Bioavailability
    qp: Optional[float] = None  # Peripheral clearance for 2-compartment
    vp: Optional[float] = None  # Peripheral volume
    
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "predicted_cl": self.predicted_cl,
            "predicted_vc": self.predicted_vc,
            "predicted_ka": self.predicted_ka,
            "predicted_f": self.predicted_f,
            "qp": self.qp,
            "vp": self.vp,
            "patient_covariates": self.covariates.to_dict(),
            "timestamp": self.timestamp,
        }


class PopPKPredictor:
    """
    Convenience class for predicting PK parameters for specific drugs.
    
    Drug-specific PopPK models can inherit from this and override methods.
    """

    def __init__(self, drug_name: str):
        """
        Initialize predictor for a drug.
        
        Args:
            drug_name: Name of the drug
        """
        self.drug_name = drug_name
        self.base_parameters = {}

    def set_base_parameters(
        self,
        cl: float,
        vc: float,
        ka: Optional[float] = None,
        f: float = 1.0,
        qp: Optional[float] = None,
        vp: Optional[float] = None
    ) -> None:
        """
        Set reference PK parameters (typically from literature or NONMEM fit).
        
        Args:
            cl: Clearance at 70 kg, normal function (L/h)
            vc: Central volume at 70 kg (L)
            ka: Absorption rate (1/h)
            f: Bioavailability
            qp: Peripheral clearance (L/h)
            vp: Peripheral volume (L)
        """
        self.base_parameters = {
            "cl": cl,
            "vc": vc,
            "ka": ka,
            "f": f,
            "qp": qp,
            "vp": vp,
        }

    def predict(
        self,
        covariates: CovariateProfile,
        renal_excretion_fraction: float = 0.0,
        hepatic_extraction_ratio: float = 0.7,
        fed_state: bool = False
    ) -> PopPKPrediction:
        """
        Predict PK parameters for a patient.
        
        Args:
            covariates: Patient covariate profile
            renal_excretion_fraction: Fraction of drug cleared by kidneys
            hepatic_extraction_ratio: Fraction cleared by liver
            fed_state: Whether patient is in fed state (affects Ka, F)
            
        Returns:
            PopPKPrediction with predicted parameters
        """
        if not self.base_parameters:
            raise ValueError(f"Base parameters not set for {self.drug_name}")

        pred_cl = PopPKModel.predict_clearance(
            self.base_parameters["cl"],
            covariates,
            renal_excretion_fraction,
            hepatic_extraction_ratio
        )

        pred_vc = PopPKModel.predict_volume(
            self.base_parameters["vc"],
            covariates
        )

        pred_ka = None
        if self.base_parameters["ka"] is not None:
            pred_ka = PopPKModel.predict_absorption_rate(
                self.base_parameters["ka"],
                covariates,
                fed_state
            )

        pred_f = PopPKModel.predict_bioavailability(
            self.base_parameters["f"],
            covariates
        )

        return PopPKPrediction(
            covariates=covariates,
            predicted_cl=pred_cl,
            predicted_vc=pred_vc,
            predicted_ka=pred_ka,
            predicted_f=pred_f,
            qp=self.base_parameters.get("qp"),
            vp=self.base_parameters.get("vp"),
        )
