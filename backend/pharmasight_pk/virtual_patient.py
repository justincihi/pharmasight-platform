"""
Virtual Patient Generator for PharmaSight.

Creates diverse patient populations for simulation and education.
Includes pre-built archetypes and random patient generation.

Use cases:
- Educational demonstrations ("flight simulator" for clinicians)
- PopPK model validation
- DDI screening across populations
- Dose individualization examples
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Any, Optional, List
import numpy as np
import logging
from enum import Enum

from .popPK import (
    CovariateProfile, Sex, RenalFunction, HepaticFunction, 
    CYPPhenotype, PopPKPredictor, PopPKPrediction
)

logger = logging.getLogger(__name__)


class PatientArchetype(Enum):
    """Pre-defined patient archetypes for common clinical scenarios."""
    HEALTHY_ADULT = "healthy_adult"
    ELDERLY = "elderly"
    RENAL_MILD = "renal_mild"
    RENAL_MODERATE = "renal_moderate"
    RENAL_SEVERE = "renal_severe"
    HEPATIC_MILD = "hepatic_mild"
    HEPATIC_MODERATE = "hepatic_moderate"
    OBESE = "obese"
    UNDERWEIGHT = "underweight"
    PREGNANT_1ST = "pregnant_1st"
    PREGNANT_2ND = "pregnant_2nd"
    PREGNANT_3RD = "pregnant_3rd"
    PEDIATRIC = "pediatric"
    CYTOPENIAS = "cytopenias"  # CYP poor metabolizer


@dataclass
class VirtualPatient:
    """
    A virtual patient instance.
    
    Attributes:
        patient_id: Unique identifier
        name: Patient display name
        covariates: PK covariate profile
        archetype: Source archetype (if applicable)
        created_at: Creation timestamp
        notes: Clinical notes/description
    """
    patient_id: str
    covariates: CovariateProfile
    name: str = ""
    archetype: Optional[PatientArchetype] = None
    created_at: str = ""
    notes: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "patient_id": self.patient_id,
            "name": self.name,
            "archetype": self.archetype.value if self.archetype else None,
            "covariates": self.covariates.to_dict(),
            "created_at": self.created_at,
            "notes": self.notes,
        }


class VirtualPatientFactory:
    """
    Factory for creating virtual patients.
    
    Provides:
    - Pre-built archetypes
    - Random patient generator
    - Population samplers
    """

    def __init__(self):
        """Initialize factory."""
        self.patient_counter = 0

    def _get_next_patient_id(self) -> str:
        """Generate unique patient ID."""
        self.patient_counter += 1
        return f"VP_{self.patient_counter:06d}"

    def create_archetype(self, archetype: PatientArchetype) -> VirtualPatient:
        """
        Create patient from predefined archetype.
        
        Args:
            archetype: PatientArchetype to instantiate
            
        Returns:
            VirtualPatient with appropriate covariates
        """
        patient_id = self._get_next_patient_id()
        
        # Define archetypes
        archetypes_db = {
            PatientArchetype.HEALTHY_ADULT: {
                "age_years": 45,
                "weight_kg": 70,
                "height_cm": 170,
                "sex": Sex.MALE,
                "egfr": 95,
                "albumin_g_dl": 4.0,
                "liver_function": HepaticFunction.NORMAL,
                "cyp3a4_phenotype": CYPPhenotype.EXTENSIVE,
                "cyp2d6_phenotype": CYPPhenotype.EXTENSIVE,
                "notes": "Healthy 45-year-old adult"
            },
            PatientArchetype.ELDERLY: {
                "age_years": 78,
                "weight_kg": 65,
                "height_cm": 165,
                "sex": Sex.FEMALE,
                "egfr": 55,
                "albumin_g_dl": 3.5,
                "liver_function": HepaticFunction.NORMAL,
                "cyp3a4_phenotype": CYPPhenotype.EXTENSIVE,
                "cyp2d6_phenotype": CYPPhenotype.EXTENSIVE,
                "notes": "Elderly patient with age-related renal decline"
            },
            PatientArchetype.RENAL_MILD: {
                "age_years": 55,
                "weight_kg": 75,
                "height_cm": 175,
                "sex": Sex.MALE,
                "egfr": 75,
                "albumin_g_dl": 3.8,
                "liver_function": HepaticFunction.NORMAL,
                "notes": "Mild renal impairment (eGFR 60-89)"
            },
            PatientArchetype.RENAL_MODERATE: {
                "age_years": 62,
                "weight_kg": 70,
                "height_cm": 170,
                "sex": Sex.FEMALE,
                "egfr": 45,
                "albumin_g_dl": 3.6,
                "liver_function": HepaticFunction.NORMAL,
                "notes": "Moderate renal impairment (eGFR 30-59)"
            },
            PatientArchetype.RENAL_SEVERE: {
                "age_years": 58,
                "weight_kg": 68,
                "height_cm": 168,
                "sex": Sex.MALE,
                "egfr": 20,
                "albumin_g_dl": 3.2,
                "liver_function": HepaticFunction.NORMAL,
                "notes": "Severe renal impairment (eGFR 15-29)"
            },
            PatientArchetype.HEPATIC_MILD: {
                "age_years": 52,
                "weight_kg": 72,
                "height_cm": 172,
                "sex": Sex.MALE,
                "egfr": 90,
                "albumin_g_dl": 3.9,
                "liver_function": HepaticFunction.MILD,
                "notes": "Mild hepatic impairment (Child-Pugh A)"
            },
            PatientArchetype.HEPATIC_MODERATE: {
                "age_years": 55,
                "weight_kg": 68,
                "height_cm": 168,
                "sex": Sex.FEMALE,
                "egfr": 85,
                "albumin_g_dl": 3.2,
                "liver_function": HepaticFunction.MODERATE,
                "notes": "Moderate hepatic impairment (Child-Pugh B)"
            },
            PatientArchetype.OBESE: {
                "age_years": 48,
                "weight_kg": 115,  # BMI ~38
                "height_cm": 173,
                "sex": Sex.FEMALE,
                "egfr": 75,
                "albumin_g_dl": 3.9,
                "liver_function": HepaticFunction.NORMAL,
                "notes": "Obese patient (BMI 30-40)"
            },
            PatientArchetype.UNDERWEIGHT: {
                "age_years": 42,
                "weight_kg": 52,  # BMI ~18
                "height_cm": 170,
                "sex": Sex.MALE,
                "egfr": 95,
                "albumin_g_dl": 3.8,
                "liver_function": HepaticFunction.NORMAL,
                "notes": "Underweight patient (BMI <20)"
            },
            PatientArchetype.PREGNANT_1ST: {
                "age_years": 32,
                "weight_kg": 68,
                "height_cm": 166,
                "sex": Sex.FEMALE,
                "egfr": 100,
                "albumin_g_dl": 3.8,
                "liver_function": HepaticFunction.NORMAL,
                "pregnancy_trimester": 1,
                "notes": "Pregnant, 1st trimester"
            },
            PatientArchetype.PREGNANT_2ND: {
                "age_years": 32,
                "weight_kg": 73,
                "height_cm": 166,
                "sex": Sex.FEMALE,
                "egfr": 105,
                "albumin_g_dl": 3.5,
                "liver_function": HepaticFunction.NORMAL,
                "pregnancy_trimester": 2,
                "notes": "Pregnant, 2nd trimester"
            },
            PatientArchetype.PREGNANT_3RD: {
                "age_years": 32,
                "weight_kg": 80,
                "height_cm": 166,
                "sex": Sex.FEMALE,
                "egfr": 110,
                "albumin_g_dl": 3.3,
                "liver_function": HepaticFunction.NORMAL,
                "pregnancy_trimester": 3,
                "notes": "Pregnant, 3rd trimester - increased CL expected"
            },
            PatientArchetype.PEDIATRIC: {
                "age_years": 8,
                "weight_kg": 25,
                "height_cm": 125,
                "sex": Sex.MALE,
                "egfr": 120,
                "albumin_g_dl": 4.1,
                "liver_function": HepaticFunction.NORMAL,
                "notes": "Pediatric patient (8 years old)"
            },
            PatientArchetype.CYTOPENIAS: {
                "age_years": 50,
                "weight_kg": 70,
                "height_cm": 170,
                "sex": Sex.MALE,
                "egfr": 90,
                "albumin_g_dl": 3.9,
                "liver_function": HepaticFunction.NORMAL,
                "cyp3a4_phenotype": CYPPhenotype.POOR,
                "cyp2d6_phenotype": CYPPhenotype.POOR,
                "cyp2c9_phenotype": CYPPhenotype.POOR,
                "notes": "CYP poor metabolizer - expected higher drug exposure"
            },
        }
        
        if archetype not in archetypes_db:
            raise ValueError(f"Unknown archetype: {archetype}")
        
        params = archetypes_db[archetype]
        notes = params.pop("notes", "")
        
        covariates = CovariateProfile(**params)
        
        patient = VirtualPatient(
            patient_id=patient_id,
            name=f"{archetype.value}_{patient_id}",
            covariates=covariates,
            archetype=archetype,
            notes=notes
        )
        
        logger.info(f"Created patient {patient.patient_id} from archetype {archetype.value}")
        return patient

    def create_random(
        self,
        age_range: tuple = (18, 85),
        weight_range: tuple = (50, 120),
        include_pregnancy: bool = False,
        include_disease_states: bool = True,
        random_seed: Optional[int] = None
    ) -> VirtualPatient:
        """
        Create random patient from parameter distributions.
        
        Args:
            age_range: (min_age, max_age) in years
            weight_range: (min_weight, max_weight) in kg
            include_pregnancy: Whether to randomly assign pregnancy
            include_disease_states: Whether to randomly assign renal/hepatic impairment
            random_seed: Random seed for reproducibility
            
        Returns:
            Randomly generated VirtualPatient
        """
        if random_seed is not None:
            np.random.seed(random_seed)
        
        patient_id = self._get_next_patient_id()
        
        # Random demographics
        age = np.random.uniform(age_range[0], age_range[1])
        weight = np.random.uniform(weight_range[0], weight_range[1])
        height = np.random.normal(170, 8)  # Normal distribution
        sex = np.random.choice([Sex.MALE, Sex.FEMALE])
        
        # Random eGFR (with age correlation)
        base_egfr = 95 - (age - 30) * 0.5  # Decline with age
        egfr = np.random.normal(base_egfr, 15)
        egfr = max(10, egfr)  # Floor at 10
        
        # Random albumin
        albumin = np.random.normal(3.8, 0.4)
        albumin = max(2.0, min(5.0, albumin))
        
        # Random liver function
        if include_disease_states and np.random.random() < 0.1:
            liver_function = np.random.choice(
                [HepaticFunction.MILD, HepaticFunction.MODERATE],
                p=[0.8, 0.2]
            )
        else:
            liver_function = HepaticFunction.NORMAL
        
        # Random CYP phenotypes
        cyp_phenotypes = [CYPPhenotype.POOR, CYPPhenotype.INTERMEDIATE, 
                         CYPPhenotype.EXTENSIVE, CYPPhenotype.ULTRA_RAPID]
        cyp_probs = [0.05, 0.20, 0.70, 0.05]  # Realistic frequencies
        
        cyp3a4 = np.random.choice(cyp_phenotypes, p=cyp_probs)
        cyp2d6 = np.random.choice(cyp_phenotypes, p=cyp_probs)
        cyp2c9 = np.random.choice(cyp_phenotypes, p=cyp_probs)
        
        # Random smoking
        smoking = np.random.random() < 0.2
        
        # Optional pregnancy
        pregnancy_trimester = None
        if include_pregnancy and sex == Sex.FEMALE and age < 50:
            if np.random.random() < 0.05:  # 5% chance
                pregnancy_trimester = np.random.randint(1, 4)
        
        covariates = CovariateProfile(
            age_years=age,
            weight_kg=weight,
            height_cm=height,
            sex=sex,
            egfr=egfr,
            albumin_g_dl=albumin,
            liver_function=liver_function,
            cyp3a4_phenotype=cyp3a4,
            cyp2d6_phenotype=cyp2d6,
            cyp2c9_phenotype=cyp2c9,
            smoking=smoking,
            pregnancy_trimester=pregnancy_trimester,
        )
        
        patient = VirtualPatient(
            patient_id=patient_id,
            name=f"Random_{patient_id}",
            covariates=covariates,
            notes="Randomly generated patient"
        )
        
        logger.info(f"Generated random patient {patient.patient_id}")
        return patient

    def create_population(
        self,
        n_patients: int,
        archetype_distribution: Optional[Dict[PatientArchetype, float]] = None,
        random_patients: bool = False,
        random_seed: Optional[int] = None
    ) -> List[VirtualPatient]:
        """
        Create a population of virtual patients.
        
        Args:
            n_patients: Number of patients to generate
            archetype_distribution: Dict of archetype → probability
            random_patients: If True, generate random patients; else use archetypes
            random_seed: Random seed for reproducibility
            
        Returns:
            List of VirtualPatient objects
        """
        if random_seed is not None:
            np.random.seed(random_seed)
        
        population = []
        
        if random_patients:
            for _ in range(n_patients):
                patient = self.create_random(random_seed=None)
                population.append(patient)
        else:
            # Sample from archetype distribution
            if archetype_distribution is None:
                # Default: uniform distribution
                all_archetypes = list(PatientArchetype)
                archetype_distribution = {a: 1.0 / len(all_archetypes) for a in all_archetypes}
            
            archetypes = list(archetype_distribution.keys())
            probabilities = list(archetype_distribution.values())
            
            for _ in range(n_patients):
                archetype = np.random.choice(archetypes, p=probabilities)
                patient = self.create_archetype(archetype)
                population.append(patient)
        
        logger.info(f"Generated population of {n_patients} patients")
        return population


def create_patient_cohort_by_archetype() -> Dict[PatientArchetype, VirtualPatient]:
    """
    Convenience function: create one patient of each archetype.
    
    Returns:
        Dict mapping archetype → VirtualPatient
    """
    factory = VirtualPatientFactory()
    cohort = {}
    
    for archetype in PatientArchetype:
        cohort[archetype] = factory.create_archetype(archetype)
    
    logger.info(f"Created reference cohort with {len(cohort)} archetypes")
    return cohort


def demonstrate_popPK_predictions(drug_name: str = "Example Drug") -> None:
    """
    Educational demonstration of PopPK predictions across archetypes.
    
    Args:
        drug_name: Name of drug to demonstrate
    """
    # Create reference cohort
    cohort = create_patient_cohort_by_archetype()
    
    # Create PopPK predictor
    predictor = PopPKPredictor(drug_name)
    predictor.set_base_parameters(
        cl=10.0,  # L/h at 70 kg
        vc=50.0,  # L
        ka=0.8,  # 1/h
        f=0.8,
    )
    
    # Predict for each archetype
    print(f"\n{'='*60}")
    print(f"Population PK Predictions for {drug_name}")
    print(f"{'='*60}\n")
    
    for archetype, patient in cohort.items():
        prediction = predictor.predict(patient.covariates)
        
        print(f"Archetype: {archetype.value}")
        print(f"  Age: {patient.covariates.age_years:.1f} yrs | Weight: {patient.covariates.weight_kg:.1f} kg")
        print(f"  eGFR: {patient.covariates.egfr:.1f} | Liver: {patient.covariates.liver_function.value}")
        print(f"  → Predicted CL: {prediction.predicted_cl:.2f} L/h")
        print(f"  → Predicted Vc: {prediction.predicted_vc:.2f} L")
        print(f"  → Predicted Ka: {prediction.predicted_ka:.3f} 1/h" if prediction.predicted_ka else "")
        print(f"  → Predicted F:  {prediction.predicted_f:.3f}\n")


if __name__ == "__main__":
    # Example usage
    demonstrate_popPK_predictions()
