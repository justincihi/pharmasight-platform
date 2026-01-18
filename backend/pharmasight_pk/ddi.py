"""
Drug-Drug Interaction (DDI) module for PharmaSight.

Provides:
- Interaction database with mechanism classification
- Quantitative DDI factor calculation
- Clinical severity assessment and recommendations
- Audit trail for all DDI calculations
- Integration with PK models for parameter adjustments
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, Any, Optional, List, Tuple, Set
import logging
from datetime import datetime
import json

logger = logging.getLogger(__name__)


class InteractionMechanism(Enum):
    """Mechanism of drug-drug interaction."""
    COMPETITIVE_INHIBITION = "competitive_inhibition"
    MECHANISM_BASED_INHIBITION = "mechanism_based_inhibition"
    TIME_DEPENDENT_INHIBITION = "time_dependent_inhibition"
    CYP_INDUCTION = "cyp_induction"
    SUBSTRATE_COMPETITION = "substrate_competition"
    PROTEIN_BINDING_DISPLACEMENT = "protein_binding_displacement"
    TRANSPORTER_INHIBITION = "transporter_inhibition"
    UNKNOWN = "unknown"


class InteractionSeverity(Enum):
    """Clinical severity of an interaction."""
    MINOR = "minor"  # Unlikely to require intervention
    MODERATE = "moderate"  # May require dosage adjustment or monitoring
    MAJOR = "major"  # Significant interaction; consider alternative
    CONTRAINDICATED = "contraindicated"  # Avoid combination unless critical


class CYPEnzyme(Enum):
    """Major human CYP450 enzymes."""
    CYP1A2 = "CYP1A2"
    CYP2B6 = "CYP2B6"
    CYP2C8 = "CYP2C8"
    CYP2C9 = "CYP2C9"
    CYP2C19 = "CYP2C19"
    CYP2D6 = "CYP2D6"
    CYP3A4 = "CYP3A4"
    CYP3A5 = "CYP3A5"


@dataclass
class InhibitionInfo:
    """Information about CYP inhibition."""
    cyp_enzyme: CYPEnzyme
    ki: Optional[float] = None  # Competitive inhibition constant (µM)
    ic50: Optional[float] = None  # Concentration for 50% inhibition (µM)
    mechanism: InteractionMechanism = InteractionMechanism.COMPETITIVE_INHIBITION
    time_dependent: bool = False  # Whether inhibition is time-dependent
    kinact: Optional[float] = None  # Inactivation rate constant (1/h) for mechanism-based
    km_inactivation: Optional[float] = None  # Km for inactivation (µM)
    
    def __post_init__(self):
        if self.ki is None and self.ic50 is None:
            raise ValueError("Either Ki or IC50 must be specified")
        if self.ki is not None and self.ki <= 0:
            raise ValueError(f"Ki must be > 0, got {self.ki}")


@dataclass
class InductionInfo:
    """Information about CYP induction."""
    cyp_enzyme: CYPEnzyme
    ec50: float  # Concentration for 50% induction (µM)
    emax: float = 2.0  # Maximum induction fold-change (default 2x)
    
    def __post_init__(self):
        if self.ec50 <= 0:
            raise ValueError(f"EC50 must be > 0, got {self.ec50}")
        if self.emax <= 0:
            raise ValueError(f"Emax must be > 0, got {self.emax}")


@dataclass
class SubstrateInfo:
    """Information about CYP substrate specificity."""
    primary_cyp: CYPEnzyme
    secondary_cyps: List[CYPEnzyme] = field(default_factory=list)
    km: Optional[float] = None  # Michaelis constant for primary pathway (µM)
    intrinsic_clearance: Optional[float] = None  # CLint (mL/min/mg protein)


@dataclass
class ClinicalRecommendation:
    """Clinical guidance for managing an interaction."""
    recommendation_type: str  # "monitor", "adjust_dose", "alternative_drug", "contraindicated"
    description: str
    monitoring_parameters: List[str] = field(default_factory=list)  # e.g., ["INR", "QTc"]
    dose_adjustment: Optional[str] = None  # e.g., "Reduce to 50% of normal dose"
    alternative_drugs: List[str] = field(default_factory=list)  # Non-interacting alternatives


@dataclass
class DrugInteractionPair:
    """Complete DDI information for two drugs."""
    drug_1: str  # Drug name/identifier
    drug_2: str  # Drug name/identifier
    severity: InteractionSeverity
    mechanism: InteractionMechanism
    
    # Inhibition details (if applicable)
    inhibitor_info: Optional[InhibitionInfo] = None  # Which drug is the inhibitor?
    inhibitor_drug: Optional[str] = None  # Name of inhibiting drug
    
    # Induction details (if applicable)
    inducer_info: Optional[InductionInfo] = None
    inducer_drug: Optional[str] = None
    
    # Substrate details (affected drug)
    substrate_info: Optional[SubstrateInfo] = None
    affected_drug: Optional[str] = None
    
    # Clinical info
    clinical_recommendation: Optional[ClinicalRecommendation] = None
    evidence_quality: str = "moderate"  # "poor", "moderate", "good", "excellent"
    references: List[str] = field(default_factory=list)  # PubMed IDs or citations
    notes: str = ""
    
    def __post_init__(self):
        # Ensure drug names are consistent
        if self.inhibitor_drug and self.inhibitor_drug not in [self.drug_1, self.drug_2]:
            raise ValueError(f"inhibitor_drug must be one of the two drugs")
        if self.inducer_drug and self.inducer_drug not in [self.drug_1, self.drug_2]:
            raise ValueError(f"inducer_drug must be one of the two drugs")
        if self.affected_drug and self.affected_drug not in [self.drug_1, self.drug_2]:
            raise ValueError(f"affected_drug must be one of the two drugs")


class DDIDatabase:
    """
    Centralized database of known drug-drug interactions.
    
    Designed to be:
    - Extensible (easy to add new interactions)
    - Queryable (search by drug, mechanism, severity)
    - AI-ready (can be populated/updated by ML models)
    - Clinically auditable (references and evidence tracking)
    """

    def __init__(self):
        """Initialize empty DDI database."""
        self.interactions: Dict[Tuple[str, str], DrugInteractionPair] = {}
        self.audit_log: List[Dict[str, Any]] = []
        logger.info("Initialized DDI Database")

    def add_interaction(self, interaction: DrugInteractionPair) -> None:
        """
        Add interaction to database.
        
        Args:
            interaction: DrugInteractionPair object
        """
        # Normalize drug order (alphabetical)
        drug_1, drug_2 = sorted([interaction.drug_1, interaction.drug_2])
        key = (drug_1, drug_2)
        
        self.interactions[key] = DrugInteractionPair(
            drug_1=drug_1,
            drug_2=drug_2,
            **{k: v for k, v in interaction.__dict__.items() if k not in ['drug_1', 'drug_2']}
        )
        
        logger.info(f"Added interaction: {drug_1} <-> {drug_2} ({interaction.severity.value})")

    def get_interaction(self, drug_1: str, drug_2: str) -> Optional[DrugInteractionPair]:
        """
        Retrieve interaction between two drugs.
        
        Args:
            drug_1: First drug name
            drug_2: Second drug name
            
        Returns:
            DrugInteractionPair if found, None otherwise
        """
        key = tuple(sorted([drug_1, drug_2]))
        return self.interactions.get(key)

    def get_interactions_for_drug(self, drug: str) -> List[DrugInteractionPair]:
        """
        Get all interactions involving a specific drug.
        
        Args:
            drug: Drug name
            
        Returns:
            List of DrugInteractionPair objects
        """
        results = []
        for (d1, d2), interaction in self.interactions.items():
            if d1.lower() == drug.lower() or d2.lower() == drug.lower():
                results.append(interaction)
        return results

    def get_interactions_by_severity(self, severity: InteractionSeverity) -> List[DrugInteractionPair]:
        """Get all interactions of a specific severity."""
        return [i for i in self.interactions.values() if i.severity == severity]

    def get_interactions_by_mechanism(self, mechanism: InteractionMechanism) -> List[DrugInteractionPair]:
        """Get all interactions with specific mechanism."""
        return [i for i in self.interactions.values() if i.mechanism == mechanism]

    def get_drugs_with_severity(self, drug_list: List[str], min_severity: InteractionSeverity) -> Dict[str, List[str]]:
        """
        Get all drug pairs above severity threshold from a drug list.
        
        Useful for screening a patient's medication list.
        
        Args:
            drug_list: List of drug names
            min_severity: Minimum severity threshold (inclusive)
            
        Returns:
            Dict mapping drugs to list of problematic interactions
        """
        severity_order = [InteractionSeverity.MINOR, InteractionSeverity.MODERATE,
                         InteractionSeverity.MAJOR, InteractionSeverity.CONTRAINDICATED]
        min_idx = severity_order.index(min_severity)
        
        result = {drug: [] for drug in drug_list}
        
        for drug_1 in drug_list:
            for drug_2 in drug_list:
                if drug_1 >= drug_2:  # Avoid duplicates
                    continue
                
                interaction = self.get_interaction(drug_1, drug_2)
                if interaction and severity_order.index(interaction.severity) >= min_idx:
                    result[drug_1].append(f"{drug_2} ({interaction.severity.value})")
                    result[drug_2].append(f"{drug_1} ({interaction.severity.value})")
        
        return result

    def to_dict(self) -> Dict[str, Any]:
        """Serialize database to dictionary."""
        return {
            "interactions": {
                f"{d1}_{d2}": {
                    "drug_1": d1,
                    "drug_2": d2,
                    "severity": i.severity.value,
                    "mechanism": i.mechanism.value,
                    "notes": i.notes,
                }
                for (d1, d2), i in self.interactions.items()
            },
            "num_interactions": len(self.interactions),
        }

    def __len__(self) -> int:
        """Return number of interactions in database."""
        return len(self.interactions)


class DDICalculator:
    """
    Calculates quantitative DDI factors for parameter adjustments.
    
    Supports:
    - Competitive inhibition: F = 1 / (1 + [I]/Ki)
    - Mechanism-based inhibition: Time-dependent models
    - Induction: E = Emax * [I] / (EC50 + [I])
    """

    def __init__(self, ddi_db: Optional[DDIDatabase] = None):
        """
        Initialize DDI calculator.
        
        Args:
            ddi_db: Optional reference to DDI database for context
        """
        self.ddi_db = ddi_db or DDIDatabase()
        self.calculation_log: List[Dict[str, Any]] = []

    def calculate_inhibition_factor(
        self,
        inhibitor_concentration: float,
        ki: float,
        mechanism: InteractionMechanism = InteractionMechanism.COMPETITIVE_INHIBITION
    ) -> float:
        """
        Calculate fraction of enzyme activity remaining under inhibition.
        
        Competitive inhibition: Fm = 1 / (1 + [I]/Ki)
        
        Args:
            inhibitor_concentration: Concentration of inhibitor (µM)
            ki: Inhibition constant (µM)
            mechanism: Type of inhibition
            
        Returns:
            Fraction of normal enzyme activity (0-1)
            
        Raises:
            ValueError: If inputs invalid
        """
        if inhibitor_concentration < 0:
            raise ValueError(f"Inhibitor concentration must be >= 0, got {inhibitor_concentration}")
        if ki <= 0:
            raise ValueError(f"Ki must be > 0, got {ki}")
        
        if mechanism == InteractionMechanism.COMPETITIVE_INHIBITION:
            fm = 1.0 / (1.0 + (inhibitor_concentration / ki))
        else:
            # For other mechanisms, default to competitive model
            logger.warning(f"Mechanism {mechanism.value} not fully implemented; using competitive model")
            fm = 1.0 / (1.0 + (inhibitor_concentration / ki))
        
        # Ensure result is in valid range
        fm = max(0.0, min(1.0, fm))
        
        self.calculation_log.append({
            "timestamp": datetime.now().isoformat(),
            "type": "inhibition",
            "inhibitor_conc": inhibitor_concentration,
            "ki": ki,
            "result": fm,
        })
        
        return fm

    def calculate_induction_factor(
        self,
        inducer_concentration: float,
        ec50: float,
        emax: float = 2.0
    ) -> float:
        """
        Calculate fold-change in enzyme activity due to induction.
        
        Sigmoid induction: E = 1 + Emax * [I] / (EC50 + [I])
        
        Args:
            inducer_concentration: Concentration of inducer (µM)
            ec50: Concentration for 50% induction (µM)
            emax: Maximum fold-change (default 2.0 = 2-fold induction)
            
        Returns:
            Enzyme activity multiplier (≥ 1.0)
            
        Raises:
            ValueError: If inputs invalid
        """
        if inducer_concentration < 0:
            raise ValueError(f"Inducer concentration must be >= 0, got {inducer_concentration}")
        if ec50 <= 0:
            raise ValueError(f"EC50 must be > 0, got {ec50}")
        if emax <= 0:
            raise ValueError(f"Emax must be > 0, got {emax}")
        
        # Sigmoid model
        e_factor = 1.0 + (emax - 1.0) * (inducer_concentration / (ec50 + inducer_concentration))
        
        self.calculation_log.append({
            "timestamp": datetime.now().isoformat(),
            "type": "induction",
            "inducer_conc": inducer_concentration,
            "ec50": ec50,
            "emax": emax,
            "result": e_factor,
        })
        
        return e_factor

    def calculate_ddi_factors_from_interaction(
        self,
        drug_1: str,
        drug_2: str,
        drug_1_concentration: float,
        drug_2_concentration: float
    ) -> Dict[str, float]:
        """
        Calculate DDI factors for a drug pair using database information.
        
        Args:
            drug_1: Name of first drug
            drug_2: Name of second drug
            drug_1_concentration: Steady-state concentration of drug 1 (µM)
            drug_2_concentration: Steady-state concentration of drug 2 (µM)
            
        Returns:
            Dict with keys like "CL" (clearance factor), "Vc" (volume factor), etc.
            
        Example:
            >>> calc = DDICalculator(ddi_db)
            >>> factors = calc.calculate_ddi_factors_from_interaction(
            ...     "midazolam", "ketoconazole", 0.01, 0.5
            ... )
            >>> # factors = {"CL": 0.25, "Vc": 1.0}  (75% CL reduction)
        """
        interaction = self.ddi_db.get_interaction(drug_1, drug_2)
        
        if interaction is None:
            logger.info(f"No known interaction between {drug_1} and {drug_2}")
            return {}
        
        ddi_factors = {}
        
        # Determine which drug is inhibitor/inducer
        if interaction.inhibitor_drug:
            inhibitor_conc = (drug_1_concentration if interaction.inhibitor_drug == drug_1
                            else drug_2_concentration)
            
            if interaction.inhibitor_info and interaction.inhibitor_info.ki:
                fm = self.calculate_inhibition_factor(
                    inhibitor_conc,
                    interaction.inhibitor_info.ki
                )
                # CL is reduced; apply to substrate drug
                ddi_factors["CL"] = fm
        
        if interaction.inducer_drug:
            inducer_conc = (drug_1_concentration if interaction.inducer_drug == drug_1
                          else drug_2_concentration)
            
            if interaction.inducer_info:
                e_factor = self.calculate_induction_factor(
                    inducer_conc,
                    interaction.inducer_info.ec50,
                    interaction.inducer_info.emax
                )
                # CL is increased
                ddi_factors["CL"] = e_factor
        
        return ddi_factors

    def get_calculation_log(self) -> List[Dict[str, Any]]:
        """Return audit log of all calculations."""
        return self.calculation_log.copy()

    def clear_log(self) -> None:
        """Clear calculation log."""
        self.calculation_log.clear()


class DDIScreener:
    """
    Clinical DDI screening tool.
    
    Analyzes medication list for:
    - Severity-based alerts
    - Drug-specific recommendations
    - Dosage adjustments
    - Alternative suggestions
    """

    def __init__(self, ddi_db: DDIDatabase):
        """
        Initialize DDI screener.
        
        Args:
            ddi_db: DDI database reference
        """
        self.ddi_db = ddi_db

    def screen_medication_list(
        self,
        drug_list: List[str],
        min_severity: InteractionSeverity = InteractionSeverity.MODERATE
    ) -> Dict[str, Any]:
        """
        Screen a patient's medication list for interactions.
        
        Args:
            drug_list: List of drug names patient is taking
            min_severity: Minimum severity to report
            
        Returns:
            Dict with alerts, recommendations, and severity breakdown
            
        Example:
            >>> screener = DDIScreener(ddi_db)
            >>> result = screener.screen_medication_list(
            ...     ["warfarin", "aspirin", "ibuprofen"],
            ...     min_severity=InteractionSeverity.MODERATE
            ... )
            >>> print(result["alerts"])
        """
        # Get all interactions above threshold
        severity_interactions = self.ddi_db.get_drugs_with_severity(drug_list, min_severity)
        
        alerts = []
        recommendations = []
        
        severity_counts = {s.value: 0 for s in InteractionSeverity}
        
        for drug in drug_list:
            for interaction_str in severity_interactions[drug]:
                # Parse interaction string: "drug_name (severity)"
                parts = interaction_str.rsplit(" (", 1)
                if len(parts) == 2:
                    other_drug = parts[0]
                    sev_str = parts[1].rstrip(")")
                    severity_counts[sev_str] += 1
                    
                    # Get full interaction details
                    full_interaction = self.ddi_db.get_interaction(drug, other_drug)
                    if full_interaction:
                        alerts.append({
                            "drug_1": drug,
                            "drug_2": other_drug,
                            "severity": full_interaction.severity.value,
                            "mechanism": full_interaction.mechanism.value,
                            "recommendation": full_interaction.clinical_recommendation.description
                            if full_interaction.clinical_recommendation else None,
                        })
        
        return {
            "medication_count": len(drug_list),
            "interaction_count": len(alerts),
            "severity_summary": severity_counts,
            "alerts": alerts,
            "timestamp": datetime.now().isoformat(),
        }

    def get_recommendations_for_pair(self, drug_1: str, drug_2: str) -> Optional[ClinicalRecommendation]:
        """Get clinical recommendations for a specific drug pair."""
        interaction = self.ddi_db.get_interaction(drug_1, drug_2)
        if interaction:
            return interaction.clinical_recommendation
        return None


def create_default_ddi_database() -> DDIDatabase:
    """
    Create a database with common clinically-relevant DDIs.
    
    This is a starting template; in production, populate from:
    - DrugBank API
    - FDA interaction checker
    - Clinical trial data
    - AI-generated predictions
    """
    db = DDIDatabase()
    
    # Example 1: Warfarin + NSAIDs (major)
    db.add_interaction(
        DrugInteractionPair(
            drug_1="warfarin",
            drug_2="ibuprofen",
            severity=InteractionSeverity.MAJOR,
            mechanism=InteractionMechanism.PROTEIN_BINDING_DISPLACEMENT,
            affected_drug="warfarin",
            clinical_recommendation=ClinicalRecommendation(
                recommendation_type="monitor",
                description="NSAIDs displace warfarin from protein binding and increase bleeding risk",
                monitoring_parameters=["INR", "bleeding signs"],
            ),
            evidence_quality="excellent",
            notes="Consider acetaminophen as alternative"
        )
    )
    
    # Example 2: Midazolam + Ketoconazole (CYP3A4 inhibition)
    db.add_interaction(
        DrugInteractionPair(
            drug_1="midazolam",
            drug_2="ketoconazole",
            severity=InteractionSeverity.MAJOR,
            mechanism=InteractionMechanism.COMPETITIVE_INHIBITION,
            inhibitor_drug="ketoconazole",
            inhibitor_info=InhibitionInfo(
                cyp_enzyme=CYPEnzyme.CYP3A4,
                ki=0.04,  # µM
                mechanism=InteractionMechanism.COMPETITIVE_INHIBITION,
            ),
            affected_drug="midazolam",
            substrate_info=SubstrateInfo(primary_cyp=CYPEnzyme.CYP3A4),
            clinical_recommendation=ClinicalRecommendation(
                recommendation_type="adjust_dose",
                description="Ketoconazole inhibits CYP3A4, reducing midazolam clearance",
                dose_adjustment="Reduce midazolam dose by 50-75%",
                monitoring_parameters=["Sedation level", "Respiratory rate"],
            ),
            evidence_quality="excellent",
            references=["12345678"],  # PubMed IDs
        )
    )
    
    # Example 3: Metoprolol + Verapamil (additive cardiac effects)
    db.add_interaction(
        DrugInteractionPair(
            drug_1="metoprolol",
            drug_2="verapamil",
            severity=InteractionSeverity.MAJOR,
            mechanism=InteractionMechanism.UNKNOWN,
            clinical_recommendation=ClinicalRecommendation(
                recommendation_type="contraindicated",
                description="Combined beta-blocker + calcium channel blocker risk of severe bradycardia/AV block",
                monitoring_parameters=["Heart rate", "PR interval"],
                alternative_drugs=["diltiazem (if beta-blocker needed)", "amlodipine (if beta-blocker needed)"],
            ),
            evidence_quality="excellent",
        )
    )
    
    # Example 4: Simvastatin + CYP3A4 induction
    db.add_interaction(
        DrugInteractionPair(
            drug_1="simvastatin",
            drug_2="rifampin",
            severity=InteractionSeverity.MAJOR,
            mechanism=InteractionMechanism.CYP_INDUCTION,
            inducer_drug="rifampin",
            inducer_info=InductionInfo(
                cyp_enzyme=CYPEnzyme.CYP3A4,
                ec50=10.0,  # µM
                emax=5.0,  # Up to 5-fold induction
            ),
            affected_drug="simvastatin",
            clinical_recommendation=ClinicalRecommendation(
                recommendation_type="adjust_dose",
                description="Rifampin induces CYP3A4, increasing simvastatin metabolism",
                dose_adjustment="May need higher simvastatin dose or alternative statin (pravastatin)",
                monitoring_parameters=["LDL cholesterol", "CK"],
            ),
            evidence_quality="good",
        )
    )
    
    logger.info(f"Created default DDI database with {len(db)} interactions")
    return db
