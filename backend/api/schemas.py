from __future__ import annotations
from pydantic import BaseModel, Field, validator
from typing import Optional, List, Dict, Any


class DrugParams(BaseModel):
    cl: float = Field(..., description="Clearance (L/h) at reference 70 kg")
    vc: float = Field(..., description="Central volume (L) at reference 70 kg")
    ka: Optional[float] = Field(None, description="Absorption rate constant (1/h)")
    f: float = Field(1.0, description="Bioavailability (0-1)")
    qp: Optional[float] = None
    vp: Optional[float] = None


class CovariatesIn(BaseModel):
    age_years: float
    weight_kg: float
    height_cm: Optional[float] = None
    sex: str
    egfr: Optional[float] = None
    albumin_g_dl: Optional[float] = None
    liver_function: Optional[str] = None
    cyp3a4_phenotype: Optional[str] = None
    cyp2d6_phenotype: Optional[str] = None
    cyp2c9_phenotype: Optional[str] = None
    smoking: Optional[bool] = False
    pregnancy_trimester: Optional[int] = None

    @validator("sex")
    def sex_lower(cls, v):
        return v.lower()


class SimulateRequest(BaseModel):
    dose_mg: float
    route: str = "oral"
    t_end_h: float = 24.0
    n_points: int = 200
    drug_name: Optional[str] = None
    drug_params: Optional[DrugParams] = None
    covariates: Optional[CovariatesIn] = None
    fed_state: Optional[bool] = False


class SimulateResult(BaseModel):
    t: List[float]
    compartments: Dict[str, List[float]]
    meta: Dict[str, Any]


class DDIScreenRequest(BaseModel):
    medications: List[str]
    min_severity: Optional[str] = "moderate"


class DDIScreenResponse(BaseModel):
    medication_count: int
    interaction_count: int
    severity_summary: Dict[str, int]
    alerts: List[Dict[str, Any]]
    timestamp: str


class DrugLookupResponse(BaseModel):
    drug_name: str
    base_parameters: DrugParams
