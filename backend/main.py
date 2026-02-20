from __future__ import annotations
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from typing import Dict, Any
import logging

from backend.api.schemas import (
    SimulateRequest, SimulateResult, DDIScreenRequest,
    DDIScreenResponse, DrugLookupResponse, DrugParams
)

# Import internal modules
from pharmasight_pk.popPK import PopPKPredictor, CovariateProfile, Sex
from pharmasight_pk.virtual_patient import VirtualPatientFactory
from pharmasight_pk.ddi import create_default_ddi_database, DDIScreener
from pharmasight_pk.models.oral_one_compartment import OralOneCompartment
from pharmasight_pk.models.one_compartment import OneCompartmentIV
from pharmasight_pk.models.base import PKParameters, AdminRoute

logger = logging.getLogger("pharmasight_api")
logging.basicConfig(level=logging.INFO)

app = FastAPI(title="PharmaSight PK API", version="0.1")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Minimal in-memory drug registry (can be extended or replaced by DB)
drug_registry: Dict[str, Dict[str, float]] = {
    "example_drug": {"cl": 10.0, "vc": 50.0, "ka": 0.8, "f": 0.8},
}

# Initialize services
vfactory = VirtualPatientFactory()
ddi_db = create_default_ddi_database()
screener = DDIScreener(ddi_db)


@app.get("/ping")
def ping():
    return {"status": "ok"}


@app.get("/pk/drug-lookup", response_model=DrugLookupResponse)
def drug_lookup(drug_name: str):
    key = drug_name.lower()
    if key not in drug_registry:
        raise HTTPException(status_code=404, detail="Drug not found in registry")
    params = drug_registry[key]
    return DrugLookupResponse(drug_name=drug_name, base_parameters=DrugParams(**params))


@app.post("/pk/simulate", response_model=SimulateResult)
def simulate(req: SimulateRequest):
    # Resolve drug parameters
    if req.drug_params:
        base = req.drug_params
    elif req.drug_name:
        key = req.drug_name.lower()
        if key not in drug_registry:
            raise HTTPException(status_code=404, detail="Drug not found in registry")
        base = DrugParams(**drug_registry[key])
    else:
        raise HTTPException(status_code=400, detail="Provide drug_name or drug_params")

    # Create covariate profile
    cov = None
    if req.covariates:
        cov_in = req.covariates
        # Map sex string to Sex enum
        sex_val = Sex.MALE if cov_in.sex.lower() == "male" else Sex.FEMALE if cov_in.sex.lower() == "female" else Sex.OTHER
        cov = CovariateProfile(
            age_years=cov_in.age_years,
            weight_kg=cov_in.weight_kg,
            height_cm=cov_in.height_cm,
            sex=sex_val,
            egfr=cov_in.egfr or 90.0,
            albumin_g_dl=cov_in.albumin_g_dl or 4.0,
            liver_function=(cov_in.liver_function or "normal"),
            cyp3a4_phenotype=(cov_in.cyp3a4_phenotype or "extensive"),
            cyp2d6_phenotype=(cov_in.cyp2d6_phenotype or "extensive"),
            cyp2c9_phenotype=(cov_in.cyp2c9_phenotype or "extensive"),
            smoking=bool(cov_in.smoking),
            pregnancy_trimester=cov_in.pregnancy_trimester,
        )

    # Use PopPKPredictor to get individualised parameters
    predictor = PopPKPredictor(req.drug_name or "custom")
    predictor.set_base_parameters(cl=base.cl, vc=base.vc, ka=base.ka, f=base.f, qp=base.qp, vp=base.vp)
    pop_pred = predictor.predict(cov if cov is not None else CovariateProfile(age_years=40, weight_kg=70, sex=Sex.MALE))

    # Build PKParameters for model
    pkparams = PKParameters(CL=pop_pred.predicted_cl, Vc=pop_pred.predicted_vc, Ka=pop_pred.predicted_ka, F=pop_pred.predicted_f)

    # Instantiate appropriate model
    route = req.route.lower()
    if route in {AdminRoute.ORAL.value, "oral", "po"}:
        model = OralOneCompartment(pkparams)
    else:
        model = OneCompartmentIV(pkparams)

    # Run simulation
    try:
        sim = model.simulate(
            dose_mg=req.dose_mg,
            route=route,
            t_end=req.t_end_h,
            n_points=req.n_points,
            covariates=(cov.to_dict() if cov else None),
        )
    except Exception as e:
        logger.exception("Simulation failed")
        raise HTTPException(status_code=500, detail=str(e))

    # Prepare response
    compartments = {}
    comp_names = sim.meta.get("compartment_names") or []
    for idx, name in enumerate(comp_names):
        compartments[name] = sim.y[idx].tolist()

    return SimulateResult(t=sim.t.tolist(), compartments=compartments, meta=sim.meta)


@app.post("/pk/screen-interactions", response_model=DDIScreenResponse)
def screen_interactions(req: DDIScreenRequest):
    # Map severity string to InteractionSeverity inside screener; screener expects enum, but its method accepts string via internal handling
    result = screener.screen_medication_list(req.medications)
    return DDIScreenResponse(**result)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("backend.main:app", host="127.0.0.1", port=8000, log_level="info")
