# DDI Analysis Fix

# In this file, we will implement the necessary fixes for the DDI analysis module.
# We will focus on providing detailed risk assessment, synergy analysis, and dosage recommendations.

# We will start by creating a more comprehensive drug interaction database.

COMPREHENSIVE_INTERACTION_DB = {
    ("sertraline", "tramadol"): {
        "risk_level": "High",
        "mechanism": "Increased risk of serotonin syndrome due to additive serotonergic effects.",
        "synergy": "Potential for dangerous synergy, leading to life-threatening serotonin toxicity.",
        "recommendation": "Avoid co-administration. If necessary, monitor closely for symptoms of serotonin syndrome (e.g., agitation, hallucinations, rapid heart rate).",
        "dosage_adjustment": "Reduce tramadol dosage by at least 50% and monitor patient response carefully."
    },
    ("fluoxetine", "tramadol"): {
        "risk_level": "High",
        "mechanism": "Fluoxetine is a potent CYP2D6 inhibitor, which can increase tramadol concentrations, raising the risk of seizures and serotonin syndrome.",
        "synergy": "Negative synergy, increasing the likelihood of adverse effects.",
        "recommendation": "Avoid this combination. Consider an alternative analgesic or antidepressant.",
        "dosage_adjustment": "If unavoidable, a significant reduction in tramadol dose is required, with close clinical and ECG monitoring."
    },
    ("alprazolam", "oxycodone"): {
        "risk_level": "High",
        "mechanism": "Additive CNS and respiratory depression.",
        "synergy": "Dangerous synergy that can lead to profound sedation, respiratory depression, coma, and death.",
        "recommendation": "Reserve concomitant prescribing for patients for whom alternative treatment options are inadequate. Limit dosages and durations to the minimum possible.",
        "dosage_adjustment": "Use the lowest effective doses. Monitor for sedation and respiratory depression."
    },
    ("diazepam", "morphine"): {
        "risk_level": "High",
        "mechanism": "Both drugs cause CNS depression, and their effects are additive.",
        "synergy": "High potential for life-threatening respiratory depression.",
        "recommendation": "Strictly avoid in most cases. In exceptional circumstances, use with extreme caution under close medical supervision.",
        "dosage_adjustment": "Significant dose reduction of both agents is necessary."
    },
    ("ketamine", "alprazolam"): {
        "risk_level": "Moderate",
        "mechanism": "Enhanced sedative and dissociative effects.",
        "synergy": "Increased risk of confusion, dizziness, and impaired motor coordination.",
        "recommendation": "Use with caution. Advise patients against operating heavy machinery or driving.",
        "dosage_adjustment": "Consider a lower dose of alprazolam."
    },
    ("mdma", "sertraline"): {
        "risk_level": "Very High",
        "mechanism": "Extreme risk of serotonin syndrome.",
        "synergy": "Life-threatening synergistic interaction.",
        "recommendation": "Absolutely contraindicated.",
        "dosage_adjustment": "N/A - Do not co-administer."
    },
    ("psilocybin", "fluoxetine"): {
        "risk_level": "Moderate",
        "mechanism": "Fluoxetine may blunt the subjective effects of psilocybin, but the interaction is not fully understood.",
        "synergy": "Unpredictable. May either reduce or alter the psychedelic experience.",
        "recommendation": "Use with caution and under supervision. Patients should be aware that the expected effects of psilocybin may be altered.",
        "dosage_adjustment": "No standard adjustment; monitor patient response."
    },
    ("aripiprazole", "fluoxetine"): {
        "risk_level": "Moderate",
        "mechanism": "Fluoxetine (a strong CYP2D6 inhibitor) can increase aripiprazole levels, increasing the risk of side effects like akathisia.",
        "synergy": "Increased likelihood of aripiprazole-related adverse effects.",
        "recommendation": "Monitor for aripiprazole side effects.",
        "dosage_adjustment": "Reduce aripiprazole dosage by 50% when co-administered with fluoxetine."
    }
}

def get_detailed_interaction_info(comp1, comp2):
    """Retrieves detailed interaction information for two compounds."""
    key1 = (comp1.lower(), comp2.lower())
    key2 = (comp2.lower(), comp1.lower())
    return COMPREHENSIVE_INTERACTION_DB.get(key1, COMPREHENSIVE_INTERACTION_DB.get(key2))


