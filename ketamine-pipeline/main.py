import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, DataStructs


# ---------------------------------------------------------------------
# Paths and basic configuration
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"

DATA_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

PATENT_JSON_PATH = PROJECT_ROOT / "patent_examples.json"
KETAMINE_SDF_PATH = DATA_DIR / "ketamine.sdf"
KETAMINE_3D_SDF_PATH = DATA_DIR / "ketamine_3d.sdf"


# ---------------------------------------------------------------------
# Utilities: patent examples and fingerprints
# ---------------------------------------------------------------------

def load_patent_examples(path: Path):
    """
    Load patent exemplified compounds and annotations from JSON.
    """
    with path.open("r") as f:
        data = json.load(f)
    return data["examples"]


def _fp(smiles: str):
    """
    Generate a Morgan fingerprint bit vector for a SMILES.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def flag_patent_like(
    analog: dict,
    patent_examples: list,
    sim_inside: float = 0.85,
    sim_near: float = 0.70,
) -> str:
    """
    Classify analog relative to the patent example set.

    analog: dict with at least:
        'smiles', 'aryl_pattern', 'n_substitution', 'prodrug_motif'
    patent_examples: list loaded from patent_examples.json

    Returns: "inside", "near-boundary", "outside", or "invalid".

    - "inside": matches a known substitution/prodrug pattern OR
                has very high structural similarity to a patent example.
    - "near-boundary": moderately high similarity but no exact pattern match.
    - "outside": below similarity thresholds and no pattern match.
    """
    a_smi = analog["smiles"]
    a_fp = _fp(a_smi)
    if a_fp is None:
        return "invalid"

    # 1) Exact pattern / motif matching (very crude Markush proxy)
    for ex in patent_examples:
        if (
            analog.get("aryl_pattern") == ex.get("aryl_pattern")
            and analog.get("n_substitution") == ex.get("n_substitution")
            and analog.get("prodrug_motif") == ex.get("prodrug_motif")
        ):
            return "inside"

    # 2) Similarity to any patent example
    max_sim = 0.0
    for ex in patent_examples:
        e_fp = _fp(ex["smiles"])
        if e_fp is None:
            continue
        sim = DataStructs.TanimotoSimilarity(a_fp, e_fp)
        if sim > max_sim:
            max_sim = sim

    if max_sim >= sim_inside:
        return "inside"
    elif max_sim >= sim_near:
        return "near-boundary"
    else:
        return "outside"


# ---------------------------------------------------------------------
# Aryl-variant generator (placeholder skeleton)
# ---------------------------------------------------------------------

def generate_aryl_variants(ketamine_smiles: str, max_variants: int = 50):
    """
    Very simplified skeleton: given ketamine SMILES, generate a small set
    of aryl-substituted analogs.

    In a production version, you would:
      - Identify the aryl ring explicitly
      - Use reaction SMARTS or atom mapping
      - Enumerate allowed substituents and positions

    For now, we keep the SMILES the same and only vary 'aryl_pattern'
    labels to get the IP-labeling logic wired up.
    """
    base = Chem.MolFromSmiles(ketamine_smiles)
    if base is None:
        raise ValueError("Bad ketamine SMILES")

    analogs = []

    # Tiny example "library" of aryl patterns to explore
    trial_patterns = [
        {"name": "2-F_4-Cl", "description": "swap F/Cl positions"},
        {"name": "3-F_4-Cl", "description": "move F to meta"},
        {"name": "2-OCF3_4-Cl", "description": "add OCF3 at ortho"},
    ]

    for pat in trial_patterns:
        # TODO: replace this placeholder with an actual aryl-ring transform.
        # For now we keep the base SMILES; only the annotation changes.
        analog_smiles = ketamine_smiles

        analogs.append(
            {
                "id": f"ARYL-{pat['name']}",
                "smiles": analog_smiles,
                "aryl_pattern": [pat["name"]],
                "n_substitution": "N-methyl",   # or parsed from the base
                "prodrug_motif": None
            }
        )

        if len(analogs) >= max_variants:
            break

    return analogs


# ---------------------------------------------------------------------
# SDF → 3D helper functions (for docking / descriptor pipelines)
# ---------------------------------------------------------------------

def load_mol_from_sdf(path: Path):
    """
    Load the first valid molecule from an SDF file.
    """
    suppl = Chem.SDMolSupplier(str(path), sanitize=True)
    for mol in suppl:
        if mol is not None:
            return mol
    raise ValueError(f"No valid molecule found in {path}")


def rebuild_from_smiles(mol):
    """
    Rebuild the molecule from a SMILES string.

    Tries to read a 'SMILES' property from the SDF first.
    If missing, falls back to a generated SMILES.
    Adds hydrogens for 3D embedding.
    """
    smi = None
    for prop_name in mol.GetPropNames():
        if prop_name.lower() == "smiles":
            smi = mol.GetProp(prop_name).strip()
            break

    if not smi:
        smi = Chem.MolToSmiles(mol)

    new_mol = Chem.MolFromSmiles(smi)
    if new_mol is None:
        raise ValueError("Could not rebuild molecule from SMILES")

    new_mol = Chem.AddHs(new_mol)
    return new_mol


def make_3d_conformer(mol):
    """
    Generate a single 3D conformer using ETKDG and optimize with UFF.
    """
    cid = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if cid == -1:
        raise RuntimeError("3D embedding failed")

    AllChem.UFFOptimizeMolecule(mol, confId=cid)
    return mol, cid


def write_sdf(mol, conf_id, out_path: Path):
    """
    Write the molecule with the given conformer ID to an SDF file.
    """
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol, confId=conf_id)
    writer.close()


# ---------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------

def main():
    # 1) Define a ketamine-like SMILES (can be updated to match your SDF)
    ketamine_smiles = "CCN(C1CCCCC1=O)c2cccc(F)c2Cl"  # arylcyclohexylamine template

    # 2) Load patent examples
    patent_examples = load_patent_examples(PATENT_JSON_PATH)

    # 3) Generate a small set of aryl variants
    analogs = generate_aryl_variants(ketamine_smiles, max_variants=20)

    # 4) Score each analog relative to the patent example space
    scored = []
    for a in analogs:
        label = flag_patent_like(a, patent_examples)
        a["ip_label"] = label
        scored.append(a)

    # 5) Print simple summary and write results to JSON for your ADMET pipeline
    print("Analog IP labels:")
    for r in scored:
        print(r["id"], r["ip_label"], r["smiles"])

    out_json = RESULTS_DIR / "ketamine_aryl_analogs_ip_labels.json"
    with out_json.open("w") as f:
        json.dump({"analogs": scored}, f, indent=2)

    print(f"\nWrote IP-labeled analogs to: {out_json}")

    # 6) Optional: if ketamine.sdf is present, generate a 3D SDF
    if KETAMINE_SDF_PATH.exists():
        print("\nFound ketamine.sdf; generating 3D conformer...")
        mol0 = load_mol_from_sdf(KETAMINE_SDF_PATH)
        mol = rebuild_from_smiles(mol0)
        mol3d, cid = make_3d_conformer(mol)
        write_sdf(mol3d, cid, KETAMINE_3D_SDF_PATH)
        print(f"Wrote 3D SDF to: {KETAMINE_3D_SDF_PATH}")
    else:
        print("\nNo ketamine.sdf found in data/, skipping 3D generation.")


if __name__ == "__main__":
    main()
