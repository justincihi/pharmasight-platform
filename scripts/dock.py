#!/usr/bin/env python3
"""
Pharmasight docking pipeline: SMILES → PDBQT → Vina → JSON result
Usage: python3 dock.py --smiles "CC..." --receptor receptors/NMDA.pdb \
       --cx 0 --cy 0 --cz 0 --sx 25 --sy 25 --sz 25 \
       --exhaustiveness 8 --num_poses 9
"""
import argparse
import json
import sys
import os
import tempfile
import subprocess
from pathlib import Path


def smiles_to_pdbqt(smiles: str, out_path: str) -> tuple[bool, str | None]:
    """Convert SMILES string to PDBQT format using RDKit + Meeko"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation
        from meeko import PDBQTWriterLegacy

        # Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, f"Invalid SMILES: {smiles}"

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)

        # Prepare molecule for docking
        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)

        # Write PDBQT
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])

        if not is_ok:
            return False, f"PDBQT writing failed: {error_msg}"

        with open(out_path, "w") as f:
            f.write(pdbqt_string)

        return True, None
    except Exception as e:
        return False, f"SMILES to PDBQT conversion failed: {str(e)}"


def pdb_to_pdbqt(pdb_path: str, out_path: str) -> tuple[bool, str | None]:
    """Convert PDB file to PDBQT format using obabel or MGLTools"""
    try:
        # Try obabel first
        result = subprocess.run(
            ["obabel", pdb_path, "-O", out_path, "-xr"],
            capture_output=True,
            text=True,
            timeout=60,
        )
        if result.returncode == 0:
            return True, None
        return False, f"obabel conversion failed: {result.stderr}"
    except FileNotFoundError:
        return False, "obabel not found - install with: sudo apt-get install openbabel"
    except Exception as e:
        return False, f"PDB to PDBQT conversion failed: {str(e)}"


def run_vina(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    cx: float,
    cy: float,
    cz: float,
    sx: float,
    sy: float,
    sz: float,
    exhaustiveness: int,
    num_poses: int,
    out_pdbqt: str,
) -> tuple[int, str, str]:
    """Run AutoDock Vina docking"""
    cmd = [
        "vina",
        "--receptor",
        receptor_pdbqt,
        "--ligand",
        ligand_pdbqt,
        "--center_x",
        str(cx),
        "--center_y",
        str(cy),
        "--center_z",
        str(cz),
        "--size_x",
        str(sx),
        "--size_y",
        str(sy),
        "--size_z",
        str(sz),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_poses),
        "--out",
        out_pdbqt,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    return result.returncode, result.stdout, result.stderr


def parse_vina_output(stdout: str) -> list[dict]:
    """Parse Vina output to extract binding affinities"""
    poses = []
    lines = stdout.split("\n")

    for line in lines:
        parts = line.split()
        if len(parts) >= 3:
            try:
                mode = int(parts[0])
                affinity = float(parts[1])
                poses.append({"mode": mode, "affinity": affinity})
            except (ValueError, IndexError):
                continue

    return poses


def main():
    parser = argparse.ArgumentParser(
        description="Pharmasight molecular docking pipeline"
    )
    parser.add_argument("--smiles", required=True, help="SMILES string of ligand")
    parser.add_argument("--receptor", required=True, help="Path to PDB receptor file")
    parser.add_argument("--cx", type=float, default=0, help="Box center X coordinate")
    parser.add_argument("--cy", type=float, default=0, help="Box center Y coordinate")
    parser.add_argument("--cz", type=float, default=0, help="Box center Z coordinate")
    parser.add_argument("--sx", type=float, default=25, help="Box size X")
    parser.add_argument("--sy", type=float, default=25, help="Box size Y")
    parser.add_argument("--sz", type=float, default=25, help="Box size Z")
    parser.add_argument(
        "--exhaustiveness", type=int, default=8, help="Vina exhaustiveness (1-32)"
    )
    parser.add_argument(
        "--num_poses", type=int, default=9, help="Number of poses to generate"
    )

    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        ligand_pdbqt = os.path.join(tmpdir, "ligand.pdbqt")
        receptor_pdbqt = args.receptor.replace(".pdb", ".pdbqt")
        out_pdbqt = os.path.join(tmpdir, "out.pdbqt")

        # Step 1: Prepare ligand (SMILES → PDBQT)
        ok, err = smiles_to_pdbqt(args.smiles, ligand_pdbqt)
        if not ok:
            print(json.dumps({"success": False, "error": f"Ligand prep failed: {err}"}))
            sys.exit(1)

        # Step 2: Prepare receptor (convert PDB → PDBQT if needed)
        if not os.path.exists(receptor_pdbqt):
            ok, err = pdb_to_pdbqt(args.receptor, receptor_pdbqt)
            if not ok:
                print(
                    json.dumps({"success": False, "error": f"Receptor prep failed: {err}"})
                )
                sys.exit(1)

        # Step 3: Run Vina
        returncode, stdout, stderr = run_vina(
            receptor_pdbqt,
            ligand_pdbqt,
            args.cx,
            args.cy,
            args.cz,
            args.sx,
            args.sy,
            args.sz,
            args.exhaustiveness,
            args.num_poses,
            out_pdbqt,
        )

        if returncode != 0:
            print(json.dumps({"success": False, "error": stderr or "Vina failed"}))
            sys.exit(1)

        # Step 4: Parse results
        poses = parse_vina_output(stdout)
        best = poses[0] if poses else None

        print(
            json.dumps(
                {
                    "success": True,
                    "binding_affinity": best["affinity"] if best else None,
                    "poses": poses,
                    "num_poses": len(poses),
                }
            )
        )


if __name__ == "__main__":
    main()
