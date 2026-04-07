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
        import shutil
        
        # If input is already PDBQT, just copy it
        if pdb_path.endswith('.pdbqt'):
            if pdb_path != out_path:  # Only copy if different paths
                shutil.copy(pdb_path, out_path)
            return True, None
        
        # Try obabel first (only for .pdb files)
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
) -> tuple[bool, str | None, list[dict]]:
    """Run AutoDock Vina docking using Python API"""
    try:
        from vina import Vina
        import sys
        from io import StringIO
        
        # Suppress Vina stdout
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        
        try:
            # Initialize Vina
            vina = Vina(sf_name='vina')
            
            # Set receptor
            vina.set_receptor(receptor_pdbqt)
            
            # Set ligand
            vina.set_ligand_from_file(ligand_pdbqt)
            
            # Set search space
            vina.compute_vina_maps(center=[cx, cy, cz], box_size=[sx, sy, sz])
            
            # Run docking
            vina.dock(exhaustiveness=exhaustiveness, n_poses=num_poses)
        finally:
            # Restore stdout
            sys.stdout = old_stdout
        
        # Extract results
        poses = []
        pose_list = vina.poses()  # Call as method
        energy_list = vina.energies()  # Call as method
        
        # Ensure we don't exceed available results
        num_results = min(len(pose_list), len(energy_list))
        for i in range(num_results):
            affinity = energy_list[i][0]  # Binding affinity is first element
            poses.append({
                "mode": i + 1,
                "affinity": float(affinity),
                "rmsd": 0.0  # RMSD not easily available from Python API
            })
        
        return True, None, poses
        
    except Exception as e:
        return False, f"Vina docking failed: {str(e)}", []


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
        
        # Handle both .pdb and .pdbqt input files
        if args.receptor.endswith('.pdbqt'):
            receptor_pdbqt = args.receptor
            # Verify PDBQT file exists
            if not os.path.exists(receptor_pdbqt):
                print(json.dumps({"success": False, "error": f"Receptor file not found: {receptor_pdbqt}"}))
                sys.exit(1)
        else:
            receptor_pdbqt = args.receptor.replace(".pdb", ".pdbqt")
            # Step 2: Prepare receptor (convert PDB → PDBQT if needed)
            if not os.path.exists(receptor_pdbqt):
                ok, err = pdb_to_pdbqt(args.receptor, receptor_pdbqt)
                if not ok:
                    print(json.dumps({"success": False, "error": f"Receptor prep failed: {err}"}))
                    sys.exit(1)
        
        # Step 1: Prepare ligand (SMILES → PDBQT)
        ok, err = smiles_to_pdbqt(args.smiles, ligand_pdbqt)
        if not ok:
            print(json.dumps({"success": False, "error": f"Ligand prep failed: {err}"}))
            sys.exit(1)

        # Step 3: Run Vina
        ok, err, poses = run_vina(
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
        )

        if not ok:
            print(json.dumps({"success": False, "error": err}))
            sys.exit(1)

        # Step 4: Return results
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
