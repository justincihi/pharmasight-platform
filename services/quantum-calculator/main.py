from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import Dict, Any
from pyscf import gto, scf, dft
import numpy as np

app = FastAPI()

class QuantumCalculationRequest(BaseModel):
    mol_geometry: str  # In Z-matrix or XYZ format
    basis: str = 'sto-3g'
    xc_functional: str = 'b3lyp'

class QuantumPropertyCalculator:
    def calculate_electronic_properties(self, req: QuantumCalculationRequest) -> Dict[str, Any]:
        """Calculate quantum mechanical properties using PySCF."""
        try:
            mol = gto.Mole()
            mol.atom = req.mol_geometry
            mol.basis = req.basis
            mol.build()
            
            # Perform DFT calculation
            mf = dft.RKS(mol)
            mf.xc = req.xc_functional
            energy = mf.kernel()
            
            homo_lumo_gap = self._calculate_gap(mf)
            dipole_moment = self._calculate_dipole(mf)
            
            return {
                'total_energy_hartree': energy,
                'homo_lumo_gap_ev': homo_lumo_gap,
                'dipole_moment_debye': dipole_moment,
                'calculation_details': {
                    'basis_set': req.basis,
                    'xc_functional': req.xc_functional,
                    'num_basis_functions': mf.mol.nbas,
                }
            }
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"PySCF calculation failed: {str(e)}")

    def _calculate_gap(self, mf) -> float:
        """Calculate HOMO-LUMO gap in eV."""
        homo_idx = mf.mo_occ.nonzero()[0][-1]
        lumo_idx = homo_idx + 1
        
        # Convert from Hartrees to eV
        gap_hartree = mf.mo_energy[lumo_idx] - mf.mo_energy[homo_idx]
        return gap_hartree * 27.2114

    def _calculate_dipole(self, mf) -> float:
        """Calculate total dipole moment in Debye."""
        dipole_components = mf.dip_moment()
        return np.linalg.norm(dipole_components)

calculator = QuantumPropertyCalculator()

@app.post("/calculate", response_model=Dict[str, Any])
async def calculate_quantum_properties(request: QuantumCalculationRequest):
    """
    Performs a quantum chemistry calculation to determine electronic properties.
    `mol_geometry` should be a string in XYZ format, e.g.,
    "O 0 0 0; H 0 1 0; H 0 0 1"
    """
    return calculator.calculate_electronic_properties(request)
