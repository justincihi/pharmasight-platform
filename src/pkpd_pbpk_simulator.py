"""
PKPD/PBPK Population Simulation Module
Simulates pharmacokinetics and pharmacodynamics across patient populations
Implements physiologically-based pharmacokinetic modeling
"""

import numpy as np
from scipy.integrate import odeint
from scipy.stats import norm, lognorm
from typing import Dict, List, Tuple, Optional
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PatientPopulation:
    """Generate virtual patient populations with physiological variability"""
    
    def __init__(self, n_patients: int = 100):
        self.n_patients = n_patients
        self.patients = []
        
    def generate_population(self, 
                          age_range: Tuple[int, int] = (18, 80),
                          weight_mean: float = 70.0,
                          weight_std: float = 15.0,
                          conditions: Optional[List[str]] = None) -> List[Dict]:
        """
        Generate virtual patient population with physiological parameters
        
        Args:
            age_range: (min_age, max_age)
            weight_mean: Mean body weight (kg)
            weight_std: Standard deviation of body weight
            conditions: List of pre-existing conditions to simulate
        """
        np.random.seed(42)  # For reproducibility
        
        patients = []
        for i in range(self.n_patients):
            # Basic demographics
            age = np.random.randint(age_range[0], age_range[1] + 1)
            weight = max(40, np.random.normal(weight_mean, weight_std))
            height = np.random.normal(170, 10)  # cm
            sex = np.random.choice(['M', 'F'])
            
            # Calculate BSA (Body Surface Area) - Mosteller formula
            bsa = np.sqrt((height * weight) / 3600)
            
            # Physiological parameters
            patient = {
                'id': f'P{i+1:03d}',
                'age': age,
                'weight': round(weight, 1),
                'height': round(height, 1),
                'bsa': round(bsa, 2),
                'sex': sex,
                
                # Organ blood flows (L/h) - scaled by weight
                'cardiac_output': round(5.0 * (weight / 70), 2),
                'liver_blood_flow': round(1.35 * (weight / 70), 2),
                'kidney_blood_flow': round(1.1 * (weight / 70), 2),
                'brain_blood_flow': round(0.7 * (weight / 70), 2),
                
                # Organ volumes (L) - scaled by weight
                'plasma_volume': round(3.0 * (weight / 70), 2),
                'liver_volume': round(1.8 * (weight / 70), 2),
                'kidney_volume': round(0.3 * (weight / 70), 2),
                'fat_volume': round(15.0 * (weight / 70), 2),
                'muscle_volume': round(28.0 * (weight / 70), 2),
                
                # Metabolic parameters (vary with age and conditions)
                'creatinine_clearance': self._calculate_crcl(age, weight, sex),
                'liver_function': self._calculate_liver_function(age, conditions),
                'cyp3a4_activity': round(np.random.lognormal(0, 0.3), 2),
                'cyp2d6_activity': round(np.random.lognormal(0, 0.4), 2),
                
                # Conditions
                'conditions': self._assign_conditions(age, conditions),
            }
            
            patients.append(patient)
        
        self.patients = patients
        return patients
    
    def _calculate_crcl(self, age: int, weight: float, sex: str) -> float:
        """Calculate creatinine clearance (Cockcroft-Gault)"""
        # Assuming serum creatinine = 1.0 mg/dL
        crcl = ((140 - age) * weight) / (72 * 1.0)
        if sex == 'F':
            crcl *= 0.85
        return round(crcl, 1)
    
    def _calculate_liver_function(self, age: int, conditions: Optional[List[str]]) -> float:
        """Calculate liver function score (0-1)"""
        function = 1.0
        
        # Age-related decline
        if age > 65:
            function *= 0.9
        if age > 75:
            function *= 0.85
        
        # Condition-related impairment
        if conditions and 'liver_disease' in conditions:
            function *= 0.6
        
        return round(function, 2)
    
    def _assign_conditions(self, age: int, requested_conditions: Optional[List[str]]) -> List[str]:
        """Assign pre-existing conditions based on age and requested conditions"""
        conditions = []
        
        if requested_conditions:
            # Assign requested conditions with probability
            for condition in requested_conditions:
                if np.random.random() < 0.3:  # 30% chance
                    conditions.append(condition)
        
        # Age-related conditions
        if age > 60 and np.random.random() < 0.2:
            conditions.append('hypertension')
        if age > 65 and np.random.random() < 0.15:
            conditions.append('diabetes')
        
        return conditions


class PKPDSimulator:
    """Pharmacokinetic and Pharmacodynamic simulator"""
    
    def __init__(self):
        self.time_points = np.linspace(0, 24, 100)  # 24 hours, 100 time points
    
    def simulate_one_compartment_pk(self, 
                                    dose: float,
                                    patient: Dict,
                                    ka: float = 1.0,
                                    cl: float = 10.0,
                                    vd: float = 50.0) -> Dict:
        """
        Simulate one-compartment PK model
        
        Args:
            dose: Dose in mg
            patient: Patient parameters
            ka: Absorption rate constant (1/h)
            cl: Clearance (L/h)
            vd: Volume of distribution (L)
        """
        # Adjust parameters based on patient characteristics
        cl_adjusted = cl * (patient['weight'] / 70) * patient['liver_function']
        vd_adjusted = vd * (patient['weight'] / 70)
        
        # Adjust for renal function if drug is renally cleared
        if patient['creatinine_clearance'] < 60:
            cl_adjusted *= (patient['creatinine_clearance'] / 100)
        
        ke = cl_adjusted / vd_adjusted  # Elimination rate constant
        
        # Analytical solution for one-compartment model with first-order absorption
        def concentration(t):
            if t == 0:
                return 0
            c = (dose * ka) / (vd_adjusted * (ka - ke)) * \
                (np.exp(-ke * t) - np.exp(-ka * t))
            return max(0, c)
        
        concentrations = [concentration(t) for t in self.time_points]
        
        # Calculate PK parameters
        cmax = max(concentrations)
        tmax = self.time_points[concentrations.index(cmax)]
        auc = np.trapz(concentrations, self.time_points)
        
        return {
            'time': self.time_points.tolist(),
            'concentration': concentrations,
            'cmax': round(cmax, 2),
            'tmax': round(tmax, 2),
            'auc': round(auc, 2),
            'clearance': round(cl_adjusted, 2),
            'volume_distribution': round(vd_adjusted, 2),
            'half_life': round(0.693 / ke, 2),
        }
    
    def simulate_two_compartment_pk(self,
                                    dose: float,
                                    patient: Dict,
                                    ka: float = 1.0,
                                    cl: float = 10.0,
                                    v1: float = 30.0,
                                    v2: float = 20.0,
                                    q: float = 5.0) -> Dict:
        """
        Simulate two-compartment PK model (central and peripheral)
        
        Args:
            dose: Dose in mg
            patient: Patient parameters
            ka: Absorption rate constant
            cl: Clearance from central compartment
            v1: Central compartment volume
            v2: Peripheral compartment volume
            q: Inter-compartmental clearance
        """
        # Adjust parameters
        cl_adj = cl * (patient['weight'] / 70) * patient['liver_function']
        v1_adj = v1 * (patient['weight'] / 70)
        v2_adj = v2 * (patient['weight'] / 70)
        q_adj = q * (patient['weight'] / 70)
        
        # Define ODE system
        def two_comp_model(y, t, ka, ke, k12, k21):
            depot, central, peripheral = y
            
            ddepot = -ka * depot
            dcentral = ka * depot - ke * central - k12 * central + k21 * peripheral
            dperipheral = k12 * central - k21 * peripheral
            
            return [ddepot, dcentral, dperipheral]
        
        # Calculate rate constants
        ke = cl_adj / v1_adj
        k12 = q_adj / v1_adj
        k21 = q_adj / v2_adj
        
        # Initial conditions
        y0 = [dose, 0, 0]  # [depot, central, peripheral]
        
        # Solve ODE
        solution = odeint(two_comp_model, y0, self.time_points, args=(ka, ke, k12, k21))
        
        # Concentration in central compartment
        concentrations = solution[:, 1] / v1_adj
        
        cmax = max(concentrations)
        tmax = self.time_points[list(concentrations).index(cmax)]
        auc = np.trapz(concentrations, self.time_points)
        
        return {
            'time': self.time_points.tolist(),
            'concentration': concentrations.tolist(),
            'cmax': round(cmax, 2),
            'tmax': round(tmax, 2),
            'auc': round(auc, 2),
            'clearance': round(cl_adj, 2),
            'vss': round(v1_adj + v2_adj, 2),  # Volume at steady state
        }
    
    def simulate_pd_effect(self, 
                          concentrations: List[float],
                          ec50: float = 10.0,
                          emax: float = 100.0,
                          hill_coefficient: float = 1.0) -> List[float]:
        """
        Simulate pharmacodynamic effect using Emax model
        
        Args:
            concentrations: Drug concentrations over time
            ec50: Concentration producing 50% of maximum effect
            emax: Maximum effect
            hill_coefficient: Hill coefficient (sigmoidicity)
        """
        effects = []
        for c in concentrations:
            effect = (emax * (c ** hill_coefficient)) / \
                    ((ec50 ** hill_coefficient) + (c ** hill_coefficient))
            effects.append(effect)
        
        return effects


class PBPKSimulator:
    """Physiologically-Based Pharmacokinetic simulator"""
    
    def __init__(self):
        self.time_points = np.linspace(0, 24, 200)
    
    def simulate_pbpk(self,
                     dose: float,
                     patient: Dict,
                     drug_properties: Dict) -> Dict:
        """
        Simulate PBPK model with multiple organs
        
        Args:
            dose: Dose in mg
            patient: Patient physiological parameters
            drug_properties: Drug-specific properties (logP, fu, etc.)
        """
        # Extract drug properties
        logp = drug_properties.get('logp', 2.0)
        fu = drug_properties.get('fraction_unbound', 0.1)
        bp_ratio = drug_properties.get('blood_plasma_ratio', 1.0)
        
        # Calculate partition coefficients (Poulin-Theil method approximation)
        kp_liver = 10 ** (0.72 * logp - 0.16)
        kp_kidney = 10 ** (0.65 * logp - 0.10)
        kp_brain = 10 ** (0.55 * logp - 0.25)
        kp_fat = 10 ** (1.2 * logp + 0.3)
        kp_muscle = 10 ** (0.45 * logp - 0.20)
        
        # Define PBPK ODE system
        def pbpk_model(y, t):
            # Compartments: gut, plasma, liver, kidney, brain, fat, muscle
            gut, plasma, liver, kidney, brain, fat, muscle = y
            
            # Absorption
            ka = 1.0  # absorption rate constant
            absorption = ka * gut
            
            # Organ concentrations
            c_plasma = plasma / patient['plasma_volume']
            c_liver = liver / patient['liver_volume']
            c_kidney = kidney / patient['kidney_volume']
            c_brain = brain / (patient['weight'] * 0.02)  # brain ~2% of body weight
            c_fat = fat / patient['fat_volume']
            c_muscle = muscle / patient['muscle_volume']
            
            # Blood flows
            q_liver = patient['liver_blood_flow']
            q_kidney = patient['kidney_blood_flow']
            q_brain = patient['brain_blood_flow']
            q_fat = patient['cardiac_output'] * 0.05
            q_muscle = patient['cardiac_output'] * 0.17
            
            # Hepatic clearance
            cl_int = 10.0 * patient['liver_function']  # Intrinsic clearance
            cl_hepatic = (q_liver * fu * cl_int) / (q_liver + fu * cl_int)
            
            # Renal clearance
            gfr = patient['creatinine_clearance'] / 60  # Convert to L/h
            cl_renal = gfr * fu
            
            # Differential equations
            dgut = -absorption
            
            dplasma = (absorption +
                      q_liver * (c_liver / kp_liver - c_plasma) +
                      q_kidney * (c_kidney / kp_kidney - c_plasma) +
                      q_brain * (c_brain / kp_brain - c_plasma) +
                      q_fat * (c_fat / kp_fat - c_plasma) +
                      q_muscle * (c_muscle / kp_muscle - c_plasma))
            
            dliver = q_liver * (c_plasma - c_liver / kp_liver) - cl_hepatic * c_liver
            dkidney = q_kidney * (c_plasma - c_kidney / kp_kidney) - cl_renal * c_kidney
            dbrain = q_brain * (c_plasma - c_brain / kp_brain)
            dfat = q_fat * (c_plasma - c_fat / kp_fat)
            dmuscle = q_muscle * (c_plasma - c_muscle / kp_muscle)
            
            return [dgut, dplasma, dliver, dkidney, dbrain, dfat, dmuscle]
        
        # Initial conditions
        y0 = [dose, 0, 0, 0, 0, 0, 0]
        
        # Solve ODE system
        try:
            solution = odeint(pbpk_model, y0, self.time_points)
            
            # Extract plasma concentrations
            plasma_conc = solution[:, 1] / patient['plasma_volume']
            liver_conc = solution[:, 2] / patient['liver_volume']
            brain_conc = solution[:, 3] / (patient['weight'] * 0.02)
            
            return {
                'time': self.time_points.tolist(),
                'plasma_concentration': plasma_conc.tolist(),
                'liver_concentration': liver_conc.tolist(),
                'brain_concentration': brain_conc.tolist(),
                'cmax_plasma': round(max(plasma_conc), 2),
                'auc_plasma': round(np.trapz(plasma_conc, self.time_points), 2),
                'brain_penetration_ratio': round(max(brain_conc) / max(plasma_conc), 3),
            }
        except Exception as e:
            logger.error(f"PBPK simulation failed: {e}")
            return {'error': str(e)}


class PopulationPKPDAnalyzer:
    """Analyze PK/PD across patient populations"""
    
    def __init__(self):
        self.pk_simulator = PKPDSimulator()
        self.pbpk_simulator = PBPKSimulator()
    
    def run_population_simulation(self,
                                 dose: float,
                                 patients: List[Dict],
                                 drug_properties: Dict,
                                 model_type: str = 'one_compartment') -> Dict:
        """
        Run PK simulation across entire patient population
        
        Args:
            dose: Dose in mg
            patients: List of patient parameters
            drug_properties: Drug properties
            model_type: 'one_compartment', 'two_compartment', or 'pbpk'
        """
        results = []
        
        for patient in patients:
            if model_type == 'one_compartment':
                pk_result = self.pk_simulator.simulate_one_compartment_pk(
                    dose, patient,
                    ka=drug_properties.get('ka', 1.0),
                    cl=drug_properties.get('cl', 10.0),
                    vd=drug_properties.get('vd', 50.0)
                )
            elif model_type == 'pbpk':
                pk_result = self.pbpk_simulator.simulate_pbpk(
                    dose, patient, drug_properties
                )
            else:
                pk_result = self.pk_simulator.simulate_one_compartment_pk(
                    dose, patient
                )
            
            pk_result['patient_id'] = patient['id']
            pk_result['patient_age'] = patient['age']
            pk_result['patient_weight'] = patient['weight']
            pk_result['patient_conditions'] = patient['conditions']
            
            results.append(pk_result)
        
        # Calculate population statistics
        cmaxs = [r['cmax'] for r in results if 'cmax' in r]
        aucs = [r['auc'] for r in results if 'auc' in r]
        
        population_stats = {
            'n_patients': len(patients),
            'cmax_mean': round(np.mean(cmaxs), 2),
            'cmax_std': round(np.std(cmaxs), 2),
            'cmax_median': round(np.median(cmaxs), 2),
            'cmax_range': (round(min(cmaxs), 2), round(max(cmaxs), 2)),
            'auc_mean': round(np.mean(aucs), 2),
            'auc_std': round(np.std(aucs), 2),
            'auc_median': round(np.median(aucs), 2),
            'auc_range': (round(min(aucs), 2), round(max(aucs), 2)),
            'variability_coefficient': round((np.std(aucs) / np.mean(aucs)) * 100, 1),
        }
        
        return {
            'individual_results': results,
            'population_statistics': population_stats,
            'dose': dose,
            'model_type': model_type,
        }
    
    def identify_high_risk_patients(self, 
                                   simulation_results: Dict,
                                   therapeutic_window: Tuple[float, float]) -> List[Dict]:
        """
        Identify patients outside therapeutic window
        
        Args:
            simulation_results: Results from population simulation
            therapeutic_window: (min_concentration, max_concentration)
        """
        min_conc, max_conc = therapeutic_window
        high_risk_patients = []
        
        for result in simulation_results['individual_results']:
            cmax = result.get('cmax', 0)
            
            risk_factors = []
            if cmax < min_conc:
                risk_factors.append('subtherapeutic')
            if cmax > max_conc:
                risk_factors.append('toxic')
            
            if risk_factors:
                high_risk_patients.append({
                    'patient_id': result['patient_id'],
                    'cmax': cmax,
                    'risk_factors': risk_factors,
                    'conditions': result.get('patient_conditions', []),
                    'recommended_dose_adjustment': self._calculate_dose_adjustment(
                        cmax, (min_conc + max_conc) / 2, result.get('dose', 100)
                    ),
                })
        
        return high_risk_patients
    
    def _calculate_dose_adjustment(self, current_cmax: float, target_cmax: float, current_dose: float) -> float:
        """Calculate recommended dose adjustment"""
        adjustment_factor = target_cmax / current_cmax
        recommended_dose = current_dose * adjustment_factor
        return round(recommended_dose, 1)


# Singleton instances
_pk_simulator = None
_pbpk_simulator = None
_population_analyzer = None

def get_pk_simulator() -> PKPDSimulator:
    global _pk_simulator
    if _pk_simulator is None:
        _pk_simulator = PKPDSimulator()
    return _pk_simulator

def get_pbpk_simulator() -> PBPKSimulator:
    global _pbpk_simulator
    if _pbpk_simulator is None:
        _pbpk_simulator = PBPKSimulator()
    return _pbpk_simulator

def get_population_analyzer() -> PopulationPKPDAnalyzer:
    global _population_analyzer
    if _population_analyzer is None:
        _population_analyzer = PopulationPKPDAnalyzer()
    return _population_analyzer

