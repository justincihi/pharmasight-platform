"""
Analog Research Integration Module
Integrates analog discoveries into the research findings system for IP protection and visibility
"""

import json
import os
from datetime import datetime
from typing import List, Dict, Any

class AnalogResearchIntegration:
    """Integrates analog discoveries into research findings system"""
    
    def __init__(self, research_findings_path: str = None):
        if research_findings_path is None:
            research_findings_path = os.path.join(
                os.path.dirname(__file__), 
                '../enhanced_research_findings.json'
            )
        self.research_findings_path = research_findings_path
        
    def load_research_findings(self) -> Dict[str, Any]:
        """Load existing research findings"""
        try:
            with open(self.research_findings_path, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            return {"findings": [], "total_findings": 0}
    
    def save_research_findings(self, data: Dict[str, Any]):
        """Save research findings"""
        with open(self.research_findings_path, 'w') as f:
            json.dump(data, f, indent=2)
    
    def convert_analog_to_research_finding(self, analog: Dict[str, Any], 
                                          parent_compound: str,
                                          session_timestamp: str) -> Dict[str, Any]:
        """Convert an analog discovery to a research finding format"""
        
        # Calculate confidence based on drug likeness and IP opportunity
        confidence = (analog['drug_likeness'] + analog['patent_opportunity_score']) / 2
        
        # Determine patent potential based on IP opportunity
        if analog['patent_opportunity_score'] >= 90:
            patent_potential = "Very High"
        elif analog['patent_opportunity_score'] >= 75:
            patent_potential = "High"
        elif analog['patent_opportunity_score'] >= 50:
            patent_potential = "Medium"
        else:
            patent_potential = "Low"
        
        # Create research finding
        finding = {
            "title": f"Novel {parent_compound} Analog Discovery: {analog['name']}",
            "description": f"Computationally generated analog of {parent_compound} using {analog['transformation_applied']} transformation strategy. Shows {analog['therapeutic_potential'].lower()} therapeutic potential with {analog['drug_likeness']}% drug-likeness score.",
            "compound_id": analog['name'].lower().replace(' ', '_'),
            "compound_name": analog['name'],
            "smiles": analog['smiles'],
            "parent_compound": parent_compound,
            "transformation_applied": analog['transformation_applied'],
            "therapeutic_area": f"{parent_compound}-based Therapeutics",
            "confidence_score": round(confidence, 1),
            "patent_potential": patent_potential,
            "ip_status": analog['patent_status'],
            "estimated_value": self._estimate_value(analog),
            "next_steps": self._determine_next_steps(analog),
            "discovery_date": session_timestamp,
            "discovery_method": "RDKit-based Computational Analog Generation",
            "molecular_properties": {
                "molecular_weight": analog['molecular_weight'],
                "logp": analog['logp'],
                "tpsa": analog['tpsa'],
                "h_bond_donors": analog['h_bond_donors'],
                "h_bond_acceptors": analog['h_bond_acceptors'],
                "rotatable_bonds": analog['rotatable_bonds'],
                "drug_likeness": analog['drug_likeness']
            },
            "patent_opportunity_score": analog['patent_opportunity_score'],
            "therapeutic_potential": analog['therapeutic_potential']
        }
        
        return finding
    
    def _estimate_value(self, analog: Dict[str, Any]) -> str:
        """Estimate commercial value based on properties"""
        if analog['drug_likeness'] >= 90 and analog['patent_opportunity_score'] >= 90:
            return "$15-25M"
        elif analog['drug_likeness'] >= 75 and analog['patent_opportunity_score'] >= 75:
            return "$8-15M"
        elif analog['drug_likeness'] >= 60:
            return "$3-8M"
        else:
            return "$1-3M"
    
    def _determine_next_steps(self, analog: Dict[str, Any]) -> str:
        """Determine next steps based on analog properties"""
        if analog['drug_likeness'] >= 85:
            return "Proceed to in vitro receptor binding assays and ADMET profiling"
        elif analog['drug_likeness'] >= 70:
            return "Conduct preliminary computational docking studies and toxicity prediction"
        else:
            return "Further structural optimization recommended before experimental validation"
    
    def integrate_analog_session(self, master_log_path: str, session_index: int = -1) -> int:
        """
        Integrate an analog discovery session into research findings
        
        Args:
            master_log_path: Path to MASTER_ANALOG_DISCOVERIES.json
            session_index: Which session to integrate (-1 for most recent)
            
        Returns:
            Number of findings added
        """
        # Load master analog discoveries
        with open(master_log_path, 'r') as f:
            master_data = json.load(f)
        
        if not master_data.get('discovery_sessions'):
            return 0
        
        # Get the specified session
        session = master_data['discovery_sessions'][session_index]
        parent_compound = session['parent_compound']
        session_timestamp = session['discovery_timestamp']
        analogs = session['analogs']
        
        # Load existing research findings
        research_data = self.load_research_findings()
        
        # Convert top analogs to research findings (top 5 based on IP opportunity)
        sorted_analogs = sorted(analogs, 
                               key=lambda x: x['patent_opportunity_score'], 
                               reverse=True)
        
        findings_added = 0
        for analog in sorted_analogs[:5]:  # Top 5 analogs
            finding = self.convert_analog_to_research_finding(
                analog, parent_compound, session_timestamp
            )
            research_data['findings'].append(finding)
            findings_added += 1
        
        # Update total count
        research_data['total_findings'] = len(research_data['findings'])
        
        # Save updated research findings
        self.save_research_findings(research_data)
        
        return findings_added
    
    def get_analog_discoveries_summary(self, master_log_path: str) -> Dict[str, Any]:
        """Get summary of all analog discoveries"""
        with open(master_log_path, 'r') as f:
            master_data = json.load(f)
        
        sessions = master_data.get('discovery_sessions', [])
        
        summary = {
            "total_sessions": len(sessions),
            "total_analogs_generated": sum(s['total_analogs'] for s in sessions),
            "parent_compounds": [s['parent_compound'] for s in sessions],
            "latest_discovery": sessions[-1]['discovery_timestamp'] if sessions else None,
            "high_potential_analogs": 0
        }
        
        # Count high potential analogs (IP opportunity >= 90)
        for session in sessions:
            for analog in session['analogs']:
                if analog['patent_opportunity_score'] >= 90:
                    summary['high_potential_analogs'] += 1
        
        return summary


if __name__ == "__main__":
    # Test integration
    integrator = AnalogResearchIntegration()
    
    master_log = "/home/ubuntu/pharmasight-latest/MASTER_ANALOG_DISCOVERIES.json"
    
    # Integrate ketamine session
    findings_added = integrator.integrate_analog_session(master_log)
    print(f"✅ Integrated {findings_added} analog discoveries into research findings")
    
    # Get summary
    summary = integrator.get_analog_discoveries_summary(master_log)
    print(f"\n📊 Discovery Summary:")
    print(f"   Total Sessions: {summary['total_sessions']}")
    print(f"   Total Analogs: {summary['total_analogs_generated']}")
    print(f"   High Potential Analogs: {summary['high_potential_analogs']}")
    print(f"   Parent Compounds: {', '.join(summary['parent_compounds'])}")

