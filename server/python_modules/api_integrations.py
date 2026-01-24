"""
API Integration Module for PharmaSight™ Platform
Provides integration with major pharmaceutical and chemical repositories.
"""

import requests
import json
import time
from typing import Dict, List, Optional, Any
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class APIIntegrationError(Exception):
    """Custom exception for API integration errors."""
    pass

class RateLimiter:
    """Rate limiter to ensure we don't exceed API limits."""
    
    def __init__(self, max_requests_per_second: float = 5.0):
        self.max_requests_per_second = max_requests_per_second
        self.min_interval = 1.0 / max_requests_per_second
        self.last_request_time = 0
    
    def wait_if_needed(self):
        """Wait if necessary to respect rate limits."""
        current_time = time.time()
        time_since_last_request = current_time - self.last_request_time
        
        if time_since_last_request < self.min_interval:
            sleep_time = self.min_interval - time_since_last_request
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()

class PubChemAPI:
    """PubChem PUG-REST API integration."""
    
    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.rate_limiter = RateLimiter(max_requests_per_second=5.0)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PharmaSight-Platform/1.0 (https://github.com/justincihi/pharmasight-platform)'
        })
    
    def get_compound_by_name(self, compound_name: str) -> Optional[Dict[str, Any]]:
        """Get compound information by name from PubChem."""
        try:
            self.rate_limiter.wait_if_needed()
            
            # First, get the CID for the compound name
            cid_url = f"{self.base_url}/compound/name/{compound_name}/cids/JSON"
            response = self.session.get(cid_url, timeout=30)
            
            if response.status_code != 200:
                logger.warning(f"Failed to get CID for {compound_name}: {response.status_code}")
                return None
            
            cid_data = response.json()
            if 'IdentifierList' not in cid_data or not cid_data['IdentifierList']['CID']:
                logger.warning(f"No CID found for {compound_name}")
                return None
            
            cid = cid_data['IdentifierList']['CID'][0]
            
            # Get detailed compound information
            return self.get_compound_by_cid(cid)
            
        except Exception as e:
            logger.error(f"Error getting compound by name {compound_name}: {str(e)}")
            return None
    
    def get_compound_by_cid(self, cid: int) -> Optional[Dict[str, Any]]:
        """Get compound information by CID from PubChem."""
        try:
            self.rate_limiter.wait_if_needed()
            
            # Get compound properties
            properties_url = f"{self.base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,InChI,InChIKey,CanonicalSMILES,IUPACName/JSON"
            response = self.session.get(properties_url, timeout=30)
            
            if response.status_code != 200:
                logger.warning(f"Failed to get properties for CID {cid}: {response.status_code}")
                return None
            
            properties_data = response.json()
            
            # Get synonyms
            synonyms = self.get_compound_synonyms(cid)
            
            # Format the response
            compound_data = {
                'source': 'PubChem',
                'cid': cid,
                'properties': properties_data['PropertyTable']['Properties'][0] if properties_data.get('PropertyTable', {}).get('Properties') else {},
                'synonyms': synonyms,
                'url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
            }
            
            return compound_data
            
        except Exception as e:
            logger.error(f"Error getting compound by CID {cid}: {str(e)}")
            return None
    
    def get_compound_synonyms(self, cid: int) -> List[str]:
        """Get synonyms for a compound by CID."""
        try:
            self.rate_limiter.wait_if_needed()
            
            synonyms_url = f"{self.base_url}/compound/cid/{cid}/synonyms/JSON"
            response = self.session.get(synonyms_url, timeout=30)
            
            if response.status_code != 200:
                return []
            
            synonyms_data = response.json()
            if 'InformationList' in synonyms_data and synonyms_data['InformationList']['Information']:
                return synonyms_data['InformationList']['Information'][0].get('Synonym', [])[:10]  # Limit to first 10 synonyms
            
            return []
            
        except Exception as e:
            logger.error(f"Error getting synonyms for CID {cid}: {str(e)}")
            return []

class FDAOpenAPI:
    """FDA openFDA API integration."""
    
    def __init__(self):
        self.base_url = "https://api.fda.gov/drug/drugsfda.json"
        self.rate_limiter = RateLimiter(max_requests_per_second=5.0)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PharmaSight-Platform/1.0 (https://github.com/justincihi/pharmasight-platform)'
        })
    
    def search_drugs_by_name(self, drug_name: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Search for drugs by name in the FDA Orange Book."""
        try:
            self.rate_limiter.wait_if_needed()
            
            # Search for drugs by brand name or active ingredient
            search_query = f'products.brand_name:"{drug_name}" OR products.active_ingredients.name:"{drug_name.upper()}"'
            params = {
                'search': search_query,
                'limit': limit
            }
            
            response = self.session.get(self.base_url, params=params, timeout=30)
            
            if response.status_code != 200:
                logger.warning(f"Failed to search FDA drugs for {drug_name}: {response.status_code}")
                return []
            
            fda_data = response.json()
            
            if 'results' not in fda_data:
                return []
            
            formatted_results = []
            for result in fda_data['results']:
                formatted_drug = self._format_fda_drug_data(result)
                if formatted_drug:
                    formatted_results.append(formatted_drug)
            
            return formatted_results
            
        except Exception as e:
            logger.error(f"Error searching FDA drugs for {drug_name}: {str(e)}")
            return []
    
    def _format_fda_drug_data(self, raw_data: Dict[str, Any]) -> Dict[str, Any]:
        """Format raw FDA API data into a standardized format."""
        try:
            products = raw_data.get('products', [])
            
            formatted_data = {
                'source': 'FDA Orange Book',
                'application_number': raw_data.get('application_number'),
                'sponsor_name': raw_data.get('sponsor_name'),
                'products': [],
                'url': f"https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&ApplNo={raw_data.get('application_number', '')}"
            }
            
            for product in products:
                formatted_product = {
                    'brand_name': product.get('brand_name', []),
                    'dosage_form': product.get('dosage_form'),
                    'route': product.get('route'),
                    'marketing_status': product.get('marketing_status'),
                    'active_ingredients': product.get('active_ingredients', []),
                    'te_code': product.get('te_code'),  # Therapeutic Equivalence code
                    'reference_drug': product.get('reference_drug'),
                    'reference_standard': product.get('reference_standard')
                }
                formatted_data['products'].append(formatted_product)
            
            return formatted_data
            
        except Exception as e:
            logger.error(f"Error formatting FDA drug data: {str(e)}")
            return None

class APIIntegrationManager:
    """Manager class for coordinating multiple API integrations."""
    
    def __init__(self):
        self.pubchem = PubChemAPI()
        self.fda = FDAOpenAPI()
        self.chembl = ChEMBLAPI()
        self.drugbank = DrugBankAPI()
        self.zinc = ZINCDatabase()
        self.opentargets = OpenTargetsAPI()
    
    def search_compound_comprehensive(self, compound_name: str) -> Dict[str, Any]:
        """Search for a compound across multiple databases."""
        results = {
            'query': compound_name,
            'sources': {},
            'summary': {
                'total_sources': 0,
                'successful_sources': 0,
                'failed_sources': []
            }
        }
        
        # Search PubChem
        try:
            pubchem_data = self.pubchem.get_compound_by_name(compound_name)
            if pubchem_data:
                results['sources']['pubchem'] = pubchem_data
                results['summary']['successful_sources'] += 1
            else:
                results['summary']['failed_sources'].append('PubChem')
        except Exception as e:
            logger.error(f"PubChem search failed: {str(e)}")
            results['summary']['failed_sources'].append('PubChem')
        
        results['summary']['total_sources'] += 1
        
        # Search FDA Orange Book
        try:
            fda_data = self.fda.search_drugs_by_name(compound_name)
            if fda_data:
                results['sources']['fda'] = fda_data
                results['summary']['successful_sources'] += 1
            else:
                results['summary']['failed_sources'].append('FDA')
        except Exception as e:
            logger.error(f"FDA search failed: {str(e)}")
            results['summary']['failed_sources'].append('FDA')
        
        results['summary']['total_sources'] += 1
        
        # Search ChEMBL
        try:
            chembl_data = self.chembl.search_compound_by_name(compound_name)
            if chembl_data:
                results['sources']['chembl'] = chembl_data
                results['summary']['successful_sources'] += 1
            else:
                results['summary']['failed_sources'].append('ChEMBL')
        except Exception as e:
            logger.error(f"ChEMBL search failed: {str(e)}")
            results['summary']['failed_sources'].append('ChEMBL')
        
        results['summary']['total_sources'] += 1
        
        # Search DrugBank
        try:
            drugbank_data = self.drugbank.search_drug_by_name(compound_name)
            if drugbank_data:
                results['sources']['drugbank'] = drugbank_data
                results['summary']['successful_sources'] += 1
            else:
                results['summary']['failed_sources'].append('DrugBank')
        except Exception as e:
            logger.error(f"DrugBank search failed: {str(e)}")
            results['summary']['failed_sources'].append('DrugBank')
        
        results['summary']['total_sources'] += 1
        
        # Search ZINC Database
        try:
            zinc_data = self.zinc.search_compound_by_name(compound_name)
            if zinc_data:
                results['sources']['zinc'] = zinc_data
                results['summary']['successful_sources'] += 1
            else:
                results['summary']['failed_sources'].append('ZINC')
        except Exception as e:
            logger.error(f"ZINC search failed: {str(e)}")
            results['summary']['failed_sources'].append('ZINC')
        
        results['summary']['total_sources'] += 1
        
        # Search OpenTargets
        try:
            opentargets_data = self.opentargets.search_drug_by_name(compound_name)
            if opentargets_data:
                results['sources']['opentargets'] = opentargets_data
                results['summary']['successful_sources'] += 1
            else:
                results['summary']['failed_sources'].append('OpenTargets')
        except Exception as e:
            logger.error(f"OpenTargets search failed: {str(e)}")
            results['summary']['failed_sources'].append('OpenTargets')
        
        results['summary']['total_sources'] += 1
        
        return results
    
    def get_compound_properties(self, compound_name: str) -> Dict[str, Any]:
        """Get comprehensive compound properties from multiple sources."""
        comprehensive_data = self.search_compound_comprehensive(compound_name)
        
        # Extract and merge properties from different sources
        merged_properties = {
            'name': compound_name,
            'molecular_formula': None,
            'molecular_weight': None,
            'inchi': None,
            'inchi_key': None,
            'smiles': None,
            'iupac_name': None,
            'synonyms': [],
            'fda_approved': False,
            'marketing_status': [],
            'therapeutic_equivalence': [],
            'bioactivity_data': [],
            'target_information': [],
            'chembl_id': None,
            'sources': list(comprehensive_data['sources'].keys())
        }
        
        # Extract PubChem data
        if 'pubchem' in comprehensive_data['sources']:
            pubchem_data = comprehensive_data['sources']['pubchem']
            properties = pubchem_data.get('properties', {})
            
            merged_properties.update({
                'molecular_formula': properties.get('MolecularFormula'),
                'molecular_weight': properties.get('MolecularWeight'),
                'inchi': properties.get('InChI'),
                'inchi_key': properties.get('InChIKey'),
                'smiles': properties.get('CanonicalSMILES'),
                'iupac_name': properties.get('IUPACName'),
                'synonyms': pubchem_data.get('synonyms', [])
            })
        
        # Extract FDA data
        if 'fda' in comprehensive_data['sources']:
            fda_data = comprehensive_data['sources']['fda']
            if fda_data:
                merged_properties['fda_approved'] = True
                for drug in fda_data:
                    for product in drug.get('products', []):
                        if product.get('marketing_status'):
                            merged_properties['marketing_status'].append(product['marketing_status'])
                        if product.get('te_code'):
                            merged_properties['therapeutic_equivalence'].append(product['te_code'])
        
        # Extract ChEMBL data
        if 'chembl' in comprehensive_data['sources']:
            chembl_data = comprehensive_data['sources']['chembl']
            if chembl_data:
                merged_properties['chembl_id'] = chembl_data.get('chembl_id')
                merged_properties['bioactivity_data'] = chembl_data.get('bioactivities', [])[:10]  # Limit to top 10
                merged_properties['target_information'] = chembl_data.get('targets', [])[:5]  # Limit to top 5
        
        # Extract DrugBank data
        if 'drugbank' in comprehensive_data['sources']:
            drugbank_data = comprehensive_data['sources']['drugbank']
            if drugbank_data:
                merged_properties['drugbank_status'] = drugbank_data.get('status')
        
        # Extract ZINC data
        if 'zinc' in comprehensive_data['sources']:
            zinc_data = comprehensive_data['sources']['zinc']
            if zinc_data:
                merged_properties['commercial_availability'] = zinc_data.get('commercial_availability', False)
                merged_properties['zinc_results_count'] = zinc_data.get('total_results', 0)
        
        # Extract OpenTargets data
        if 'opentargets' in comprehensive_data['sources']:
            opentargets_data = comprehensive_data['sources']['opentargets']
            if opentargets_data:
                merged_properties['opentargets_hits'] = len(opentargets_data.get('hits', []))
        
        return merged_properties

# This will be moved to the end of the file


class ChEMBLAPI:
    """ChEMBL API integration for bioactivity data."""
    
    def __init__(self):
        self.base_url = "https://www.ebi.ac.uk/chembl/api/data"
        self.rate_limiter = RateLimiter(max_requests_per_second=10.0)  # ChEMBL has no strict rate limits
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PharmaSight-Platform/1.0 (https://github.com/justincihi/pharmasight-platform)',
            'Accept': 'application/json'
        })
    
    def search_compound_by_name(self, compound_name: str) -> Optional[Dict[str, Any]]:
        """Search for compound in ChEMBL by name."""
        try:
            self.rate_limiter.wait_if_needed()
            
            # Search for molecules by name
            search_url = f"{self.base_url}/molecule/search.json"
            params = {
                'q': compound_name,
                'format': 'json'
            }
            
            response = self.session.get(search_url, params=params, timeout=30)
            
            if response.status_code != 200:
                logger.warning(f"Failed to search ChEMBL for {compound_name}: {response.status_code}")
                return None
            
            search_data = response.json()
            
            if not search_data.get('molecules') or len(search_data['molecules']) == 0:
                return None
            
            # Get the first match
            molecule = search_data['molecules'][0]
            chembl_id = molecule.get('molecule_chembl_id')
            
            if chembl_id:
                return self.get_compound_by_chembl_id(chembl_id)
            
            return None
            
        except Exception as e:
            logger.error(f"Error searching ChEMBL for {compound_name}: {str(e)}")
            return None
    
    def get_compound_by_chembl_id(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        """Get detailed compound information by ChEMBL ID."""
        try:
            self.rate_limiter.wait_if_needed()
            
            # Get molecule details
            molecule_url = f"{self.base_url}/molecule/{chembl_id}.json"
            response = self.session.get(molecule_url, timeout=30)
            
            if response.status_code != 200:
                return None
            
            molecule_data = response.json()
            
            # Get bioactivity data
            bioactivities = self.get_bioactivities(chembl_id)
            
            # Get target information
            targets = self.get_targets_for_compound(chembl_id)
            
            formatted_data = {
                'source': 'ChEMBL',
                'chembl_id': chembl_id,
                'molecule_data': molecule_data,
                'bioactivities': bioactivities,
                'targets': targets,
                'url': f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/"
            }
            
            return formatted_data
            
        except Exception as e:
            logger.error(f"Error getting ChEMBL compound {chembl_id}: {str(e)}")
            return None
    
    def get_bioactivities(self, chembl_id: str, limit: int = 20) -> List[Dict[str, Any]]:
        """Get bioactivity data for a compound."""
        try:
            self.rate_limiter.wait_if_needed()
            
            bioactivity_url = f"{self.base_url}/activity.json"
            params = {
                'molecule_chembl_id': chembl_id,
                'limit': limit,
                'format': 'json'
            }
            
            response = self.session.get(bioactivity_url, params=params, timeout=30)
            
            if response.status_code != 200:
                return []
            
            bioactivity_data = response.json()
            return bioactivity_data.get('activities', [])
            
        except Exception as e:
            logger.error(f"Error getting bioactivities for {chembl_id}: {str(e)}")
            return []
    
    def get_targets_for_compound(self, chembl_id: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Get target information for a compound."""
        try:
            self.rate_limiter.wait_if_needed()
            
            # Get targets through bioactivities
            bioactivities = self.get_bioactivities(chembl_id, limit=50)
            
            target_ids = set()
            for activity in bioactivities:
                if activity.get('target_chembl_id'):
                    target_ids.add(activity['target_chembl_id'])
            
            targets = []
            for target_id in list(target_ids)[:limit]:
                target_data = self.get_target_details(target_id)
                if target_data:
                    targets.append(target_data)
            
            return targets
            
        except Exception as e:
            logger.error(f"Error getting targets for {chembl_id}: {str(e)}")
            return []
    
    def get_target_details(self, target_chembl_id: str) -> Optional[Dict[str, Any]]:
        """Get detailed target information."""
        try:
            self.rate_limiter.wait_if_needed()
            
            target_url = f"{self.base_url}/target/{target_chembl_id}.json"
            response = self.session.get(target_url, timeout=30)
            
            if response.status_code != 200:
                return None
            
            return response.json()
            
        except Exception as e:
            logger.error(f"Error getting target details for {target_chembl_id}: {str(e)}")
            return None

class DrugBankAPI:
    """DrugBank API integration for comprehensive drug information."""
    
    def __init__(self):
        # Note: DrugBank requires API key for full access
        # For now, we'll use the public search functionality
        self.base_url = "https://go.drugbank.com"
        self.rate_limiter = RateLimiter(max_requests_per_second=2.0)  # Conservative rate limiting
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PharmaSight-Platform/1.0 (https://github.com/justincihi/pharmasight-platform)'
        })
    
    def search_drug_by_name(self, drug_name: str) -> Optional[Dict[str, Any]]:
        """Search for drug information by name using DrugBank's public interface."""
        try:
            self.rate_limiter.wait_if_needed()
            
            # Note: This is a simplified implementation
            # Full DrugBank API access requires authentication
            search_data = {
                'source': 'DrugBank',
                'drug_name': drug_name,
                'status': 'Limited access - API key required for full data',
                'url': f"https://go.drugbank.com/drugs?q={drug_name.replace(' ', '+')}"
            }
            
            # For demonstration, we'll return basic structure
            # In production, this would use the authenticated API
            return search_data
            
        except Exception as e:
            logger.error(f"Error searching DrugBank for {drug_name}: {str(e)}")
            return None

class ZINCDatabase:
    """ZINC Database integration for commercially available compounds."""
    
    def __init__(self):
        self.base_url = "https://zinc.docking.org"
        self.api_base = "https://zinc.docking.org/api"
        self.rate_limiter = RateLimiter(max_requests_per_second=3.0)  # Conservative rate limiting
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PharmaSight-Platform/1.0 (https://github.com/justincihi/pharmasight-platform)',
            'Accept': 'application/json'
        })
    
    def search_compound_by_name(self, compound_name: str) -> Optional[Dict[str, Any]]:
        """Search for compounds in ZINC database."""
        try:
            self.rate_limiter.wait_if_needed()
            
            # ZINC has a REST API for substance search
            search_url = f"{self.api_base}/substances/search"
            params = {
                'q': compound_name,
                'format': 'json',
                'limit': 10
            }
            
            response = self.session.get(search_url, params=params, timeout=30)
            
            if response.status_code != 200:
                logger.warning(f"ZINC search failed for {compound_name}: {response.status_code}")
                # Return basic info even if API fails
                return {
                    'source': 'ZINC',
                    'compound_name': compound_name,
                    'status': 'API request failed - using fallback',
                    'url': f"https://zinc.docking.org/substances/search/?q={compound_name.replace(' ', '+')}"
                }
            
            try:
                search_data = response.json()
                
                formatted_results = {
                    'source': 'ZINC',
                    'compound_name': compound_name,
                    'results': search_data.get('results', [])[:5],  # Limit to top 5 results
                    'total_results': search_data.get('total', 0),
                    'commercial_availability': True if search_data.get('results') else False,
                    'url': f"https://zinc.docking.org/substances/search/?q={compound_name.replace(' ', '+')}"
                }
                
                return formatted_results
                
            except json.JSONDecodeError:
                # Handle non-JSON response
                return {
                    'source': 'ZINC',
                    'compound_name': compound_name,
                    'status': 'Response parsing failed',
                    'url': f"https://zinc.docking.org/substances/search/?q={compound_name.replace(' ', '+')}"
                }
            
        except Exception as e:
            logger.error(f"Error searching ZINC for {compound_name}: {str(e)}")
            # Return basic info even on error
            return {
                'source': 'ZINC',
                'compound_name': compound_name,
                'status': f'Search error: {str(e)[:100]}',
                'url': f"https://zinc.docking.org/substances/search/?q={compound_name.replace(' ', '+')}"
            }
    
    def get_compound_by_zinc_id(self, zinc_id: str) -> Optional[Dict[str, Any]]:
        """Get detailed compound information by ZINC ID."""
        try:
            self.rate_limiter.wait_if_needed()
            
            compound_url = f"{self.api_base}/substances/{zinc_id}"
            response = self.session.get(compound_url, timeout=30)
            
            if response.status_code != 200:
                return None
            
            compound_data = response.json()
            
            return {
                'source': 'ZINC',
                'zinc_id': zinc_id,
                'compound_data': compound_data,
                'url': f"https://zinc.docking.org/substances/{zinc_id}/"
            }
            
        except Exception as e:
            logger.error(f"Error getting ZINC compound {zinc_id}: {str(e)}")
            return None

class OpenTargetsAPI:
    """OpenTargets platform integration for target-disease associations."""
    
    def __init__(self):
        self.base_url = "https://api.platform.opentargets.org/api/v4/graphql"
        self.rate_limiter = RateLimiter(max_requests_per_second=10.0)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PharmaSight-Platform/1.0 (https://github.com/justincihi/pharmasight-platform)',
            'Content-Type': 'application/json'
        })
    
    def search_drug_by_name(self, drug_name: str) -> Optional[Dict[str, Any]]:
        """Search for drug information in OpenTargets."""
        try:
            self.rate_limiter.wait_if_needed()
            
            # GraphQL query for drug search
            query = """
            query searchDrug($queryString: String!) {
              search(queryString: $queryString, entityNames: ["drug"]) {
                hits {
                  id
                  name
                  description
                  entity
                }
              }
            }
            """
            
            variables = {"queryString": drug_name}
            
            response = self.session.post(
                self.base_url,
                json={"query": query, "variables": variables},
                timeout=30
            )
            
            if response.status_code != 200:
                logger.warning(f"Failed to search OpenTargets for {drug_name}: {response.status_code}")
                return None
            
            data = response.json()
            
            if data.get('data', {}).get('search', {}).get('hits'):
                hits = data['data']['search']['hits']
                if hits:
                    return {
                        'source': 'OpenTargets',
                        'drug_name': drug_name,
                        'hits': hits,
                        'url': f"https://platform.opentargets.org/search?q={drug_name.replace(' ', '+')}"
                    }
            
            return None
            
        except Exception as e:
            logger.error(f"Error searching OpenTargets for {drug_name}: {str(e)}")
            return None


# Global instance for use in the main application
api_manager = APIIntegrationManager()

# Test function
if __name__ == "__main__":
    print("Testing Enhanced API Integration...")
    result = api_manager.get_compound_properties("aspirin")
    print("Test result sources:", result.get('sources', []))
    print("ChEMBL ID:", result.get('chembl_id'))
    print("Bioactivity data count:", len(result.get('bioactivity_data', [])))
    print("Target information count:", len(result.get('target_information', [])))
