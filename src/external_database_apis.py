"""
External Database API Integration Module
Connects to 6 major pharmaceutical databases for compound discovery and validation
"""

import requests
import time
import json
from typing import Dict, List, Optional
from functools import lru_cache
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ExternalDatabaseAPI:
    """Base class for external database API connections"""
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PharmaSight/3.0 (Research Platform; +https://pharmasight.com)'
        })
        self.rate_limit_delay = 0.5  # seconds between requests
        self.last_request_time = 0
        
    def _rate_limit(self):
        """Implement rate limiting"""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.rate_limit_delay:
            time.sleep(self.rate_limit_delay - elapsed)
        self.last_request_time = time.time()
    
    def _make_request(self, url: str, params: dict = None, timeout: int = 10) -> Optional[dict]:
        """Make rate-limited API request with error handling"""
        self._rate_limit()
        try:
            response = self.session.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json() if response.content else None
        except requests.exceptions.RequestException as e:
            logger.error(f"API request failed for {url}: {e}")
            return None


class PubChemAPI(ExternalDatabaseAPI):
    """PubChem API integration for compound data"""
    
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    @lru_cache(maxsize=1000)
    def search_compound(self, name: str) -> Optional[Dict]:
        """Search for compound by name"""
        url = f"{self.BASE_URL}/compound/name/{name}/JSON"
        data = self._make_request(url)
        
        if not data or 'PC_Compounds' not in data:
            return None
        
        compound = data['PC_Compounds'][0]
        return self._parse_pubchem_compound(compound)
    
    @lru_cache(maxsize=1000)
    def get_compound_by_smiles(self, smiles: str) -> Optional[Dict]:
        """Get compound data by SMILES"""
        url = f"{self.BASE_URL}/compound/smiles/{smiles}/JSON"
        data = self._make_request(url)
        
        if not data or 'PC_Compounds' not in data:
            return None
        
        return self._parse_pubchem_compound(data['PC_Compounds'][0])
    
    def search_similar_compounds(self, smiles: str, threshold: float = 0.9) -> List[Dict]:
        """Find structurally similar compounds"""
        url = f"{self.BASE_URL}/compound/fastsimilarity_2d/smiles/{smiles}/JSON"
        params = {'Threshold': int(threshold * 100)}
        
        data = self._make_request(url, params=params)
        if not data or 'IdentifierList' not in data:
            return []
        
        # Get details for similar compounds (limit to top 10)
        cids = data['IdentifierList']['CID'][:10]
        return [self.get_compound_by_cid(cid) for cid in cids if self.get_compound_by_cid(cid)]
    
    @lru_cache(maxsize=1000)
    def get_compound_by_cid(self, cid: int) -> Optional[Dict]:
        """Get compound by PubChem CID"""
        url = f"{self.BASE_URL}/compound/cid/{cid}/JSON"
        data = self._make_request(url)
        
        if not data or 'PC_Compounds' not in data:
            return None
        
        return self._parse_pubchem_compound(data['PC_Compounds'][0])
    
    def _parse_pubchem_compound(self, compound: dict) -> Dict:
        """Parse PubChem compound data"""
        props = compound.get('props', [])
        
        def get_prop(label: str):
            for prop in props:
                if prop.get('urn', {}).get('label') == label:
                    return prop.get('value', {}).get('sval') or prop.get('value', {}).get('fval')
            return None
        
        return {
            'source': 'PubChem',
            'cid': compound.get('id', {}).get('id', {}).get('cid'),
            'molecular_formula': get_prop('Molecular Formula'),
            'molecular_weight': get_prop('Molecular Weight'),
            'canonical_smiles': get_prop('SMILES') or get_prop('Canonical SMILES'),
            'iupac_name': get_prop('IUPAC Name'),
            'inchi': get_prop('InChI'),
            'inchi_key': get_prop('InChIKey'),
        }


class ChEMBLAPI(ExternalDatabaseAPI):
    """ChEMBL API integration for bioactivity data"""
    
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
    
    @lru_cache(maxsize=1000)
    def search_compound(self, name: str) -> Optional[Dict]:
        """Search for compound by name"""
        url = f"{self.BASE_URL}/molecule.json"
        params = {'molecule_synonyms__molecule_synonym__iexact': name}
        
        data = self._make_request(url, params=params)
        if not data or not data.get('molecules'):
            return None
        
        return self._parse_chembl_molecule(data['molecules'][0])
    
    @lru_cache(maxsize=1000)
    def get_bioactivity_data(self, chembl_id: str) -> List[Dict]:
        """Get bioactivity data for a compound"""
        url = f"{self.BASE_URL}/activity.json"
        params = {
            'molecule_chembl_id': chembl_id,
            'limit': 100
        }
        
        data = self._make_request(url, params=params)
        if not data or not data.get('activities'):
            return []
        
        return [self._parse_bioactivity(act) for act in data['activities']]
    
    @lru_cache(maxsize=1000)
    def search_by_smiles(self, smiles: str, similarity: int = 90) -> List[Dict]:
        """Search for similar compounds by SMILES"""
        url = f"{self.BASE_URL}/similarity/{smiles}/{similarity}.json"
        
        data = self._make_request(url)
        if not data or not data.get('molecules'):
            return []
        
        return [self._parse_chembl_molecule(mol) for mol in data['molecules'][:10]]
    
    def _parse_chembl_molecule(self, molecule: dict) -> Dict:
        """Parse ChEMBL molecule data"""
        return {
            'source': 'ChEMBL',
            'chembl_id': molecule.get('molecule_chembl_id'),
            'pref_name': molecule.get('pref_name'),
            'max_phase': molecule.get('max_phase'),
            'molecular_weight': molecule.get('molecule_properties', {}).get('full_mwt'),
            'alogp': molecule.get('molecule_properties', {}).get('alogp'),
            'canonical_smiles': molecule.get('molecule_structures', {}).get('canonical_smiles'),
            'therapeutic_flag': molecule.get('therapeutic_flag'),
        }
    
    def _parse_bioactivity(self, activity: dict) -> Dict:
        """Parse bioactivity data"""
        return {
            'assay_chembl_id': activity.get('assay_chembl_id'),
            'target_chembl_id': activity.get('target_chembl_id'),
            'target_pref_name': activity.get('target_pref_name'),
            'standard_type': activity.get('standard_type'),
            'standard_value': activity.get('standard_value'),
            'standard_units': activity.get('standard_units'),
            'activity_comment': activity.get('activity_comment'),
        }


class DrugBankAPI(ExternalDatabaseAPI):
    """DrugBank API integration (requires API key for full access)"""
    
    # Note: DrugBank requires paid API access for programmatic queries
    # This implementation uses the open data dumps when available
    
    def __init__(self, api_key: Optional[str] = None):
        super().__init__()
        self.api_key = api_key
        if api_key:
            self.session.headers.update({'Authorization': f'Bearer {api_key}'})
    
    def search_compound(self, name: str) -> Optional[Dict]:
        """Search DrugBank (limited without API key)"""
        # Without API key, return placeholder
        # With API key, implement full search
        if not self.api_key:
            logger.warning("DrugBank API key not provided. Using limited data.")
            return None
        
        # Implement API search when key is available
        return None


class ZINCAPI(ExternalDatabaseAPI):
    """ZINC database API for purchasable compounds"""
    
    BASE_URL = "https://zinc.docking.org/api"
    
    def search_compound(self, smiles: str) -> Optional[Dict]:
        """Search ZINC by SMILES"""
        # ZINC API implementation
        # Note: ZINC has specific query formats
        logger.info("ZINC API search not yet fully implemented")
        return None
    
    def search_analogs(self, smiles: str, similarity: float = 0.8) -> List[Dict]:
        """Find purchasable analogs"""
        logger.info("ZINC analog search not yet fully implemented")
        return []


class OpenTargetsAPI(ExternalDatabaseAPI):
    """Open Targets Platform API for target-disease associations"""
    
    BASE_URL = "https://api.platform.opentargets.org/api/v4/graphql"
    
    def search_target(self, gene_symbol: str) -> Optional[Dict]:
        """Search for target information"""
        query = """
        query target($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            id
            approvedSymbol
            approvedName
            biotype
            functionDescriptions
          }
        }
        """
        
        # GraphQL query implementation
        logger.info("OpenTargets API search not yet fully implemented")
        return None
    
    def get_disease_associations(self, target_id: str) -> List[Dict]:
        """Get disease associations for a target"""
        logger.info("OpenTargets disease associations not yet fully implemented")
        return []


class FDAOrangeBookAPI(ExternalDatabaseAPI):
    """FDA Orange Book API for patent and exclusivity data"""
    
    BASE_URL = "https://api.fda.gov/drug"
    
    @lru_cache(maxsize=1000)
    def search_drug(self, drug_name: str) -> Optional[Dict]:
        """Search FDA Orange Book"""
        url = f"{self.BASE_URL}/label.json"
        params = {
            'search': f'openfda.brand_name:"{drug_name}"',
            'limit': 1
        }
        
        data = self._make_request(url, params=params)
        if not data or not data.get('results'):
            return None
        
        return self._parse_fda_drug(data['results'][0])
    
    def _parse_fda_drug(self, drug: dict) -> Dict:
        """Parse FDA drug data"""
        openfda = drug.get('openfda', {})
        return {
            'source': 'FDA',
            'brand_name': openfda.get('brand_name', [None])[0],
            'generic_name': openfda.get('generic_name', [None])[0],
            'manufacturer_name': openfda.get('manufacturer_name', [None])[0],
            'product_type': openfda.get('product_type', [None])[0],
            'route': openfda.get('route', []),
            'substance_name': openfda.get('substance_name', []),
        }


class UnifiedDatabaseSearch:
    """Unified interface to search across all databases"""
    
    def __init__(self, drugbank_api_key: Optional[str] = None):
        self.pubchem = PubChemAPI()
        self.chembl = ChEMBLAPI()
        self.drugbank = DrugBankAPI(api_key=drugbank_api_key)
        self.zinc = ZINCAPI()
        self.opentargets = OpenTargetsAPI()
        self.fda = FDAOrangeBookAPI()
    
    def search_compound_all_databases(self, name: str) -> Dict[str, Optional[Dict]]:
        """Search for compound across all databases"""
        results = {
            'pubchem': None,
            'chembl': None,
            'drugbank': None,
            'fda': None,
        }
        
        # Search PubChem
        try:
            results['pubchem'] = self.pubchem.search_compound(name)
        except Exception as e:
            logger.error(f"PubChem search failed: {e}")
        
        # Search ChEMBL
        try:
            results['chembl'] = self.chembl.search_compound(name)
        except Exception as e:
            logger.error(f"ChEMBL search failed: {e}")
        
        # Search DrugBank
        try:
            results['drugbank'] = self.drugbank.search_compound(name)
        except Exception as e:
            logger.error(f"DrugBank search failed: {e}")
        
        # Search FDA
        try:
            results['fda'] = self.fda.search_drug(name)
        except Exception as e:
            logger.error(f"FDA search failed: {e}")
        
        return results
    
    def find_novel_analogs(self, smiles: str, similarity_threshold: float = 0.85) -> List[Dict]:
        """
        Find novel analogs not in major databases
        Returns compounds with high similarity but not exact matches
        """
        all_analogs = []
        
        # Search PubChem for similar compounds
        try:
            pubchem_similar = self.pubchem.search_similar_compounds(smiles, similarity_threshold)
            all_analogs.extend(pubchem_similar)
        except Exception as e:
            logger.error(f"PubChem similarity search failed: {e}")
        
        # Search ChEMBL for similar compounds
        try:
            chembl_similar = self.chembl.search_by_smiles(smiles, int(similarity_threshold * 100))
            all_analogs.extend(chembl_similar)
        except Exception as e:
            logger.error(f"ChEMBL similarity search failed: {e}")
        
        # Filter for truly novel compounds (not exact matches)
        novel_analogs = []
        for analog in all_analogs:
            analog_smiles = analog.get('canonical_smiles', '')
            if analog_smiles and analog_smiles != smiles:
                # Check if compound has limited data (indicator of novelty)
                if self._is_novel_compound(analog):
                    novel_analogs.append(analog)
        
        return novel_analogs
    
    def _is_novel_compound(self, compound: Dict) -> bool:
        """Determine if compound is novel based on available data"""
        # Novel compounds typically have:
        # - Limited bioactivity data
        # - No FDA approval
        # - Limited literature references
        
        source = compound.get('source')
        
        if source == 'PubChem':
            # Check if it has limited data
            return compound.get('iupac_name') is None or compound.get('molecular_formula') is None
        
        elif source == 'ChEMBL':
            # Check development phase
            max_phase = compound.get('max_phase', 0)
            return max_phase == 0 or max_phase is None
        
        return True
    
    def get_comprehensive_compound_data(self, name: str, smiles: Optional[str] = None) -> Dict:
        """Get comprehensive data from all sources"""
        all_data = self.search_compound_all_databases(name)
        
        # If SMILES provided, get bioactivity data
        if smiles:
            try:
                # Get ChEMBL ID from search results
                chembl_data = all_data.get('chembl')
                if chembl_data and chembl_data.get('chembl_id'):
                    bioactivity = self.chembl.get_bioactivity_data(chembl_data['chembl_id'])
                    all_data['bioactivity'] = bioactivity
            except Exception as e:
                logger.error(f"Bioactivity data retrieval failed: {e}")
        
        return all_data


# Singleton instance
_unified_search = None

def get_unified_search(drugbank_api_key: Optional[str] = None) -> UnifiedDatabaseSearch:
    """Get or create unified search instance"""
    global _unified_search
    if _unified_search is None:
        _unified_search = UnifiedDatabaseSearch(drugbank_api_key=drugbank_api_key)
    return _unified_search

