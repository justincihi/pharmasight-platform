"""
Unified Search Engine for PharmaSight™ Platform
Provides intelligent search and data aggregation across all integrated pharmaceutical databases.
"""

import logging
import json
from typing import Dict, List, Any, Optional, Set, Tuple
from dataclasses import dataclass, field
from datetime import datetime, timedelta
import hashlib
import pickle
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

from api_integrations import APIIntegrationManager

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class SearchResult:
    """Standardized search result structure."""
    compound_name: str
    source: str
    confidence_score: float
    data: Dict[str, Any]
    timestamp: datetime = field(default_factory=datetime.now)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'compound_name': self.compound_name,
            'source': self.source,
            'confidence_score': self.confidence_score,
            'data': self.data,
            'timestamp': self.timestamp.isoformat()
        }

@dataclass
class AggregatedResult:
    """Aggregated search result from multiple sources."""
    compound_name: str
    primary_identifier: Optional[str]
    molecular_formula: Optional[str]
    molecular_weight: Optional[float]
    sources: List[str]
    confidence_score: float
    comprehensive_data: Dict[str, Any]
    search_metadata: Dict[str, Any]
    timestamp: datetime = field(default_factory=datetime.now)

class IntelligentCache:
    """Intelligent caching system with TTL and smart invalidation."""
    
    def __init__(self, cache_dir: str = "/tmp/pharmasight_cache", ttl_hours: int = 24):
        self.cache_dir = cache_dir
        self.ttl = timedelta(hours=ttl_hours)
        os.makedirs(cache_dir, exist_ok=True)
    
    def _get_cache_key(self, compound_name: str) -> str:
        """Generate cache key for compound."""
        return hashlib.md5(compound_name.lower().encode()).hexdigest()
    
    def _get_cache_path(self, cache_key: str) -> str:
        """Get cache file path."""
        return os.path.join(self.cache_dir, f"{cache_key}.pkl")
    
    def get(self, compound_name: str) -> Optional[AggregatedResult]:
        """Get cached result if valid."""
        try:
            cache_key = self._get_cache_key(compound_name)
            cache_path = self._get_cache_path(cache_key)
            
            if not os.path.exists(cache_path):
                return None
            
            with open(cache_path, 'rb') as f:
                cached_result = pickle.load(f)
            
            # Check TTL
            if datetime.now() - cached_result.timestamp > self.ttl:
                os.remove(cache_path)
                return None
            
            return cached_result
            
        except Exception as e:
            logger.error(f"Cache get error: {str(e)}")
            return None
    
    def set(self, compound_name: str, result: AggregatedResult) -> None:
        """Cache result."""
        try:
            cache_key = self._get_cache_key(compound_name)
            cache_path = self._get_cache_path(cache_key)
            
            with open(cache_path, 'wb') as f:
                pickle.dump(result, f)
                
        except Exception as e:
            logger.error(f"Cache set error: {str(e)}")

class DataAggregator:
    """Intelligent data aggregation and fusion system."""
    
    @staticmethod
    def calculate_confidence_score(sources: List[str], data_quality: Dict[str, float]) -> float:
        """Calculate overall confidence score based on sources and data quality."""
        source_weights = {
            'pubchem': 0.25,
            'chembl': 0.25,
            'fda': 0.20,
            'drugbank': 0.15,
            'opentargets': 0.10,
            'zinc': 0.05
        }
        
        weighted_score = 0.0
        total_weight = 0.0
        
        for source in sources:
            if source in source_weights:
                weight = source_weights[source]
                quality = data_quality.get(source, 0.5)
                weighted_score += weight * quality
                total_weight += weight
        
        return min(weighted_score / total_weight if total_weight > 0 else 0.0, 1.0)
    
    @staticmethod
    def merge_molecular_data(source_data: Dict[str, Any]) -> Dict[str, Any]:
        """Merge molecular data from multiple sources with conflict resolution."""
        merged = {
            'molecular_formula': None,
            'molecular_weight': None,
            'smiles': None,
            'inchi': None,
            'inchi_key': None,
            'iupac_name': None,
            'synonyms': set(),
            'identifiers': {}
        }
        
        # Priority order for molecular data
        priority_sources = ['pubchem', 'chembl', 'drugbank', 'fda']
        
        for source in priority_sources:
            if source not in source_data:
                continue
                
            data = source_data[source]
            
            # Handle different data structures (dict vs list)
            if isinstance(data, list):
                # For sources that return lists (like FDA), skip molecular data extraction
                continue
            elif not isinstance(data, dict):
                continue
            
            # Extract molecular formula
            if not merged['molecular_formula']:
                formula = data.get('molecular_formula') or data.get('properties', {}).get('MolecularFormula')
                if formula:
                    merged['molecular_formula'] = formula
            
            # Extract molecular weight
            if not merged['molecular_weight']:
                weight = (data.get('molecular_weight') or 
                         data.get('properties', {}).get('MolecularWeight'))
                if weight:
                    try:
                        merged['molecular_weight'] = float(weight)
                    except (ValueError, TypeError):
                        pass
            
            # Extract SMILES
            if not merged['smiles']:
                smiles = (data.get('smiles') or 
                         data.get('properties', {}).get('CanonicalSMILES'))
                if smiles:
                    merged['smiles'] = smiles
            
            # Extract InChI
            if not merged['inchi']:
                inchi = (data.get('inchi') or 
                        data.get('properties', {}).get('InChI'))
                if inchi:
                    merged['inchi'] = inchi
            
            # Extract InChI Key
            if not merged['inchi_key']:
                inchi_key = (data.get('inchi_key') or 
                           data.get('properties', {}).get('InChIKey'))
                if inchi_key:
                    merged['inchi_key'] = inchi_key
            
            # Extract IUPAC name
            if not merged['iupac_name']:
                iupac = (data.get('iupac_name') or 
                        data.get('properties', {}).get('IUPACName'))
                if iupac:
                    merged['iupac_name'] = iupac
            
            # Collect synonyms
            synonyms = data.get('synonyms')
            if synonyms:
                if isinstance(synonyms, list):
                    merged['synonyms'].update(synonyms)
                elif isinstance(synonyms, str):
                    merged['synonyms'].add(synonyms)
            
            # Collect identifiers
            if source == 'pubchem':
                cid = data.get('cid') or data.get('properties', {}).get('CID')
                if cid:
                    merged['identifiers']['pubchem_cid'] = str(cid)
            elif source == 'chembl':
                chembl_id = data.get('chembl_id')
                if chembl_id:
                    merged['identifiers']['chembl_id'] = chembl_id
        
        # Convert synonyms set back to list
        merged['synonyms'] = list(merged['synonyms'])[:20]  # Limit to 20 synonyms
        
        return merged
    
    @staticmethod
    def aggregate_bioactivity_data(source_data: Dict[str, Any]) -> Dict[str, Any]:
        """Aggregate bioactivity data from multiple sources."""
        bioactivity = {
            'total_bioactivities': 0,
            'unique_targets': set(),
            'activity_types': set(),
            'bioactivity_summary': []
        }
        
        # Process ChEMBL bioactivity data
        if 'chembl' in source_data and source_data['chembl'].get('bioactivity_data'):
            activities = source_data['chembl']['bioactivity_data']
            bioactivity['total_bioactivities'] += len(activities)
            
            for activity in activities[:10]:  # Limit to top 10
                if activity.get('target_chembl_id'):
                    bioactivity['unique_targets'].add(activity['target_chembl_id'])
                if activity.get('standard_type'):
                    bioactivity['activity_types'].add(activity['standard_type'])
                
                bioactivity['bioactivity_summary'].append({
                    'target': activity.get('target_pref_name', 'Unknown'),
                    'activity_type': activity.get('standard_type'),
                    'value': activity.get('standard_value'),
                    'units': activity.get('standard_units'),
                    'relation': activity.get('standard_relation')
                })
        
        # Convert sets to lists for JSON serialization
        bioactivity['unique_targets'] = list(bioactivity['unique_targets'])
        bioactivity['activity_types'] = list(bioactivity['activity_types'])
        
        return bioactivity

class UnifiedSearchEngine:
    """Main unified search engine for pharmaceutical data."""
    
    def __init__(self):
        self.api_manager = APIIntegrationManager()
        self.cache = IntelligentCache()
        self.aggregator = DataAggregator()
        
    def search_compound(self, compound_name: str, use_cache: bool = True, 
                       parallel_search: bool = True) -> AggregatedResult:
        """
        Perform unified search across all databases with intelligent aggregation.
        
        Args:
            compound_name: Name of compound to search
            use_cache: Whether to use cached results
            parallel_search: Whether to search databases in parallel
        
        Returns:
            AggregatedResult with comprehensive data
        """
        # Check cache first
        if use_cache:
            cached_result = self.cache.get(compound_name)
            if cached_result:
                logger.info(f"Returning cached result for {compound_name}")
                return cached_result
        
        # Perform comprehensive search
        if parallel_search:
            comprehensive_data = self._parallel_search(compound_name)
        else:
            comprehensive_data = self.api_manager.search_compound_comprehensive(compound_name)
        
        # Aggregate and process results
        aggregated_result = self._aggregate_results(compound_name, comprehensive_data)
        
        # Cache result
        if use_cache:
            self.cache.set(compound_name, aggregated_result)
        
        return aggregated_result
    
    def _parallel_search(self, compound_name: str) -> Dict[str, Any]:
        """Search databases in parallel for faster results."""
        search_functions = [
            ('pubchem', lambda: self.api_manager.pubchem.get_compound_by_name(compound_name)),
            ('fda', lambda: self.api_manager.fda.search_drugs_by_name(compound_name)),
            ('chembl', lambda: self.api_manager.chembl.search_compound_by_name(compound_name)),
            ('drugbank', lambda: self.api_manager.drugbank.search_drug_by_name(compound_name)),
            ('zinc', lambda: self.api_manager.zinc.search_compound_by_name(compound_name)),
            ('opentargets', lambda: self.api_manager.opentargets.search_drug_by_name(compound_name))
        ]
        
        results = {
            'query': compound_name,
            'sources': {},
            'summary': {
                'total_sources': len(search_functions),
                'successful_sources': 0,
                'failed_sources': []
            }
        }
        
        with ThreadPoolExecutor(max_workers=6) as executor:
            # Submit all searches
            future_to_source = {
                executor.submit(func): source 
                for source, func in search_functions
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_source, timeout=60):
                source = future_to_source[future]
                try:
                    result = future.result(timeout=30)
                    if result:
                        results['sources'][source] = result
                        results['summary']['successful_sources'] += 1
                    else:
                        results['summary']['failed_sources'].append(source)
                except Exception as e:
                    logger.error(f"Parallel search failed for {source}: {str(e)}")
                    results['summary']['failed_sources'].append(source)
        
        return results
    
    def _aggregate_results(self, compound_name: str, comprehensive_data: Dict[str, Any]) -> AggregatedResult:
        """Aggregate results from multiple sources into unified format."""
        sources = list(comprehensive_data.get('sources', {}).keys())
        source_data = comprehensive_data.get('sources', {})
        
        # Calculate data quality scores
        data_quality = {}
        for source, data in source_data.items():
            quality_score = self._calculate_data_quality(source, data)
            data_quality[source] = quality_score
        
        # Calculate overall confidence
        confidence_score = self.aggregator.calculate_confidence_score(sources, data_quality)
        
        # Merge molecular data
        molecular_data = self.aggregator.merge_molecular_data(source_data)
        
        # Aggregate bioactivity data
        bioactivity_data = self.aggregator.aggregate_bioactivity_data(source_data)
        
        # Extract regulatory information
        regulatory_info = self._extract_regulatory_info(source_data)
        
        # Extract commercial availability
        commercial_info = self._extract_commercial_info(source_data)
        
        # Create comprehensive data structure
        comprehensive_data_dict = {
            'molecular_data': molecular_data,
            'bioactivity_data': bioactivity_data,
            'regulatory_info': regulatory_info,
            'commercial_info': commercial_info,
            'raw_source_data': source_data
        }
        
        # Create search metadata
        search_metadata = {
            'search_time': datetime.now().isoformat(),
            'sources_queried': comprehensive_data.get('summary', {}).get('total_sources', 0),
            'sources_successful': len(sources),
            'data_quality_scores': data_quality,
            'search_strategy': 'parallel' if len(sources) > 1 else 'sequential'
        }
        
        return AggregatedResult(
            compound_name=compound_name,
            primary_identifier=molecular_data.get('identifiers', {}).get('chembl_id') or 
                             molecular_data.get('identifiers', {}).get('pubchem_cid'),
            molecular_formula=molecular_data.get('molecular_formula'),
            molecular_weight=molecular_data.get('molecular_weight'),
            sources=sources,
            confidence_score=confidence_score,
            comprehensive_data=comprehensive_data_dict,
            search_metadata=search_metadata
        )
    
    def _calculate_data_quality(self, source: str, data: Any) -> float:
        """Calculate data quality score for a source."""
        if not data:
            return 0.0
        
        quality_score = 0.5  # Base score
        
        if source == 'pubchem':
            if isinstance(data, dict):
                if data.get('properties', {}).get('MolecularFormula'):
                    quality_score += 0.2
                if data.get('properties', {}).get('CanonicalSMILES'):
                    quality_score += 0.2
                if data.get('synonyms'):
                    quality_score += 0.1
        
        elif source == 'chembl':
            if isinstance(data, dict):
                if data.get('chembl_id'):
                    quality_score += 0.2
                if data.get('bioactivities'):
                    quality_score += 0.2
                if data.get('targets'):
                    quality_score += 0.1
        
        elif source == 'fda':
            if isinstance(data, list) and len(data) > 0:
                quality_score += 0.3
        
        return min(quality_score, 1.0)
    
    def _extract_regulatory_info(self, source_data: Dict[str, Any]) -> Dict[str, Any]:
        """Extract regulatory information from source data."""
        regulatory = {
            'fda_approved': False,
            'marketing_status': [],
            'therapeutic_equivalence': [],
            'approval_date': None
        }
        
        if 'fda' in source_data and source_data['fda']:
            regulatory['fda_approved'] = True
            fda_data = source_data['fda']
            
            if isinstance(fda_data, list):
                for drug in fda_data:
                    for product in drug.get('products', []):
                        if product.get('marketing_status'):
                            regulatory['marketing_status'].append(product['marketing_status'])
                        if product.get('te_code'):
                            regulatory['therapeutic_equivalence'].append(product['te_code'])
        
        return regulatory
    
    def _extract_commercial_info(self, source_data: Dict[str, Any]) -> Dict[str, Any]:
        """Extract commercial availability information."""
        commercial = {
            'commercially_available': False,
            'zinc_results': 0,
            'suppliers': []
        }
        
        if 'zinc' in source_data and source_data['zinc']:
            zinc_data = source_data['zinc']
            commercial['commercially_available'] = zinc_data.get('commercial_availability', False)
            commercial['zinc_results'] = zinc_data.get('total_results', 0)
        
        return commercial
    
    def batch_search(self, compound_names: List[str], max_workers: int = 4) -> Dict[str, AggregatedResult]:
        """Perform batch search for multiple compounds."""
        results = {}
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_compound = {
                executor.submit(self.search_compound, compound): compound 
                for compound in compound_names
            }
            
            for future in as_completed(future_to_compound):
                compound = future_to_compound[future]
                try:
                    result = future.result()
                    results[compound] = result
                except Exception as e:
                    logger.error(f"Batch search failed for {compound}: {str(e)}")
                    # Create minimal result for failed searches
                    results[compound] = AggregatedResult(
                        compound_name=compound,
                        primary_identifier=None,
                        molecular_formula=None,
                        molecular_weight=None,
                        sources=[],
                        confidence_score=0.0,
                        comprehensive_data={'error': str(e)},
                        search_metadata={'error': True}
                    )
        
        return results

# Global instance
unified_search = UnifiedSearchEngine()

# Test function
if __name__ == "__main__":
    print("Testing Unified Search Engine...")
    
    # Test single compound search
    result = unified_search.search_compound("aspirin")
    print(f"Search completed for aspirin:")
    print(f"- Sources: {result.sources}")
    print(f"- Confidence: {result.confidence_score:.2f}")
    print(f"- Molecular Formula: {result.molecular_formula}")
    print(f"- Molecular Weight: {result.molecular_weight}")
    print(f"- Primary ID: {result.primary_identifier}")
    
    # Test batch search
    print("\nTesting batch search...")
    batch_results = unified_search.batch_search(["aspirin", "ibuprofen"])
    for compound, result in batch_results.items():
        print(f"{compound}: {len(result.sources)} sources, confidence {result.confidence_score:.2f}")
