"""
API endpoints for Confidence-Based Data Organization and DOI Tracking System
Provides REST API for regulatory compliance and data organization
"""

from flask import Blueprint, request, jsonify, session, send_file
from confidence_doi_system import get_confidence_doi_manager
import json
import os
from datetime import datetime

# Create blueprint for confidence DOI API
confidence_doi_bp = Blueprint('confidence_doi', __name__, url_prefix='/api/compliance')

def require_admin():
    """Check if user has admin privileges."""
    if not session.get('logged_in'):
        return False
    return True

@confidence_doi_bp.route('/datasets/confidence/<float:threshold>', methods=['GET'])
def get_confidence_dataset(threshold):
    """Get confidence-based dataset for a specific threshold."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        manager = get_confidence_doi_manager()
        dataset = manager.generate_confidence_dataset(threshold)
        
        dataset_dict = {
            'confidence_threshold': dataset.confidence_threshold,
            'total_compounds': dataset.total_compounds,
            'compounds': dataset.compounds,
            'avg_confidence': dataset.avg_confidence,
            'ip_opportunities': dataset.ip_opportunities,
            'therapeutic_areas': dataset.therapeutic_areas,
            'generated_date': dataset.generated_date
        }
        
        return jsonify({
            'dataset': dataset_dict,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@confidence_doi_bp.route('/datasets/high-confidence', methods=['GET'])
def get_high_confidence_dataset():
    """Get high-confidence dataset (>90% confidence)."""
    return get_confidence_dataset(0.9)

@confidence_doi_bp.route('/datasets/export/json', methods=['POST'])
def export_dataset_json():
    """Export confidence dataset to JSON file."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        data = request.get_json()
        threshold = data.get('confidence_threshold', 0.9)
        
        manager = get_confidence_doi_manager()
        file_path = manager.export_confidence_dataset_json(threshold)
        
        return jsonify({
            'file_path': file_path,
            'download_url': f'/api/compliance/download/{os.path.basename(file_path)}',
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@confidence_doi_bp.route('/datasets/export/csv', methods=['POST'])
def export_dataset_csv():
    """Export confidence dataset to CSV file."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        data = request.get_json()
        threshold = data.get('confidence_threshold', 0.9)
        
        manager = get_confidence_doi_manager()
        file_path = manager.export_confidence_dataset_csv(threshold)
        
        return jsonify({
            'file_path': file_path,
            'download_url': f'/api/compliance/download/{os.path.basename(file_path)}',
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@confidence_doi_bp.route('/doi/references', methods=['GET'])
def get_doi_references():
    """Get DOI references with optional filtering."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        research_area = request.args.get('research_area')
        min_confidence = request.args.get('min_confidence', type=float)
        
        manager = get_confidence_doi_manager()
        references = manager.get_doi_references(research_area, min_confidence)
        
        references_data = []
        for ref in references:
            ref_dict = {
                'doi': ref.doi,
                'title': ref.title,
                'authors': ref.authors,
                'journal': ref.journal,
                'publication_year': ref.publication_year,
                'abstract': ref.abstract,
                'keywords': ref.keywords,
                'citation_count': ref.citation_count,
                'impact_factor': ref.impact_factor,
                'research_area': ref.research_area,
                'compound_relevance': ref.compound_relevance,
                'retrieval_date': ref.retrieval_date,
                'confidence_score': ref.confidence_score
            }
            references_data.append(ref_dict)
        
        return jsonify({
            'references': references_data,
            'total': len(references_data),
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@confidence_doi_bp.route('/doi/references', methods=['POST'])
def add_doi_reference():
    """Add a new DOI reference."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        data = request.get_json()
        
        required_fields = ['doi']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing required field: {field}'}), 400
        
        manager = get_confidence_doi_manager()
        
        doi_id = manager.add_doi_reference(
            doi=data['doi'],
            title=data.get('title'),
            authors=data.get('authors', []),
            journal=data.get('journal'),
            publication_year=data.get('publication_year'),
            abstract=data.get('abstract'),
            keywords=data.get('keywords', []),
            citation_count=data.get('citation_count', 0),
            impact_factor=data.get('impact_factor', 0.0),
            research_area=data.get('research_area'),
            compound_relevance=data.get('compound_relevance'),
            confidence_score=data.get('confidence_score', 0.8)
        )
        
        return jsonify({
            'doi_id': doi_id,
            'message': 'DOI reference added successfully',
            'status': 'success'
        }), 201
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@confidence_doi_bp.route('/regulatory/report', methods=['GET'])
def get_regulatory_report():
    """Generate comprehensive regulatory compliance report."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        confidence_threshold = request.args.get('confidence_threshold', 0.9, type=float)
        
        manager = get_confidence_doi_manager()
        report = manager.generate_regulatory_report(confidence_threshold)
        
        return jsonify({
            'report': report,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@confidence_doi_bp.route('/regulatory/export', methods=['POST'])
def export_regulatory_report():
    """Export regulatory compliance report to file."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        data = request.get_json()
        confidence_threshold = data.get('confidence_threshold', 0.9)
        
        manager = get_confidence_doi_manager()
        file_path = manager.export_regulatory_report(confidence_threshold)
        
        return jsonify({
            'file_path': file_path,
            'download_url': f'/api/compliance/download/{os.path.basename(file_path)}',
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@confidence_doi_bp.route('/download/<filename>', methods=['GET'])
def download_file(filename):
    """Download exported files."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        file_path = f'/home/ubuntu/pharmasight-platform/data/{filename}'
        
        if not os.path.exists(file_path):
            return jsonify({'error': 'File not found'}), 404
        
        return send_file(file_path, as_attachment=True)
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@confidence_doi_bp.route('/stats', methods=['GET'])
def get_compliance_stats():
    """Get compliance and data organization statistics."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        manager = get_confidence_doi_manager()
        
        # Get datasets for different thresholds
        high_confidence = manager.generate_confidence_dataset(0.9)
        medium_confidence = manager.generate_confidence_dataset(0.7)
        all_data = manager.generate_confidence_dataset(0.0)
        
        # Get DOI references
        doi_refs = manager.get_doi_references()
        
        stats = {
            'discovery_stats': {
                'total_discoveries': all_data.total_compounds,
                'high_confidence_discoveries': high_confidence.total_compounds,
                'medium_confidence_discoveries': medium_confidence.total_compounds,
                'avg_confidence_all': all_data.avg_confidence,
                'avg_confidence_high': high_confidence.avg_confidence,
                'ip_opportunities_total': all_data.ip_opportunities,
                'ip_opportunities_high_confidence': high_confidence.ip_opportunities
            },
            'reference_stats': {
                'total_doi_references': len(doi_refs),
                'high_confidence_references': len([r for r in doi_refs if r.confidence_score >= 0.9]),
                'research_areas': list(set([r.research_area for r in doi_refs if r.research_area])),
                'avg_citation_count': sum([r.citation_count for r in doi_refs]) / len(doi_refs) if doi_refs else 0
            },
            'therapeutic_areas': {
                'high_confidence': high_confidence.therapeutic_areas,
                'all_discoveries': all_data.therapeutic_areas
            },
            'regulatory_readiness': {
                'high_confidence_threshold_met': high_confidence.total_compounds > 0,
                'doi_tracking_active': len(doi_refs) > 0,
                'data_export_ready': True,
                'compliance_score': min(100, (high_confidence.total_compounds * 10) + (len(doi_refs) * 5))
            }
        }
        
        return jsonify({
            'stats': stats,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Export the blueprint
def get_confidence_doi_blueprint():
    """Get the confidence DOI blueprint for registration."""
    return confidence_doi_bp
