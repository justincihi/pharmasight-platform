"""
API endpoints for Compound Discovery Logging System
Provides REST API for logging and querying discovery data
"""

from flask import Blueprint, request, jsonify, session
from discovery_logging_system import get_discovery_logger
import json
from datetime import datetime

# Create blueprint for discovery logging API
discovery_logging_bp = Blueprint("discovery_logging", __name__, url_prefix="/api/logging")

def require_admin():
    """Check if user has admin privileges."""
    if not session.get("logged_in"):
        return False
    return True

@discovery_logging_bp.route("/discoveries", methods=["GET"])
def get_discoveries():
    """Get discovery logs with optional filtering."""
    if not require_admin():
        return jsonify({"error": "Admin access required"}), 403
    
    try:
        logger = get_discovery_logger()
        
        confidence_threshold = request.args.get("confidence_threshold", type=float)
        ip_status = request.args.get("ip_status")
        limit = request.args.get("limit", 100, type=int)
        offset = request.args.get("offset", 0, type=int)
        
        discoveries = logger.get_discoveries(confidence_threshold, ip_status, limit, offset)
        
        discoveries_data = []
        for discovery in discoveries:
            discovery_dict = {
                "id": discovery.id,
                "compound_name": discovery.compound_name,
                "smiles": discovery.smiles,
                "category": discovery.category,
                "confidence_score": discovery.confidence_score,
                "ip_status": discovery.ip_status,
                "discovery_date": discovery.discovery_date,
                "source_research_topic_id": discovery.source_research_topic_id,
                "generated_by": discovery.generated_by,
                "molecular_weight": discovery.molecular_weight,
                "logp": discovery.logp,
                "receptor_activity": discovery.receptor_activity,
                "binding_affinity": discovery.binding_affinity,
                "notes": discovery.notes
            }
            discoveries_data.append(discovery_dict)
        
        return jsonify({
            "discoveries": discoveries_data,
            "total": len(discoveries_data),
            "status": "success"
        })
    
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@discovery_logging_bp.route("/discoveries/<int:discovery_id>/references", methods=["GET"])
def get_references(discovery_id):
    """Get all research references for a specific discovery."""
    if not require_admin():
        return jsonify({"error": "Admin access required"}), 403
    
    try:
        logger = get_discovery_logger()
        references = logger.get_references_for_discovery(discovery_id)
        
        references_data = []
        for ref in references:
            ref_dict = {
                "id": ref.id,
                "doi": ref.doi,
                "title": ref.title,
                "authors": ref.authors,
                "journal": ref.journal,
                "publication_year": ref.publication_year,
                "retrieval_date": ref.retrieval_date
            }
            references_data.append(ref_dict)
        
        return jsonify({
            "references": references_data,
            "total": len(references_data),
            "status": "success"
        })
    
    except Exception as e:
        return jsonify({"error": str(e)}), 500

# Export the blueprint
def get_discovery_logging_blueprint():
    """Get the discovery logging blueprint for registration."""
    return discovery_logging_bp
