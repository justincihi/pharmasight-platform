"""
API endpoints for Administrator Research Topic Management
Provides REST API for managing research topics dynamically
"""

from flask import Blueprint, request, jsonify, session
from admin_research_manager import get_research_topic_manager
import json
from datetime import datetime

# Create blueprint for admin research API
admin_research_bp = Blueprint('admin_research', __name__, url_prefix='/api/admin/research')

def require_admin():
    """Check if user has admin privileges."""
    # For now, any logged-in user can manage research topics
    # In production, implement proper role-based access control
    if not session.get('logged_in'):
        return False
    return True

@admin_research_bp.route('/topics', methods=['GET'])
def get_topics():
    """Get all research topics with optional filtering."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        manager = get_research_topic_manager()
        status = request.args.get('status')
        search = request.args.get('search')
        
        if search:
            topics = manager.search_topics(search)
        else:
            topics = manager.get_all_topics(status)
        
        # Convert to dict format for JSON response
        topics_data = []
        for topic in topics:
            topic_dict = {
                'id': topic.id,
                'title': topic.title,
                'description': topic.description,
                'keywords': topic.keywords,
                'priority': topic.priority,
                'status': topic.status,
                'created_by': topic.created_by,
                'created_date': topic.created_date,
                'last_modified': topic.last_modified,
                'target_compounds': topic.target_compounds,
                'confidence_threshold': topic.confidence_threshold,
                'ip_focus': topic.ip_focus,
                'therapeutic_areas': topic.therapeutic_areas
            }
            topics_data.append(topic_dict)
        
        return jsonify({
            'topics': topics_data,
            'total': len(topics_data),
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@admin_research_bp.route('/topics', methods=['POST'])
def create_topic():
    """Create a new research topic."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        data = request.get_json()
        
        # Validate required fields
        required_fields = ['title', 'description', 'keywords']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing required field: {field}'}), 400
        
        manager = get_research_topic_manager()
        
        # Parse keywords if they're provided as a string
        keywords = data['keywords']
        if isinstance(keywords, str):
            # Split by comma and clean up
            keywords = [k.strip() for k in keywords.split(',') if k.strip()]
        
        # Parse target compounds
        target_compounds = data.get('target_compounds', [])
        if isinstance(target_compounds, str):
            target_compounds = [c.strip() for c in target_compounds.split(',') if c.strip()]
        
        # Parse therapeutic areas
        therapeutic_areas = data.get('therapeutic_areas', [])
        if isinstance(therapeutic_areas, str):
            therapeutic_areas = [a.strip() for a in therapeutic_areas.split(',') if a.strip()]
        
        topic_id = manager.create_topic(
            title=data['title'],
            description=data['description'],
            keywords=keywords,
            priority=data.get('priority', 3),
            created_by=session.get('user', 'admin'),
            target_compounds=target_compounds,
            confidence_threshold=data.get('confidence_threshold', 0.8),
            ip_focus=data.get('ip_focus', True),
            therapeutic_areas=therapeutic_areas
        )
        
        return jsonify({
            'topic_id': topic_id,
            'message': 'Research topic created successfully',
            'status': 'success'
        }), 201
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@admin_research_bp.route('/topics/<int:topic_id>', methods=['GET'])
def get_topic(topic_id):
    """Get a specific research topic by ID."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        manager = get_research_topic_manager()
        topic = manager.get_topic_by_id(topic_id)
        
        if not topic:
            return jsonify({'error': 'Topic not found'}), 404
        
        topic_dict = {
            'id': topic.id,
            'title': topic.title,
            'description': topic.description,
            'keywords': topic.keywords,
            'priority': topic.priority,
            'status': topic.status,
            'created_by': topic.created_by,
            'created_date': topic.created_date,
            'last_modified': topic.last_modified,
            'target_compounds': topic.target_compounds,
            'confidence_threshold': topic.confidence_threshold,
            'ip_focus': topic.ip_focus,
            'therapeutic_areas': topic.therapeutic_areas
        }
        
        return jsonify({
            'topic': topic_dict,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@admin_research_bp.route('/topics/<int:topic_id>', methods=['PUT'])
def update_topic(topic_id):
    """Update a research topic."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        data = request.get_json()
        manager = get_research_topic_manager()
        
        # Parse keywords if they're provided as a string
        if 'keywords' in data and isinstance(data['keywords'], str):
            data['keywords'] = [k.strip() for k in data['keywords'].split(',') if k.strip()]
        
        # Parse target compounds
        if 'target_compounds' in data and isinstance(data['target_compounds'], str):
            data['target_compounds'] = [c.strip() for c in data['target_compounds'].split(',') if c.strip()]
        
        # Parse therapeutic areas
        if 'therapeutic_areas' in data and isinstance(data['therapeutic_areas'], str):
            data['therapeutic_areas'] = [a.strip() for a in data['therapeutic_areas'].split(',') if a.strip()]
        
        # Add modified_by
        data['modified_by'] = session.get('user', 'admin')
        
        success = manager.update_topic(topic_id, **data)
        
        if success:
            return jsonify({
                'message': 'Research topic updated successfully',
                'status': 'success'
            })
        else:
            return jsonify({'error': 'Topic not found'}), 404
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@admin_research_bp.route('/topics/<int:topic_id>', methods=['DELETE'])
def delete_topic(topic_id):
    """Delete a research topic."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        manager = get_research_topic_manager()
        success = manager.delete_topic(topic_id, session.get('user', 'admin'))
        
        if success:
            return jsonify({
                'message': 'Research topic deleted successfully',
                'status': 'success'
            })
        else:
            return jsonify({'error': 'Topic not found'}), 404
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@admin_research_bp.route('/topics/<int:topic_id>/history', methods=['GET'])
def get_topic_history(topic_id):
    """Get the history of changes for a topic."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        manager = get_research_topic_manager()
        history = manager.get_topic_history(topic_id)
        
        return jsonify({
            'history': history,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@admin_research_bp.route('/goals/active', methods=['GET'])
def get_active_goals():
    """Get active research goals formatted for the research engine."""
    try:
        manager = get_research_topic_manager()
        goals = manager.get_active_research_goals()
        
        return jsonify({
            'goals': goals,
            'total': len(goals),
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@admin_research_bp.route('/stats', methods=['GET'])
def get_research_stats():
    """Get statistics about research topics."""
    if not require_admin():
        return jsonify({'error': 'Admin access required'}), 403
    
    try:
        manager = get_research_topic_manager()
        
        all_topics = manager.get_all_topics()
        active_topics = manager.get_all_topics('active')
        completed_topics = manager.get_all_topics('completed')
        paused_topics = manager.get_all_topics('paused')
        
        # Calculate priority distribution
        priority_dist = {}
        for topic in active_topics:
            priority_dist[topic.priority] = priority_dist.get(topic.priority, 0) + 1
        
        # Calculate therapeutic area distribution
        therapeutic_areas = {}
        for topic in active_topics:
            for area in topic.therapeutic_areas:
                therapeutic_areas[area] = therapeutic_areas.get(area, 0) + 1
        
        stats = {
            'total_topics': len(all_topics),
            'active_topics': len(active_topics),
            'completed_topics': len(completed_topics),
            'paused_topics': len(paused_topics),
            'priority_distribution': priority_dist,
            'therapeutic_area_distribution': therapeutic_areas,
            'avg_confidence_threshold': sum(t.confidence_threshold for t in active_topics) / len(active_topics) if active_topics else 0,
            'ip_focused_topics': sum(1 for t in active_topics if t.ip_focus)
        }
        
        return jsonify({
            'stats': stats,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Export the blueprint
def get_admin_research_blueprint():
    """Get the admin research blueprint for registration."""
    return admin_research_bp
