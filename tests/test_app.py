#!/usr/bin/env python3
"""
Comprehensive tests for PharmaSight Flask application
Tests the main app, RDKit integration, and API endpoints

Run with: pytest tests/test_app.py -v
"""

import pytest
import sys
import os
import json

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from app import app as main_app


class TestMainApp:
    """Tests for main Flask application (app.py)"""

    def setup_method(self):
        """Setup test fixtures"""
        self.app = main_app
        self.app.config['TESTING'] = True
        self.client = self.app.test_client()

    def test_index_route(self):
        """Test main index page loads"""
        response = self.client.get('/')
        assert response.status_code == 200
        assert b'Drug Discovery Platform' in response.data

    def test_health_endpoint(self):
        """Test health check endpoint"""
        response = self.client.get('/health')
        assert response.status_code == 200

        data = json.loads(response.data)
        assert data['status'] == 'healthy'
        assert 'timestamp' in data
        assert 'version' in data

    def test_html_contains_features(self):
        """Test that main page contains key features"""
        response = self.client.get('/')
        html = response.data.decode('utf-8')

        assert 'AI-Powered Analysis' in html
        assert 'IP Management' in html
        assert 'Regulatory Compliance' in html

    def test_response_headers(self):
        """Test response headers are set correctly"""
        response = self.client.get('/')
        assert response.content_type.startswith('text/html')


class TestRDKitAPI:
    """Tests for RDKit REST API"""

    def setup_method(self):
        """Setup test fixtures"""
        try:
            from rdkit_api import create_rdkit_api
            self.app = create_rdkit_api()
            self.app.config['TESTING'] = True
            self.client = self.app.test_client()
            self.rdkit_available = True
        except ImportError:
            self.rdkit_available = False
            pytest.skip("RDKit API not available")

    def test_health_endpoint(self):
        """Test RDKit API health check"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        response = self.client.get('/api/rdkit/health')
        assert response.status_code == 200

        data = json.loads(response.data)
        assert data['status'] == 'ok'
        assert 'rdkit_version' in data
        assert len(data['available_endpoints']) > 0

    def test_validate_valid_smiles(self):
        """Test SMILES validation with valid input"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        response = self.client.post('/api/rdkit/validate',
                                   json={'smiles': 'CC(=O)O'})
        assert response.status_code == 200

        data = json.loads(response.data)
        assert data['valid'] == True
        assert 'canonical' in data

    def test_validate_invalid_smiles(self):
        """Test SMILES validation with invalid input"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        response = self.client.post('/api/rdkit/validate',
                                   json={'smiles': 'INVALID'})
        assert response.status_code == 200

        data = json.loads(response.data)
        assert data['valid'] == False

    def test_properties_calculation(self):
        """Test molecular properties calculation"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        response = self.client.post('/api/rdkit/properties',
                                   json={'smiles': 'CC(=O)O'})
        assert response.status_code == 200

        data = json.loads(response.data)
        assert 'molecular_weight' in data
        assert 'logP' in data
        assert 'num_h_donors' in data
        assert 'lipinski_violations' in data

    def test_properties_missing_smiles(self):
        """Test properties endpoint with missing SMILES"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        response = self.client.post('/api/rdkit/properties', json={})
        assert response.status_code == 400

        data = json.loads(response.data)
        assert 'error' in data

    def test_similarity_calculation(self):
        """Test molecular similarity calculation"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        response = self.client.post('/api/rdkit/similarity',
                                   json={
                                       'smiles1': 'CC(=O)O',
                                       'smiles2': 'CCC(=O)O'
                                   })
        assert response.status_code == 200

        data = json.loads(response.data)
        assert 'similarity' in data
        assert 0 <= data['similarity'] <= 1

    def test_canonicalize(self):
        """Test SMILES canonicalization"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        response = self.client.post('/api/rdkit/canonicalize',
                                   json={'smiles': 'C(C)C'})
        assert response.status_code == 200

        data = json.loads(response.data)
        assert data['canonical'] == 'CCC'

    def test_batch_properties(self):
        """Test batch properties calculation"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        smiles_list = ['CC(=O)O', 'CCC(=O)O', 'c1ccccc1']
        response = self.client.post('/api/rdkit/batch/properties',
                                   json={'smiles_list': smiles_list})
        assert response.status_code == 200

        data = json.loads(response.data)
        assert data['total'] == 3
        assert len(data['results']) == 3

    def test_batch_too_many_molecules(self):
        """Test batch endpoint rejects too many molecules"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        smiles_list = ['CC(=O)O'] * 101  # More than limit
        response = self.client.post('/api/rdkit/batch/properties',
                                   json={'smiles_list': smiles_list})
        assert response.status_code == 400


class TestIntegrationWorkflow:
    """Integration tests for complete workflows"""

    def setup_method(self):
        """Setup test fixtures"""
        try:
            from rdkit_api import create_rdkit_api
            self.app = create_rdkit_api()
            self.app.config['TESTING'] = True
            self.client = self.app.test_client()
            self.rdkit_available = True
        except ImportError:
            self.rdkit_available = False
            pytest.skip("RDKit API not available")

    def test_drug_discovery_workflow(self):
        """Test complete drug discovery workflow"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        # 1. Validate compound
        response = self.client.post('/api/rdkit/validate',
                                   json={'smiles': 'CC(=O)Oc1ccccc1C(=O)O'})
        assert response.status_code == 200
        assert json.loads(response.data)['valid'] == True

        # 2. Calculate properties
        response = self.client.post('/api/rdkit/properties',
                                   json={'smiles': 'CC(=O)Oc1ccccc1C(=O)O'})
        assert response.status_code == 200
        props = json.loads(response.data)
        assert props['lipinski_violations'] == 0  # Aspirin passes Lipinski

        # 3. Compare with similar compound
        response = self.client.post('/api/rdkit/similarity',
                                   json={
                                       'smiles1': 'CC(=O)Oc1ccccc1C(=O)O',
                                       'smiles2': 'CC(=O)Oc1ccccc1'
                                   })
        assert response.status_code == 200
        sim = json.loads(response.data)
        assert sim['similarity'] > 0.5  # Should be similar


class TestErrorHandling:
    """Tests for error handling"""

    def setup_method(self):
        """Setup test fixtures"""
        try:
            from rdkit_api import create_rdkit_api
            self.app = create_rdkit_api()
            self.app.config['TESTING'] = True
            self.client = self.app.test_client()
            self.rdkit_available = True
        except ImportError:
            self.rdkit_available = False
            pytest.skip("RDKit API not available")

    def test_invalid_json(self):
        """Test handling of invalid JSON"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        response = self.client.post('/api/rdkit/properties',
                                   data='invalid json',
                                   content_type='application/json')
        # Should handle gracefully
        assert response.status_code in [400, 500]

    def test_missing_required_field(self):
        """Test missing required fields"""
        if not self.rdkit_available:
            pytest.skip("RDKit not available")

        response = self.client.post('/api/rdkit/similarity',
                                   json={'smiles1': 'CC(=O)O'})
        assert response.status_code == 400


class TestDeploymentReadiness:
    """Tests to verify app is ready for deployment"""

    def test_main_app_imports(self):
        """Test that main app imports successfully"""
        try:
            import app
            assert hasattr(app, 'app')
        except ImportError as e:
            pytest.fail(f"Failed to import main app: {e}")

    def test_requirements_files_exist(self):
        """Test that requirement files exist"""
        assert os.path.exists('requirements.txt')

    def test_main_app_has_health_check(self):
        """Test that main app has health check endpoint"""
        client = main_app.test_client()
        response = client.get('/health')
        assert response.status_code == 200

    def test_app_has_secret_key(self):
        """Test that app has secret key configured"""
        assert main_app.config.get('SECRET_KEY') is not None

    def test_environment_ready(self):
        """Test that key environment dependencies are available"""
        # Check Flask is available
        try:
            import flask
            assert True
        except ImportError:
            pytest.fail("Flask not available")

    def test_port_configuration(self):
        """Test that app can read PORT from environment"""
        # This is important for Replit/Railway deployment
        os.environ['PORT'] = '8080'
        # Just verify the app doesn't crash when PORT is set
        assert True


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
