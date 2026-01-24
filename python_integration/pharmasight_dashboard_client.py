#!/usr/bin/env python3
"""
PharmaSight Dashboard API Client
Python integration for autonomous research engine

Usage:
    from pharmasight_dashboard_client import PharmaSightDashboard
    
    dashboard = PharmaSightDashboard(
        api_url="https://your-dashboard.manus.space",
        api_key="your-api-key-here"
    )
    
    # Import discoveries
    discoveries = [
        {
            "compoundId": "KETAMINE-NEW-001",
            "compoundName": "Novel Ketamine Analog",
            "smiles": "CNC1(c2ccccc2F)CCCCC1=O",
            "parentCompound": "Ketamine",
            "mechanismOfAction": "NMDA antagonist",
            "keyDifferences": "Fluorinated derivative",
            "confidence": 90,
            "similarity": 88,
            "safetyScore": 85,
            "efficacyScore": 90,
            "drugLikenessScore": 95,
            "patentStatus": "patent-free",
            "marketValue": "$50M",
            "discoveryMethod": "ai-generation"
        }
    ]
    
    result = dashboard.import_discoveries(discoveries)
    print(f"Imported: {result['imported']}, Skipped: {result['skipped']}")
"""

import requests
import json
from typing import List, Dict, Optional
from datetime import datetime

class PharmaSightDashboard:
    """Client for interacting with PharmaSight Admin Dashboard API"""
    
    def __init__(self, api_url: str, api_key: str):
        """
        Initialize dashboard client
        
        Args:
            api_url: Base URL of the dashboard (e.g., "https://your-dashboard.manus.space")
            api_key: API key for authentication (PLATFORM_API_KEY)
        """
        self.api_url = api_url.rstrip('/')
        self.api_key = api_key
        self.session = requests.Session()
    
    def health_check(self) -> Dict:
        """Check if the dashboard API is available"""
        try:
            response = self.session.get(f"{self.api_url}/api/platform/health", timeout=10)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            return {"status": "error", "message": str(e)}
    
    def import_discoveries(self, discoveries: List[Dict]) -> Dict:
        """
        Import analog discoveries to the dashboard
        
        Args:
            discoveries: List of discovery dictionaries with required fields:
                - compoundId: Unique identifier
                - compoundName: Human-readable name
                - smiles: SMILES notation
                - parentCompound: Parent compound name
                - mechanismOfAction: Description of mechanism
                - keyDifferences: What makes this analog unique
                - confidence: Confidence score (0-100)
                - similarity: Similarity to parent (0-100)
                - safetyScore: Safety score (0-100)
                - efficacyScore: Efficacy score (0-100)
                - drugLikenessScore: Drug-likeness score (0-100)
                - patentStatus: "patent-free", "patent-opportunity", etc.
                - marketValue: Estimated market value (e.g., "$50M")
                - discoveryMethod: How it was discovered
        
        Returns:
            Dict with import results: {imported, skipped, errors, details}
        """
        payload = {
            "apiKey": self.api_key,
            "discoveries": discoveries
        }
        
        try:
            response = self.session.post(
                f"{self.api_url}/api/platform/discoveries/import",
                json=payload,
                timeout=60
            )
            response.raise_for_status()
            return response.json()
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "imported": 0,
                "skipped": 0,
                "errors": len(discoveries)
            }
    
    def get_recent_discoveries(self, limit: int = 10) -> Dict:
        """
        Get recent analog discoveries from the dashboard
        
        Args:
            limit: Maximum number of discoveries to return
        
        Returns:
            Dict with discoveries list
        """
        try:
            response = self.session.get(
                f"{self.api_url}/api/platform/discoveries/recent",
                params={"apiKey": self.api_key, "limit": limit},
                timeout=30
            )
            response.raise_for_status()
            return response.json()
        except Exception as e:
            return {"success": False, "error": str(e), "discoveries": []}
    
    def get_analog(self, compound_id: str) -> Optional[Dict]:
        """
        Get specific analog by compound ID
        
        Args:
            compound_id: Compound identifier
        
        Returns:
            Analog data or None if not found
        """
        try:
            response = self.session.get(
                f"{self.api_url}/api/platform/analogs/{compound_id}",
                params={"apiKey": self.api_key},
                timeout=30
            )
            response.raise_for_status()
            result = response.json()
            return result.get("analog")
        except Exception as e:
            print(f"Error fetching analog {compound_id}: {e}")
            return None
    
    def update_analog(self, compound_id: str, update_data: Dict) -> bool:
        """
        Update analog data
        
        Args:
            compound_id: Compound identifier
            update_data: Fields to update
        
        Returns:
            True if successful, False otherwise
        """
        payload = {
            "apiKey": self.api_key,
            **update_data
        }
        
        try:
            response = self.session.put(
                f"{self.api_url}/api/platform/analogs/{compound_id}",
                json=payload,
                timeout=30
            )
            response.raise_for_status()
            return True
        except Exception as e:
            print(f"Error updating analog {compound_id}: {e}")
            return False


# Example usage for autonomous research engine
if __name__ == "__main__":
    import os
    
    # Configuration
    API_URL = os.getenv("PHARMASIGHT_DASHBOARD_URL", "http://localhost:3000")
    API_KEY = os.getenv("PLATFORM_API_KEY", "")
    
    if not API_KEY:
        print("Error: PLATFORM_API_KEY environment variable not set")
        exit(1)
    
    # Initialize client
    dashboard = PharmaSightDashboard(API_URL, API_KEY)
    
    # Health check
    health = dashboard.health_check()
    print(f"Dashboard status: {health.get('status')}")
    
    if health.get('status') != 'ok':
        print(f"Dashboard not available: {health.get('message')}")
        exit(1)
    
    # Example: Import a discovery
    test_discovery = {
        "compoundId": f"TEST-{datetime.now().strftime('%Y%m%d-%H%M%S')}",
        "compoundName": "Test Analog",
        "smiles": "CC(C)C",
        "parentCompound": "Test Parent",
        "mechanismOfAction": "Test mechanism",
        "keyDifferences": "Test differences",
        "confidence": 85,
        "similarity": 90,
        "safetyScore": 80,
        "efficacyScore": 85,
        "drugLikenessScore": 90,
        "patentStatus": "patent-free",
        "marketValue": "$10M",
        "discoveryMethod": "test-script"
    }
    
    print("\nImporting test discovery...")
    result = dashboard.import_discoveries([test_discovery])
    print(f"Result: {json.dumps(result, indent=2)}")
    
    # Get recent discoveries
    print("\nFetching recent discoveries...")
    recent = dashboard.get_recent_discoveries(limit=5)
    print(f"Found {recent.get('count', 0)} recent discoveries")
