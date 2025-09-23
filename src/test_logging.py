import os
from discovery_logging_system import DiscoveryLogger

print('Testing Discovery Logging System...')
logger = DiscoveryLogger()
discoveries = logger.get_discoveries()
print(f'Total discoveries: {len(discoveries)}')

if discoveries:
    for discovery in discoveries[:2]:
        print(f'- {discovery.compound_name} (SMILES: {discovery.smiles})')
else:
    # Log a test discovery if none exist
    logger.log_discovery(
        compound_name="TestCompound-001",
        smiles="C1=CC=C(C=C1)C(=O)O",
        category="Test-2024-A1",
        confidence_score=0.95,
        ip_status="Potential Opportunity",
        source_research_topic_id=1,
        molecular_weight=122.12,
        logp=1.87
    )
    print("Logged a test discovery.")
    discoveries = logger.get_discoveries()
    print(f'Total discoveries now: {len(discoveries)}')


print('Discovery Logging System working successfully!')

