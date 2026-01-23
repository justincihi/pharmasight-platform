"""
BioTransformer Service Configuration
"""

import os

# Service configuration
SERVICE_NAME = "biotransformer"
VERSION = "1.0.0"
PORT = int(os.getenv('PORT', 8000))
DEBUG = os.getenv('DEBUG', 'False').lower() == 'true'

# BioTransformer configuration
BIOTRANSFORMER_JAR = os.getenv('BIOTRANSFORMER_JAR', '/opt/BioTransformer3.0.jar')
JAVA_OPTS = os.getenv('JAVA_OPTS', '-Xmx4g')  # 4GB heap size

# Timeout settings
PREDICTION_TIMEOUT = int(os.getenv('PREDICTION_TIMEOUT', 300))  # 5 minutes
BATCH_TIMEOUT = int(os.getenv('BATCH_TIMEOUT', 1800))  # 30 minutes

# Default settings
DEFAULT_METABOLISM_TYPE = 'human'
DEFAULT_STEPS = 1
MAX_STEPS = 3
MAX_BATCH_SIZE = int(os.getenv('MAX_BATCH_SIZE', 100))

# Logging
LOG_LEVEL = os.getenv('LOG_LEVEL', 'INFO')
