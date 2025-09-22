import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

from pharmasight_complete import app

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5008, debug=False)
