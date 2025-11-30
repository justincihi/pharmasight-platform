# Quick Fixes for PharmaSight Platform

## Issue 1: Analog Generation API Mismatch

**File**: `src/analog_generation_fix.py`

**Current Code**:
```python
def resolve_compound_name(compound_name):
    # Returns string SMILES
    return smiles_string
```

**Fixed Code**:
```python
def resolve_compound_name(compound_name):
    # Returns dict with name and SMILES
    smiles = get_smiles(compound_name)
    return {
        "name": compound_name,
        "smiles": smiles,
        "source": "database" if smiles else "not_found"
    }
```

---

## Issue 2: Research Findings Function Signature

**File**: `src/research_findings_fix.py`

**Current Code**:
```python
def get_research_findings_with_hypotheses():
    # No parameters
    return all_findings
```

**Fixed Code**:
```python
def get_research_findings_with_hypotheses(compound_name=None):
    # Optional compound filter
    all_findings = load_findings()
    if compound_name:
        return [f for f in all_findings if f['compound'] == compound_name]
    return all_findings
```

---

## Issue 3: Frontend JavaScript Login Button

**File**: `src/pharmasight_complete.py` (in HTML template)

**Add JavaScript**:
```javascript
document.addEventListener('DOMContentLoaded', function() {
    const loginBtn = document.getElementById('loginBtn');
    if (loginBtn) {
        loginBtn.addEventListener('click', function(e) {
            e.preventDefault();
            const username = document.getElementById('username').value;
            const password = document.getElementById('password').value;
            
            // Simple validation
            if (username && password) {
                // Hide login, show main content
                document.getElementById('loginPortal').style.display = 'none';
                document.getElementById('mainContent').style.display = 'block';
            } else {
                alert('Please enter username and password');
            }
        });
    }
});
```

---

## Issue 4: Missing Static Resources

**Create**: `static/` directory structure

```bash
mkdir -p static/css static/js static/images
```

**Add to Flask app**:
```python
app = Flask(__name__, static_folder='static', static_url_path='/static')
```

---

## Testing Commands

```bash
# Test RDKit integration
python3.11 -c "from rdkit import Chem; print('✅ RDKit OK')"

# Test Flask app
curl http://localhost:5000/health

# Run diagnostic suite
python3.11 diagnostic_test.py
```
