# PharmaSight Platform - Deployment Guide

Complete guide for deploying PharmaSight to Replit, Railway, and other platforms.

---

## 🚀 Quick Deploy Options

### **Option 1: Replit (Easiest - 2 minutes)**

[![Run on Replit](https://replit.com/badge/github/yourusername/pharmasight-platform)](https://replit.com)

**Steps:**

1. **Go to [Replit.com](https://replit.com)**

2. **Click "Create Repl"** → "Import from GitHub"

3. **Paste your repository URL:**
   ```
   https://github.com/justincihi/pharmasight-platform
   ```

4. **Replit will auto-detect** the configuration from `.replit` file

5. **Click "Run"** - That's it! Your app will be live.

6. **Access your app:**
   - Development: `https://your-repl-name.your-username.repl.co`
   - After deployment: You'll get a permanent URL

**Configuration:**
- ✅ `.replit` file included
- ✅ `replit.nix` for dependencies
- ✅ Runs on port 5000 automatically
- ✅ Environment variables set

---

### **Option 2: Railway (Fast - 3 minutes)**

[![Deploy on Railway](https://railway.app/button.svg)](https://railway.app)

**Steps:**

1. **Go to [Railway.app](https://railway.app)**

2. **Click "New Project"** → "Deploy from GitHub repo"

3. **Connect your GitHub account** and select `pharmasight-platform`

4. **Railway auto-detects** configuration from `railway.json`

5. **Add environment variables** (optional):
   ```
   PORT=5000
   ```

6. **Deploy!** Railway will:
   - Install dependencies from `requirements.txt`
   - Run `python app.py` automatically
   - Provide HTTPS URL

7. **Access your app:**
   ```
   https://pharmasight-production.up.railway.app
   ```

**Configuration:**
- ✅ `railway.json` included
- ✅ `Procfile` for production
- ✅ Auto-scaling enabled
- ✅ HTTPS by default

---

### **Option 3: Render**

**Steps:**

1. **Go to [Render.com](https://render.com)**

2. **Click "New +"** → "Web Service"

3. **Connect GitHub** and select repository

4. **Configure:**
   - **Build Command:** `pip install -r requirements.txt`
   - **Start Command:** `gunicorn app:app`
   - **Environment:** Python 3.11

5. **Deploy!**

---

### **Option 4: Heroku**

**Steps:**

1. **Install Heroku CLI:**
   ```bash
   brew install heroku/brew/heroku  # Mac
   # or download from heroku.com
   ```

2. **Login and create app:**
   ```bash
   heroku login
   heroku create pharmasight-platform
   ```

3. **Deploy:**
   ```bash
   git push heroku main
   ```

4. **Open app:**
   ```bash
   heroku open
   ```

**Configuration:**
- ✅ `Procfile` included
- ✅ `runtime.txt` specifies Python version
- ✅ `requirements.txt` has all dependencies

---

## 🌐 Local Preview (Test Before Deploy)

### **Method 1: Simple Python Server**

```bash
cd /home/user/pharmasight-platform
python app.py
```

Open: `http://localhost:5000`

### **Method 2: Production-like with Gunicorn**

```bash
gunicorn app:app --bind 0.0.0.0:5000
```

Open: `http://localhost:5000`

### **Method 3: Test with ngrok (Public URL)**

```bash
# Install ngrok
curl -s https://ngrok-agent.s3.amazonaws.com/ngrok.asc | sudo tee /etc/apt/trusted.gpg.d/ngrok.asc
echo "deb https://ngrok-agent.s3.amazonaws.com buster main" | sudo tee /etc/apt/sources.list.d/ngrok.list
sudo apt update && sudo apt install ngrok

# Run your app
python app.py &

# Expose to internet
ngrok http 5000
```

You'll get a public URL like: `https://abcd-12-34-56-78.ngrok-free.app`

---

## 📋 Pre-Deployment Checklist

Before deploying, verify everything works:

### **1. Run Tests**

```bash
# Activate conda environment
conda activate pharmasight

# Run all tests
pytest tests/test_app.py -v
pytest tests/test_rdkit_integration.py -v

# Expected output:
# 30+ tests passing
```

### **2. Test Main App**

```bash
python app.py
# Visit http://localhost:5000
# Verify the UI loads
```

### **3. Test RDKit API**

```bash
cd src
python rdkit_api.py
# Visit http://localhost:5000/api/rdkit/health
# Should see: {"status": "ok", "rdkit_version": "2025.09.1"}
```

### **4. Check Dependencies**

```bash
pip install -r requirements.txt
# Should install without errors
```

---

## ⚙️ Configuration Files Reference

### **requirements.txt**
```
Flask==3.1.1
gunicorn==21.2.0
flask-cors==4.0.0
```

**For RDKit features, add:**
```bash
# Note: RDKit is best installed via conda
# conda install -c conda-forge rdkit
```

### **Procfile** (for Heroku/Railway)
```
web: gunicorn app:app
```

### **.replit** (for Replit)
```
run = "python app.py"
entrypoint = "app.py"
```

### **railway.json** (for Railway)
```json
{
  "build": { "builder": "NIXPACKS" },
  "deploy": { "startCommand": "python app.py" }
}
```

---

## 🔧 Environment Variables

Set these in your deployment platform:

| Variable | Value | Required | Description |
|----------|-------|----------|-------------|
| `PORT` | `5000` | No | Port to run on (auto-set by most platforms) |
| `FLASK_ENV` | `production` | Recommended | Flask environment mode |
| `SECRET_KEY` | `your-secret-key` | Recommended | Flask secret key for sessions |

**Setting in Replit:**
- Secrets tab → Add `SECRET_KEY`

**Setting in Railway:**
- Variables tab → Add variables

**Setting in Heroku:**
```bash
heroku config:set SECRET_KEY=your-secret-key
```

---

## 🧪 Testing Deployed App

Once deployed, test these endpoints:

### **1. Health Check**
```bash
curl https://your-app.com/health

# Expected:
# {"status": "healthy", "timestamp": "...", "version": "..."}
```

### **2. Main Page**
```bash
curl https://your-app.com/

# Should return HTML with "Drug Discovery Platform"
```

### **3. RDKit API (if integrated)**
```bash
curl https://your-app.com/api/rdkit/health

# Expected:
# {"status": "ok", "rdkit_version": "2025.09.1"}
```

---

## 🐛 Troubleshooting

### **App won't start**

**Check logs:**
```bash
# Replit: View console output
# Railway: Deployments → View logs
# Heroku: heroku logs --tail
```

**Common issues:**
- Missing dependencies: Check `requirements.txt`
- Port binding: Use `PORT` from environment
- Python version: Verify `runtime.txt` matches available versions

### **Import errors**

**Solution:**
```bash
# Ensure all imports are available
pip install -r requirements.txt

# For RDKit features, note that RDKit needs conda
# Most platforms don't support conda, so RDKit features
# may not work in basic deployments
```

### **Static files not loading**

**For deployment, ensure Flask serves static files:**
```python
# In app.py
app = Flask(__name__, static_folder='molecular_images')
```

---

## 📦 Deployment Modes

### **Mode 1: Basic Flask App (No RDKit)**

**What works:**
- ✅ Main UI
- ✅ Health check
- ✅ Static content
- ❌ RDKit molecular analysis

**Platforms:** Replit, Railway, Render, Heroku

**Setup:** Just deploy - works out of the box!

### **Mode 2: Full RDKit Integration**

**What works:**
- ✅ Everything in Mode 1
- ✅ RDKit API endpoints
- ✅ Molecular visualization
- ✅ ADMET predictions

**Platforms:** Requires conda support or Docker

**Setup:**
1. Use Docker deployment
2. Or use platforms supporting conda (limited)

---

## 🐳 Docker Deployment (Advanced)

For full RDKit support, use Docker:

**Create `Dockerfile`:**
```dockerfile
FROM continuumio/miniconda3:latest

WORKDIR /app
COPY . /app

# Install dependencies
RUN conda env create -f environment.yml
RUN echo "source activate pharmasight" > ~/.bashrc

# Expose port
EXPOSE 5000

# Run app
CMD ["conda", "run", "-n", "pharmasight", "python", "app.py"]
```

**Deploy to:**
- Google Cloud Run
- AWS ECS
- Railway (with Dockerfile)
- Render (with Docker)

---

## 🎯 Recommended Deployment Strategy

### **For Quick Demo:**
→ **Use Replit** (2 minutes, free tier, public URL)

### **For Production:**
→ **Use Railway** (auto-scaling, HTTPS, $5/month)

### **For Full RDKit Features:**
→ **Use Docker on Cloud Run** (conda support, pay-per-use)

### **For Development:**
→ **Local with ngrok** (test with public URL)

---

## 📊 Platform Comparison

| Feature | Replit | Railway | Render | Heroku |
|---------|--------|---------|--------|--------|
| Free Tier | ✅ Yes | ✅ Limited | ✅ Yes | ✅ Limited |
| Auto HTTPS | ✅ | ✅ | ✅ | ✅ |
| Custom Domain | ✅ Paid | ✅ | ✅ | ✅ |
| RDKit Support | ❌ | ⚠️ Docker | ⚠️ Docker | ⚠️ Docker |
| Deploy Time | 1-2 min | 2-3 min | 3-5 min | 2-4 min |
| Pricing | Free/$7/mo | $5/mo | Free/$7/mo | $7/mo |

---

## ✅ Success Checklist

After deployment:

- [ ] App loads at public URL
- [ ] Health endpoint returns 200
- [ ] Main page displays correctly
- [ ] No console errors
- [ ] Mobile responsive
- [ ] HTTPS working
- [ ] Logs show no errors

---

## 🆘 Support

**Issues?**
1. Check deployment logs
2. Verify `requirements.txt` is complete
3. Test locally first
4. Check platform-specific docs

**Need help?**
- Replit: https://docs.replit.com
- Railway: https://docs.railway.app
- Render: https://render.com/docs

---

## 🎉 Next Steps

Once deployed:

1. **Share your URL** with users
2. **Monitor logs** for errors
3. **Set up custom domain** (optional)
4. **Add analytics** (Google Analytics, etc.)
5. **Enable error tracking** (Sentry, etc.)

---

**Your app is ready to deploy!** Choose a platform and go live in minutes! 🚀
