# PharmaSight™ Deployment Guide

## 🚀 Permanent Deployment Options

This guide covers deploying PharmaSight™ to production cloud platforms.

---

## Option 1: Railway (Recommended)

Railway offers the easiest deployment with automatic GitHub integration.

### Steps:

1. **Sign up at Railway**
   - Go to https://railway.app
   - Sign in with your GitHub account

2. **Create New Project**
   - Click "New Project"
   - Select "Deploy from GitHub repo"
   - Choose `justincihi/pharmasight-platform`
   - Select branch: `analog-discoveries-ip-protected`

3. **Configure Environment**
   - Railway will auto-detect the configuration from `railway.json`
   - No additional environment variables needed
   - The app will automatically use the `$PORT` variable

4. **Deploy**
   - Railway will automatically build and deploy
   - You'll get a public URL like: `https://pharmasight-production.up.railway.app`

5. **Custom Domain (Optional)**
   - Go to Settings → Domains
   - Add your custom domain
   - Update DNS records as instructed

### Configuration Files:
- ✅ `railway.json` - Railway-specific configuration
- ✅ `Procfile` - Process configuration
- ✅ `requirements.txt` - Python dependencies
- ✅ `runtime.txt` - Python version

---

## Option 2: Render

Render provides a simple deployment with free tier.

### Steps:

1. **Sign up at Render**
   - Go to https://render.com
   - Sign in with GitHub

2. **Create Web Service**
   - Click "New +" → "Web Service"
   - Connect your GitHub repository
   - Select `justincihi/pharmasight-platform`
   - Branch: `analog-discoveries-ip-protected`

3. **Configure Service**
   - **Name:** pharmasight-platform
   - **Environment:** Python 3
   - **Build Command:** `pip install -r requirements.txt`
   - **Start Command:** `gunicorn -w 4 -b 0.0.0.0:$PORT src.pharmasight_complete:app --timeout 120`

4. **Deploy**
   - Click "Create Web Service"
   - Render will build and deploy automatically
   - You'll get a URL like: `https://pharmasight-platform.onrender.com`

---

## Option 3: Fly.io

Fly.io offers global deployment with edge computing.

### Steps:

1. **Install Fly CLI**
   ```bash
   curl -L https://fly.io/install.sh | sh
   ```

2. **Login to Fly**
   ```bash
   fly auth login
   ```

3. **Launch App**
   ```bash
   cd /home/ubuntu/pharmasight-latest
   fly launch --name pharmasight-platform
   ```

4. **Deploy**
   ```bash
   fly deploy
   ```

5. **Open App**
   ```bash
   fly open
   ```

---

## Option 4: Heroku

Classic platform with extensive documentation.

### Steps:

1. **Install Heroku CLI**
   ```bash
   curl https://cli-assets.heroku.com/install.sh | sh
   ```

2. **Login to Heroku**
   ```bash
   heroku login
   ```

3. **Create App**
   ```bash
   cd /home/ubuntu/pharmasight-latest
   heroku create pharmasight-platform
   ```

4. **Push to Heroku**
   ```bash
   git push heroku analog-discoveries-ip-protected:main
   ```

5. **Open App**
   ```bash
   heroku open
   ```

---

## 🔧 Configuration Details

### Required Files

All deployment files are already configured:

1. **requirements.txt** - Python dependencies
   - Flask 3.1.2
   - RDKit 2025.9.1
   - Gunicorn 23.0.0
   - All scientific libraries

2. **Procfile** - Process configuration
   ```
   web: gunicorn -w 4 -b 0.0.0.0:$PORT src.pharmasight_complete:app --timeout 120
   ```

3. **runtime.txt** - Python version
   ```
   python-3.11.14
   ```

4. **railway.json** - Railway-specific config
   - Auto-restart on failure
   - Nixpacks builder
   - Optimized start command

### Environment Variables

No environment variables are required for basic deployment. The application uses:
- JSON file storage (no database needed)
- Free PubMed API (no API key needed)
- Local RDKit processing (no external services)

### Port Configuration

The application automatically uses the `$PORT` environment variable provided by the hosting platform.

---

## 📊 Resource Requirements

### Minimum Requirements:
- **Memory:** 512 MB (recommended: 1 GB)
- **CPU:** 1 core (recommended: 2 cores)
- **Disk:** 1 GB
- **Python:** 3.11+

### Platform Recommendations:

| Platform | Free Tier | Recommended Plan | Cost |
|----------|-----------|------------------|------|
| Railway | 500 hours/month | Hobby ($5/month) | $5-10/month |
| Render | 750 hours/month | Starter ($7/month) | $7-15/month |
| Fly.io | 3 VMs free | Hobby ($1.94/month) | $2-10/month |
| Heroku | Eco ($5/month) | Basic ($7/month) | $7-15/month |

---

## 🔄 Automatic Deployments

All platforms support automatic deployments from GitHub:

1. **Enable Auto-Deploy**
   - Connect your GitHub repository
   - Select the branch: `analog-discoveries-ip-protected`
   - Enable automatic deployments

2. **Every Git Push Triggers Deployment**
   - Push to GitHub → Automatic build → Automatic deploy
   - No manual intervention needed

3. **Rollback Support**
   - All platforms support instant rollback to previous versions
   - Keep deployment history for debugging

---

## 🔒 Security Considerations

### HTTPS
- All platforms provide automatic HTTPS
- SSL certificates are managed automatically
- No configuration needed

### Data Protection
- All analog discoveries are already in GitHub (public registry)
- Research articles database is public
- No sensitive user data stored

### API Rate Limiting
- Autonomous research engine limited to 50 calls/day
- Built-in rate limiting (0.4s between calls)
- Cost-controlled operation

---

## 📝 Post-Deployment Checklist

After deployment, verify:

- [ ] Platform is accessible via public URL
- [ ] All pages load correctly
- [ ] Analog generation works
- [ ] Research article database is accessible
- [ ] Export functions work (PDF, CSV, Excel)
- [ ] PKPD simulations run
- [ ] DDI analysis functions
- [ ] Public registry is viewable

---

## 🆘 Troubleshooting

### Build Failures

**Issue:** RDKit installation fails
**Solution:** Most platforms support RDKit via pip. If issues occur, contact platform support.

**Issue:** Memory errors during build
**Solution:** Increase memory allocation in platform settings.

### Runtime Errors

**Issue:** Application won't start
**Solution:** Check logs for import errors. Verify all files are committed to GitHub.

**Issue:** Port binding errors
**Solution:** Ensure Procfile uses `$PORT` variable, not hardcoded port.

### Performance Issues

**Issue:** Slow response times
**Solution:** Increase worker count in Procfile or upgrade to higher tier.

**Issue:** Timeout errors
**Solution:** Increase timeout in Procfile (currently 120s).

---

## 📞 Support

For deployment issues:
- **Railway:** https://railway.app/help
- **Render:** https://render.com/docs
- **Fly.io:** https://fly.io/docs
- **Heroku:** https://devcenter.heroku.com

For PharmaSight™ issues:
- Submit feedback at https://help.manus.im
- Check GitHub Issues: https://github.com/justincihi/pharmasight-platform/issues

---

## ✅ Recommended: Railway Deployment

**Why Railway?**
- Easiest setup (one-click GitHub integration)
- Automatic HTTPS and domain management
- Built-in monitoring and logs
- Generous free tier (500 hours/month)
- Excellent Python/Flask support
- Auto-detects configuration from `railway.json`

**Quick Start:**
1. Go to https://railway.app
2. Sign in with GitHub
3. New Project → Deploy from GitHub
4. Select `justincihi/pharmasight-platform`
5. Done! Your app is live in ~2 minutes

---

**Last Updated:** 2025-11-11  
**Status:** Production Ready ✅

