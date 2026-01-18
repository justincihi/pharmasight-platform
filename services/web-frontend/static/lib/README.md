# JavaScript Libraries

This directory should contain:

1. **three.min.js** - Download from https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.min.js
2. **gsap.min.js** - Download from https://cdn.jsdelivr.net/npm/gsap@3.12.5/dist/gsap.min.js

## Quick Setup:

```bash
cd services/web-frontend/static/lib

# Download Three.js
curl -o three.min.js https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.min.js

# Download GSAP
curl -o gsap.min.js https://cdn.jsdelivr.net/npm/gsap@3.12.5/dist/gsap.min.js
```

Alternatively, you can use CDN links directly in the HTML (update index.html):

```html
<script src="https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/gsap@3.12.5/dist/gsap.min.js"></script>
```
