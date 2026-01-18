# PharmaSight™ Web Frontend

**Futuristic, AI-Powered Drug Discovery Platform Interface**

A next-generation web interface featuring glassmorphism design, animated molecular visualizations, and comprehensive integrations showcase.

## Features

### 🎨 Design System
- **Glassmorphism UI**: Translucent cards with blur effects
- **Gradient Accents**: Cyan-to-purple color scheme
- **Smooth Animations**: 60fps CSS and GSAP animations
- **Responsive**: Mobile-first, works on all devices

### 🧬 Molecular Visualizations
- **Particle Field Background**: Interactive molecular particle system
- **3D Molecule Viewer**: Three.js-powered molecular models
- **Animated Dissociation**: Molecules breaking apart and reassembling
- **Protein-Ligand Docking**: Simulated binding animations
- **Neurotransmitter Synapse**: SERT transporter and serotonin visualization
- **Receptor Cascades**: G-protein and β-arrestin pathway animations

### 📚 Educational Content
- **Loading Screens**: 20+ rotating educational facts
- **Processing Animations**: Stage-by-stage computation feedback
- **Interactive Diagrams**: Pharmacology visualizations
- **Integrated Documentation**: Complete integration list

### 🔌 Integrations Showcase
- **Current Tools**: RDKit, PubChem, ChEMBL, AutoDock Vina, PySCF
- **AI Models**: OpenAI GPT-4, Google Gemini, Anthropic Claude
- **Planned Features**: CRISPR tools, IBM RXN, bioNEMO, flow chemistry

## Quick Start

### Local Development

```bash
# Navigate to frontend directory
cd services/web-frontend

# Install dependencies
pip install -r requirements.txt

# Run development server
python main.py

# Access at http://localhost:8090
```

### Docker Deployment

```bash
# From project root
docker-compose up web-frontend

# Access at http://localhost:8090
```

## API Endpoints

### Static Content
- `GET /` - Landing page
- `GET /static/*` - CSS, JS, images

### File Upload
- `POST /api/upload/image` - Upload single image
- `POST /api/upload/bulk` - Upload multiple images (max 10)
- `GET /api/images` - List all uploaded images
- `DELETE /api/images/{filename}` - Delete image

### Data APIs
- `GET /api/integrations` - Current and planned integrations
- `GET /health` - Health check endpoint

## Directory Structure

```
web-frontend/
├── main.py                 # FastAPI application
├── Dockerfile              # Container configuration
├── requirements.txt        # Python dependencies
├── templates/
│   └── index.html         # Main landing page
├── static/
│   ├── css/
│   │   ├── main.css       # Design system
│   │   └── animations.css # Animation library
│   ├── js/
│   │   ├── app.js         # Main app logic
│   │   ├── particles.js   # Particle background
│   │   ├── molecular-animations.js  # Three.js animations
│   │   └── loading-screens.js       # Educational slides
│   └── assets/
│       ├── images/        # Uploaded graphics
│       ├── videos/        # Hero background videos
│       └── educational/   # Diagram slides
└── INTEGRATIONS.md        # Integration documentation
```

## Animations

### Molecular Loader
Rotating atoms orbiting a central nucleus with pulsing effects.

### Particle Field
100 particles with:
- Physics-based movement
- Mouse interaction (repulsion)
- Connection lines within 150px
- Fade effects

### Three.js Molecular Animations

#### 1. Dissociation
- Molecule breaks into individual atoms
- Atoms drift apart
- Smooth transition to reassembly

#### 2. Protein-Ligand Docking
- Ligand approaches protein
- Rotation and positioning
- Binding visualization

#### 3. Neurotransmitter Synapse
- Presynaptic/postsynaptic neurons
- Serotonin molecules falling through cleft
- SERT transporter reuptake animation
- Pulsing receptors

#### 4. Receptor Cascade
- Ligand binding to receptor
- Conformational change
- G-protein activation
- Signal propagation

## Customization

### Colors
Edit CSS variables in `static/css/main.css`:

```css
:root {
    --primary-cyan: #00d4ff;
    --primary-purple: #9333ea;
    --primary-blue: #0066ff;
}
```

### Educational Slides
Add to `static/js/loading-screens.js`:

```javascript
this.educationalSlides = [
    "Your educational fact here...",
    // Add more facts
];
```

### Integrations
Update `main.py` `/api/integrations` endpoint.

## Upload Graphics

### Via API

```bash
# Upload single image
curl -X POST http://localhost:8090/api/upload/image \
  -F "file=@my-image.png"

# Upload multiple images
curl -X POST http://localhost:8090/api/upload/bulk \
  -F "files=@image1.png" \
  -F "files=@image2.jpg"
```

### Via Web Interface

1. Navigate to `/admin/upload` (coming soon)
2. Drag and drop images
3. Images appear in `/static/assets/images/`

## Performance

- **Initial Load**: <2s (with CDN caching)
- **Animation FPS**: 60fps (hardware accelerated)
- **Image Optimization**: Auto-compression on upload
- **Caching**: Browser cache + service worker

## Browser Support

- Chrome 90+
- Firefox 88+
- Safari 14+
- Edge 90+

## Dependencies

### Python
- FastAPI 0.109.0
- Uvicorn 0.27.0
- Jinja2 3.1.3
- Python-multipart 0.0.6

### JavaScript (CDN)
- Three.js 0.160.0
- GSAP 3.12.5

## Environment Variables

```env
PORT=8090
API_GATEWAY_URL=http://localhost:8080
```

## Production Deployment

### Optimization Checklist
- [ ] Minify CSS/JS
- [ ] Compress images (WebP format)
- [ ] Enable CDN for static assets
- [ ] Configure nginx reverse proxy
- [ ] Set up SSL/TLS
- [ ] Enable HTTP/2
- [ ] Configure caching headers

### Performance Targets
- Lighthouse Score: 95+
- Core Web Vitals: All green
- Time to Interactive: <3s
- First Contentful Paint: <1.5s

## Development

### Adding New Pages

1. Create template in `templates/`
2. Add route in `main.py`:

```python
@app.get("/new-page", response_class=HTMLResponse)
async def new_page(request: Request):
    return templates.TemplateResponse("new-page.html", {"request": request})
```

### Adding Animations

1. Create CSS keyframes in `animations.css`
2. Add JavaScript logic in appropriate file
3. Trigger via scroll observer or user interaction

## Troubleshooting

### Issue: Animations not working
**Solution**: Check browser console for JavaScript errors. Ensure Three.js and GSAP loaded from CDN.

### Issue: Images not uploading
**Solution**: Check file size (<10MB) and format (PNG, JPG, SVG, WEBP only).

### Issue: Slow loading
**Solution**: Enable CDN caching, compress images, check network tab in dev tools.

## Contributing

1. Fork the repository
2. Create feature branch
3. Make changes
4. Test thoroughly
5. Submit pull request

## License

Proprietary - PharmaSight™ 2026

## Contact

- **GitHub**: https://github.com/justincihi/pharmasight-platform
- **Issues**: https://github.com/justincihi/pharmasight-platform/issues
- **Email**: support@pharmasight.com

---

**Built with ❤️ for accelerating drug discovery**
