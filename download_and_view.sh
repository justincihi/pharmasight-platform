#!/bin/bash
# Script to help you download and view the molecular visualizations

echo "=========================================="
echo "  PharmaSight - Download Instructions"
echo "=========================================="
echo ""
echo "Choose an option:"
echo ""
echo "1. EASIEST: Using VS Code Remote"
echo "   - Right-click 'molecular_viewer.html'"
echo "   - Select 'Download...'"
echo "   - Also download the 'molecular_images/' folder"
echo "   - Open molecular_viewer.html in your browser"
echo ""
echo "2. Using SCP (from your local machine):"
echo "   scp -r user@server:/home/user/pharmasight-platform/molecular_images ./"
echo "   scp user@server:/home/user/pharmasight-platform/molecular_viewer.html ./"
echo "   # Then open molecular_viewer.html locally"
echo ""
echo "3. Using SSH Port Forwarding:"
echo "   ssh -L 8000:localhost:8000 user@server"
echo "   # Then open http://localhost:8000/molecular_viewer.html"
echo ""
echo "4. Create a ZIP file to download:"
echo ""

# Create a zip file with everything needed
cd /home/user/pharmasight-platform
zip -r pharmasight_viewer.zip molecular_viewer.html molecular_images/ >/dev/null 2>&1

if [ $? -eq 0 ]; then
    echo "   ✅ Created: pharmasight_viewer.zip"
    echo "   Download this file and extract it locally"
    echo "   Size: $(du -h pharmasight_viewer.zip | cut -f1)"
    echo ""
    echo "   Location: /home/user/pharmasight-platform/pharmasight_viewer.zip"
else
    echo "   ℹ️  ZIP not available (install with: apt-get install zip)"
fi

echo ""
echo "=========================================="
