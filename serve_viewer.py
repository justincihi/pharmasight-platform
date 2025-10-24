#!/usr/bin/env python3
"""
Simple HTTP server to view molecular visualizations in browser

Usage:
    python serve_viewer.py

Then open: http://localhost:8000
"""

import http.server
import socketserver
import os
import webbrowser
from threading import Timer

# Change to the pharmasight-platform directory
os.chdir('/home/user/pharmasight-platform')

PORT = 8000

class MyHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Add CORS headers
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate')
        super().end_headers()

    def log_message(self, format, *args):
        # Custom logging
        print(f"[{self.log_date_time_string()}] {format % args}")

def open_browser():
    """Open browser after a short delay"""
    webbrowser.open(f'http://localhost:{PORT}/molecular_viewer.html')

if __name__ == '__main__':
    print("=" * 70)
    print("  PharmaSight Molecular Viewer")
    print("=" * 70)
    print(f"\nStarting HTTP server on port {PORT}...")
    print(f"\n🌐 Open in browser: http://localhost:{PORT}/molecular_viewer.html")
    print("\nAvailable files:")
    print("  • molecular_viewer.html - Main viewer interface")
    print("  • molecular_images/     - Generated molecule images")
    print("\nPress Ctrl+C to stop the server")
    print("=" * 70 + "\n")

    # Optionally open browser automatically after 1 second
    # Timer(1.0, open_browser).start()

    with socketserver.TCPServer(("", PORT), MyHTTPRequestHandler) as httpd:
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n\nShutting down server...")
            httpd.shutdown()
