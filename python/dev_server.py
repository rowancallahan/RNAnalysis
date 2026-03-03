#!/usr/bin/env python3
"""
Development HTTP server that wraps the existing Python handlers.
Serves the frontend from src/ and exposes the analysis API over HTTP.
This replaces Electron's IPC for browser-based development.
"""

import sys
import os
import json
import traceback
from http.server import HTTPServer, SimpleHTTPRequestHandler
from urllib.parse import urlparse, parse_qs
import mimetypes

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from main_handler import state, COMMANDS

# Project root (parent of python/)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(PROJECT_ROOT, 'src')
NODE_MODULES = os.path.join(PROJECT_ROOT, 'node_modules')

PORT = 3000


class DevHandler(SimpleHTTPRequestHandler):
    """HTTP handler that serves static files and API endpoints."""

    def do_POST(self):
        parsed = urlparse(self.path)

        if parsed.path == '/api/command':
            self._handle_api_command()
        elif parsed.path == '/api/upload':
            self._handle_file_upload()
        else:
            self.send_error(404, 'Not found')

    def do_GET(self):
        parsed = urlparse(self.path)

        if parsed.path == '/api/status':
            self._send_json({'ready': True})
            return

        # Serve static files from src/
        self._serve_static(parsed.path)

    def _handle_api_command(self):
        """Route API commands to existing Python handlers."""
        try:
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length)
            request = json.loads(body)

            command = request.get('command')
            params = request.get('params', {})

            if command not in COMMANDS:
                self._send_json({'error': f'Unknown command: {command}'}, 400)
                return

            result = COMMANDS[command](params)
            self._send_json({'data': result})

        except Exception as e:
            traceback.print_exc()
            self._send_json({'error': str(e)}, 500)

    def _handle_file_upload(self):
        """Handle file uploads for counts/metadata loading."""
        try:
            content_type = self.headers.get('Content-Type', '')

            if 'multipart/form-data' in content_type:
                # Parse multipart form data
                import cgi
                form = cgi.FieldStorage(
                    fp=self.rfile,
                    headers=self.headers,
                    environ={
                        'REQUEST_METHOD': 'POST',
                        'CONTENT_TYPE': content_type,
                    }
                )

                file_field = form['file']
                command = form.getvalue('command', 'load_counts')

                # Save uploaded file to temp location
                import tempfile
                suffix = os.path.splitext(file_field.filename)[1] if file_field.filename else '.csv'
                with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                    tmp.write(file_field.file.read())
                    tmp_path = tmp.name

                # Run the load command with the temp file path
                if command in COMMANDS:
                    result = COMMANDS[command]({'file_path': tmp_path})
                    self._send_json({'data': result})
                else:
                    self._send_json({'error': f'Unknown command: {command}'}, 400)

                # Clean up temp file
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass
            else:
                self._send_json({'error': 'Expected multipart/form-data'}, 400)

        except Exception as e:
            traceback.print_exc()
            self._send_json({'error': str(e)}, 500)

    def _serve_static(self, path):
        """Serve static files from src/ directory."""
        if path == '/' or path == '':
            path = '/index.html'

        # Try src/ first
        file_path = os.path.join(SRC_DIR, path.lstrip('/'))

        if not os.path.isfile(file_path):
            # Try node_modules for CDN fallback
            if path.startswith('/node_modules/'):
                file_path = os.path.join(PROJECT_ROOT, path.lstrip('/'))
            else:
                self.send_error(404, f'File not found: {path}')
                return

        if not os.path.isfile(file_path):
            self.send_error(404, f'File not found: {path}')
            return

        # Determine content type
        content_type, _ = mimetypes.guess_type(file_path)
        if content_type is None:
            content_type = 'application/octet-stream'

        try:
            with open(file_path, 'rb') as f:
                content = f.read()

            self.send_response(200)
            self.send_header('Content-Type', content_type)
            self.send_header('Content-Length', len(content))
            self.send_header('Cache-Control', 'no-cache')
            self.send_header('Access-Control-Allow-Origin', '*')
            self.end_headers()
            self.wfile.write(content)
        except IOError:
            self.send_error(500, 'Error reading file')

    def _send_json(self, data, status=200):
        """Send a JSON response."""
        body = json.dumps(data).encode('utf-8')
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', len(body))
        self.send_header('Access-Control-Allow-Origin', '*')
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        """Handle CORS preflight."""
        self.send_response(200)
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.end_headers()

    def log_message(self, format, *args):
        """Custom log format."""
        if '/api/' in (args[0] if args else ''):
            print(f"[API] {args[0]}")
        elif not any(ext in (args[0] if args else '') for ext in ['.js', '.css', '.png', '.ico']):
            print(f"[DEV] {args[0]}")


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else PORT

    server = HTTPServer(('localhost', port), DevHandler)
    print(f"Dev server running at http://localhost:{port}")
    print(f"Serving frontend from: {SRC_DIR}")
    print(f"Python API ready at: http://localhost:{port}/api/command")
    print()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down dev server...")
        server.server_close()


if __name__ == '__main__':
    main()
