"""
dev_server.py -- preview the website locally.

    python scripts/dev_server.py            # then open http://localhost:8000/

Serves the docs/ folder. The calculator normally runs the Python modules in
the browser (Pyodide). For quick local work, open

    http://localhost:8000/app.html?engine=http

and the page sends its requests to this server instead, which answers them
with the same web_api.handle() running in ordinary Python.
"""

import http.server
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import web_api                                                     # noqa: E402


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kw):
        super().__init__(*args, directory=os.path.join(ROOT, "docs"), **kw)

    def do_POST(self):
        if self.path != "/api":
            self.send_error(404)
            return
        body = self.rfile.read(int(self.headers.get("Content-Length", 0))).decode("utf-8")
        out = web_api.handle(body).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(out)))
        self.end_headers()
        self.wfile.write(out)

    def log_message(self, *args):
        pass


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    print(f"Serving docs/ at http://localhost:{port}/  (calculator with the local engine: "
          f"http://localhost:{port}/app.html?engine=http)")
    http.server.ThreadingHTTPServer(("127.0.0.1", port), Handler).serve_forever()
