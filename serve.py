#!/usr/bin/env python3
"""Serve the Solenoid Review Bench on this machine and keep each reviewer's labels in a JSON file.

The page reads its review data from DATA_DIR under /data/ and saves every label through PUT /api/labels/<id> to
DATA_DIR/labels/<reviewer>.json, which holds {label id: label}. Each save rewrites the file whole (a temporary file
renamed over the old one), so the file on disk is always complete. The other files in that folder are other
reviewers' labels, which the page shows read-only. Python 3.9 or newer, standard library only.

usage: python3 serve.py [--data DIR] [--reviewer NAME | --labels FILE] [--port N] [--open]
"""
import argparse
import getpass
import json
import mimetypes
import os
import shutil
import tempfile
import threading
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent
MAX_LABEL_BYTES = 64 * 1024
TYPES = {".html": "text/html; charset=utf-8", ".json": "application/json", ".txt": "text/plain; charset=utf-8",
         ".gz": "application/gzip", ".png": "image/png", ".jpg": "image/jpeg"}


class Labels:
    def __init__(self, path):
        self.path = path
        self.lock = threading.Lock()
        self.labels = json.loads(path.read_text()) if path.exists() else {}

    def snapshot(self):
        with self.lock:
            return dict(self.labels)

    def put(self, label_id, label):
        with self.lock:
            self.labels[label_id] = label
            self.path.parent.mkdir(parents=True, exist_ok=True)
            handle, temporary = tempfile.mkstemp(dir=self.path.parent, prefix=f".{self.path.name}.")
            with os.fdopen(handle, "w") as out:
                json.dump(self.labels, out, indent=1, sort_keys=True)
                out.write("\n")
                out.flush()
                os.fsync(out.fileno())
            os.replace(temporary, self.path)


def other_labels(own):
    """{reviewer: labels} for every other labels file beside ``own``."""
    others = {}
    for path in sorted(own.parent.glob("*.json")) if own.parent.is_dir() else []:
        if path != own and not path.name.startswith("."):
            try:
                others[path.stem] = json.loads(path.read_text())
            except ValueError as error:
                print(f"skipping {path}: {error}", flush=True)
    return others


class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        path = unquote(urlsplit(self.path).path)
        if path == "/api/labels":
            labels = self.server.labels
            self.send_json({"file": self.server.labels_shown, "reviewer": labels.path.stem,
                            "labels": labels.snapshot(), "others": other_labels(labels.path)})
        elif path in ("/", "/index.html"):
            self.send_file(ROOT / "index.html")
        elif path.startswith("/data/"):
            # Path parts are checked rather than resolved, so the data folder may hold symbolic links.
            parts = path[len("/data/"):].split("/")
            target = self.server.data_dir.joinpath(*parts)
            if all(part not in ("", ".", "..") for part in parts) and target.is_file():
                self.send_file(target)
            else:
                self.send_error(HTTPStatus.NOT_FOUND)
        else:
            self.send_error(HTTPStatus.NOT_FOUND)

    def do_PUT(self):
        path = unquote(urlsplit(self.path).path)
        length = int(self.headers.get("Content-Length") or 0)
        label_id = path[len("/api/labels/"):] if path.startswith("/api/labels/") else ""
        if not label_id or length > MAX_LABEL_BYTES:
            self.send_error(HTTPStatus.BAD_REQUEST)
            return
        try:
            label = json.loads(self.rfile.read(length))
        except ValueError:
            label = None
        if not isinstance(label, dict):
            self.send_error(HTTPStatus.BAD_REQUEST, "a label is a JSON object")
            return
        self.server.labels.put(label_id, label)
        self.send_json({"saved": label_id})

    def send_json(self, payload):
        body = json.dumps(payload).encode()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def send_file(self, path):
        kind = TYPES.get(path.suffix) or mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(path.stat().st_size))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        with open(path, "rb") as source:
            shutil.copyfileobj(source, self.wfile)

    def log_request(self, code="-", size="-"):
        if isinstance(code, int) and code >= 400:
            super().log_request(code, size)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data", type=Path, default=ROOT / "data",
                        help="review data folder (default: data/ beside this script)")
    parser.add_argument("--reviewer", default=getpass.getuser(),
                        help="saves labels to DATA/labels/<reviewer>.json (default: your user name)")
    parser.add_argument("--labels", type=Path, help="labels file to use instead of DATA/labels/<reviewer>.json")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--open", action="store_true", help="open the bench in the default browser")
    args = parser.parse_args()

    data_dir = args.data.resolve()
    if not (data_dir / "manifest.json").is_file() and not (data_dir / "giants" / "index.json").is_file():
        parser.error(f"{data_dir} holds neither manifest.json nor giants/index.json (see docs/data-format.md)")
    labels_path = (args.labels or data_dir / "labels" / f"{args.reviewer}.json").resolve()
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    server.data_dir = data_dir
    server.labels = Labels(labels_path)
    server.labels_shown = str(labels_path.relative_to(ROOT) if labels_path.is_relative_to(ROOT) else labels_path)
    url = f"http://{args.host}:{args.port}/"
    others = other_labels(labels_path)
    print(f"Solenoid Review Bench at {url}\n  data: {data_dir}\n  labels: {labels_path} "
          f"({len(server.labels.labels)} saved)"
          + "".join(f"\n  {name}: {len(labels)} labels, shown read-only" for name, labels in others.items())
          + "\nCtrl-C stops the server.", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
