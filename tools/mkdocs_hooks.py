# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""MkDocs hook: cache busting for the documentation's own stylesheets and scripts.

Every local file of extra_css and extra_javascript is linked as path?v=h, with h the first 10 hex digits of
the SHA-1 of its content, so that a browser fetches it again whenever it changes (the web server lets browsers
keep these files for hours, and only Material's own files carry a hash in their names).
"""
import hashlib
import os


def _versioned(path, docs_dir):
    """path?v=<SHA-1 prefix of the file> for a file in docs_dir; other entries (URLs) unchanged."""
    name = str(path)
    full = os.path.join(docs_dir, name)
    if "://" in name or "?" in name or not os.path.isfile(full):
        return path
    with open(full, "rb") as f:
        return f"{name}?v={hashlib.sha1(f.read()).hexdigest()[:10]}"


def on_config(config):
    """Replaces the local extra_css and extra_javascript entries by their versioned URLs."""
    config["extra_css"] = [_versioned(p, config["docs_dir"]) for p in config["extra_css"]]
    for i, script in enumerate(config["extra_javascript"]):
        path = getattr(script, "path", script)
        new = _versioned(path, config["docs_dir"])
        if new != path:
            if hasattr(script, "path"):
                script.path = new
            else:
                config["extra_javascript"][i] = new
    return config
