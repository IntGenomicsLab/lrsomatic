"""Rewrite repo-relative links that point outside docs/ to GitHub URLs pinned to the built ref."""

import os
import re

REPO = "IntGenomicsLab/lrsomatic"
REF = os.environ.get("DOCS_REF", "dev")
# ../path (relative to docs/) and /path (repo-root absolute, as GitHub renders README links)
LINK = re.compile(r"(!?)\[([^\]]*)\]\((?:\.\./|/)([^)\s#]+)(#[^)\s]*)?\)")
BADGE = re.compile(r"^\s*\[!\[")


def on_page_markdown(markdown, page, config, files):
    docs_dir = config["docs_dir"]

    def fix(m):
        bang, text, path, anchor = m.group(1), m.group(2), m.group(3), m.group(4) or ""
        if path.startswith("docs/") and os.path.exists(os.path.join(docs_dir, path[5:])):
            path = path[5:]
        if os.path.exists(os.path.join(docs_dir, path)):
            return f"{bang}[{text}]({path}{anchor})"
        if bang:
            return f"![{text}](https://raw.githubusercontent.com/{REPO}/{REF}/{path})"
        return f"[{text}](https://github.com/{REPO}/blob/{REF}/{path}{anchor})"

    markdown = LINK.sub(fix, markdown)
    if page.file.src_uri == "index.md":
        # README title, badge wall and the "documentation is here" line are for GitHub readers
        markdown = "\n".join(
            line
            for line in markdown.splitlines()
            if not BADGE.match(line)
            and line.strip() != "# IntGenomicsLab/lrsomatic"
            and "intgenomicslab.github.io/lrsomatic" not in line
        )
    return markdown
