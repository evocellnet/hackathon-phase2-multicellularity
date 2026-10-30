#!/usr/bin/env python
"""Render report/REPORT.md (+ supplement) into a self-contained HTML page with embedded figures.
Usage: conda_envs/annot/bin/python scripts/report/render_html.py report/REPORT.md report/report.html [--title T]"""
import base64, re, sys, html, pathlib, markdown

src, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
title = sys.argv[sys.argv.index("--title") + 1] if "--title" in sys.argv else "Multicellularity Origins Report"
text = src.read_text()

def embed(m):
    alt, path = m.group(1), m.group(2)
    p = (src.parent / path)
    if not p.exists():
        return m.group(0)
    b64 = base64.b64encode(p.read_bytes()).decode()
    return f'<figure><img alt="{html.escape(alt)}" src="data:image/png;base64,{b64}"></figure>'

text = re.sub(r"!\[([^\]]*)\]\(([^)]+\.png)\)", embed, text)
body = markdown.markdown(text, extensions=["tables", "fenced_code", "toc", "sane_lists"])
# wrap tables for horizontal scroll
body = body.replace("<table>", '<div class="tw"><table>').replace("</table>", "</table></div>")
toc_items = re.findall(r'<h2 id="([^"]+)">(.*?)</h2>', body)
toc = "".join(f'<li><a href="#{i}">{re.sub("<[^>]+>", "", t)}</a></li>' for i, t in toc_items)

css = """
/* Layout: one reading column (~70ch) with a compact contents list; journal-like serif body, sans labels and captions. */
:root{--bg:#f6f7f4;--surface:#ffffff;--fg:#1d2420;--muted:#5b665f;--rule:#d8ddd6;--accent:#2f6b4f;
--display:"Source Serif 4",Georgia,serif;--body:"Source Serif 4",Georgia,serif;--ui:"IBM Plex Sans",system-ui,sans-serif;--mono:"IBM Plex Mono",ui-monospace,monospace}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#151a17;--surface:#1c2320;--fg:#e4e9e5;--muted:#9eaaa2;--rule:#33403a;--accent:#7cc3a0;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#151a17;--surface:#1c2320;--fg:#e4e9e5;--muted:#9eaaa2;--rule:#33403a;--accent:#7cc3a0;color-scheme:dark}
body{background:var(--bg);color:var(--fg);font-family:var(--body);font-size:1.05rem;line-height:1.62}
.wrap{max-width:46rem;margin:0 auto;padding-inline:16px;padding-block:2.5rem 4rem}
h1{font-family:var(--display);font-weight:600;font-size:clamp(1.6rem,4.2vw,2.25rem);line-height:1.18;text-wrap:balance;margin:.2rem 0 1.2rem}
h2{font-family:var(--ui);font-size:1.02rem;letter-spacing:.04em;text-transform:uppercase;color:var(--accent);border-bottom:1px solid var(--rule);padding-bottom:.35rem;margin:2.6rem 0 .8rem;text-wrap:balance}
h3{font-family:var(--display);font-weight:600;font-size:1.15rem;margin:1.8rem 0 .5rem;text-wrap:balance}
p,li{margin:.55rem 0} a{color:var(--accent)}
code{font-family:var(--mono);font-size:.84em;background:var(--surface);border:1px solid var(--rule);border-radius:3px;padding:0 .25em;overflow-wrap:anywhere}
figure{margin:1.6rem 0;background:#fff;border:1px solid var(--rule);border-radius:6px;padding:.6rem}
figure img{display:block;width:100%;height:auto}
.tw{overflow-x:auto;margin:1rem 0}
table{border-collapse:collapse;font-family:var(--ui);font-size:.86rem;font-variant-numeric:tabular-nums;min-width:100%}
th,td{border-bottom:1px solid var(--rule);padding:.35rem .55rem;text-align:left;vertical-align:top}
th{color:var(--muted);font-weight:600}
blockquote{margin:1rem 0;padding:.6rem 1rem;background:var(--surface);border:1px solid var(--rule);border-radius:6px}
nav.toc{font-family:var(--ui);font-size:.86rem;background:var(--surface);border:1px solid var(--rule);border-radius:6px;padding:.8rem 1.1rem;margin:1.2rem 0 1.8rem}
nav.toc ol{columns:2 16rem;margin:.3rem 0 0;padding-left:1.2rem} nav.toc li{margin:.15rem 0}
.eyebrow{font-family:var(--ui);font-size:.75rem;letter-spacing:.09em;text-transform:uppercase;color:var(--accent);font-weight:600}
"""
head = (f'<title>{html.escape(title)}</title>\n<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono&display=swap">\n<style>{css}</style>\n')
# put contents after the abstract heading's preceding block: insert before first h2
first_h2 = body.find("<h2")
nav = f'<nav class="toc"><strong>Contents</strong><ol>{toc}</ol></nav>'
body = body[:first_h2] + nav + body[first_h2:]
out.write_text(head + '<div class="wrap"><div class="eyebrow">Agentic comparative genomics · evolution of multicellularity</div>\n' + body + "</div>\n")
print(out, out.stat().st_size)
