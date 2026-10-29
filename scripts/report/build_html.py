#!/usr/bin/env python
"""Render report/REPORT.md as a self-contained HTML page (report/report.html) for sharing.

Figures referenced as ../figures/*.png are embedded as data: URIs; HTML source comments (<!-- src: ... -->) are kept, so every
number keeps its provenance in the page source. Headings get ids and a generated table of contents.
Run with an environment that has the `markdown` package (conda_envs/annot).
"""
import base64
import html
import os
import re

import markdown

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
src = open(f"{P}/report/REPORT.md").read()

title_md = src.splitlines()[0].lstrip("# ").strip()
body_md = "\n".join(src.splitlines()[1:])

md = markdown.Markdown(extensions=["tables", "toc", "sane_lists", "attr_list"],
                       extension_configs={"toc": {"toc_depth": "2-3", "permalink": False}})
body = md.convert(body_md)
toc = md.toc


def embed(m):
    path = os.path.normpath(os.path.join(f"{P}/report", m.group(1)))
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode()
    return f'src="data:image/png;base64,{b64}"'


body = re.sub(r'src="(\.\./figures/[^"]+\.png)"', embed, body)
# wrap tables so only they scroll horizontally; figures get a figure wrapper
body = re.sub(r"(<table>.*?</table>)", r'<div class="tablewrap">\1</div>', body, flags=re.S)
body = re.sub(r'<p>(<img [^>]*>)</p>', r'<figure>\1</figure>', body)
# mark the human-written disclaimer block
body = body.replace("<p>(Human-written text starts)", '<div class="disclaimer"><p class="dlabel">Human-written text</p><p>', 1)
body = re.sub(r"\(Human-written text ends\)</p>", "</p></div>", body, count=1)

CSS = r"""
/* Layout: one reading column (~68ch) with a sticky contents rail on wide screens; figures may run wider than the text. */
:root{
  --bg:#f7f8f6; --surface:#ffffff; --fg:#1d2321; --muted:#5b6662; --rule:#d8ddd9; --accent:#2f6f62; --accent-soft:#e3efeb;
  --code:#eef1ef;
  --serif:"Source Serif 4","Georgia","Times New Roman",serif;
  --sans:"IBM Plex Sans","Helvetica Neue",Arial,sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,"SFMono-Regular",Menlo,monospace;
}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){
  --bg:#141917; --surface:#1b2220; --fg:#e4e9e6; --muted:#9aa6a1; --rule:#2e3835; --accent:#7cc3b2; --accent-soft:#1f302b; --code:#222b28;
  color-scheme:dark}}
:root[data-theme="dark"]{
  --bg:#141917; --surface:#1b2220; --fg:#e4e9e6; --muted:#9aa6a1; --rule:#2e3835; --accent:#7cc3b2; --accent-soft:#1f302b; --code:#222b28;
  color-scheme:dark}
html{scroll-behavior:smooth}
@media (prefers-reduced-motion: reduce){html{scroll-behavior:auto}}
body{background:var(--bg);color:var(--fg);font-family:var(--serif);font-size:17px;line-height:1.62;padding-inline:16px;padding-block:24px 64px}
.page{max-width:1180px;margin:0 auto;display:grid;grid-template-columns:minmax(0,1fr);gap:40px}
@media (min-width:1080px){.page{grid-template-columns:240px minmax(0,1fr)}}
nav.toc{font-family:var(--sans);font-size:13.5px;line-height:1.45}
@media (min-width:1080px){nav.toc{position:sticky;top:calc(env(safe-area-inset-top,0px) + 16px);align-self:start;max-height:calc(100vh - 32px);overflow:auto}}
nav.toc .label{text-transform:uppercase;letter-spacing:.08em;font-size:11.5px;color:var(--muted);margin:0 0 8px}
nav.toc ul{list-style:none;margin:0;padding:0}
nav.toc ul ul{padding-left:12px;margin:2px 0 6px}
nav.toc a{color:var(--fg);text-decoration:none;display:block;padding:3px 0;border-left:2px solid transparent;padding-left:8px}
nav.toc a:hover,nav.toc a:focus-visible{color:var(--accent);border-left-color:var(--accent);outline:none}
main{min-width:0;max-width:72ch}
header.top{margin-bottom:28px}
.kicker{font-family:var(--sans);text-transform:uppercase;letter-spacing:.09em;font-size:12px;color:var(--accent);margin:0 0 10px}
h1{font-family:var(--sans);font-weight:600;font-size:clamp(26px,3.6vw,36px);line-height:1.18;letter-spacing:-.01em;margin:0;text-wrap:balance}
h2{font-family:var(--sans);font-weight:600;font-size:23px;line-height:1.25;margin:2.2em 0 .6em;padding-top:.6em;border-top:1px solid var(--rule);text-wrap:balance}
h3{font-family:var(--sans);font-weight:600;font-size:18px;line-height:1.3;margin:1.8em 0 .4em;text-wrap:balance}
p,li{margin:0 0 .9em}
a{color:var(--accent)}
strong{font-weight:650}
hr{display:none}
code{font-family:var(--mono);font-size:.86em;background:var(--code);padding:.08em .32em;border-radius:3px}
blockquote{margin:1em 0;padding:.2em 1em;border-left:3px solid var(--accent);color:var(--muted)}
.disclaimer{background:var(--accent-soft);border:1px solid var(--rule);border-radius:6px;padding:14px 18px;font-family:var(--sans);font-size:14.5px;line-height:1.55;margin:0 0 28px}
.disclaimer p{margin:0 0 .7em}
.disclaimer .dlabel{text-transform:uppercase;letter-spacing:.08em;font-size:11.5px;color:var(--accent);margin-bottom:6px}
.tablewrap{overflow-x:auto;margin:1em 0 1.4em;border-top:1px solid var(--rule);border-bottom:1px solid var(--rule)}
table{border-collapse:collapse;font-family:var(--sans);font-size:13.5px;line-height:1.4;font-variant-numeric:tabular-nums;width:100%}
th,td{text-align:left;vertical-align:top;padding:7px 10px;border-bottom:1px solid var(--rule)}
th{font-weight:600;color:var(--muted);font-size:12.5px;text-transform:uppercase;letter-spacing:.04em}
tr:last-child td{border-bottom:none}
figure{margin:1.6em 0;width:min(100%, 1000px)}
@media (min-width:1080px){figure{width:calc(100% + 18vw);max-width:1000px}}
figure img{display:block;width:100%;height:auto;border-radius:4px;background:#ffffff;border:1px solid var(--rule)}
figure + p em, figure + p{font-size:.94em}
main > p:first-of-type{font-size:1.03em}
::selection{background:var(--accent-soft)}
"""

page = f"""<title>Seven Origins of Multicellularity</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400&family=IBM+Plex+Sans:wght@400;600&family=Source+Serif+4:opsz,wght@8..60,400;8..60,650&display=swap">
<style>{CSS}</style>
<div class="page">
<nav class="toc" aria-label="Contents"><p class="label">Contents</p>{toc}</nav>
<main>
<header class="top"><p class="kicker">Agent-run comparative genomics · Ancestra · 187 eukaryotic proteomes</p>
<h1>{html.escape(title_md)}</h1></header>
{body}
</main>
</div>
"""
os.makedirs(f"{P}/report", exist_ok=True)
with open(f"{P}/report/report.html", "w") as f:
    f.write(page)
print(f"wrote report/report.html ({os.path.getsize(f'{P}/report/report.html') / 1e6:.1f} MB)")
