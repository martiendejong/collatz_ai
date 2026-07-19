# R3292: generate self-contained HTML blog package with embedded images.
import os, re, shutil, sys
sys.stdout.reconfigure(encoding="utf-8")

BLOG = r"E:\projects\collatz\research\blog"
OUT = os.path.join(BLOG, "html")
IMG = os.path.join(OUT, "img")
VIZ = r"E:\projects\collatz\research\viz"
os.makedirs(IMG, exist_ok=True)

POSTS = ["01-the-experiment", "02-the-machine", "03-the-fuel-economy",
         "04-the-fair-coin", "05-the-comma-calendar", "06-what-is-new"]
TITLES = {}

# images per post: {post: [(after_heading_fragment, [(file, caption)]), ...]}
IMAGES = {
    "01-the-experiment": [
        ("What came out", [("ca_26623.png",
          "A Collatz orbit as a base-6 cellular automaton: the record-holder 26623. "
          "Black triangles are fuel burning; the noise is a provably fair coin.")]),
    ],
    "02-the-machine": [
        ("The spacetime pictures", [
            ("ca_repunit20.png", "2^20 − 1: a pure fuel tank burning. The giant black triangle "
             "is twenty binary ones being converted; the rest is the fair-coin regime."),
            ("ca_27.png", "The famous 27: no big tank, just a long fever of small flare-ups."),
            ("ca_837799.png", "The champion 837799 (524 steps)."),
        ]),
    ],
    "03-the-fuel-economy": [
        ("The Audit", [("phase_ca_837799.png",
          "The champion in phase time: one row per full climb+fall cycle. Orange marks the "
          "fuel about to burn; red marks the big packets.")]),
        ("The decompression principle", [("verify_repunit20.png",
          "Verification plate: the full automaton (left) and the phase machine (right), "
          "linked row by row. The green bars are whole binary→ternary conversions done in one step.")]),
    ],
    "04-the-fair-coin": [
        ("No base helps", [("bases_repunit20.png",
          "The same orbit in six bases. Structure is visible exactly in the 2-3-smooth bases; "
          "base 7 is provably blind.")]),
    ],
    "05-the-comma-calendar": [
        ("The control experiment: 3n+5", [("phasebases_repunit20.png",
          "The phase machine of 2^20 − 1 rendered in six bases.")]),
    ],
    "06-what-is-new": [
        ("The three candidates that survived the literature check", [("bases_837799.png",
          "The champion orbit in six bases — the kind of picture this project leaves behind.")]),
    ],
}

CSS = """
body{font-family:Georgia,serif;max-width:760px;margin:2em auto;padding:0 1em;line-height:1.65;color:#222}
h1{font-size:1.55em;line-height:1.3} h2{margin-top:1.6em;border-bottom:2px solid #f4632e;padding-bottom:.2em}
a{color:#b33d10} img{max-width:100%;border:1px solid #eee;display:block;margin:1em auto}
.cap{font-size:.85em;color:#666;text-align:center;margin-top:-0.6em;margin-bottom:1.4em}
blockquote{background:#faf6ef;border-left:4px solid #f4632e;padding:.6em 1em;margin:1em 0;font-style:italic}
.nav{display:flex;justify-content:space-between;margin:2.5em 0 1em;font-size:.95em;border-top:1px solid #ddd;padding-top:1em}
.series{font-size:.85em;color:#888;margin-bottom:1.5em}
code{background:#f4f4f4;padding:1px 5px}
"""

def md_to_html(md):
    lines = md.split("\n")
    out = []
    para = []
    def flush():
        if para:
            out.append("<p>" + " ".join(para) + "</p>")
            para.clear()
    in_list = None
    for ln in lines:
        s = ln.rstrip()
        if s.startswith("# "):
            flush(); out.append(f"<h1>{s[2:]}</h1>"); continue
        if s.startswith("## "):
            flush(); out.append(f"<h2>{s[3:]}</h2>"); continue
        if s.startswith("> "):
            flush(); out.append(f"<blockquote>{s[2:]}</blockquote>"); continue
        m = re.match(r"^(\d+)\.\s+(.*)", s)
        if m:
            if in_list != "ol": flush(); out.append("<ol>"); in_list = "ol"
            out.append(f"<li>{m.group(2)}</li>"); continue
        if s.startswith("- "):
            if in_list != "ul": flush(); out.append("<ul>"); in_list = "ul"
            out.append(f"<li>{s[2:]}</li>"); continue
        if in_list and not s:
            out.append(f"</{in_list}>"); in_list = None; continue
        if not s:
            flush(); continue
        para.append(s)
    flush()
    if in_list: out.append(f"</{in_list}>")
    html = "\n".join(out)
    html = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", html)
    html = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", html)
    html = re.sub(r"\[(.+?)\]\((.+?)\.md\)", r'<a href="\2.html">\1</a>', html)
    html = re.sub(r"\[(.+?)\]\((https?://.+?)\)", r'<a href="\2">\1</a>', html)
    html = re.sub(r"`(.+?)`", r"<code>\1</code>", html)
    return html

used = set()
for i, post in enumerate(POSTS):
    md = open(os.path.join(BLOG, post + ".md"), encoding="utf-8").read()
    title = md.split("\n")[0].lstrip("# ").strip()
    TITLES[post] = title
    body = md_to_html(md)
    # inject images after their headings
    for frag, imgs in IMAGES.get(post, []):
        block = ""
        for f, cap in imgs:
            src = os.path.join(VIZ, f)
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(IMG, f))
                used.add(f)
                block += f'<img src="img/{f}" alt="{cap}"><div class="cap">{cap}</div>\n'
        # insert after the first heading containing frag, else append
        pat = re.compile(r"(<h2>[^<]*" + re.escape(frag[:20]) + r"[^<]*</h2>)")
        if pat.search(body):
            body = pat.sub(r"\1\n" + block, body, count=1)
        else:
            body += "\n" + block
    prev = f'<a href="{POSTS[i-1]}.html">&larr; previous</a>' if i > 0 else '<a href="index.html">&larr; series index</a>'
    nxt = f'<a href="{POSTS[i+1]}.html">next &rarr;</a>' if i < len(POSTS) - 1 else '<a href="index.html">series index &rarr;</a>'
    page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>{CSS}</style></head><body>
<div class="series">The Collatz Sessions &mdash; part {i+1} of 6 &middot; <a href="index.html">series index</a></div>
{body}
<div class="nav">{prev}{nxt}</div>
</body></html>"""
    open(os.path.join(OUT, post + ".html"), "w", encoding="utf-8").write(page)
    print(f"{post}.html written")

# index page
items = "\n".join(f'<li><a href="{p}.html">{TITLES[p].split("&mdash;")[-1] if "&mdash;" in TITLES[p] else TITLES[p]}</a></li>' for p in POSTS)
idx = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>The Collatz Sessions</title><style>{CSS}</style></head><body>
<h1>The Collatz Sessions</h1>
<p><i>A six-part series on a months-long human&ndash;AI research collaboration attacking the
Collatz conjecture &mdash; everything verified, everything honestly labeled, everything public:
<a href="https://github.com/martiendejong/collatz_ai">github.com/martiendejong/collatz_ai</a>.</i></p>
<ol>{items}</ol>
<img src="img/ca_repunit20.png" alt="A fuel tank burning"><div class="cap">2^20 &minus; 1: a pure fuel tank burning through the Collatz automaton.</div>
</body></html>"""
open(os.path.join(OUT, "index.html"), "w", encoding="utf-8").write(idx)
print(f"index.html written; images copied: {sorted(used)}")
