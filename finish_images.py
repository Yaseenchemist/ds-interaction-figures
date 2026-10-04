#!/usr/bin/env python3
"""
finish_images.py  --  final touch for the Discovery Studio images (run after ds_batch_interactions.pl)

  2D  *_2D.svg  ->  final/<name>_2D.png   (600 dpi, white background, tight crop)
        * ligand bonds take the 3D carbon colour (same as the 3D picture)
        * clean readable font
        * van der Waals-only residues (pale green circles) removed  (--keep-vdw to keep them)
        * hydrogens removed from the ligand                          (--keep-h to keep them)
  3D  *_3D.png  ->  final/<name>_3D.png   (white margins cropped, 600 dpi tag)

Usage:   python3 ~/Documents/ds_automation/finish_images.py
Options: --carbon 0,190,190   --keep-vdw   --keep-h   --width-in 8   --dpi 600   --folder PATH
"""
import argparse, glob, os, re, math
import cairosvg
from PIL import Image, ImageChops

ap = argparse.ArgumentParser()
ap.add_argument("--folder", default=os.path.expanduser("~/Documents/ds_automation/out"))
ap.add_argument("--carbon", default="0,190,190", help="ligand carbon colour R,G,B (must equal the 3D one)")
ap.add_argument("--keep-vdw", action="store_true")
ap.add_argument("--keep-h", action="store_true")
ap.add_argument("--inner-width", type=float, default=0.12, help="width of the inner aromatic-ring bonds (same stick look as outer bonds)")
ap.add_argument("--no-legend", action="store_true")
ap.add_argument("--no-dist", action="store_true")
ap.add_argument("--keep-labels", action="store_true", help="keep element letters on the ligand")
ap.add_argument("--width-in", type=float, default=8.0)
ap.add_argument("--dpi", type=int, default=600)
ap.add_argument("--font", default="DejaVu Sans")
ap.add_argument("--upscale-3d", type=float, default=1.0, help="enlarge the 3D picture (Lanczos) after cropping (1 = none, the DS picture is already large)")
ap.add_argument("--dec", type=int, default=1, help="decimals of the distances (1 or 2)")
ap.add_argument("--res-font", type=float, default=0.0165, help="residue label size as a fraction of the picture width")
ap.add_argument("--dist-font", type=float, default=0.0140, help="distance label size as a fraction of the picture width")
ap.add_argument("--line-w", type=float, default=0.0022, help="interaction line width as a fraction of the picture width")
ap.add_argument("--ds-lines", action="store_true", help="DS drew its own interaction lines (SHOW_DS_LINES=1) -> do not redraw them")
ap.add_argument("--res-fmt", default="3+num", help="3+num -> 'ASP 1046'; cap -> 'Asp1046'; raw -> 'ASP1046'")
ap.add_argument("--only", default="", help="process only this name (e.g. lig1)")
ap.add_argument("--metal-cut", type=float, default=2.45, help="metal coordination cut-off (A) for N/O donors (S/Cl use +0.35)")
ap.add_argument("--no-metal", action="store_true", help="do not add metal (Zn ...) coordination contacts")
a = ap.parse_args()

R, G, B = [int(v) for v in a.carbon.split(",")]
OUT = os.path.join(a.folder, "final")
os.makedirs(OUT, exist_ok=True)


def fnum(tag, attr):
    m = re.search(r'\b%s="([-\d.eE]+)"' % attr, tag)
    return float(m.group(1)) if m else None


def autocrop(im, margin_frac=0.04, bg=(255, 255, 255)):
    im = im.convert("RGB")
    diff = ImageChops.difference(im, Image.new("RGB", im.size, bg))
    diff = diff.point(lambda p: 255 if p > 12 else 0)
    box = diff.getbbox()
    if not box:
        return im
    x0, y0, x1, y1 = box
    m = int(max(x1 - x0, y1 - y0) * margin_frac)
    return im.crop((max(0, x0 - m), max(0, y0 - m), min(im.width, x1 + m), min(im.height, y1 + m)))



LEGEND = {   # DS 2D interaction colour -> name   (edit the names here if your DS version differs)
    "rgb(0,255,0)": "Conventional Hydrogen Bond",
    "rgb(220,255,220)": "Carbon Hydrogen Bond",
    "rgb(255,200,255)": "Alkyl / Pi-Alkyl",
    "rgb(255,100,200)": "Pi-Pi T-shaped / Stacked",
    "rgb(200,100,255)": "Pi-Sigma",
    "rgb(255,150,0)": "Pi-Cation / Pi-Anion",
    "rgb(255,200,0)": "Pi-Sulfur",
    "rgb(0,255,255)": "Halogen",
    "rgb(170,255,0)": "Pi-Lone Pair",
    "rgb(130,90,180)": "Metal-Acceptor",
}
DS_TYPE_COL = [   # interaction type text (lower case, contains) -> colour DS uses in its 2D diagram
    ("metal", "rgb(130,90,180)"), ("carbon hydrogen", "rgb(220,255,220)"), ("pi-donor", "rgb(220,255,220)"),
    ("hydrogen", "rgb(0,255,0)"), ("halogen", "rgb(0,255,255)"), ("lone pair", "rgb(170,255,0)"),
    ("pi-pi", "rgb(255,100,200)"), ("pi pi", "rgb(255,100,200)"), ("sigma", "rgb(200,100,255)"),
    ("alkyl", "rgb(255,200,255)"), ("cation", "rgb(255,150,0)"), ("anion", "rgb(255,150,0)"), ("charge", "rgb(255,150,0)"),
    ("salt", "rgb(255,150,0)"), ("sulfur", "rgb(255,200,0)"), ("sulphur", "rgb(255,200,0)"),
]
MCOL = "rgb(130,90,180)"            # colour of metal-coordination lines (2D legend / 2D lines)
LEGEND_KEYS = {   # words searched in the interaction type names exported by DS
    "rgb(0,255,0)": ["hydrogen"], "rgb(220,255,220)": ["carbon hydrogen", "hydrogen"],
    "rgb(255,200,255)": ["alkyl"], "rgb(255,100,200)": ["pi-pi", "pi pi", "stack", "t-shaped"],
    "rgb(200,100,255)": ["sigma"], "rgb(255,150,0)": ["cation", "anion", "charge"], "rgb(255,200,0)": ["sulfur", "sulphur"],
}


def legend_names(name):
    """colour -> readable interaction names that really occur in this complex (from the DS table + metal contacts)"""
    out = {}
    types = [r.get("type", "") for r in read_tsv(os.path.join(a.folder, name + "_nonbonds.tsv"))]
    if not a.no_metal and metal_contacts(name)[0]:
        types.append("Metal-Acceptor")
    for t in types:
        tl = t.lower()
        for k, col in DS_TYPE_COL:
            if k in tl:
                nm = re.sub(r"\s*\(.*?\)", "", t).strip()
                if nm and nm not in out.setdefault(col, []):
                    out[col].append(nm)
                break
    return out


def add_legend(src, name=""):
    used = []
    for m in re.finditer(r'<line\b[^>]*stroke:(rgb\([0-9,]+\)); stroke-width:0.05[^>]*/>', src):
        if m.group(1) in LEGEND and m.group(1) not in used:
            used.append(m.group(1))
    if not used:
        return src
    names = legend_names(name) if name else {}
    def lname(col):
        n = names.get(col)
        txt = " / ".join(n[:3]) if n else LEGEND[col]
        return txt if len(txt) <= 46 else txt[:44] + "..."
    vb = re.search(r'viewBox="([-\d.eE]+)\s+([-\d.eE]+)\s+([-\d.eE]+)\s+([-\d.eE]+)"', src)
    x0, y0, w, h = [float(v) for v in vb.groups()]
    row = 0.62
    rows_n = (len(used) + 1) // 2 if len(used) > 3 else len(used)
    colw = 9.6
    items = ['<text x="%.3f" y="%.3f" font-family="%s" font-size="0.52" font-weight="bold" fill="rgb(0,0,0)">Interactions</text>' % (x0 + 0.3, y0 + h + 0.75, a.font)]
    for i, col in enumerate(used):
        cidx, ridx = (i // rows_n, i % rows_n)
        y = y0 + h + 1.05 + ridx * row
        xx = x0 + 0.3 + cidx * colw
        items.append('<rect x="%.3f" y="%.3f" width="0.42" height="0.42" fill="%s" stroke="rgb(90,90,90)" stroke-width="0.02"/>' % (xx, y, col))
        items.append('<text x="%.3f" y="%.3f" font-family="%s" font-size="0.46" fill="rgb(0,0,0)">%s</text>' % (xx + 0.65, y + 0.37, a.font, lname(col)))
    newh = h + 1.3 + rows_n * row
    src = src.replace(vb.group(0), 'viewBox="%s %s %s %.4f"' % (vb.group(1), vb.group(2), vb.group(3), newh), 1)
    return src.replace("</svg>", "".join(items) + "</svg>")


def add_distances(src, name):
    """put the 3D distance (A) on every coloured dashed line of the 2D diagram (needs <name>_nonbonds.tsv written by DS)"""
    tsv = os.path.join(a.folder, name + "_nonbonds.tsv")
    if not os.path.isfile(tsv):
        return src
    rows = []
    pdbfile = os.path.join(os.path.dirname(os.path.abspath(a.folder)), "pdb", name + ".pdb")
    atoms_, lig_keys = [], set()
    if os.path.isfile(pdbfile):
        for ln in open(pdbfile, errors="ignore"):
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    rn, num = ln[17:20].strip().upper(), ln[22:26].strip()
                    xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                except ValueError:
                    continue
                het = ln.startswith("HETATM") and rn not in ("HOH", "WAT", "ZN", "MG", "CA", "NA", "K", "CL", "FE", "MN", "CU", "CO", "NI", "CD", "HG")
                atoms_.append((rn + num, xyz, het))
    def res_at(xyz):
        """nearest non-ligand residue to a point, or None when the point belongs to the ligand"""
        if not atoms_:
            return None
        best = min(atoms_, key=lambda t: (t[1][0] - xyz[0]) ** 2 + (t[1][1] - xyz[1]) ** 2 + (t[1][2] - xyz[2]) ** 2)
        d = math.dist(best[1], xyz)
        return None if best[2] or d > 3.3 else best[0]
    lines_ = open(tsv, encoding="utf-8", errors="ignore").read().splitlines()
    head = lines_[0].split("\t") if lines_ else []
    for ln in lines_[1:]:
        p = ln.split("\t")
        try:
            dist = float(p[1])
        except (ValueError, IndexError):
            continue
        keys = []
        if "from_res" in head and len(p) >= 6:
            for rtxt, xtxt in ((p[2], p[3]), (p[4], p[5])):
                names = [r.strip().upper() for r in rtxt.split(";") if r.strip() and not r.strip().startswith("*")]
                names = [r for r in names if not re.fullmatch(r"[A-Z0-9]{3}\d{4}", r) or True]
                if names and atoms_:
                    nonlig = [r for r in names if any(t[0] == r and not t[2] for t in atoms_)]
                    names = nonlig
                if names:
                    keys += names
                elif xtxt.count(",") == 2:
                    try:
                        r = res_at(tuple(float(v) for v in xtxt.split(",")))
                    except ValueError:
                        r = None
                    if r:
                        keys.append(r)
        elif len(p) >= 3:
            keys = [r.strip().upper() for r in p[2].split(",") if r.strip()]
        if keys:
            rows.append((p[0].lower(), dist, sorted(set(keys))))
    if not rows:
        return src
    # residue circles (coloured ones) and their names
    circ = [(fnum(m.group(0), "cx"), fnum(m.group(0), "cy"), re.search(r'url\(#(\w+)\)', m.group(0)).group(1))
            for m in re.finditer(r"<circle\b[^>]*fill=\"url\(#\w+\)\"[^>]*/>", src)]
    txt = [(fnum(m.group(0), "x"), fnum(m.group(0), "y"), m.group(1).upper() + m.group(2))
           for m in re.finditer(r'<text\b[^>]*>\s*([A-Za-z0-9]{1,4})\s*</text>\s*<text[^>]*>\s*[A-Za-z0-9]*:(-?\d+[A-Za-z]?)\s*</text>', src)]
    named = []
    for cx, cy, g in circ:
        if txt:
            t = min(txt, key=lambda t: abs(t[0] - cx) + abs(t[1] - cy))
            if abs(t[0] - cx) + abs(t[1] - cy) < 1.0:
                named.append((cx, cy, t[2]))
    lines = []
    for m in re.finditer(r'<line\b[^>]*stroke:(rgb\([0-9,]+\)); stroke-width:0.05[^>]*/>', src):
        ln = m.group(0)
        if m.group(1) == MCOL:
            continue
        x1, y1, x2, y2 = fnum(ln, "x1"), fnum(ln, "y1"), fnum(ln, "x2"), fnum(ln, "y2")
        best = None
        for cx, cy, rn in named:
            for (ex, ey) in ((x1, y1), (x2, y2)):
                d = math.hypot(ex - cx, ey - cy)
                if best is None or d < best[0]:
                    best = (d, rn)
        if best and best[0] < 1.2:
            lines.append((m.group(1), best[1], x1, y1, x2, y2))
    used = {}
    out = []
    for col, rn, x1, y1, x2, y2 in sorted(lines, key=lambda l: math.hypot(l[4] - l[2], l[5] - l[3])):
        cand = [r for r in rows if rn in r[2]]
        keys = LEGEND_KEYS.get(col, [])
        cand_t = [r for r in cand if any(k in r[0] for k in keys)] or cand
        if not cand_t:
            continue
        cand_t = sorted(cand_t, key=lambda r: r[1])
        k = used.get((rn, col), 0)
        used[(rn, col)] = k + 1
        dist = cand_t[min(k, len(cand_t) - 1)][1]
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        base = ('<text x="%.4f" y="%.4f" font-family="%s" font-size="0.44" font-weight="bold" text-anchor="middle" ' % (mx, my + 0.15, a.font))
        out.append(base + 'fill="rgb(255,255,255)" stroke="rgb(255,255,255)" stroke-width="0.14">%.*f</text>' % (a.dec, dist))   # white halo
        out.append(base + 'fill="rgb(200,0,0)">%.*f</text>' % (a.dec, dist))
    return src.replace("</svg>", "".join(out) + "</svg>")


def kekule_inner(src):
    """DS draws the aromatic ring as a dashed inner ring (one dashed line per bond).
    Keep every OTHER inner line (alternating single / double bonds, Kekule style) and draw it as a solid stick."""
    inner = [m for m in re.finditer(r"<line\b[^>]*/>", src)
             if "oid=" not in m.group(0) and "stroke-dasharray" in m.group(0) and "rgb(146,146,146)" in m.group(0)]
    if not inner:
        return src
    P = []
    for m in inner:
        t = m.group(0)
        P.append(((fnum(t, "x1"), fnum(t, "y1")), (fnum(t, "x2"), fnum(t, "y2"))))
    n = len(P)
    def near(i, j):
        return min(math.hypot(p[0] - q[0], p[1] - q[1]) for p in P[i] for q in P[j]) < 0.2
    adj = {i: [j for j in range(n) if j != i and near(i, j)] for i in range(n)}
    seen, keep = set(), set()
    for st in range(n):
        if st in seen:
            continue
        order, stack = [], [st]
        # walk the ring: follow neighbours, preferring unseen
        cur, prev = st, None
        while cur is not None and cur not in seen:
            seen.add(cur); order.append(cur)
            nxt = [j for j in adj[cur] if j not in seen]
            cur = nxt[0] if nxt else None
        for k, idx in enumerate(order):
            if k % 2 == 0 and not (len(order) % 2 == 1 and k == len(order) - 1):
                keep.add(idx)
    c = (R, G, B)
    out, last = [], 0
    for i, m in enumerate(inner):
        out.append(src[last:m.start()])
        if i in keep:
            t = re.sub(r'\s*stroke-dasharray="[^"]*"|\s*stroke-dashoffset="[^"]*"', "", m.group(0))
            t = t.replace("stroke:rgb(146,146,146)", "stroke:rgb(%d,%d,%d)" % c).replace(
                "stroke-width:0.025", "stroke-width:%.3f; stroke-linecap:round" % a.inner_width)
            out.append(t)
        last = m.end()
    out.append(src[last:])
    return "".join(out)



# =====================================================================================================
#  METAL COORDINATION (Zn, Mg ...): DS does not list it as a non-bond, so it is computed here from the PDB
# =====================================================================================================
METALS = {"ZN", "MG", "MN", "FE", "FE2", "CU", "CU1", "CO", "NI", "CD", "HG", "CA", "MO"}
WATER = {"HOH", "WAT", "H2O", "TIP", "SOL", "DOD"}
COV = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07, "F": 0.57, "CL": 1.02, "BR": 1.20, "I": 1.39, "H": 0.31}
_PDB_CACHE = {}

def load_pdb(name):
    if name in _PDB_CACHE:
        return _PDB_CACHE[name]
    f = os.path.join(os.path.dirname(os.path.abspath(a.folder)), "pdb", name + ".pdb")
    at = []
    if os.path.isfile(f):
        for ln in open(f, errors="ignore"):
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            try:
                xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                continue
            nm = ln[12:16].strip().upper()
            rn = ln[17:20].strip().upper()
            el = (ln[76:78].strip() or re.sub(r"[^A-Z]", "", nm)[:1]).upper()
            if ln.startswith("ATOM") and el in ("CA",):
                el = "C"
            at.append({"het": ln.startswith("HETATM"), "name": nm, "rn": rn, "res": rn + ln[22:26].strip(), "el": el, "xyz": xyz,
                       "alt": ln[16:17].strip()})
    _PDB_CACHE[name] = at
    return at


def shown_residues(name):
    return {r["text"].strip().upper() for r in read_tsv(os.path.join(a.folder, name + "_labels.tsv")) if r.get("kind") == "RES"}


def wl_hashes(elems, adj, rounds=3):
    h = {i: elems[i] for i in elems}
    for _ in range(rounds):
        h = {i: elems[i] + "(" + ",".join(sorted(h[j] for j in adj[i])) + ")" for i in elems}
    return h


def metal_contacts(name):
    """-> (rows, ligand_graph).  rows: dict(type, distance, from_xyz (metal), to_xyz (donor), metal, donor_res, donor_el,
    in_lig, wl (graph hash of the donor inside the ligand))"""
    if a.no_metal:
        return [], None
    atoms = [x for x in load_pdb(name) if x["alt"] in ("", "A")]
    metals = [x for x in atoms if x["het"] and x["rn"] in METALS]
    if not metals:
        return [], None
    lig = [x for x in atoms if x["het"] and x["rn"] not in METALS and x["rn"] not in WATER and x["el"] != "H"]
    ligres = {x["res"] for x in lig}
    shown = shown_residues(name)
    # ligand covalent graph
    idx = {id(x): i for i, x in enumerate(lig)}
    adj = {i: [] for i in range(len(lig))}
    for i, x in enumerate(lig):
        for j in range(i + 1, len(lig)):
            y = lig[j]
            if x["res"] != y["res"]:
                continue
            if math.dist(x["xyz"], y["xyz"]) < COV.get(x["el"], .8) + COV.get(y["el"], .8) + 0.4:
                adj[i].append(j); adj[j].append(i)
    elems = {i: x["el"] for i, x in enumerate(lig)}
    wl = wl_hashes(elems, adj) if lig else {}
    rows = []
    for m in metals:
        for x in atoms:
            if x is m or x["el"] not in ("N", "O", "S", "CL") or (x["het"] and x["rn"] in METALS):
                continue
            cut = a.metal_cut + (0.35 if x["el"] in ("S", "CL") else 0)
            d = math.dist(m["xyz"], x["xyz"])
            if d > cut:
                continue
            in_lig = x["res"] in ligres
            if not in_lig and x["res"] not in shown:
                continue
            rows.append({"type": "Metal-Acceptor", "distance": "%.2f" % d, "from_xyz": "%.3f,%.3f,%.3f" % m["xyz"],
                         "to_xyz": "%.3f,%.3f,%.3f" % x["xyz"], "metal": m["res"], "metal_name": m["rn"], "donor_res": x["res"],
                         "donor_el": x["el"], "in_lig": in_lig, "wl": wl.get(idx.get(id(x)), "") if in_lig else "",
                         "metal_xyz": m["xyz"]})
    return rows, lig


def metal_2d(src, name):
    rows, _ = metal_contacts(name)
    if not rows:
        return src
    vb = re.search(r'viewBox="([-\d.eE]+)\s+([-\d.eE]+)\s+([-\d.eE]+)\s+([-\d.eE]+)"', src)
    vx, vy, vw, vh = [float(v) for v in vb.groups()]
    ELE = {(R, G, B): "C", (38, 64, 242): "N", (230, 25, 25): "O", (242, 204, 25): "S", (25, 191, 25): "CL", (240, 240, 240): "H"}
    # --- ligand atoms/bonds from the (already restyled) ligand half-bonds
    halves = []
    for m in re.finditer(r"<line\b[^>]*/>", src):
        ln = m.group(0)
        if "oid=" not in ln:
            continue
        col = re.search(r"stroke:rgb\((\d+),(\d+),(\d+)\)", ln)
        if not col:
            continue
        c = tuple(int(v) for v in col.groups())
        w = re.search(r"stroke-width:([0-9.]+)", ln)
        w = float(w.group(1)) if w else 0
        x1, y1, x2, y2 = fnum(ln, "x1"), fnum(ln, "y1"), fnum(ln, "x2"), fnum(ln, "y2")
        if None in (x1, y1, x2, y2):
            continue
        if c in ELE and w > 0.1:
            el, ax, ay = ELE[c], x2, y2
        elif w < 0.05 and c == (218, 187, 255):                 # DS colour of Zn (metal half-bonds)
            L = math.hypot(x2 - x1, y2 - y1) or 1
            el, ax, ay = "M", x1 + (x2 - x1) / L * 0.5, y1 + (y2 - y1) / L * 0.5
        else:
            continue
        halves.append((round(x1, 2), round(y1, 2), el, ax, ay, ln))
    nodes = []          # [x, y, el, [count]]
    def node_of(x, y, el):
        for n in nodes:
            if math.hypot(n[0] - x, n[1] - y) < 0.3 and (n[2] == el or "M" in (n[2], el)):
                return nodes.index(n)
        nodes.append([x, y, el]); return len(nodes) - 1
    bonds = {}
    for h in halves:
        bonds.setdefault((h[0], h[1]), []).append(h)
    edges = set()
    for k, hs in bonds.items():
        ids = [node_of(h[3], h[4], h[2]) for h in hs]
        if len(ids) == 2 and ids[0] != ids[1]:
            edges.add(tuple(sorted(ids)))
    mnodes = [i for i, n in enumerate(nodes) if n[2] == "M"]
    heavy = [i for i, n in enumerate(nodes) if n[2] not in ("M", "H")]
    adj = {i: [] for i in heavy}
    for p, q in edges:
        if p in adj and q in adj:
            adj[p].append(q); adj[q].append(p)
    elems = {i: nodes[i][2] for i in heavy}
    wl2 = wl_hashes(elems, adj) if heavy else {}
    # --- residue circles (names) for protein donors
    circ = [(fnum(m.group(0), "cx"), fnum(m.group(0), "cy")) for m in re.finditer(r"<circle\b[^>]*fill=\"url\(#\w+\)\"[^>]*/>", src)]
    txt = [(fnum(m.group(0), "x"), fnum(m.group(0), "y"), m.group(1).upper() + m.group(2))
           for m in re.finditer(r'<text\b[^>]*>\s*([A-Za-z0-9]{1,4})\s*</text>\s*<text[^>]*>\s*[A-Za-z0-9]*:(-?\d+[A-Za-z]?)\s*</text>', src)]
    named = {}
    for cx, cy in circ:
        if txt:
            t = min(txt, key=lambda t: abs(t[0] - cx) + abs(t[1] - cy))
            if abs(t[0] - cx) + abs(t[1] - cy) < 1.0:
                named[t[2]] = (cx, cy)
    out = []
    removed_m = False
    used = set()
    by_metal = {}
    for r in rows:
        by_metal.setdefault(r["metal"], []).append(r)
    lig_cx = sum(nodes[i][0] for i in heavy) / len(heavy) if heavy else vx + vw / 2
    lig_cy = sum(nodes[i][1] for i in heavy) / len(heavy) if heavy else vy + vh / 2
    for mres, rs in by_metal.items():
        # donors in the 2D ligand
        donors = []
        for r in rs:
            if not r["in_lig"]:
                continue
            cand = [i for i in heavy if i not in used and nodes[i][2] == r["donor_el"]]
            exact = [i for i in cand if wl2.get(i) == r["wl"]]
            near_m = [i for i in (exact or cand) if any(((i, mm) in edges or (mm, i) in edges) for mm in mnodes)]
            pick = (near_m or exact or cand or [None])[0]
            if pick is not None:
                used.add(pick)
            donors.append((r, pick))
        # position of the metal circle
        zpos = None
        new_circle = False
        mt = [(fnum(m.group(0), "x"), fnum(m.group(0), "y")) for m in re.finditer(r"<text\b[^>]*oid=[^>]*>\s*%s\s*</text>" % re.escape(rs[0]["metal_name"]), src, re.I)]
        if mt and len(by_metal) == 1:
            zpos = (mt[0][0], mt[0][1] + 0.05)
        elif mnodes and len(by_metal) == 1:
            zpos = (nodes[mnodes[0]][0], nodes[mnodes[0]][1])
        metal_in_2d = zpos is not None
        if zpos is None:
            pts = [(nodes[p][0], nodes[p][1]) for r, p in donors if p is not None]
            if pts:
                dx = sum(p[0] for p in pts) / len(pts); dy = sum(p[1] for p in pts) / len(pts)
                ux, uy = dx - lig_cx, dy - lig_cy
                n_ = math.hypot(ux, uy) or 1.0
                base = math.atan2(uy / n_, ux / n_)
                others = [v for v in named.values()] + [(nodes[i][0], nodes[i][1]) for i in heavy]
                best = None
                for da in (0, 25, -25, 50, -50, 80, -80, 110, -110, 150, -150, 180):
                    ang = base + math.radians(da)
                    cand = (dx + math.cos(ang) * 1.55, dy + math.sin(ang) * 1.55)
                    free = min(math.hypot(cand[0] - o[0], cand[1] - o[1]) for o in others) if others else 9
                    if best is None or free > best[0]:
                        best = (free, cand)
                    if free > 1.25:
                        break
                zpos = best[1]
                new_circle = True
        if zpos is None:
            continue
        # drop DS's own pale Zn bonds (re-drawn below as dashed coordination lines)
        for h in halves:
            if metal_in_2d and h[2] == "M" and h[5] in src:
                src = src.replace(h[5], "", 1); removed_m = True
        if new_circle:
            nm = rs[0]["metal_name"].capitalize()
            out.append('<circle cx="%.4f" cy="%.4f" r="0.46" fill="rgb(218,187,255)" stroke="rgb(60,0,100)" stroke-width="0.02"/>' % zpos)
            out.append('<text x="%.4f" y="%.4f" font-family="%s" font-size="0.30" font-weight="bold" text-anchor="middle" fill="rgb(0,0,0)">%s</text>' % (zpos[0], zpos[1] - 0.02, a.font, nm.upper()))
            out.append('<text x="%.4f" y="%.4f" font-family="%s" font-size="0.24" text-anchor="middle" fill="rgb(0,0,0)">A:%s</text>' % (zpos[0], zpos[1] + 0.25, a.font, re.sub(r"\D", "", mres)))
        targets = []
        for r, p in donors:
            if p is not None:
                targets.append((r, (nodes[p][0], nodes[p][1]), 0.22))
        for r in rs:
            if not r["in_lig"]:
                key = r["donor_res"]
                if key in named:
                    targets.append((r, named[key], 0.65))
        for r, (tx, ty), trim in targets:
            L = math.hypot(tx - zpos[0], ty - zpos[1]) or 1
            ux, uy = (tx - zpos[0]) / L, (ty - zpos[1]) / L
            x0, y0 = zpos[0] + ux * 0.42, zpos[1] + uy * 0.42
            x1_, y1_ = tx - ux * trim, ty - uy * trim
            out.append('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" style="stroke:%s; stroke-width:0.05;" stroke-dasharray="0.16,0.10" stroke-linecap="round"/>' % (x0, y0, x1_, y1_, MCOL))
            Lseg = math.hypot(x1_ - x0, y1_ - y0) or 1
            px_, py_ = -(y1_ - y0) / Lseg, (x1_ - x0) / Lseg
            if (zpos[0] + px_ * 0.4 - lig_cx) ** 2 + (zpos[1] + py_ * 0.4 - lig_cy) ** 2 < (zpos[0] - px_ * 0.4 - lig_cx) ** 2 + (zpos[1] - py_ * 0.4 - lig_cy) ** 2:
                px_, py_ = -px_, -py_                      # put the number on the side away from the ligand
            off = 0.34 if Lseg < 2.0 else 0.0
            mx, my = (x0 + x1_) / 2 + px_ * off, (y0 + y1_) / 2 + py_ * off
            base = ('<text x="%.4f" y="%.4f" font-family="%s" font-size="0.44" font-weight="bold" text-anchor="middle" ' % (mx, my + 0.15, a.font))
            out.append(base + 'fill="rgb(255,255,255)" stroke="rgb(255,255,255)" stroke-width="0.14">%.*f</text>' % (a.dec, float(r["distance"])))
            out.append(base + 'fill="rgb(90,40,150)">%.*f</text>' % (a.dec, float(r["distance"])))
        # enlarge the drawing area if the metal circle sits outside
        x_lo, x_hi, y_lo, y_hi = vx, vx + vw, vy, vy + vh
        x_lo, x_hi = min(x_lo, zpos[0] - 1.0), max(x_hi, zpos[0] + 1.0)
        y_lo, y_hi = min(y_lo, zpos[1] - 1.0), max(y_hi, zpos[1] + 1.0)
        vx, vy, vw, vh = x_lo, y_lo, x_hi - x_lo, y_hi - y_lo
    src = src.replace(vb.group(0), 'viewBox="%.5f %.5f %.5f %.5f"' % (vx, vy, vw, vh), 1)
    return src.replace("</svg>", "".join(out) + "</svg>")


def clean_svg(src, name=''):
    src = src.replace("Luxi Sans", a.font)

    # --- remove van der Waals-only residues (circles filled with rad_grad0) ---
    if not a.keep_vdw:
        texts = [(m.group(0), fnum(m.group(0), "x"), fnum(m.group(0), "y"))
                 for m in re.finditer(r"<text\b[^>]*>[^<]*</text>", src)]
        for m in list(re.finditer(r"<circle\b[^>]*/>", src)):
            tag = m.group(0)
            if "rad_grad0" not in tag:
                continue
            cx, cy = fnum(tag, "cx"), fnum(tag, "cy")
            src = src.replace(tag, "", 1)
            for t, x, y in texts:
                if x is not None and abs(x - cx) < 0.06 and abs(y - cy) < 0.7:
                    src = src.replace(t, "", 1)

    # --- remove hydrogens (labels + C-H / N-H bond halves) ---
    if not a.keep_h:
        lines = re.findall(r"<line\b[^>]*/>", src)
        for ln in lines:
            if "rgb(255,255,255)" in ln:                     # white half-bond = goes to an H atom
                x1, y1 = fnum(ln, "x1"), fnum(ln, "y1")
                oid = re.search(r'oid="([^"]*)"', ln)
                src = src.replace(ln, "", 1)
                for other in lines:
                    if other is ln or other not in src:
                        continue
                    if oid and oid.group(0) in other:
                        ox, oy = fnum(other, "x1"), fnum(other, "y1")
                        if abs(ox - x1) < 1e-3 and abs(oy - y1) < 1e-3:
                            src = src.replace(other, "", 1)
        src = re.sub(r"<text\b[^>]*>\s*H\s*</text>\s*", "", src)

    # --- ligand: ball & stick, same colours as the 3D picture ---
    PAL = {"rgb(146,146,146)": (R, G, B), "rgb(143,143,255)": (38, 64, 242), "rgb(240,0,0)": (230, 25, 25),
           "rgb(255,200,50)": (242, 204, 25), "rgb(32,240,32)": (25, 191, 25), "rgb(255,255,255)": (240, 240, 240),
           "rgb(255,200,0)": (242, 204, 25), "rgb(255,150,0)": (255, 128, 0)}
    src = kekule_inner(src)
    balls = []
    def fix_line(m):
        ln = m.group(0)
        if "oid=" not in ln or "stroke-width:0.025" not in ln:
            return ln
        col = re.search(r"stroke:(rgb\([0-9,]+\))", ln)
        if not col or col.group(1) not in PAL:
            return ln
        x1, y1, x2, y2 = fnum(ln, "x1"), fnum(ln, "y1"), fnum(ln, "x2"), fnum(ln, "y2")
        L = math.hypot(x2 - x1, y2 - y1) or 1
        ax, ay = x1 + (x2 - x1) / L * 0.5, y1 + (y2 - y1) / L * 0.5      # real atom position (bond length 1.0)
        c = PAL[col.group(1)]
        balls.append((ax, ay, c))
        new = ln.replace(col.group(0), "stroke:rgb(%d,%d,%d)" % c).replace("stroke-width:0.025", "stroke-width:0.16; stroke-linecap:round")
        new = re.sub(r"\bx2=\"[-\d.eE]+\"", 'x2="%.5f"' % ax, new)
        new = re.sub(r"\by2=\"[-\d.eE]+\"", 'y2="%.5f"' % ay, new)
        return new
    src = re.sub(r"<line\b[^>]*/>", fix_line, src)
    if not a.keep_labels:
        src = re.sub(r"<text\b[^>]*oid=[^>]*>\s*(?:C|N|O|S|P|F|Cl|Br|I|H)\s*</text>\s*", "", src)
    merged = []
    for x, y, c in balls:
        for m_ in merged:
            if math.hypot(m_[0] - x, m_[1] - y) < 0.3:
                m_[3].append((x, y)); break
        else:
            merged.append([x, y, c, [(x, y)]])
    balls = [(sum(p[0] for p in m_[3]) / len(m_[3]), sum(p[1] for p in m_[3]) / len(m_[3]), m_[2]) for m_ in merged]
    dots = "".join('<circle cx="%.5f" cy="%.5f" r="0.215" fill="rgb(%d,%d,%d)" stroke="rgb(%d,%d,%d)" stroke-width="0.03"/>' %
                   (x, y, c[0], c[1], c[2], int(c[0] * .55), int(c[1] * .55), int(c[2] * .55)) for x, y, c in balls)
    src = src.replace("</svg>", dots + "</svg>")
    if not a.no_dist:
        src = add_distances(src, name)
    if not a.no_metal:
        try:
            src = metal_2d(src, name)
        except Exception as e:          # never stop the batch because of the metal step
            print("   (metal contacts skipped in 2D for %s: %s)" % (name, e))
    if not a.no_legend:
        src = add_legend(src, name)
    src = src.replace("stroke-width:0.05;", "stroke-width:0.07;")      # slightly thicker interaction dashes
    return src


# =====================================================================================================
#  3D:  DS geometry picture (no text)  +  crisp labels drawn here
# =====================================================================================================
import numpy as np
from PIL import ImageDraw, ImageFont

FONT_FILES = ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
              "/usr/share/fonts/TTF/DejaVuSans-Bold.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
              "/usr/share/fonts/truetype/freefont/FreeSansBold.ttf"]
def get_font(px):
    cands = list(FONT_FILES)
    try:
        import matplotlib.font_manager as fm
        cands.insert(0, fm.findfont(fm.FontProperties(family="DejaVu Sans", weight="bold")))
    except Exception:
        pass
    for f in cands:
        if os.path.isfile(f):
            return ImageFont.truetype(f, int(px))
    return ImageFont.load_default()

# interaction type -> (line/number colour)
TYPE_COL = [("unfavorable", (220, 0, 0)), ("metal", (130, 90, 180)), ("lone pair", (120, 180, 0)), ("carbon hydrogen", (60, 160, 110)), ("hydrogen", (0, 160, 0)), ("halogen", (0, 150, 150)),
            ("pi-pi", (225, 40, 150)), ("pi pi", (225, 40, 150)), ("stack", (225, 40, 150)), ("t-shaped", (225, 40, 150)),
            ("sigma", (140, 50, 210)), ("alkyl", (190, 60, 190)), ("cation", (240, 120, 0)), ("anion", (240, 120, 0)),
            ("charge", (240, 120, 0)), ("electrostatic", (240, 120, 0)), ("sulfur", (200, 160, 0)), ("sulphur", (200, 160, 0)),
            ("metal", (110, 110, 110)), ("donor", (0, 160, 0))]
def type_col(t):
    t = t.lower()
    for k, c in TYPE_COL:
        if k in t:
            return c
    return (200, 0, 0)

HUES = [(255, 0, 0), (255, 255, 0), (0, 255, 0), (0, 255, 255), (0, 0, 255), (255, 0, 255)]
HUE_DEG = [0, 60, 120, 180, 240, 300]

def find_blob(img_hsv, hue_idx):
    h, s, v = img_hsv
    d = np.abs(((h - HUE_DEG[hue_idx] / 360.0) + 0.5) % 1.0 - 0.5)
    m = (d < 0.045) & (s > 0.65) & (v > 0.25)
    ys, xs = np.nonzero(m)
    if len(xs) < 30:
        return None
    cx, cy = float(np.median(xs)), float(np.median(ys))
    for _ in range(4):                     # densest cluster (ignore stray pixels of the same hue)
        k = (xs - cx) ** 2 + (ys - cy) ** 2 < 70 ** 2
        if k.sum() < 20:
            break
        cx, cy = float(xs[k].mean()), float(ys[k].mean())
    k = (xs - cx) ** 2 + (ys - cy) ** 2 < 70 ** 2
    return (cx, cy) if k.sum() >= 30 else None

def fit_projection(X, uv):
    """3x4 camera matrix from 3D points and pixels (DLT with outlier rejection). returns (P, rms, n)"""
    X = np.asarray(X, float); uv = np.asarray(uv, float)
    keep = np.ones(len(X), bool)
    for _ in range(6):
        if keep.sum() < 6:
            return None, 1e9, int(keep.sum())
        Xk, uk = X[keep], uv[keep]
        mx, sx = Xk.mean(0), Xk.std() or 1.0
        mu, su = uk.mean(0), uk.std() or 1.0
        Xn, un = (Xk - mx) / sx, (uk - mu) / su
        A = []
        for (x, y, z), (u, v) in zip(Xn, un):
            A.append([x, y, z, 1, 0, 0, 0, 0, -u * x, -u * y, -u * z, -u])
            A.append([0, 0, 0, 0, x, y, z, 1, -v * x, -v * y, -v * z, -v])
        _, _, Vt = np.linalg.svd(np.array(A))
        Pn = Vt[-1].reshape(3, 4)
        T3 = np.array([[1 / sx, 0, 0, -mx[0] / sx], [0, 1 / sx, 0, -mx[1] / sx], [0, 0, 1 / sx, -mx[2] / sx], [0, 0, 0, 1]])
        T2i = np.array([[su, 0, mu[0]], [0, su, mu[1]], [0, 0, 1]])
        P = T2i @ Pn @ T3
        def proj(Xs):
            h = np.c_[Xs, np.ones(len(Xs))] @ P.T
            return h[:, :2] / h[:, 2:3]
        err = np.linalg.norm(proj(X) - uv, axis=1)
        rms = float(np.sqrt((err[keep] ** 2).mean()))
        bad = keep & (err > max(6.0, 3 * rms))
        if not bad.any():
            return P, rms, int(keep.sum())
        keep[np.argmax(np.where(keep, err, -1))] = False
    return P, rms, int(keep.sum())

def read_tsv(path):
    if not os.path.isfile(path):
        return []
    L = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    if not L:
        return []
    head = L[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in L[1:] if ln.strip()]

def xyz3(t):
    try:
        v = [float(q) for q in t.split(",")]
        return v if len(v) == 3 else None
    except Exception:
        return None

def res_text(raw):
    m = re.fullmatch(r"([A-Za-z]+)(-?\d+[A-Za-z]?)", raw.strip())
    if not m or a.res_fmt == "raw":
        return raw.strip().upper()
    if a.res_fmt == "cap":
        return m.group(1).capitalize() + m.group(2)
    return m.group(1).upper() + " " + m.group(2)

def calibrate(name, tag):
    pts = read_tsv(os.path.join(a.folder, name + "_calib.tsv"))
    X, uv = [], []
    for p in range(2):
        mp = os.path.join(a.folder, "calib", "%s_%s_m%d.png" % (name, tag, p))
        if not os.path.isfile(mp):
            continue
        im = Image.open(mp).convert("RGB")
        hsv = np.asarray(im.convert("HSV"), float) / 255.0
        hsv = (hsv[..., 0], hsv[..., 1], hsv[..., 2])
        for r in pts:
            if int(r["pass"]) != p:
                continue
            b = find_blob(hsv, int(r["hue"]))
            if b:
                X.append([float(r["x"]), float(r["y"]), float(r["z"])]); uv.append(b)
    if len(X) < 6:
        return None, 0, len(X)
    return fit_projection(X, uv)

def dashed(draw, p0, p1, col, w, dash, gap):
    x0, y0 = p0; x1, y1 = p1
    L = math.hypot(x1 - x0, y1 - y0)
    if L < 1:
        return
    ux, uy = (x1 - x0) / L, (y1 - y0) / L
    t = 0.0
    while t < L:
        e = min(t + dash, L)
        draw.line([(x0 + ux * t, y0 + uy * t), (x0 + ux * e, y0 + uy * e)], fill=col, width=int(w))
        r = w / 2
        for q in (t, e):
            draw.ellipse([x0 + ux * q - r, y0 + uy * q - r, x0 + ux * q + r, y0 + uy * q + r], fill=col)
        t += dash + gap

def rect_overlap(r1, r2, pad=4):
    return not (r1[2] + pad < r2[0] or r2[2] + pad < r1[0] or r1[3] + pad < r2[1] or r2[3] + pad < r1[1])

def compose(name, tag, geo, out):
    P, rms, n = calibrate(name, tag)
    img = Image.open(geo).convert("RGB")
    W0, H0 = img.size
    if P is None or rms > 0.01 * W0:
        print("   !! %s %s: camera calibration failed (points=%d, rms=%.1f px) -> picture saved WITHOUT labels" % (name, tag, n, rms))
        autocrop(img, 0.05).save(out, dpi=(a.dpi, a.dpi))
        return False
    M = int(0.11 * W0)                                    # white margin for labels outside the molecule
    canvas = Image.new("RGB", (W0 + 2 * M, H0 + 2 * M), (255, 255, 255))
    canvas.paste(img, (M, M))
    W, H = canvas.size
    def proj(p):
        h = P @ np.array([p[0], p[1], p[2], 1.0])
        return (h[0] / h[2] + M, h[1] / h[2] + M)
    # obstacle map (non-white pixels) as an integral image
    arr = np.asarray(canvas, np.int16)
    occ = (np.abs(arr - 255).sum(2) > 60).astype(np.int32)
    ii = np.pad(occ.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    def occupancy(r):
        x0, y0, x1, y1 = [int(round(v)) for v in r]
        x0, y0, x1, y1 = max(0, x0), max(0, y0), min(W, x1), min(H, y1)
        if x1 <= x0 or y1 <= y0:
            return 1.0
        return (ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0]) / float((x1 - x0) * (y1 - y0))
    d = ImageDraw.Draw(canvas)
    fr = get_font(a.res_font * W0); fd = get_font(a.dist_font * W0)
    lw = max(3, a.line_w * W0)
    nbs = read_tsv(os.path.join(a.folder, name + "_nonbonds.tsv"))
    rows = read_tsv(os.path.join(a.folder, name + "_labels.tsv"))
    mrows, _lg = metal_contacts(name)
    for r in mrows:
        nbs.append(r)
        if not any(x.get("kind") == "RES" and x["text"].strip().upper() == r["metal"] for x in rows):
            mx_, my_, mz_ = r["metal_xyz"]
            rows.append({"kind": "RES", "text": r["metal"], "x": str(mx_), "y": str(my_), "z": str(mz_)})
    lig = [r for r in rows if r["kind"] == "LIGAND"]
    res = [r for r in rows if r["kind"] == "RES"]
    lig_uv = proj([float(lig[0]["x"]), float(lig[0]["y"]), float(lig[0]["z"])]) if lig else (W / 2, H / 2)

    # 1) interaction lines (strong colours, dashed)
    segs = []
    if not a.ds_lines:
        for r in nbs:
            p0, p1 = xyz3(r.get("from_xyz", "")), xyz3(r.get("to_xyz", ""))
            if p0 and p1:
                dashed(d, proj(p0), proj(p1), type_col(r["type"]), lw, lw * 3.2, lw * 2.2)
    for r in nbs:
        p0, p1 = xyz3(r.get("from_xyz", "")), xyz3(r.get("to_xyz", ""))
        if p0 and p1:
            segs.append((r, proj(p0), proj(p1)))

    placed = []
    def tsize(txt, font):
        b = d.textbbox((0, 0), txt, font=font, stroke_width=5)
        return b[2] - b[0], b[3] - b[1], b[0], b[1]

    # 2) residue names: pushed outward from the ligand into the emptiest nearby spot
    jobs = []
    for r in res:
        A = proj([float(r["x"]), float(r["y"]), float(r["z"])])
        jobs.append((math.hypot(A[0] - lig_uv[0], A[1] - lig_uv[1]), r, A))
    for _, r, A in sorted(jobs, key=lambda j: -j[0]):
        txt = res_text(r["text"])
        tw, th, ox, oy = tsize(txt, fr)
        dx, dy = A[0] - lig_uv[0], A[1] - lig_uv[1]
        n_ = math.hypot(dx, dy) or 1.0
        base = math.atan2(dy / n_, dx / n_)
        best = None
        for rad in [0.03, 0.05, 0.075, 0.10, 0.13, 0.17, 0.22]:
            for da in [0, 20, -20, 40, -40, 65, -65, 90, -90, 120, -120, 150, -150, 180]:
                ang = base + math.radians(da)
                cx, cy = A[0] + math.cos(ang) * rad * W0, A[1] + math.sin(ang) * rad * W0
                box = (cx - tw / 2 - 6, cy - th / 2 - 6, cx + tw / 2 + 6, cy + th / 2 + 6)
                if box[0] < 8 or box[1] < 8 or box[2] > W - 8 or box[3] > H - 8:
                    continue
                if any(rect_overlap(box, p) for p in placed):
                    continue
                sc = occupancy(box) * 12 + rad * 3 + abs(da) / 180.0 * 0.6
                if best is None or sc < best[0]:
                    best = (sc, box, cx, cy)
        if best is None:
            continue
        _, box, cx, cy = best
        placed.append(box)
        # thin leader to the residue when the label is not touching it
        if math.hypot(cx - A[0], cy - A[1]) > th * 1.3:
            ex, ey = cx, cy
            t = max(0.0, 1.0 - (th * 0.75) / (math.hypot(cx - A[0], cy - A[1]) or 1))
            d.line([A, (A[0] + (ex - A[0]) * t, A[1] + (ey - A[1]) * t)], fill=(120, 120, 120), width=max(2, int(lw * 0.45)))
        d.text((cx - tw / 2 - ox, cy - th / 2 - oy), txt, font=fr, fill=(0, 0, 0), stroke_width=5, stroke_fill=(255, 255, 255))
        d.text((cx - tw / 2 - ox, cy - th / 2 - oy), txt, font=fr, fill=(0, 0, 0))

    # 3) distances: at the middle of each line, coloured like the line, white halo
    for r, p0, p1 in segs:
        try:
            dist = float(r["distance"])
        except (ValueError, KeyError):
            continue
        txt = "%.*f" % (a.dec, dist)
        tw, th, ox, oy = tsize(txt, fd)
        mx, my = (p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2
        L = math.hypot(p1[0] - p0[0], p1[1] - p0[1]) or 1
        px, py = -(p1[1] - p0[1]) / L, (p1[0] - p0[0]) / L
        best = None
        for t in (0.5, 0.42, 0.58, 0.35, 0.65):
            for off in (0, 1, -1, 2, -2):
                cx = p0[0] + (p1[0] - p0[0]) * t + px * off * th * 0.9
                cy = p0[1] + (p1[1] - p0[1]) * t + py * off * th * 0.9
                box = (cx - tw / 2 - 5, cy - th / 2 - 5, cx + tw / 2 + 5, cy + th / 2 + 5)
                if box[0] < 4 or box[1] < 4 or box[2] > W - 4 or box[3] > H - 4:
                    continue
                if any(rect_overlap(box, p, 2) for p in placed):
                    continue
                sc = abs(off) * 0.8 + abs(t - 0.5) * 2 + occupancy(box) * 0.5
                if best is None or sc < best[0]:
                    best = (sc, box, cx, cy)
        if best is None:
            continue
        _, box, cx, cy = best
        placed.append(box)
        col = type_col(r["type"])
        dark = tuple(int(c * 0.78) for c in col)
        d.text((cx - tw / 2 - ox, cy - th / 2 - oy), txt, font=fd, fill=dark, stroke_width=5, stroke_fill=(255, 255, 255))
        d.text((cx - tw / 2 - ox, cy - th / 2 - oy), txt, font=fd, fill=dark)

    im = autocrop(canvas, 0.02)
    if a.upscale_3d and a.upscale_3d != 1:
        im = im.resize((int(im.width * a.upscale_3d), int(im.height * a.upscale_3d)), Image.LANCZOS)
    im.save(out, dpi=(a.dpi, a.dpi))
    print("   labelled %s %s  (camera rms %.1f px, %d points, %d residues, %d distances)" % (name, tag, rms, n, len(res), len(segs)))
    return True


for svg in sorted(glob.glob(os.path.join(a.folder, "*_2D.svg"))):
    name = os.path.basename(svg)[:-4]
    if a.only and name[:-3] != a.only:
        continue
    tmp = os.path.join(OUT, name + "_tmp.png")
    src = clean_svg(open(svg, encoding="utf-8").read(), name[:-3])
    cairosvg.svg2png(bytestring=src.encode("utf-8"), write_to=tmp,
                     output_width=int(a.width_in * a.dpi * 1.4), background_color="white")
    im = autocrop(Image.open(tmp))
    im.save(os.path.join(OUT, name + ".png"), dpi=(a.dpi, a.dpi))
    os.remove(tmp)
    print("2D OK", name, im.size)



jobs3d = []
for geo in sorted(glob.glob(os.path.join(a.folder, "*_3D_geo.png"))):
    nm = os.path.basename(geo)[:-len("_3D_geo.png")]
    jobs3d.append((nm, "v0", geo, os.path.join(OUT, nm + "_3D.png")))
vdir_out = os.path.join(OUT, "views"); os.makedirs(vdir_out, exist_ok=True)
for geo in sorted(glob.glob(os.path.join(a.folder, "views", "*_view*_geo.png"))):
    m = re.match(r"(.+)_view(\d+)_geo\.png$", os.path.basename(geo))
    if m:
        jobs3d.append((m.group(1), "v" + m.group(2), geo, os.path.join(vdir_out, "%s_view%s.png" % (m.group(1), m.group(2)))))
for nm, tag, geo, out in jobs3d:
    if a.only and nm != a.only:
        continue
    compose(nm, tag, geo, out)
    print("3D OK", os.path.basename(out))

print("\nFinal images are in:", OUT)
