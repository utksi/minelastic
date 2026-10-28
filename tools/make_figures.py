# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
"""Figures and symbol of the README and the documentation: what the minimal strain sets do.

The landing figure is drawn twice, on a transparent background, with the colours of light and of dark pages.
Left: the usual approach, six strain patterns, each one strain component. Right: a cubic crystal under the
single pattern u = xx + yz of minelastic, with the stress components that give its three constants,
sigma_xx = c11 delta, sigma_yy = sigma_zz = c12 delta and sigma_yz = c44 delta (Voigt notation, engineering
shear strain). Bottom: the number of strain patterns K of the minimal set for each crystal system, from
minelastic (one representative space group each), out of the six of the usual approach.

Cells are parallelepipeds with edge vectors a_i' = F a_i, F = I + e (symmetric strain tensor e, drawn
exaggerated), in the oblique projection (x, y, z) -> (x + k y cos t, z + k y sin t), k = 0.5, t = 30 deg.
A face is drawn if its outward normal n has n . d > 0, with d = (k cos t, -1, k sin t) the direction
towards the viewer (P d = 0). The stress arrows are the tractions t = sigma n on the visible faces.

The symbol (docs/assets/logo.svg) is a unit cube under the shear F = I + e, e_xz = e_zx = 0.17, in the
isometric projection (x, y, z) -> ((x - y) cos 30, z - (x + y) sin 30), white on a rounded blue tile, with a
teal stress arrow normal to its +x face. The wordmark is the symbol followed by the name in bold Roboto, with
the letters converted to outlines.

Text is set in Roboto, the font of the documentation: regular and bold from the paths in MINELASTIC_FONTS
(colon-separated .ttf files) or DEFAULT_FONTS, medium (with Greek letters) for the formulas. Without them,
DejaVu Sans is used.

    python tools/make_figures.py    # writes docs/assets: hero-light.svg, hero-dark.svg, social-preview.png,
                                    # logo.svg, wordmark-light.svg, wordmark-dark.svg

The icon of the Claude Code plugin is the symbol rendered with librsvg:

    convert -background none -density 1536 docs/assets/logo.svg -resize 1024x1024 -depth 8 .claude-plugin/icon.png
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from ase.spacegroup import crystal  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402

from minelastic import ElasticSetup  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "assets")
_R = "/cvmfs/software.eessi.io/versions/2025.06/software/linux/x86_64/amd/zen4/software"
DEFAULT_FONTS = [f"{_R}/R-bundle-CRAN/2025.11-foss-2025b/gfonts/assets/fonts/roboto/roboto-v20-latin-regular.ttf",
                 f"{_R}/R-bundle-CRAN/2025.11-foss-2025b/gfonts/assets/fonts/roboto/roboto-v20-latin-700.ttf",
                 f"{_R}/R/4.4.2-gfbf-2024a/lib/R/library/grDevices/fonts/Roboto/Roboto-Medium.ttf"]


def use_fonts():
    """Registers the Roboto files of MINELASTIC_FONTS or DEFAULT_FONTS with matplotlib and sets the text and
    formula fonts (formulas upright, in Roboto Medium); returns the family name used."""
    paths = [p for p in os.environ.get("MINELASTIC_FONTS", ":".join(DEFAULT_FONTS)).split(":") if os.path.exists(p)]
    for p in paths:
        font_manager.fontManager.addfont(p)
    family = "Roboto" if paths else "DejaVu Sans"
    plt.rcParams.update({"font.family": family, "svg.fonttype": "path", "mathtext.fontset": "custom",
                         "mathtext.rm": f"{family}:medium", "mathtext.it": f"{family}:medium",
                         "mathtext.bf": f"{family}:bold", "mathtext.fallback": "stixsans"})
    return family


def projection(k, t):
    """(P, d): the oblique projection P = [[1, k cos t, 0], [0, k sin t, 1]] (t in degrees) and the direction
    d = (k cos t, -1, k sin t) towards the viewer, P d = 0."""
    c, s = k * np.cos(np.radians(t)), k * np.sin(np.radians(t))
    return np.array([[1.0, c, 0.0], [0.0, s, 1.0]]), np.array([c, -1.0, s])


MAIN, ICON = projection(0.5, 30), projection(0.75, 42)       # the cells, and the crystal-system icons
P = MAIN[0]
LIGHT = np.array([0.3, -0.5, 1.0]) / np.linalg.norm([0.3, -0.5, 1.0])

# colours for light and dark pages (transparent background); ink2 is the secondary text, C the active set
THEMES = {"light": dict(ink="#16202e", ink2="#3d5a80", card_left="#f1f3f7", card_right="#e9f1fc", rule="#dfe4ec",
                        empty="#b7c0cf", shadow="#1b1f27", gray_face="#dfe2e8", gray_edge="#7d8492",
                        blue="#2a78d6", blue_text="#1d5aa6", blue_face="#bcd5f5", blue_edge="#1d5aa6",
                        atom="#1d5aa6", stress="#0f9e7b", name="#142c52"),
          "dark": dict(ink="#eef1f6", ink2="#a9c3e8", card_left="#262b35", card_right="#1e2a3e", rule="#3a404c",
                       empty="#5f6878", shadow="#000000", gray_face="#626a78", gray_edge="#b9c1ce",
                       blue="#4a90e8", blue_text="#86b6f2", blue_face="#5b8fd8", blue_edge="#c3dafa",
                       atom="#9cc3f5", stress="#2fe0b0", name="#ffffff")}
C = dict(THEMES["light"])
SURFACE = "#fbfcfe"                                   # opaque background of the social preview
LOGO = dict(tile="#2a78d6", stroke="#1d5aa6", top="#ffffff", right="#d6e6fb", left="#a3c4ef", arrow="#2fe0b0")

# crystal systems: name, representative space group, cell (a, b, c in A; angles in degrees)
SYSTEMS = [("triclinic", 2, [4.1, 4.7, 5.3, 81, 97, 103]),
           ("monoclinic", 12, [12.2, 3.0, 5.8, 90, 103.8, 90]),
           ("orthorhombic", 70, [8.3, 4.8, 8.6, 90, 90, 90]),
           ("tetragonal", 136, [4.6, 4.6, 2.96, 90, 90, 90]),
           ("trigonal", 154, [4.9, 4.9, 5.4, 90, 90, 120]),
           ("hexagonal", 194, [3.2, 3.2, 5.2, 90, 90, 120]),
           ("cubic", 225, [4.04, 4.04, 4.04, 90, 90, 90])]
# edge vectors (rows a, b, c) of the cell icons of the crystal systems
ICONS = {"triclinic": [(1.0, 0, 0), (0.25, 0.85, 0.08), (0.3, 0.18, 1.0)],
         "monoclinic": [(1.05, 0, 0), (0, 0.85, 0), (0.38, 0, 1.0)],
         "orthorhombic": [(1.25, 0, 0), (0, 0.65, 0), (0, 0, 0.95)],
         "tetragonal": [(0.8, 0, 0), (0, 0.8, 0), (0, 0, 1.35)],
         "trigonal": [(0.9, 0, 0), (0.308, 0.846, 0), (0.308, 0.216, 0.818)],      # rhombohedron, 70 deg
         "cubic": [(1, 0, 0), (0, 1, 0), (0, 0, 1)]}


def n_patterns():
    """{system: K}, the number of strain patterns of the minimal set of a representative crystal."""
    out = {}
    for name, sg, cellpar in SYSTEMS:
        atoms = crystal(["Si", "O"], basis=[(0.1123, 0.2234, 0.3345), (0.4111, 0.0789, 0.2222)],
                        spacegroup=sg, cellpar=cellpar)
        out[name] = len(ElasticSetup(atoms).patterns)
    return out


def shade(color, normal):
    """RGB of color (hex) for a face of unit normal n: c (0.66 + 0.34 t) + (1 - c) 0.45 t^2, t = max(0, n . LIGHT),
    a Lambert factor with the lit faces lifted towards white."""
    rgb = np.array(matplotlib.colors.to_rgb(color))
    t = max(0.0, float(normal @ LIGHT))
    return np.clip(rgb * (0.66 + 0.34 * t) + (1 - rgb) * 0.45 * t ** 2, 0, 1)


def hull(points):
    """Convex hull (counter-clockwise vertices) of 2D points, by the monotone-chain algorithm."""
    pts = sorted(map(tuple, np.round(points, 9)))
    cross = lambda o, a, b: (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])      # noqa: E731
    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return np.array(lower[:-1] + upper[:-1])


def shadow(ax, vertices, offset=(0.22, -0.16), scale=1.0, proj=MAIN):
    """Soft contact shadow: the hull of the footprint (all vertices dropped to the floor z = min z), drawn ten
    times with opacity 0.02, shifted on the page by offset times scale times 0.2, 0.4, ..., 2."""
    floor = vertices.copy()
    floor[:, 2] = vertices[:, 2].min()
    base = hull((proj[0] @ floor.T).T)
    for k in np.linspace(0.2, 2.0, 10):
        ax.add_patch(Polygon(base + scale * k * np.asarray(offset), closed=True, facecolor=C["shadow"],
                             edgecolor="none", alpha=0.02, zorder=1))


def atoms(ax, vertices, radius, color, z=5):
    """Spheres at the vertices of a parallelepiped (all but vertex 2, (0, 1, 0), which faces away from the
    viewer): a disc with a darker rim and an offset highlight."""
    rgb = np.array(matplotlib.colors.to_rgb(color))
    for k, v in enumerate(vertices):
        if k == 2:
            continue
        c = P @ v
        ax.add_patch(matplotlib.patches.Circle(c, radius, facecolor=rgb, edgecolor=rgb * 0.6, lw=0.8, zorder=z))
        ax.add_patch(matplotlib.patches.Circle(c + radius * np.array([-0.3, 0.32]), 0.38 * radius,
                                               facecolor=(1, 1, 1, 0.55), edgecolor="none", zorder=z))


def faces_of(vertices, faces, proj=MAIN):
    """[(page polygon (n, 2), unit outward normal (3,))] of the faces (vertex index lists) of a convex body."""
    center = vertices.mean(axis=0)
    out = []
    for f in faces:
        v = vertices[f]
        n = np.cross(v[1] - v[0], v[2] - v[0])
        n /= np.linalg.norm(n)
        if n @ (v.mean(axis=0) - center) < 0:
            n = -n
        out.append(((proj[0] @ v.T).T, n))
    return out


def parallelepiped(origin, edges):
    """Vertices (8, 3) and faces of the cell origin + i a + j b + k c, i, j, k in {0, 1}; edges = rows a, b, c."""
    a, b, c = np.asarray(edges, float)
    v = np.array([origin + i * a + j * b + k * c for i in (0, 1) for j in (0, 1) for k in (0, 1)])
    faces = [[0, 1, 3, 2], [4, 5, 7, 6], [0, 1, 5, 4], [2, 3, 7, 6], [0, 2, 6, 4], [1, 3, 7, 5]]
    return v, faces


def prism(origin, radius, height, n=6):
    """Vertices and faces of a regular n-gonal prism of circumradius radius about the z axis (hexagonal cell)."""
    phi = np.radians(np.arange(n) * 360 / n + 180 / n)
    ring = np.c_[radius * np.cos(phi), radius * np.sin(phi), np.zeros(n)]
    v = np.vstack([origin + ring, origin + ring + [0, 0, height]])
    faces = [list(range(n)), list(range(n, 2 * n))] + [[i, (i + 1) % n, n + (i + 1) % n, n + i] for i in range(n)]
    return v, faces


def draw_body(ax, vertices, faces, color, edge, lw=1.3, z=3, proj=MAIN):
    """Visible faces (n . d > 0, d towards the viewer) of a convex body, shaded, with edges."""
    for poly, n in faces_of(vertices, faces, proj):
        if n @ proj[1] > 1e-9:
            ax.add_patch(Polygon(poly, closed=True, facecolor=shade(color, n), edgecolor=edge, lw=lw,
                                 joinstyle="round", zorder=z))


def ghost(ax, vertices, faces, z=2):
    """Dashed outline of the visible faces of the undeformed body."""
    for poly, n in faces_of(vertices, faces):
        if n @ MAIN[1] > 1e-9:
            ax.add_patch(Polygon(poly, closed=True, fill=False, edgecolor=C["ink2"], lw=0.9, ls=(0, (3, 2.5)),
                                 zorder=z, alpha=0.7))


def strained_cube(ax, origin, size, eps, color, edge, lw=1.3, atom=None, atom_color=None):
    """Cube of edge size at origin, deformed by F = I + eps about origin, over the dashed undeformed cube,
    with atoms of radius atom at its corners (none if atom is None). Returns F."""
    F = np.eye(3) + eps
    v0, faces = parallelepiped(np.asarray(origin, float), size * np.eye(3))
    v, faces = parallelepiped(np.asarray(origin, float), (F @ (size * np.eye(3))).T)
    shadow(ax, v, scale=size / 2)
    ghost(ax, v0, faces)
    draw_body(ax, v, faces, color, edge, lw=lw)
    if atom:
        atoms(ax, v, atom, atom_color or edge)
    return F


def arrow(ax, start3, vec3, label=None, offset=(0, 0), ha="left"):
    """Arrow from the 3D point start3 along vec3 (both projected), with an optional label at its tip."""
    s, e = P @ start3, P @ (start3 + vec3)
    ax.add_patch(FancyArrowPatch(s, e, arrowstyle="-|>", mutation_scale=16, lw=2.2, color=C["stress"], zorder=6,
                                 shrinkA=0, shrinkB=0))
    if label:
        ax.text(e[0] + offset[0], e[1] + offset[1], label, fontsize=13.5, color=C["ink"], ha=ha, va="center",
                zorder=7)


def iso(p):
    """Isometric page coordinates (SVG, y down) of a 3D point, viewed along (1, 1, 1)."""
    x, y, z = p
    return np.array([(x - y) * np.cos(np.radians(30)), -(z - (x + y) * np.sin(np.radians(30)))])


def symbol_geometry(size=64.0):
    """Shapes of the symbol in a size x size box (SVG coordinates, y down): the three visible faces of the
    sheared unit cube [(polygon, colour)] and the arrow (shaft start, shaft end, head triangle), with the cube
    and arrow tip fitted into the tile with a margin of 11/64 of its side."""
    F = np.eye(3) + np.array([[0, 0, 0.17], [0, 0, 0], [0.17, 0, 0]])
    v = lambda i, j, k: F @ np.array([i, j, k], float)                                   # noqa: E731
    faces = [([v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)], LOGO["top"]),             # +z
             ([v(1, 0, 0), v(1, 1, 0), v(1, 1, 1), v(1, 0, 1)], LOGO["right"]),           # +x
             ([v(0, 1, 0), v(1, 1, 0), v(1, 1, 1), v(0, 1, 1)], LOGO["left"])]            # +y
    center = np.mean(faces[1][0], axis=0)
    start, tip = center + np.array([0.04, 0, 0]), center + np.array([0.82, 0, 0])
    pts = np.array([iso(q) for f, _ in faces for q in f] + [iso(tip)])
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    scale = (size - 2 * 11 * size / 64) / (hi - lo).max()
    shift = (size - scale * (hi - lo)) / 2 - scale * lo
    m = lambda q: scale * iso(q) + shift                                                   # noqa: E731
    a, b = m(start), m(tip)
    d = (b - a) / np.linalg.norm(b - a)
    n = np.array([-d[1], d[0]])
    base = b - 7.5 * size / 64 * d
    head = [b, base + 5.2 * size / 64 * n, base - 5.2 * size / 64 * n]
    return [([m(q) for q in f], col) for f, col in faces], (a, base, head)


def logo_elements(size=64.0, x0=0.0, y0=0.0):
    """SVG elements of the symbol in the box [x0, x0 + size] x [y0, y0 + size]: tile, faces, arrow with a white
    edge."""
    faces, (a, base, head) = symbol_geometry(size)
    k = size / 64
    sh = np.array([x0, y0])
    fmt = lambda pts: " ".join(f"{x:.2f},{y:.2f}" for x, y in np.asarray(pts) + sh)      # noqa: E731
    out = [f'<rect x="{x0 + k:.2f}" y="{y0 + k:.2f}" width="{size - 2 * k:.2f}" height="{size - 2 * k:.2f}" '
           f'rx="{0.22 * size:.2f}" fill="{LOGO["tile"]}"/>']
    for poly, col in faces:
        out.append(f'<polygon points="{fmt(poly)}" fill="{col}" stroke="{LOGO["stroke"]}" '
                   f'stroke-width="{1.3 * k:.2f}" stroke-linejoin="round"/>')
    (ax_, ay), (bx, by) = a + sh, base + sh
    for col, extra in (("#ffffff", 1.0), (LOGO["arrow"], 0.0)):                         # white edge, then fill
        out.append(f'<line x1="{ax_:.2f}" y1="{ay:.2f}" x2="{bx:.2f}" y2="{by:.2f}" stroke="{col}" '
                   f'stroke-width="{(4.2 + 2 * extra) * k:.2f}" stroke-linecap="round"/>')
        out.append(f'<polygon points="{fmt(head)}" fill="{col}" stroke="{col}" '
                   f'stroke-width="{(2 * extra + 0.6) * k:.2f}" stroke-linejoin="round"/>')
    return out


def logo_svg(size=64):
    """SVG text of the symbol alone."""
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {size} {size}">\n'
            + "\n".join(logo_elements(size)) + "\n</svg>\n")


def path_d(path, scale, x0, y0):
    """SVG path data of a matplotlib Path (font units, y up), scaled by scale and placed with its origin at
    (x0, y0) of an SVG box (y down)."""
    out, verts, codes, i = [], path.vertices, path.codes, 0
    tr = lambda p: (x0 + scale * p[0], y0 - scale * p[1])                               # noqa: E731
    while i < len(codes):
        c = codes[i]
        if c == Path.MOVETO:
            out.append("M%.2f %.2f" % tr(verts[i]))
            i += 1
        elif c == Path.LINETO:
            out.append("L%.2f %.2f" % tr(verts[i]))
            i += 1
        elif c == Path.CURVE3:
            out.append("Q%.2f %.2f %.2f %.2f" % (tr(verts[i]) + tr(verts[i + 1])))
            i += 2
        elif c == Path.CURVE4:
            out.append("C%.2f %.2f %.2f %.2f %.2f %.2f" % (tr(verts[i]) + tr(verts[i + 1]) + tr(verts[i + 2])))
            i += 3
        else:                                                                             # CLOSEPOLY, STOP
            out.append("Z")
            i += 1
    return " ".join(out)


def wordmark_svg(color, height=64):
    """SVG text of the wordmark: the symbol (height x height) and the name in bold, as outlines, in color; the
    cap height of the name is 0.46 of the symbol height, centred on it."""
    use_fonts()
    prop = font_manager.FontProperties(family=plt.rcParams["font.family"], weight="bold")
    text = TextPath((0, 0), "minelastic", size=100, prop=prop)
    cap = TextPath((0, 0), "H", size=100, prop=prop).get_extents().height
    scale = 0.46 * height / cap
    ext = text.get_extents()
    x0 = height + 0.22 * height - scale * ext.x0
    y0 = height / 2 + 0.46 * height / 2
    width = x0 + scale * ext.x1 + 2
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.1f} {height}">\n'
            + "\n".join(logo_elements(height))
            + f'\n<path d="{path_d(text, scale, x0, y0)}" fill="{color}"/>\n</svg>\n')


def draw_symbol(ax, x, y, size):
    """The symbol drawn with matplotlib, its tile spanning [x, x + size] x [y, y + size] in data units."""
    faces, (a, base, head) = symbol_geometry(64.0)
    to = lambda pts: np.c_[x + np.asarray(pts)[:, 0] * size / 64, y + size - np.asarray(pts)[:, 1] * size / 64]  # noqa
    ax.add_patch(FancyBboxPatch((x + size / 64, y + size / 64), size * 62 / 64, size * 62 / 64,
                                boxstyle=f"round,pad=0,rounding_size={0.22 * size:.3f}", facecolor=LOGO["tile"],
                                edgecolor="none", zorder=8))
    lw = 1.3 * size / 64 * 72 / 1.0 / 2.67 / 1.0
    for poly, col in faces:
        ax.add_patch(Polygon(to(poly), closed=True, facecolor=col, edgecolor=LOGO["stroke"], lw=lw * 0.6,
                             joinstyle="round", zorder=9))
    shaft = to([a, base])
    for col, w in (("#ffffff", 1.6), (LOGO["arrow"], 1.0)):
        ax.plot(shaft[:, 0], shaft[:, 1], color=col, lw=lw * 2.2 * w, solid_capstyle="round", zorder=10)
        ax.add_patch(Polygon(to(head), closed=True, facecolor=col, edgecolor=col, lw=lw * 0.8 * w, zorder=10))


def draw(path, title=False, dpi=100, theme="light"):
    """The landing figure, written to path (format from the extension), in the colours of theme ('light' or
    'dark') on a transparent background; title adds the wordmark and tagline on an opaque light background.
    Returns {system: K}."""
    use_fonts()
    C.update(THEMES[theme])
    counts = n_patterns()
    fig = plt.figure(figsize=(12, 6), facecolor=SURFACE if title else "none")
    ax = fig.add_axes([0, 0, 1, 1])
    if title:                                    # a title row above the same drawing, at the same aspect
        ax.set_xlim(-1.2, 33.2)
        ax.set_ylim(0, 17.2)
        draw_symbol(ax, 0.4, 15.75, 1.3)
        ax.text(2.0, 16.4, "minelastic", fontsize=24, weight="bold", color=C["name"], va="center")
        ax.text(8.6, 16.36, "elastic constants from the fewest stress calculations", fontsize=15, color=C["ink2"],
                va="center")
    else:
        ax.set_xlim(0, 32)
        ax.set_ylim(0, 16)
    ax.set_aspect("equal")
    ax.axis("off")
    card = dict(boxstyle="round,pad=0,rounding_size=0.45", edgecolor="none", zorder=0)
    ax.add_patch(FancyBboxPatch((0.35, 4.85), 12.0, 10.75, facecolor=C["card_left"], **card))
    ax.add_patch(FancyBboxPatch((16.15, 4.85), 15.5, 10.75, facecolor=C["card_right"], **card))

    # left: the six single strains
    ax.text(1.0, 14.85, "Usual approach", fontsize=15, weight="bold", color=C["ink"], va="center")
    ax.text(1.0, 14.05, "six strain patterns, one strain component each", fontsize=12, color=C["ink2"],
            va="center")
    a, s = 0.24, 1.9
    singles = [("xx", 0, 0, a), ("yy", 1, 1, a), ("zz", 2, 2, a),
               ("yz", 1, 2, a / 1.4), ("xz", 0, 2, a / 1.4), ("xy", 0, 1, a / 1.4)]
    for k, (name, i, j, amp) in enumerate(singles):
        eps = np.zeros((3, 3))
        eps[i, j] = eps[j, i] = amp
        x0, y0 = 1.3 + 3.75 * (k % 3), 9.95 - 3.95 * (k // 3)
        strained_cube(ax, (x0, 0, y0), s, eps, C["gray_face"], C["gray_edge"], lw=1.1)
        ax.text(x0 + 1.25, y0 - 0.62, rf"$\varepsilon_{{{name}}}$", fontsize=13, color=C["ink"], ha="center",
                va="center")

    # middle: six to one
    y_mid = 9.9
    ax.add_patch(FancyArrowPatch((12.75, y_mid), (15.85, y_mid), mutation_scale=18, color=C["blue"], zorder=5,
                                 arrowstyle="simple,head_width=1.6,head_length=1.2,tail_width=0.55"))
    ax.text(14.3, y_mid + 1.15, r"$6 \rightarrow 1$", fontsize=22, color=C["blue"], ha="center", va="center")
    ax.text(14.3, y_mid - 1.05, "cubic crystal", fontsize=11, color=C["ink2"], ha="center", va="center")

    # right: one cell under xx + yz, and the stresses that give c11, c12 and c44
    ax.text(16.8, 14.85, "minelastic", fontsize=15, weight="bold", color=C["blue_text"], va="center")
    ax.text(16.8, 14.05, "one strain pattern; each constant from one stress component", fontsize=12,
            color=C["ink2"], va="center")
    S, origin = 4.2, np.array([19.9, 0.0, 6.3])
    F = strained_cube(ax, origin, S, np.array([[0.16, 0, 0], [0, 0, 0.11], [0, 0.11, 0]]), C["blue_face"],
                      C["blue_edge"], lw=1.6, atom=0.27, atom_color=C["atom"])
    cx = origin + F @ np.array([S, S / 2, S / 2])                  # centre of the +x face
    cz = origin + F @ np.array([S / 2, S / 2, S])                  # centre of the +z face
    arrow(ax, cx, np.array([2.1, 0, 0]), r"$\sigma_{xx} = c_{11}\,\delta$", offset=(0.25, 0))
    arrow(ax, cz, np.array([0, 0, 1.6]), r"$\sigma_{zz} = c_{12}\,\delta$", offset=(0.25, 0.05))
    yy = origin + F @ np.array([S * 0.12, 0.0, S * 0.3])             # on the front (-y) face, near its edge
    arrow(ax, yy, np.array([0, -4.6, 0]))
    tip = P @ (yy + np.array([0, -4.6, 0]))
    ax.text(tip[0] - 0.4, tip[1] - 0.55, r"$\sigma_{yy} = c_{12}\,\delta$", fontsize=13.5, color=C["ink"],
            ha="center", va="center", zorder=7)
    start, shear = origin + F @ np.array([S * 0.18, S * 0.12, S]), F @ np.array([0, S * 0.62, 0])
    arrow(ax, start, shear)
    lab = P @ (start + shear)
    ax.text(lab[0] - 1.05, lab[1] + 0.8, r"$\sigma_{yz} = c_{44}\,\delta$", fontsize=13.5, color=C["ink"],
            ha="right", va="center", zorder=7)
    ax.text(31.2, 6.0, "pattern", fontsize=11, color=C["ink2"], ha="right", va="center")
    ax.text(31.2, 5.35, r"$\varepsilon_{xx} + \varepsilon_{yz}$, drawn enlarged", fontsize=11, color=C["ink2"],
            ha="right", va="center")

    # bottom: strain patterns per crystal system
    ax.text(1.0, 4.1, "Strain patterns per crystal system", fontsize=12, weight="bold", color=C["ink"],
            va="center")
    ax.text(31.0, 4.1, "filled: minelastic, all six: usual approach", fontsize=11, color=C["ink2"], va="center",
            ha="right")
    for k, (name, _, _) in enumerate(SYSTEMS):
        xc, yb = 2.6 + 4.55 * k, 1.75
        if name == "hexagonal":
            v, f = prism(np.array([xc - 0.6, 0.75, yb]), 0.75, 1.15)
        else:
            v, f = parallelepiped(np.array([xc - 0.8, 0.0, yb]), np.array(ICONS[name], float))
        shadow(ax, v, scale=0.5, proj=ICON)
        draw_body(ax, v, f, C["blue_face"], C["blue_edge"], lw=1.0, proj=ICON)
        ax.text(xc, 1.2, name, fontsize=11.5, color=C["ink"], ha="center", va="center")
        for i in range(6):
            filled = i < counts[name]
            ax.add_patch(FancyBboxPatch((xc - 1.17 + 0.4 * i, 0.38), 0.3, 0.3,
                                        boxstyle="round,pad=0,rounding_size=0.07", lw=1.0,
                                        facecolor=C["blue"] if filled else "none",
                                        edgecolor=C["blue"] if filled else C["empty"], zorder=4))
    fig.savefig(path, dpi=dpi, facecolor=SURFACE if title else "none", transparent=not title)
    plt.close(fig)
    return counts


def main():
    """Writes the symbol, the wordmarks, the landing figure (for light and dark pages) and the 1280 x 640
    social preview of the repository."""
    os.makedirs(OUT, exist_ok=True)
    print("font:", use_fonts())
    for name, text in (("logo.svg", logo_svg()), ("wordmark-light.svg", wordmark_svg(THEMES["light"]["name"])),
                       ("wordmark-dark.svg", wordmark_svg(THEMES["dark"]["name"]))):
        with open(os.path.join(OUT, name), "w") as f:
            f.write(text)
    counts = draw(os.path.join(OUT, "hero-light.svg"))
    draw(os.path.join(OUT, "hero-dark.svg"), theme="dark")
    draw(os.path.join(OUT, "social-preview.png"), title=True, dpi=1280 / 12)
    print(counts)


if __name__ == "__main__":
    main()
