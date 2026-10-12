"""Schematic illustrations for the 'How it works' page. Each returns an inline <svg>. Colours come from the page's CSS classes, so they follow light and dark mode."""
import math

def _svg(vb, label, body):
    w, h = vb
    return (f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}"><title>{label}</title>{body}</svg>')

def cut(F, R, noisy=0, label=""):
    """An amplicon on a 280-unit scale, with the forward read from the left and the reverse read from the right."""
    L, x0 = 280, 20
    fx0, fx1 = x0, x0 + F; rx0, rx1 = x0 + L - R, x0 + L
    b = [f'<rect class="card" x="{x0}" y="12" width="{L}" height="16" rx="4"/>', f'<text class="tm" x="{x0 + 6}" y="24">amplicon</text>']
    ov0, ov1 = rx0, fx1
    if ov1 > ov0:
        b.append(f'<rect class="band" x="{ov0}" y="36" width="{ov1 - ov0}" height="52"/>')
    else:
        b.append(f'<rect class="warns" x="{fx1}" y="36" width="{rx0 - fx1}" height="52" stroke-dasharray="4 3"/>')
    b.append(f'<rect class="acc" x="{fx0}" y="44" width="{F}" height="12"/><polygon class="acc" points="{fx1},44 {fx1 + 9},50 {fx1},56"/>')
    b.append(f'<rect class="rev" x="{rx0}" y="68" width="{R}" height="12"/><polygon class="rev" points="{rx0},68 {rx0 - 9},74 {rx0},80"/>')
    if noisy:
        b.append(f'<rect class="warn" x="{fx1 - noisy}" y="44" width="{noisy}" height="12"/><rect class="warn" x="{rx0}" y="68" width="{noisy}" height="12"/>')
    b.append(f'<text class="{"tw" if not ov1 > ov0 else "ta"}" x="160" y="108" text-anchor="middle">{label}</text>')
    return _svg((320, 118), f"Schematic: {label}", "".join(b))

def sample():
    b = []
    for r in range(5):
        for c in range(16):
            pick = (r * 16 + c) % 4 == 0
            b.append(f'<circle class="{"acc" if pick else "dim"}" cx="{28 + c * 17.5:.1f}" cy="{20 + r * 18}" r="{5 if pick else 4}"/>')
    b.append('<text class="tm" x="160" y="118" text-anchor="middle">every Nth pair, across all samples</text>')
    return _svg((320, 128), "Schematic: a regular sample of read pairs, highlighted among all pairs", "".join(b))

def _bars(n=22):
    return [max(0.04, math.exp(-((i - 7) ** 2) / 26.0) + 0.55 * math.exp(-((i - 13) ** 2) / 14.0) * (i > 7) * 0.4 + (0.03 if i > 18 else 0)) for i in range(n)]

def hist(floor=False, cap=0.0):
    hs = _bars(); m = max(hs); b = []
    bw = 270 / len(hs); cutoff = int(len(hs) * 0.88)
    for i, v in enumerate(hs):
        h = 70 * v / m; x = 24 + i * bw
        late = floor and i >= cutoff
        b.append(f'<rect class="{"dim" if late else "acc"}" x="{x:.1f}" y="{104 - h:.1f}" width="{bw - 2:.1f}" height="{h:.1f}" rx="1.5"/>')
    b.append('<line class="ink" x1="22" y1="104" x2="298" y2="104"/>')
    if floor:
        xc = 24 + cutoff * bw - 1; xf = xc + 28
        b.append(f'<line class="warns" x1="{xc:.1f}" y1="22" x2="{xc:.1f}" y2="104" stroke-dasharray="4 3"/>')
        b.append(f'<text class="tw" x="{xc - 4:.1f}" y="14" text-anchor="end">longest {cap} set aside</text>')
        b.append(f'<line class="accs" x1="{xf:.1f}" y1="22" x2="{xf:.1f}" y2="104" stroke-width="2"/><text class="ta" x="{xf + 4:.1f}" y="14">floor</text>')
        b.append(f'<line class="ink" x1="{xc:.1f}" y1="114" x2="{xf:.1f}" y2="114"/><line class="ink" x1="{xc:.1f}" y1="110" x2="{xc:.1f}" y2="118"/><line class="ink" x1="{xf:.1f}" y1="110" x2="{xf:.1f}" y2="118"/>')
        b.append(f'<text class="tm" x="{(xc + xf) / 2:.1f}" y="129" text-anchor="middle">+ 12</text><text class="tm" x="24" y="122">amplicon length</text>')
    else:
        b.append('<text class="tm" x="160" y="122" text-anchor="middle">the lengths of your merged reads = your amplicon lengths</text>')
    return _svg((320, 134), "Schematic: a histogram of amplicon lengths" + (", with the longest ones set aside and the floor marked" if floor else ""), "".join(b))

def search():
    b = ['<text class="tm" x="20" y="12">quality table</text><text class="tm" x="176" y="12">yield by total length</text>']
    n, s = 8, 14
    for i in range(n):          # rows: reverse length, columns: forward length
        for j in range(n):
            tot = i + j; v = max(0.08, 1.0 - 0.07 * tot + (0.1 if tot > 11 else 0))
            b.append(f'<rect class="acc" x="{20 + j * s}" y="{20 + i * s}" width="{s - 1}" height="{s - 1}" opacity="{v:.2f}"/>')
            if tot == 7: b.append(f'<rect class="ink" x="{20 + j * s}" y="{20 + i * s}" width="{s - 1}" height="{s - 1}"/>')
            if (i, j) == (3, 4): b.append(f'<circle class="warn" cx="{20 + j * s + 6}" cy="{20 + i * s + 6}" r="3.5"/>')
    b.append('<text class="tm" x="20" y="146">forward length \u2192</text><text class="tm" x="20" y="159">outline = one total length, dot = its best split</text>')
    pts = [(184, 112), (204, 86), (224, 62), (244, 50), (264, 47), (284, 46)]
    b.append('<line class="ink" x1="180" y1="118" x2="296" y2="118"/><line class="ink" x1="180" y1="118" x2="180" y2="34"/>')
    b.append('<polyline class="accs" fill="none" stroke-width="2" points="' + " ".join(f"{x},{y}" for x, y in pts) + '"/>')
    for x, y in pts: b.append(f'<circle class="acc" cx="{x}" cy="{y}" r="3"/>')
    b.append('<line class="warns" x1="264" y1="32" x2="264" y2="118" stroke-dasharray="4 3"/><text class="tw" x="264" y="29" text-anchor="middle">stops here</text>')
    b.append('<text class="tm" x="238" y="132" text-anchor="middle">longer total \u2192</text>')
    return _svg((320, 166), "Schematic: a table of pass shares for every cut, and a curve of yield that rises and then levels off", "".join(b))

def pick():
    ys = [("shorter", 0.62, False), ("", 0.84, False), ("", 0.915, True), ("", 0.92, True), ("longer", 0.91, True)]
    b = ['<text class="tm" x="72" y="12">bars = yield</text>', '<rect class="band" x="14" y="68" width="292" height="78" rx="6"/>',
         '<text class="ta" x="298" y="92" text-anchor="end">within 1 point</text><text class="ta" x="298" y="106" text-anchor="end">of the best</text>']
    for i, (lab, v, near) in enumerate(ys):
        y = 22 + i * 26
        b.append(f'<text class="tm" x="18" y="{y + 11}">{lab}</text>')
        b.append(f'<rect class="{"acc" if near else "dim"}" x="72" y="{y}" width="{120 * v:.1f}" height="14" rx="3"/>')
    b.append(f'<text class="ok" x="192" y="{22 + 2 * 26 + 13}" style="font:700 15px system-ui,sans-serif">\u2713</text>')
    b.append('<text class="tm" x="160" y="160" text-anchor="middle">the tie goes to the cut with the least overlap</text>')
    return _svg((320, 168), "Schematic: five cuts ranked by yield; the cuts within one point of the best are tied, and the shortest of them is chosen", "".join(b))

def ee():
    xs = [i / 200.0 for i in range(201)]; ys = [(x ** 1.4) * math.exp(-x * 5.2) for x in xs]; ymax = max(ys); tot = sum(ys)
    cum, p75 = 0.0, 0.5
    for x, y in zip(xs, ys):
        cum += y
        if cum / tot >= 0.75: p75 = x; break                                           # the shaded part really is 75% of the curve
    px = lambda x: 24 + 270 * x; py = lambda y: 96 - 62 * y / ymax
    pts = [(px(x), py(y)) for x, y in zip(xs, ys)]
    area = "M24,96 " + " ".join(f"L{a:.1f},{c:.1f}" for a, c in pts) + " L294,96 Z"
    shade = "M24,96 " + " ".join(f"L{a:.1f},{c:.1f}" for (a, c), x in zip(pts, xs) if x <= p75) + f" L{px(p75):.1f},96 Z"
    b = [f'<path class="dim" d="{area}"/><path class="acc" d="{shade}" opacity="0.6"/><line class="ink" x1="22" y1="96" x2="298" y2="96"/>']
    b.append(f'<line class="warns" x1="{px(p75):.1f}" y1="30" x2="{px(p75):.1f}" y2="96" stroke-dasharray="4 3"/><text class="tw" x="{px(p75):.1f}" y="24" text-anchor="middle">75% of reads</text>')
    b.append('<text class="tm" x="24" y="110">expected errors \u2192</text>')
    x2, x5, xl, xe = px(0.22), px(0.80), px(p75 + 0.12), px(0.97)
    b.append(f'<line class="accs" x1="{x2:.1f}" y1="128" x2="{xe:.1f}" y2="128" stroke-width="3"/><circle class="acc" cx="{x2:.1f}" cy="128" r="4"/><polygon class="acc" points="{xe:.1f},123 {xe + 9:.1f},128 {xe:.1f},133"/>')
    b.append(f'<line class="warns" x1="{x5:.1f}" y1="121" x2="{x5:.1f}" y2="135" stroke-dasharray="3 2"/>')
    b.append(f'<polygon class="warn" points="{xl - 5:.1f},116 {xl + 5:.1f},116 {xl:.1f},124"/><text class="tw" x="{xl:.1f}" y="112" text-anchor="middle">limit</text>')
    b.append(f'<text class="ta" x="{x2:.1f}" y="146" text-anchor="middle">2</text><text class="tw" x="{x5:.1f}" y="146" text-anchor="middle">5: a note</text>')
    b.append('<text class="tm" x="160" y="162" text-anchor="middle">limit = 75th percentile + 1, from 2 upward, no upper limit</text>')
    return _svg((320, 170), "Schematic: the distribution of expected errors, with the 75th percentile marked, a limit that starts at 2 and has no upper limit, and a note above 5", "".join(b))

def report():
    def doc(x, name, glyph):
        return (f'<rect class="card" x="{x}" y="10" width="74" height="86" rx="6"/><polygon class="dim" points="{x + 56},10 {x + 74},28 {x + 56},28"/>{glyph}'
                f'<text class="t" x="{x + 37}" y="114" text-anchor="middle" style="font:600 10.5px system-ui,sans-serif">{name}</text>')
    md = "".join(f'<line class="grid" x1="{14 + 12}" y1="{44 + i * 10}" x2="{14 + 12 + 46}" y2="{44 + i * 10}"/>' for i in range(4))
    js = '<text class="ta" x="160" y="66" text-anchor="middle" style="font:700 22px ui-monospace,Menlo,Consolas,monospace">{ }</text>'
    png = '<rect class="acc" x="244" y="56" width="8" height="26"/><rect class="acc" x="256" y="42" width="8" height="40" opacity=".7"/><rect class="acc" x="268" y="50" width="8" height="32" opacity=".5"/><line class="ink" x1="240" y1="82" x2="282" y2="82"/>'
    return _svg((320, 122), "Schematic: the three output files", doc(14, "REPORT.md", md) + doc(123, "summary.json", js) + doc(232, "summary.png", png))

def its():
    b = ['<text class="tm" x="20" y="16">one read, from the start</text>']
    b.append('<rect class="acc" x="20" y="26" width="170" height="16" rx="3"/><rect class="dim" x="190" y="26" width="100" height="16" rx="3"/>')
    b.append('<line class="warns" x1="190" y1="20" x2="190" y2="48" stroke-width="2"/><text class="tw" x="190" y="62" text-anchor="middle">first base of quality 2 or below</text>')
    b.append('<text class="t" x="20" y="92">kept if the cut read is</text>')
    for i, w in enumerate(["long enough", "free of N", "within maxEE"]):
        x = 4 + i * 106
        b.append(f'<rect class="card" x="{x}" y="100" width="100" height="24" rx="12"/><text class="ok" x="{x + 9}" y="116" style="font:700 12px system-ui,sans-serif">\u2713</text><text class="tm" x="{x + 24}" y="116">{w}</text>')
    return _svg((320, 134), "Schematic: for ITS, each read is cut at its first low-quality base and kept only if it is long enough, free of N and within the error limit", "".join(b))

