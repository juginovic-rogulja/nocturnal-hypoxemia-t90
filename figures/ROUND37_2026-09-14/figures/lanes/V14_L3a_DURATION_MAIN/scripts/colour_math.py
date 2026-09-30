"""sRGB <-> CIELAB (D65, 2 degree), CIEDE2000, and WCAG 2 contrast. Pure Python, no packages."""
import math

def hex_to_rgb(h):
    h = h.lstrip("#"); return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))

def rgb_to_hex(rgb):
    return "#%02x%02x%02x" % tuple(int(round(max(0.0, min(1.0, v)) * 255)) for v in rgb)

def _lin(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

def _gam(c):
    return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055

XN, YN, ZN = 0.95047, 1.0, 1.08883

def rgb_to_xyz(rgb):
    r, g, b = (_lin(v) for v in rgb)
    return (0.4124564 * r + 0.3575761 * g + 0.1804375 * b,
            0.2126729 * r + 0.7151522 * g + 0.0721750 * b,
            0.0193339 * r + 0.1191920 * g + 0.9503041 * b)

def xyz_to_rgb(xyz):
    x, y, z = xyz
    r = 3.2404542 * x - 1.5371385 * y - 0.4985314 * z
    g = -0.9692660 * x + 1.8760108 * y + 0.0415560 * z
    b = 0.0556434 * x - 0.2040259 * y + 1.0572252 * z
    return tuple(_gam(v) for v in (r, g, b))

def _f(t):
    return t ** (1 / 3) if t > (6 / 29) ** 3 else t / (3 * (6 / 29) ** 2) + 4 / 29

def _finv(t):
    return t ** 3 if t > 6 / 29 else 3 * (6 / 29) ** 2 * (t - 4 / 29)

def rgb_to_lab(rgb):
    x, y, z = rgb_to_xyz(rgb)
    fx, fy, fz = _f(x / XN), _f(y / YN), _f(z / ZN)
    return (116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz))

def lab_to_rgb(lab):
    L, a, b = lab
    fy = (L + 16) / 116; fx = fy + a / 500; fz = fy - b / 200
    return xyz_to_rgb((XN * _finv(fx), YN * _finv(fy), ZN * _finv(fz)))

def in_gamut(rgb, tol=1e-6):
    return all(-tol <= v <= 1 + tol for v in rgb)

def ciede2000(lab1, lab2):
    L1, a1, b1 = lab1; L2, a2, b2 = lab2
    C1 = math.hypot(a1, b1); C2 = math.hypot(a2, b2); Cm = (C1 + C2) / 2
    G = 0.5 * (1 - math.sqrt(Cm ** 7 / (Cm ** 7 + 25 ** 7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = math.hypot(a1p, b1), math.hypot(a2p, b2)
    def hp(a, b):
        if a == 0 and b == 0: return 0.0
        h = math.degrees(math.atan2(b, a)); return h + 360 if h < 0 else h
    h1p, h2p = hp(a1p, b1), hp(a2p, b2)
    dLp = L2 - L1; dCp = C2p - C1p
    if C1p * C2p == 0: dhp = 0.0
    elif abs(h2p - h1p) <= 180: dhp = h2p - h1p
    elif h2p - h1p > 180: dhp = h2p - h1p - 360
    else: dhp = h2p - h1p + 360
    dHp = 2 * math.sqrt(C1p * C2p) * math.sin(math.radians(dhp / 2))
    Lm = (L1 + L2) / 2; Cmp = (C1p + C2p) / 2
    if C1p * C2p == 0: hmp = h1p + h2p
    elif abs(h1p - h2p) <= 180: hmp = (h1p + h2p) / 2
    elif h1p + h2p < 360: hmp = (h1p + h2p + 360) / 2
    else: hmp = (h1p + h2p - 360) / 2
    T = (1 - 0.17 * math.cos(math.radians(hmp - 30)) + 0.24 * math.cos(math.radians(2 * hmp))
         + 0.32 * math.cos(math.radians(3 * hmp + 6)) - 0.20 * math.cos(math.radians(4 * hmp - 63)))
    dth = 30 * math.exp(-((hmp - 275) / 25) ** 2)
    RC = 2 * math.sqrt(Cmp ** 7 / (Cmp ** 7 + 25 ** 7))
    SL = 1 + 0.015 * (Lm - 50) ** 2 / math.sqrt(20 + (Lm - 50) ** 2)
    SC = 1 + 0.045 * Cmp; SH = 1 + 0.015 * Cmp * T
    RT = -math.sin(math.radians(2 * dth)) * RC
    return math.sqrt((dLp / SL) ** 2 + (dCp / SC) ** 2 + (dHp / SH) ** 2 + RT * (dCp / SC) * (dHp / SH))

def luminance(rgb):
    r, g, b = (_lin(v) for v in rgb)   # WCAG uses the 0.03928 knee; identical to 4 decimals
    return 0.2126 * r + 0.7152 * g + 0.0722 * b

def contrast(rgb1, rgb2):
    l1, l2 = luminance(rgb1), luminance(rgb2)
    hi, lo = max(l1, l2), min(l1, l2)
    return (hi + 0.05) / (lo + 0.05)

def ladder_lab_line(top_hex, n_rungs):
    """Equal CIELAB lightness steps from white toward top_hex: rung k (1..n) sits at fraction
    k/n along the straight Lab line from white to the top, so the top rung is exactly top_hex."""
    white = rgb_to_lab((1.0, 1.0, 1.0)); top = rgb_to_lab(hex_to_rgb(top_hex))
    out = []
    for k in range(1, n_rungs + 1):
        t = k / n_rungs
        lab = tuple(w + t * (p - w) for w, p in zip(white, top))
        rgb = lab_to_rgb(lab)
        out.append((rgb_to_hex(rgb) if k < n_rungs else top_hex.lower(), lab, in_gamut(rgb)))
    return out

if __name__ == "__main__":
    TOP = "#39c445"; INK = "#1a1d20"
    print("top rung", TOP, "Lab", tuple(round(v, 2) for v in rgb_to_lab(hex_to_rgb(TOP))))
    print("\nA. straight CIELAB line white -> top, 4 rungs:")
    lad = ladder_lab_line(TOP, 4)
    for i, (h, lab, ok) in enumerate(lad, 1):
        rr = rgb_to_lab(hex_to_rgb(h))
        print(f"  rung {i}  {h}  target Lab ({lab[0]:.2f},{lab[1]:.2f},{lab[2]:.2f})  realised Lab ({rr[0]:.2f},{rr[1]:.2f},{rr[2]:.2f})  in gamut {ok}"
              f"   ink {contrast(hex_to_rgb(h), hex_to_rgb(INK)):.2f}:1  white {contrast(hex_to_rgb(h), (1,1,1)):.2f}:1")
    hexes = ["#ffffff"] + [h for h, _, _ in lad]
    for a, b in zip(hexes, hexes[1:]):
        print(f"  dE00 {a} -> {b}: {ciede2000(rgb_to_lab(hex_to_rgb(a)), rgb_to_lab(hex_to_rgb(b))):.2f}")
    print("\nB. sRGB white mix, mix fraction solved so L matches the same equal steps:")
    white = (1.0, 1.0, 1.0); top = hex_to_rgb(TOP); Lt = rgb_to_lab(top)[0]
    hexesB = ["#ffffff"]
    for k in (1, 2, 3):
        Ltarget = 100 + k / 4 * (Lt - 100)
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2; rgb = tuple(w + mid * (p - w) for w, p in zip(white, top))
            if rgb_to_lab(rgb)[0] > Ltarget: lo = mid
            else: hi = mid
        rgb = tuple(w + lo * (p - w) for w, p in zip(white, top)); h = rgb_to_hex(rgb)
        lab = rgb_to_lab(hex_to_rgb(h)); hexesB.append(h)
        print(f"  rung {k}  {h}  mix {lo*100:.1f}% of top  Lab ({lab[0]:.2f},{lab[1]:.2f},{lab[2]:.2f})  ink {contrast(hex_to_rgb(h), hex_to_rgb(INK)):.2f}:1")
    hexesB.append(TOP)
    for a, b in zip(hexesB, hexesB[1:]):
        print(f"  dE00 {a} -> {b}: {ciede2000(rgb_to_lab(hex_to_rgb(a)), rgb_to_lab(hex_to_rgb(b))):.2f}")
    print("\nold ladder for reference:")
    old = ["#ffffff", "#c8eecb", "#92df99", "#5dcf66", "#298d32"]
    for h in old[1:]:
        lab = rgb_to_lab(hex_to_rgb(h)); print(f"  {h} Lab ({lab[0]:.2f},{lab[1]:.2f},{lab[2]:.2f})  ink {contrast(hex_to_rgb(h), hex_to_rgb(INK)):.2f}:1  white {contrast(hex_to_rgb(h), (1,1,1)):.2f}:1")
    for a, b in zip(old, old[1:]):
        print(f"  dE00 {a} -> {b}: {ciede2000(rgb_to_lab(hex_to_rgb(a)), rgb_to_lab(hex_to_rgb(b))):.2f}")
    print("\nblue headers, for the record: #1a1d21 on #0288d1 = %.2f:1, on #3f9fd8 = %.2f:1, on #7cc0e9 = %.2f:1, on #b3dcf2 = %.2f:1" % tuple(
        contrast(hex_to_rgb("#1a1d21"), hex_to_rgb(h)) for h in ("#0288d1", "#3f9fd8", "#7cc0e9", "#b3dcf2")))
    print("sanity: dE00 of the Sharma pair (50,2.6772,-79.7751)-(50,0,-82.7485) = %.4f (expected 2.0425)" % ciede2000((50,2.6772,-79.7751),(50,0,-82.7485)))
    print("sanity: dE00 pair (50,2.5,0)-(73,25,-18) = %.4f (expected 27.1492)" % ciede2000((50,2.5,0),(73,25,-18)))
