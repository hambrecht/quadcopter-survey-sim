"""Colour-vision check of the design palette (port of the OKLab / Machado-2009 CVD check).
Prints the worst pairwise OKLab Delta E x100 under normal, protan, deutan and tritan vision,
and contrast against white. Targets: CVD >= 8 (6-8 only with secondary encoding), normal >= 15."""
import itertools, math, sys
M = {"protan": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
     "deutan": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
     "tritan": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]]}
def lin(h):
    h = h.lstrip("#"); c = [int(h[i:i+2], 16) / 255 for i in (0, 2, 4)]
    return [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in c]
def oklab(rgb):
    r, g, b = rgb
    l = (0.4122214708*r + 0.5363325363*g + 0.0514459929*b) ** (1/3)
    m = (0.2119034982*r + 0.6806995451*g + 0.1073969566*b) ** (1/3)
    s = (0.0883024619*r + 0.2817188376*g + 0.6299787005*b) ** (1/3)
    return [0.2104542553*l + 0.7936177850*m - 0.0040720468*s, 1.9779984951*l - 2.4285922050*m + 0.4505937099*s,
            0.0259040371*l + 0.7827717662*m - 0.8086757660*s]
def sim(h, k):
    r = lin(h); return [min(1, max(0, sum(M[k][i][j] * r[j] for j in range(3)))) for i in range(3)]
def dE(a, b, k=None):
    A = oklab(sim(a, k) if k else lin(a)); B = oklab(sim(b, k) if k else lin(b))
    return 100 * math.dist(A, B)
def contrast(h, bg="#FFFFFF"):
    Y = lambda x: sum(w * c for w, c in zip((0.2126, 0.7152, 0.0722), lin(x)))
    hi, lo = sorted([Y(h), Y(bg)], reverse=True); return (hi + 0.05) / (lo + 0.05)
def report(pal, names):
    for k in (None, "protan", "deutan", "tritan"):
        worst = min(((dE(pal[i], pal[j], k), names[i], names[j]) for i, j in itertools.combinations(range(len(pal)), 2)))
        print(f"{k or 'normal':7s} worst dE {worst[0]:5.1f}  ({worst[1]} vs {worst[2]})")
    for n, h in zip(names, pal): print(f"  {n:7s} {h}  L={oklab(lin(h))[0]:.2f}  contrast vs white {contrast(h):.2f}:1")
if __name__ == "__main__":
    from style import DESIGN_COLORS
    report(list(DESIGN_COLORS.values()), list(DESIGN_COLORS))
