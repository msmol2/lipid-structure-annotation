"""Assign one plausible molecular structure (SMILES) to each lipid annotation.

Rule-based: split the sum composition (e.g. PC 34:1) into common mammalian fatty
acyl chains, place double bonds at their usual positions (cis), and assemble the
SMILES from headgroup templates. Also computes formula / monoisotopic mass with a
small SMILES parser and compares the expected adduct m/z with the reported m/z.
"""
import csv, itertools, re, sys
from collections import Counter

# ---------------- fatty acids: weight (abundance prior) + double-bond positions -------------
FA = {  # "n:k": (weight, [delta positions, all cis])
 "10:0":(0.5,[]), "12:0":(1,[]), "13:0":(0.3,[]), "14:0":(4,[]), "14:1":(1,[9]), "15:0":(3,[]),
 "16:0":(10,[]), "16:1":(5,[9]), "16:2":(0.5,[9,12]), "17:0":(3,[]), "17:1":(1,[10]),
 "18:0":(8,[]), "18:1":(9,[9]), "18:2":(9,[9,12]), "18:3":(4,[9,12,15]), "18:4":(1,[6,9,12,15]),
 "19:0":(1,[]), "20:0":(2,[]), "20:1":(3,[11]), "20:2":(2,[11,14]), "20:3":(4,[8,11,14]),
 "20:4":(7,[5,8,11,14]), "20:5":(4,[5,8,11,14,17]), "21:0":(0.5,[]), "22:0":(2,[]), "22:1":(1,[13]),
 "22:2":(1,[13,16]), "22:3":(1,[13,16,19]), "22:4":(3,[7,10,13,16]), "22:5":(4,[7,10,13,16,19]),
 "22:6":(6,[4,7,10,13,16,19]), "23:0":(1,[]), "24:0":(2,[]), "24:1":(2,[15]), "25:0":(0.5,[]),
 "26:0":(0.5,[]), "26:1":(0.3,[17]),
}
# sphingolipid N-acyl priors (very-long saturated chains dominate in plasma SM)
SM_ACYL = {"16:0":10,"24:1":8,"24:0":8,"22:0":6,"23:0":5,"18:0":5,"20:0":3,"14:0":2,"15:0":2,
           "17:0":2,"25:0":1,"26:0":1,"26:1":1,"22:1":2,"18:1":2,"20:1":1,"23:1":1,"25:1":1}
SPHINGOID = {"d18:1":10,"d18:2":2,"d18:0":2,"d16:1":1,"d17:1":1}

def nk(s):
    n, k = s.split(":"); return int(n), int(k)

def positions(n, k):
    key = f"{n}:{k}"
    if key in FA: return FA[key][1]
    if k == 0: return []
    # generic: methylene-interrupted cis bonds counted back from the omega-3 end
    p = [n - 3 - 3 * i for i in range(k)]
    return sorted(p) if min(p) >= 3 else None

def weight(n, k, table=None):
    if table is not None:
        if f"{n}:{k}" in table: return table[f"{n}:{k}"]
    if f"{n}:{k}" in FA: return FA[f"{n}:{k}"][0]
    return 0.01 if positions(n, k) is not None and n >= 2 else 0  # fallback: odd/rare chain

def chain_atoms(n, starts, first=2):
    """SMILES for carbons first..n, cis double bonds starting at the given carbons."""
    out = []
    for i in range(first, n + 1):
        if i in starts: out.append("/C=")
        elif i - 1 in starts: out.append("C\\")
        else: out.append("C")
    return "".join(out)

def acyl_fwd(n, k):          # C2..Cn, written from the carbonyl outward
    return chain_atoms(n, positions(n, k))
def acyl_rev(n, k):          # written methyl-end first, ending next to the carbonyl
    return acyl_fwd(n, k)[::-1]
def alkenyl_rev(n, k):       # plasmalogen sn-1 O-(1Z)-alkenyl, k excludes the vinyl bond
    extra = positions(n, k) or []
    return chain_atoms(n, [1] + extra, first=1)[::-1]

def best_split(total, nchains, table=None, sat_first=True, forbid=()):
    n_tot, k_tot = nk(total)
    names = list(FA) if table is None else sorted(set(FA) | set(table))
    cands = [(nk(x)) for x in names]
    best = None
    def consider(combo):
        nonlocal best
        w = 1.0
        for n, k in combo: w *= weight(n, k, table)
        if w > 0 and (best is None or w > best[0]): best = (w, combo)
    for combo in itertools.combinations_with_replacement(cands, nchains - 1):
        n_last = n_tot - sum(c[0] for c in combo); k_last = k_tot - sum(c[1] for c in combo)
        if n_last < 2 or k_last < 0 or k_last > (n_last - 4) // 3 + 1: continue
        if positions(n_last, k_last) is None: continue
        consider(tuple(combo) + ((n_last, k_last),))
    if best is None: return None
    chains = sorted(best[1], key=lambda c: (c[1], c[0]))  # most saturated first
    return chains

def fa_str(c): return f"{c[0]}:{c[1]}"

GPC = "COP(=O)([O-])OCC[N+](C)(C)C"
HEAD = {"PC": GPC, "PE": "COP(=O)(O)OCCN", "PS": "COP(=O)(O)OC[C@H](N)C(=O)O",
        "PI": "COP(=O)(O)OC1C(O)C(O)C(O)C(O)C1O", "PG": "COP(=O)(O)OCC(O)CO", "PA": "COP(=O)(O)O"}
CHOL = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"

def sphingoid_tail(sp):
    n = int(sp[1:3]); db = sp.split(":")[1]
    if db == "0": return "C" * (n - 3)
    if db == "1": return "/C=C/" + "C" * (n - 5)
    return "/C=C/" + "C" * 8 + "/C=C\\" + "C" * (n - 15)   # d18:2 (4E,14Z)

def build(cls, species):
    tot = species.split(" ")[1]
    n_tot, k_tot = nk(tot)
    if cls in ("PC", "PE", "PS", "PI", "PG", "PA", "DAG"):
        c1, c2 = best_split(tot, 2)
        head = "CO" if cls == "DAG" else HEAD[cls]
        smi = f"{acyl_rev(*c1)}C(=O)OC[C@H]({head})OC(=O){acyl_fwd(*c2)}"
        label = f"{cls if cls!='DAG' else 'DG'} {fa_str(c1)}/{fa_str(c2)}/0:0" if cls == "DAG" else f"{cls} {fa_str(c1)}/{fa_str(c2)}"
        return label, smi, [c1, c2]
    if cls == "TAG":
        ch = best_split(tot, 3)
        sat, mid, last = ch[0], ch[-1], ch[1]  # most unsaturated chain at sn-2
        smi = f"{acyl_rev(*sat)}C(=O)OCC(OC(=O){acyl_fwd(*mid)})COC(=O){acyl_fwd(*last)}"
        return f"TG {fa_str(sat)}/{fa_str(mid)}/{fa_str(last)}", smi, [sat, mid, last]
    if cls in ("LPC", "LPE"):
        c = (n_tot, k_tot); head = GPC if cls == "LPC" else HEAD["PE"]
        return f"{cls} {tot}/0:0", f"{acyl_rev(*c)}C(=O)OC[C@H]({head})O", [c]
    if cls == "MAG":
        c = (n_tot, k_tot)
        return f"MG {tot}/0:0/0:0", f"{acyl_rev(*c)}C(=O)OCC(O)CO", [c]
    if cls == "CE":
        c = (n_tot, k_tot)
        return f"CE {tot}", CHOL.replace("C4)O)C)C", f"C4)OC(=O){acyl_fwd(*c)})C)C"), [c]
    if cls == "FFA":
        c = (n_tot, k_tot)
        return f"FA {tot}", f"OC(=O){acyl_fwd(*c)}", [c]
    if cls == "pPE":
        best = None
        for p in ("16:0", "18:0", "18:1"):
            pn, pk = nk(p); rn, rk = n_tot - pn, k_tot - pk
            if rk < 0 or positions(rn, rk) is None: continue
            w = {"16:0": 4, "18:0": 5, "18:1": 3}[p] * weight(rn, rk)
            if best is None or w > best[0]: best = (w, (pn, pk), (rn, rk))
        _, a, c2 = best
        smi = f"{alkenyl_rev(*a)}OC[C@H]({HEAD['PE']})OC(=O){acyl_fwd(*c2)}"
        return f"PE P-{fa_str(a)}/{fa_str(c2)}", smi, [a, c2]
    if cls == "pPC":   # annotated 'pPC 18:0': lyso-plasmalogen PC
        a = (n_tot, k_tot)
        return f"LPC P-{tot}", f"{alkenyl_rev(*a)}OC[C@H]({GPC})O", [a]
    if cls in ("SM", "CerP"):
        best = None
        for sp, ws in SPHINGOID.items():
            sn, sk = int(sp[1:3]), int(sp[-1]); rn, rk = n_tot - sn, k_tot - sk
            if rn < 2 or rk < 0 or positions(rn, rk) is None: continue
            w = ws * (SM_ACYL.get(f"{rn}:{rk}", 0) or weight(rn, rk) * 0.1)
            if w > 0 and (best is None or w > best[0]): best = (w, sp, (rn, rk))
        _, sp, c = best
        head = GPC if cls == "SM" else "COP(=O)(O)O"
        smi = f"{acyl_rev(*c)}C(=O)N[C@@H]({head})[C@@H]({sphingoid_tail(sp)})O"
        lab = "SM" if cls == "SM" else "CerP"
        return f"{lab} {sp}/{fa_str(c)}", smi, [c]
    raise ValueError(cls)

# ---------------- tiny SMILES -> formula / mass (independent of RDKit) ----------------------
MASS = {"C": 12.0, "H": 1.00782503223, "N": 14.00307400443, "O": 15.99491461957, "P": 30.97376199842}
VAL = {"C": [4], "N": [3], "O": [2], "P": [3, 5]}
def formula(smi):
    tok = re.findall(r"\[[^\]]+\]|Cl|Br|[CNOP]|\(|\)|=|#|/|\\|\d", smi)
    atoms, bonds, stack, prev, order, rings = [], [], [], None, 1, {}
    for t in tok:
        if t == "(": stack.append(prev); continue
        if t == ")": prev = stack.pop(); continue
        if t in "/\\": continue
        if t == "=": order = 2; continue
        if t.isdigit():
            if t in rings:
                a, o = rings.pop(t); bonds[a] += max(o, order); bonds[prev] += max(o, order)
            else: rings[t] = (prev, order)
            order = 1; continue
        if t.startswith("["):
            m = re.match(r"\[([A-Z][a-z]?)[@]*(H\d?)?([+-])?\]", t)
            el, h, ch = m.group(1), m.group(2), m.group(3)
            atoms.append([el, int(h[1:] or 1) if h else 0, ch]); bonds.append(0)
        else:
            atoms.append([t, None, None]); bonds.append(0)
        i = len(atoms) - 1
        if prev is not None: bonds[prev] += order; bonds[i] += order
        prev, order = i, 1
    cnt, charge = Counter(), 0
    for (el, h, ch), b in zip(atoms, bonds):
        cnt[el] += 1
        if h is None:
            v = next((v for v in VAL[el] if v >= b), None)
            if v is None: raise ValueError(f"valence error on {el} in {smi}")
            h = v - b
        if ch: charge += 1 if ch == "+" else -1
        cnt["H"] += h
    assert not rings and not stack, "unclosed ring/branch"
    mono = sum(MASS[e] * n for e, n in cnt.items())
    f = "".join(f"{e}{cnt[e] if cnt[e] > 1 else ''}" for e in ["C", "H", "N", "O", "P"] if cnt[e])
    return f, mono, charge

ADDUCT = {"[M+H]": 1.007276, "[M]": 1.007276, "[M+Na]": 22.989218, "[M+NH4]": 18.033823,
          "[M+Li]": 7.015455, "[M-H]": -1.007276, "[M-Ac-H]": 59.013305, "": -1.007276}

if __name__ == "__main__":
    src, out = sys.argv[1], sys.argv[2]
    rows = list(csv.DictReader(open(src)))
    species = {}
    for r in rows:
        species.setdefault((r["lipid_class"], r["lipid_species"]), []).append(r)
    res = []
    for (cls, sp), rs in species.items():
        label, smi, chains = build(cls, sp)
        f, mono, charge = formula(smi)
        n_tot, k_tot = nk(sp.split(" ")[1])
        if cls in ("SM", "CerP"):
            spb = label.split(" ")[1].split("/")[0]; nsum = int(spb[1:3]) + chains[0][0]; ksum = int(spb[-1]) + chains[0][1]
        else:
            nsum, ksum = sum(c[0] for c in chains), sum(c[1] for c in chains)
        ppms = []
        for r in rs:
            calc = mono + ADDUCT[r["adduct"]]
            ppms.append((r["adduct"] or "[M-H]", float(r["mz"]), calc, (float(r["mz"]) - calc) / calc * 1e6))
        worst = max(ppms, key=lambda x: abs(x[3]))
        res.append({
            "lipid_annotation": sp, "lipid_class": cls, "plausible_molecular_species": label,
            "chain_sum_check": f"{nsum}:{ksum} {'OK' if (nsum, ksum) == (n_tot, k_tot) else 'MISMATCH'}",
            "smiles": smi, "formula": f, "monoisotopic_mass": round(mono, 4),
            "adducts_reported": "; ".join(f"{a} {m:.4f}" for a, m, _, _ in ppms),
            "max_abs_mass_error_ppm": round(abs(worst[3]), 1),
            "structure_source": "Rule-based construction using LIPID MAPS shorthand; common mammalian acyl chains, cis double bonds at typical positions",
            "note": "One plausible molecular structure consistent with the lipid annotation; not the uniquely identified structure",
        })
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
    print(len(rows), "rows ->", len(res), "unique annotations")
    print("chain-sum mismatches:", sum("MISMATCH" in r["chain_sum_check"] for r in res))
    bad = [r for r in res if r["max_abs_mass_error_ppm"] > 10]
    print("mass error >10 ppm:", len(bad))
    for r in bad: print("  ", r["lipid_annotation"], r["plausible_molecular_species"], r["adducts_reported"], r["max_abs_mass_error_ppm"])
