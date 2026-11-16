"""
build_compound_db.py -- build data/compounds.json from published data tables.

The compound database is generated, not typed by hand, so every number can be
traced to its source. To add compounds, add rows to data/compound_list.csv
(id, name, formula, CAS number, family) and run this script again:

    git clone --depth 1 https://github.com/CalebBell/chemicals
    git clone --depth 1 https://github.com/CalebBell/thermo
    python scripts/build_compound_db.py --chemicals chemicals --thermo thermo

Sources (all read from the data files shipped with the open-source
`chemicals` and `thermo` libraries by Caleb Bell and contributors):

  psrk4     Tc, Pc, acentric factor.  Horstmann, Jabłoniec, Krafczyk, Fischer,
            Gmehling, "PSRK group contribution equation of state: comprehensive
            revision and extension IV, including critical constants and
            alpha-function parameters for 1000 components", Fluid Phase
            Equilibria 227 (2005) 157-164.   [chemicals: Critical Properties/
            Appendix to PSRK Revision 4.tsv]
  poling    Antoine constants and their temperature range.  Poling, Prausnitz,
            O'Connell, "The Properties of Gases and Liquids", 5th ed.,
            McGraw-Hill, 2001, Appendix A; converted by `chemicals` to
            log10(P/Pa) = A - B/(T/K + C).   [chemicals: Vapor Pressure/
            Antoine Collection Poling.tsv]
  chemsep_nrtl   NRTL binary parameters ("DECHEMA NRTL at P = 1 atm").
            ChemSep interaction parameter library, H. Kooijman and R. Taylor,
            distributed under the Artistic License 2.0; `thermo` stores them
            as b_ij in kelvin with tau_ij = b_ij / T.   [thermo: Interaction
            Parameters/ChemSep/nrtl.json]
  chemsep_pr     Peng-Robinson binary interaction parameters k_ij ("DECHEMA
            Peng-Robinson parameters"), same library.   [thermo: Interaction
            Parameters/ChemSep/pr.json]
  iupac_aw  Molar mass from the formula and IUPAC abridged standard atomic
            weights (2021).
  derived   Normal boiling point: the temperature at which the Antoine
            equation gives 101 325 Pa.
  eppr78    Group-interaction parameters A_kl, B_kl (MPa) of the E-PPR78
            model, which predicts a temperature-dependent Peng-Robinson
            k_ij(T) from the groups in each molecule.  Jaubert, Qian, Lasala,
            Privat, Fluid Phase Equilibria 560 (2022) 113456.   [thermo:
            group_contribution/ppr78.py, table EPPR78_INTERACTIONS]. The
            groups of each compound are given in data/compound_list.csv.

A compound needs critical constants; Antoine constants are optional (a gas
such as CO2 has none in Poling's table and is then used with the equations of
state only).
"""

import argparse
import csv
import datetime
import json
import math
import os
import re
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ATOMIC_WEIGHT = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
                 "S": 32.06, "Cl": 35.45, "Br": 79.904}

SOURCES = {
    "psrk4": "Horstmann, S., Jabłoniec, A., Krafczyk, J., Fischer, K., Gmehling, J. "
             "PSRK group contribution equation of state: comprehensive revision and "
             "extension IV, including critical constants and alpha-function parameters "
             "for 1000 components. Fluid Phase Equilibria 227 (2005) 157-164. "
             "Read from the chemicals library data files.",
    "poling": "Poling, B. E., Prausnitz, J. M., O'Connell, J. P. The Properties of Gases "
              "and Liquids, 5th ed., McGraw-Hill, 2001, Appendix A. Read from the "
              "chemicals library data files, which store log10(P/Pa) = A - B/(T/K + C).",
    "chemsep_nrtl": "ChemSep interaction parameter data, file 'DECHEMA NRTL at P = 1 atm' "
                    "(copyright 2008 H. Kooijman and R. Taylor, Artistic License 2.0). "
                    "The file gives A12 and A21 in cal/mol with tau_ij = A_ij / RT; read "
                    "from the thermo library data files and converted to J/mol.",
    "chemsep_pr": "ChemSep interaction parameter data, file 'DECHEMA Peng-Robinson "
                  "Parameters' (copyright 2009 H. Kooijman and R. Taylor, Artistic "
                  "License 2.0). Read from the thermo library data files.",
    "iupac_aw": "Molar mass from the molecular formula and IUPAC abridged standard "
                "atomic weights (2021).",
    "derived": "Calculated here: temperature at which the Antoine equation gives 101 325 Pa.",
    "eppr78": "Jaubert, J.-N., Qian, J.-W., Lasala, S., Privat, R. The impressive impact of "
              "including enthalpy and heat capacity of mixing data when parameterising "
              "equations of state. Application to the development of the E-PPR78 "
              "(Enhanced-Predictive-Peng-Robinson-78) model. Fluid Phase Equilibria 560 "
              "(2022) 113456. Group-interaction parameters A_kl and B_kl (MPa) read from "
              "the thermo library (group_contribution/ppr78.py); the groups of each "
              "molecule are assigned in data/compound_list.csv.",
    "libraries": "Bell, C. and contributors. chemicals: Chemical properties component of "
                 "Chemical Engineering Design Library (ChEDL), and thermo: Chemical "
                 "properties component of ChEDL. https://github.com/CalebBell/chemicals, "
                 "https://github.com/CalebBell/thermo (MIT License).",
}


def tsv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def molar_mass(formula):
    total = 0.0
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        total += ATOMIC_WEIGHT[el] * (int(n) if n else 1)
    return round(total, 3)


def commit(path):
    try:
        return subprocess.check_output(["git", "-C", path, "log", "-1", "--format=%H %ad",
                                        "--date=short"], text=True).strip()
    except Exception:
        return "unknown"


def eppr78_table(thermo_path):
    """{(group_k, group_l): (A_kl, B_kl)} in MPa from the thermo source file."""
    src = open(os.path.join(thermo_path, "thermo", "group_contribution", "ppr78.py"),
               encoding="utf-8").read()
    pat = re.compile(r'^EPPR78_INTERACTIONS\[\(EPPR78_GROUPS\["([^"]+)"\], '
                     r'EPPR78_GROUPS\["([^"]+)"\]\)\] = \(([-0-9.eE]+), ([-0-9.eE]+)\)', re.M)
    out = {}
    for g1, g2, A, B in pat.findall(src):
        out[(g1, g2)] = (float(A), float(B))
    return out


def parse_groups(text):
    groups = {}
    for part in filter(None, (t.strip() for t in (text or "").split(";"))):
        g, n = part.split(":")
        groups[g.strip()] = int(n)
    return groups


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--chemicals", required=True, help="path to a checkout of CalebBell/chemicals")
    ap.add_argument("--thermo", required=True, help="path to a checkout of CalebBell/thermo")
    a = ap.parse_args()

    antoine = {r["CAS"]: r for r in tsv(os.path.join(
        a.chemicals, "chemicals", "Vapor Pressure", "Antoine Collection Poling.tsv"))}
    psrk = {r["CAS"]: r for r in tsv(os.path.join(
        a.chemicals, "chemicals", "Critical Properties", "Appendix to PSRK Revision 4.tsv"))}
    ip = os.path.join(a.thermo, "thermo", "Interaction Parameters", "ChemSep")
    nrtl = json.load(open(os.path.join(ip, "nrtl.json")))["data"]
    pr = json.load(open(os.path.join(ip, "pr.json")))["data"]

    with open(os.path.join(ROOT, "data", "compound_list.csv"), encoding="utf-8") as f:
        wanted = list(csv.DictReader(f))

    compounds, skipped = [], []
    for w in wanted:
        cas = w["cas"]
        if cas not in psrk:
            skipped.append(f"{w['id']} (no critical constants)")
            continue
        c, an = psrk[cas], antoine.get(cas)
        rec = {
            "id": w["id"], "name": w["name"], "formula": w["formula"], "cas": cas,
            "family": w["family"],
            "M": {"value": molar_mass(w["formula"]), "unit": "g/mol", "source": "iupac_aw"},
            "Tc": {"value": float(c["Tc"]), "unit": "K", "source": "psrk4"},
            "Pc": {"value": float(c["Pc"]), "unit": "Pa", "source": "psrk4"},
            "omega": {"value": float(c["omega"]), "unit": "", "source": "psrk4"},
            "antoine": None, "Tb": None,
        }
        if an is not None:
            A, B, C = float(an["A"]), float(an["B"]), float(an["C"])
            Tb = B / (A - math.log10(101325.0)) - C
            rec["antoine"] = {"A": A, "B": B, "C": C, "Tmin": float(an["Tmin"]),
                              "Tmax": float(an["Tmax"]), "P_unit": "Pa", "T_unit": "K",
                              "form": "log10(P/Pa) = A - B/(T/K + C)", "source": "poling"}
            rec["Tb"] = {"value": round(Tb, 2), "unit": "K", "source": "derived"}
        groups = parse_groups(w.get("ppr78_groups"))
        if groups:
            rec["ppr78"] = {"groups": groups, "source": "eppr78"}
        compounds.append(rec)

    # E-PPR78 group interactions, only for the groups the compounds use
    table = eppr78_table(a.thermo)
    used = sorted({g for c in compounds for g in c.get("ppr78", {}).get("groups", {})})
    unknown = [g for g in used if not any(g in k for k in table)]
    if unknown:
        raise SystemExit(f"E-PPR78 groups not in the table: {unknown}")
    ppr78 = {"groups": used, "unit": "MPa", "source": "eppr78",
             "form": "kij(T) = [-1/2 sum_k sum_l (a_ik - a_jk)(a_il - a_jl) A_kl (298.15/T)^(B_kl/A_kl - 1)"
                     " - (sqrt(a_i)/b_i - sqrt(a_j)/b_j)^2] / [2 sqrt(a_i a_j)/(b_i b_j)]",
             "AB": {}, "missing": []}
    for i, g1 in enumerate(used):
        for g2 in used[i + 1:]:
            v = table.get((g1, g2)) or table.get((g2, g1))
            if v is None:                      # not fitted by the authors
                ppr78["missing"].append(f"{g1}|{g2}")
                v = (0.0, 0.0)
            ppr78["AB"][f"{g1}|{g2}"] = list(v)

    pairs = []
    for i, ci in enumerate(compounds):
        for cj in compounds[i + 1:]:
            entry = {"a": ci["id"], "b": cj["id"]}
            kab, kba = f"{ci['cas']} {cj['cas']}", f"{cj['cas']} {ci['cas']}"
            if kab in nrtl and kba in nrtl:
                # b_ab: tau_ab = b_ab / T, with tau_ab = (g_ab - g_bb) / RT
                entry["nrtl"] = {"b_ab": nrtl[kab]["bij"], "b_ba": nrtl[kba]["bij"],
                                 "alpha": nrtl[kab]["alphaij"], "unit": "K",
                                 "source": "chemsep_nrtl"}
            k = pr.get(kab) or pr.get(kba)
            if k is not None:
                entry["kij_pr"] = {"value": k["kij"], "source": "chemsep_pr"}
            if len(entry) > 2:
                pairs.append(entry)

    out = {
        "meta": {
            "generated": datetime.date.today().isoformat(),
            "generator": "scripts/build_compound_db.py",
            "upstream": {"chemicals": commit(a.chemicals), "thermo": commit(a.thermo)},
            "sources": SOURCES,
            "notes": "Values are reproduced from the cited sources for teaching use. "
                     "Check them against the original publication before relying on them.",
        },
        "compounds": compounds,
        "pairs": pairs,
        "ppr78": ppr78,
    }
    path = os.path.join(ROOT, "data", "compounds.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    n_nrtl = sum(1 for p in pairs if "nrtl" in p)
    n_pr = sum(1 for p in pairs if "kij_pr" in p)
    print(f"wrote {path}: {len(compounds)} compounds, {n_nrtl} NRTL pairs, {n_pr} PR kij pairs")
    if skipped:
        print("skipped:", "; ".join(skipped))


if __name__ == "__main__":
    main()
