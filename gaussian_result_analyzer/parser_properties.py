"""Charges, dipole, spin, dispersion, post-HF energies and NBO from a Gaussian log."""

import re
import logging

from .parser_structure import _route_text

#: Debye per atomic unit of dipole moment
AU_TO_DEBYE = 2.541746473

_ATOM_CHARGE = re.compile(
    r"^\s*(\d+)\s+([A-Za-z]{1,3})\s+(-?\d+\.\d+)(?:\s+(-?\d+\.\d+))?"
)
_DIPOLE = re.compile(
    r"X=\s*(-?\d+\.\d+)\s+Y=\s*(-?\d+\.\d+)\s+Z=\s*(-?\d+\.\d+)\s+Tot=\s*(-?\d+\.\d+)"
)
_S2 = re.compile(r"^\s*<Sx>=.*<S\*\*2>=\s*(-?\d+\.\d+)")
_DISPERSION = re.compile(r"Dispersion energy=\s*(-?\d+\.\d+)\s+Hartrees")
_NBO_ROW = re.compile(
    r"^\s*(\d+)\.\s+((?:BD|CR|LP|LV|RY|3C|3Cn)\*?\s*\(\s*\d+\))\s*(.+?)\s+"
    r"(\d+\.\d+)\s+(-?\d+\.\d+)"
)
_E2_ROW = re.compile(
    r"^\s*(\d+)\.\s+(.+?)\s*/\s*(\d+)\.\s+(.+?)\s+"
    r"(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$"
)
# NBO 7 prints the same E(2) table without the "/" between donor and acceptor.
_E2_ROW_NBO7 = re.compile(
    r"^\s*(\d+)\.\s+(\S.*?\S)\s+(\d+)\.\s+(\S.*?\S)\s+"
    r"(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
)
_NBO_TOTAL = re.compile(
    r"^\s*(Total Lewis|Valence non-Lewis|Rydberg non-Lewis|Total non-Lewis)"
    r"\s+(\d+\.\d+)\s+\(\s*(\d+\.\d+)%"
)
_NPA_ROW = re.compile(
    r"^\s*([A-Z][a-z]?)\s+(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+"
    r"(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$"
)
_CONFIG_ORBITAL = re.compile(r"(\d+[a-z][a-z0-9]*)\(\s*(\d+\.\d+)\)", re.IGNORECASE)
_NAO_ROW = re.compile(
    r"^\s*(\d+)\s+([A-Z][a-z]?)\s+(\d+)\s+(\S+)\s+"
    r"(Cor|Val|Ryd)\(\s*([^)]+)\)\s+(\d+\.\d+)(?:\s+(-?\d+\.\d+))?\s*$"
)
_BENDING_ROW = re.compile(
    r"^\s*(\d+)\.\s+(.+?)\s{2,}"
    r"((?:--|-?\d+\.\d+)(?:\s+(?:--|-?\d+\.\d+)){7})\s*$"
)
_D_FLOAT = r"(-?\d*\.\d+(?:[DdEe][-+]?\d+)?)"


def _float_d(text):
    """Float from Fortran notation, which writes the exponent with a D."""
    return float(text.replace("D", "E").replace("d", "e"))


class _PropertyParsingMixin:
    """Parsing of scalar and per-atom properties."""

    def _last_index(self, predicate):
        """Index of the last line satisfying *predicate*, or None."""
        for i in range(len(self.lines) - 1, -1, -1):
            if predicate(self.lines[i]):
                return i
        return None

    def parse_dipole(self):
        """The last dipole moment, printed in Debye in the standard orientation."""
        self.data["dipoles"] = None
        self.data["dipole"] = None
        i = self._last_index(lambda l: "Dipole moment (field-independent basis" in l)
        if i is None or i + 1 >= len(self.lines):
            return
        match = _DIPOLE.search(self.lines[i + 1])
        if not match:
            return
        x, y, z, total = (float(match.group(k)) for k in (1, 2, 3, 4))
        self.data["dipoles"] = {
            "vector_au": (x / AU_TO_DEBYE, y / AU_TO_DEBYE, z / AU_TO_DEBYE),
            "vector_debye": (x, y, z),
            "magnitude_au": total / AU_TO_DEBYE,
            "magnitude_debye": total,
            "vector": (x, y, z),
            "magnitude": total,
        }
        self.data["dipole"] = self.data["dipoles"]

    def parse_spin_contamination(self):
        """<S**2> of the final wavefunction against the ideal S(S+1)."""
        self.data["spin_s2"] = None
        actual = None
        for line in self.lines:
            match = _S2.search(line)
            if match:
                actual = float(match.group(1))
        if actual is None:
            return
        s = (self.data.get("mult", 1) - 1) / 2.0
        ideal = s * (s + 1)
        self.data["spin_s2"] = {
            "actual": actual,
            "ideal": ideal,
            "contamination": actual - ideal,
        }

    def parse_dispersion(self):
        """The empirical dispersion energy, when one was requested."""
        self.data["dispersion"] = None
        for line in self.lines:
            match = _DISPERSION.search(line)
            if match:
                self.data["dispersion"] = float(match.group(1))

    def parse_energy_components(self):
        """Post-HF energies, in the order the methods build on each other."""
        found = {}
        patterns = (
            ("MP2 Correlation", re.compile(r"\bE2\s*=\s*" + _D_FLOAT), False),
            ("MP2 Total", re.compile(r"EUMP2\s*=\s*" + _D_FLOAT), False),
            (
                "Double-Hybrid Correlation",
                re.compile(r"\bE2\(\w+\)\s*=\s*" + _D_FLOAT),
                False,
            ),
            (
                "Double-Hybrid Total",
                re.compile(r"\bE2\(\w+\)\s*=.*\bE\(\w+\)\s*=\s*" + _D_FLOAT),
                False,
            ),
            (
                "CCSD Total",
                re.compile(r"Wavefunction amplitudes converged\. E\(Corr\)=\s*" + _D_FLOAT),
                False,
            ),
            ("CCSD(T) Total", re.compile(r"CCSD\(T\)=\s*" + _D_FLOAT), False),
            ("T1 Diagnostic", re.compile(r"T1 Diagnostic\s*=\s*" + _D_FLOAT), True),
            (
                "Counterpoise Corrected",
                re.compile(r"Counterpoise corrected energy\s*=\s*" + _D_FLOAT),
                False,
            ),
            ("BSSE", re.compile(r"BSSE energy\s*=\s*" + _D_FLOAT), False),
        )
        for line in self.lines:
            for label, pattern, dimensionless in patterns:
                match = pattern.search(line)
                if match:
                    found[label] = (_float_d(match.group(1)), dimensionless)
        self.data["energy_components"] = [
            {"label": label, "value": found[label][0], "dimensionless": found[label][1]}
            for label, _p, _d in patterns
            if label in found
        ]

    def parse_nbo_orbitals(self):
        """The Natural Bond Orbitals summary of the last NBO analysis.

        Both NBO generations Gaussian can drive are supported: the built-in
        NBO 3.1 ("Natural Bond Orbitals (Summary):") and a linked NBO 7
        ("NATURAL BOND ORBITALS (Summary):"), whose Rydberg rows are printed
        as RY rather than RY*. If a job prints the summary more than once
        (spin sections, optimization steps), the last one wins, matching the
        ORCA analyzer.
        """
        orbitals = []
        start = self._last_index(
            lambda line: "natural bond orbitals (summary)" in line.lower()
        )
        if start is not None:
            for row in self.lines[start + 1 : start + 5000]:
                if row.strip().startswith("Total Lewis") or (
                    "NBO analysis completed" in row
                ):
                    break
                match = _NBO_ROW.match(row)
                if match:
                    atoms = match.group(3).strip()
                    orbitals.append(
                        {
                            "index": int(match.group(1)),
                            "type": " ".join(match.group(2).split()),
                            "atoms": atoms,
                            "atom_indices": [
                                int(n) - 1 for n in re.findall(r"\b\d+\b", atoms)
                            ],
                            "occupancy": float(match.group(4)),
                            "energy": float(match.group(5)),
                        }
                    )
        self.data["nbo_orbitals"] = orbitals

    def parse_nbo_hybrids(self):
        """Augment nbo_orbitals with per-atom hybridization (%s / %p / %d).

        From the detailed "(Occupancy) Bond orbital / Coefficients / Hybrids"
        block, e.g.:
            1. (1.99996) CR ( 1) O  1            s(100.00%)
            2. (1.99698) LP ( 1) O  1            s(  0.00%)p 1.00( 99.92%) ...
            4. (1.99864) BD ( 1) O  1- H  2
                       ( 72.35%)   0.8506* O  1 s( 22.56%)p 3.43( 77.29%) ...
                       ( 27.65%)   0.5258* H  2 s( 99.76%)p 0.00(  0.24%)
        Adds a "hybrids" list to each matching NBO: one component per atom
        with {atom_sym, atom_idx, s_pct, p_pct, d_pct, weight_pct?, label}.

        If the block appears more than once, the last occurrence wins (it
        corresponds to the same NBO analysis run as nbo_orbitals /
        nbo_perturbation, which both already pick their last occurrence).
        Ported from the ORCA Result Analyzer.
        """
        orbitals = self.data.get("nbo_orbitals", [])
        if not orbitals:
            return

        start = self._last_index(
            lambda line: "Bond orbital / Coefficients / Hybrids" in line
        )
        if start is None:
            return

        s_re = re.compile(r"s\(\s*([\d.]+)%\)")
        p_re = re.compile(r"p\s*[\d.]*\(\s*([\d.]+)%\)")
        d_re = re.compile(r"d\s*[\d.]*\(\s*([\d.]+)%\)")
        header_re = re.compile(
            r"^\s*(\d+)\.\s*\(\s*[\d.]+\)\s*\S+?\s*\(\s*\d+\)\s*(.*)$"
        )
        pol_re = re.compile(
            r"^\s*\(\s*([\d.]+)%\)\s+[-\d.]+\*?\s*([A-Za-z]+)\s+(\d+)\s+(.*)$"
        )

        def _label(s_pct, p_pct, d_pct):
            # Below ~5% s the sp^n ratio explodes and is meaningless (typical
            # of diffuse Rydberg orbitals); fall back to the dominant character.
            if s_pct < 5.0:
                return "d" if d_pct > p_pct else "p"
            if p_pct <= 1.0 and d_pct <= 1.0:
                return "s"
            base = f"sp{p_pct / s_pct:.2f}" if p_pct > 1.0 else "s"
            # Annotate d only when chemically meaningful (> ~1%).
            if d_pct > 1.0:
                base += f"d{d_pct / s_pct:.2f}"
            return base

        def _extract(text, atom_sym=None, atom_num=None, weight=None):
            sm = s_re.search(text)
            if not sm:
                return None
            pm = p_re.search(text)
            dm = d_re.search(text)
            s_pct = float(sm.group(1))
            p_pct = float(pm.group(1)) if pm else 0.0
            d_pct = float(dm.group(1)) if dm else 0.0
            comp = {
                "atom_sym": atom_sym,
                "atom_idx": (atom_num - 1) if atom_num is not None else None,
                "s_pct": s_pct,
                "p_pct": p_pct,
                "d_pct": d_pct,
                "label": _label(s_pct, p_pct, d_pct),
                "raw": text[sm.start() :].rstrip(),
            }
            if weight is not None:
                comp["weight_pct"] = weight
            return comp

        hybrids = {}  # nbo index -> list of components
        current_idx = None
        end_markers = (
            "SECOND ORDER PERTURBATION",
            "NATURAL BOND ORBITALS (SUMMARY)",
            "NBO ANALYSIS COMPLETED",
        )
        for row in self.lines[start + 1 :]:
            if any(marker in row.upper() for marker in end_markers):
                break

            hm = header_re.match(row)
            if hm:
                current_idx = int(hm.group(1))
                rest = hm.group(2)
                comp = _extract(rest)
                if comp is not None:
                    am = re.match(r"\s*([A-Za-z]+)\s+(\d+)", rest)
                    if am:
                        comp["atom_sym"] = am.group(1)
                        comp["atom_idx"] = int(am.group(2)) - 1
                    hybrids.setdefault(current_idx, []).append(comp)
                continue

            pm = pol_re.match(row)
            if pm and current_idx is not None:
                comp = _extract(
                    pm.group(4), pm.group(2), int(pm.group(3)), float(pm.group(1))
                )
                if comp is not None:
                    hybrids.setdefault(current_idx, []).append(comp)

        for orb in orbitals:
            orb["hybrids"] = hybrids.get(orb["index"], [])

    def parse_nbo_perturbation(self):
        """Donor/acceptor rows of the second-order perturbation table.

        NBO 3.1 separates donor and acceptor with "/" under a title-case
        header; NBO 7 prints an all-caps header and no separator. The last
        table in the file wins, matching parse_nbo_orbitals.
        """
        interactions = []
        start = self._last_index(
            lambda line: "second order perturbation theory analysis" in line.lower()
        )
        if start is not None:
            for row in self.lines[start + 1 : start + 20000]:
                if "natural bond orbitals (summary)" in row.lower() or (
                    "NBO analysis completed" in row
                ):
                    break
                match = _E2_ROW.match(row)
                if match is None:
                    match = _E2_ROW_NBO7.match(row)
                    if match and (
                        "/" in match.group(2) or "/" in match.group(4)
                    ):
                        # A "/" means this is an NBO 3.1 row the strict
                        # pattern rejected; do not guess at its fields.
                        match = None
                if match:
                    interactions.append(
                        {
                            "donor": f"{match.group(1)}. {' '.join(match.group(2).split())}",
                            "acceptor": f"{match.group(3)}. {' '.join(match.group(4).split())}",
                            "e2_kcal": float(match.group(5)),
                            "e_diff": float(match.group(6)),
                            "fock": float(match.group(7)),
                        }
                    )
        self.data["nbo_perturbation"] = interactions

    def parse_nbo_totals(self):
        """Lewis / non-Lewis occupancy totals of the last NBO summary.

        The footer under the summary table, e.g.:
            Total Lewis    9.98992  ( 99.8992%)
            Valence non-Lewis    0.00027  (  0.0027%)
            Rydberg non-Lewis    0.00982  (  0.0982%)
        Stored as {"lewis", "lewis_pct", "valence_nonlewis",
        "valence_nonlewis_pct", "rydberg_nonlewis", "rydberg_nonlewis_pct"}
        (occupancies in electrons, percents of the total), or None when the
        job printed no NBO summary. The same lines also appear in the NPA
        summary block, and open-shell jobs print one footer per spin set;
        the first occurrence of each label wins (the alpha set in
        open-shell jobs), consistent with the NPA charge reader.
        """
        keys = {
            "Total Lewis": "lewis",
            "Valence non-Lewis": "valence_nonlewis",
            "Rydberg non-Lewis": "rydberg_nonlewis",
            "Total non-Lewis": "total_nonlewis",
        }
        found = {}
        for line in self.lines:
            match = _NBO_TOTAL.match(line)
            if match and keys[match.group(1)] not in found:
                found[keys[match.group(1)]] = (
                    float(match.group(2)),
                    float(match.group(3)),
                )
        totals = None
        if "lewis" in found:
            totals = {}
            for key, (occupancy, pct) in found.items():
                totals[key] = occupancy
                totals[f"{key}_pct"] = pct
        self.data["nbo_totals"] = totals

    def parse_nbo_population(self):
        """Per-atom NPA breakdown: charge, Core, Valence, Rydberg, Total.

        Reads the first "Summary of Natural Population Analysis" table —
        the same table parse_charges takes the NBO charges from, keeping
        its remaining columns. The first table is the spin-summed one in
        open-shell jobs (per-spin tables follow it). Each entry is
        {"index", "element", "charge", "core", "valence", "rydberg",
        "total"} with electron counts; [] when no table was printed.
        """
        rows = []
        start = next(
            (
                i
                for i, line in enumerate(self.lines)
                if "Summary of Natural Population Analysis" in line
            ),
            None,
        )
        if start is not None:
            for line in self.lines[start : start + 60]:
                match = _NPA_ROW.match(line)
                if match:
                    rows.append(
                        {
                            "index": int(match.group(2)),
                            "element": match.group(1),
                            "charge": float(match.group(3)),
                            "core": float(match.group(4)),
                            "valence": float(match.group(5)),
                            "rydberg": float(match.group(6)),
                            "total": float(match.group(7)),
                        }
                    )
                elif rows and "===" in line:
                    break
        self.data["nbo_population"] = rows

    def parse_nbo_configurations(self):
        """Natural electron configuration of each atom.

        The "Natural Electron Configuration" table, e.g.
        ``O  1  [core]2s( 1.75)2p( 5.12)3p( 0.01)`` becomes
        {"index", "element", "core": True,
         "orbitals": {"2s": 1.75, "2p": 5.12, "3p": 0.01},
         "raw": "[core]2s( 1.75)2p( 5.12)3p( 0.01)"}.
        First table wins, as for parse_nbo_population; [] when absent.
        """
        configs = []
        start = next(
            (
                i
                for i, line in enumerate(self.lines)
                if "Natural Electron Configuration" in line
            ),
            None,
        )
        if start is not None:
            for line in self.lines[start + 1 : start + 80]:
                match = re.match(r"^\s*([A-Z][a-z]?)\s+(\d+)\s+(\S.*)$", line)
                if match:
                    orbitals = {
                        label.lower(): float(occ)
                        for label, occ in _CONFIG_ORBITAL.findall(match.group(3))
                    }
                    if orbitals:
                        configs.append(
                            {
                                "index": int(match.group(2)),
                                "element": match.group(1),
                                "core": "[core]" in match.group(3),
                                "orbitals": orbitals,
                                "raw": match.group(3).strip(),
                            }
                        )
                        continue
                if configs:
                    break
        self.data["nbo_configurations"] = configs

    def parse_nbo_nao(self):
        """Natural atomic orbital (NAO) occupancies and energies.

        Reads the first "Natural atomic orbital occupancies" table.
        Each entry is {"index", "element", "atom", "lang", "kind",
        "ao", "occupancy", "energy"} where kind is Cor/Val/Ryd and ao
        the orbital label inside the parentheses (e.g. "2p"). The first
        table is the spin-summed one in open-shell jobs; NBO 3.1 prints
        no Energy column there, so energy is None for its rows.
        [] when no table was printed.
        """
        rows = []
        start = next(
            (
                i
                for i, line in enumerate(self.lines)
                if "natural atomic orbital occupancies" in line.lower()
            ),
            None,
        )
        if start is not None:
            gap = 0
            for line in self.lines[start + 1 : start + 5000]:
                if rows and "natural atomic orbital occupancies" in line.lower():
                    break  # the next (per-spin) NAO table begins
                match = _NAO_ROW.match(line)
                if match:
                    gap = 0
                    rows.append(
                        {
                            "index": int(match.group(1)),
                            "element": match.group(2),
                            "atom": int(match.group(3)),
                            "lang": match.group(4),
                            "kind": match.group(5),
                            "ao": match.group(6).strip(),
                            "occupancy": float(match.group(7)),
                            "energy": (
                                float(match.group(8)) if match.group(8) else None
                            ),
                        }
                    )
                elif rows:
                    # Atom blocks inside the table are separated by
                    # blank lines; tolerate a few before giving up.
                    gap += 1
                    if gap > 10:
                        break
        self.data["nbo_nao"] = rows

    def parse_nbo_bending(self):
        """NHO directionality / bond-bending deviations.

        Reads the last "NHO Directionality and Bond Bending" table that
        prints any rows (it belongs to the final NBO analysis, like the
        summary; open-shell jobs print one table per spin set, and a set
        with no above-threshold rows prints an empty table). Each
        entry is {"index", "type", "atoms", "loc_theta", "loc_phi",
        "h1_theta", "h1_phi", "h1_dev", "h2_theta", "h2_phi",
        "h2_dev"} in degrees, with None where the table prints "--"
        (e.g. a lone pair has no line of centers). [] when no table
        printed any rows — NBO only prints rows above its thresholds.
        """
        entries = []
        starts = [
            i
            for i, line in enumerate(self.lines)
            if "nho directionality" in line.lower()
        ]
        for start in reversed(starts):
            entries = []
            i = start
            while i < len(self.lines) and "====" not in self.lines[i]:
                i += 1
            for line in self.lines[i + 1 : i + 200]:
                match = _BENDING_ROW.match(line)
                if not match:
                    if entries:
                        break
                    continue
                descriptor = match.group(2).strip()
                head = re.match(
                    r"^(CR|LP|LV|BD\*|BD|RY\*|RY)\s*\(\s*(\d+)\)\s*(.*)$",
                    descriptor,
                )
                if not head:
                    if entries:
                        break
                    continue
                values = [
                    None if token == "--" else float(token)
                    for token in match.group(3).split()
                ]
                entries.append(
                    {
                        "index": int(match.group(1)),
                        "type": head.group(1),
                        "atoms": [
                            (element, int(number))
                            for element, number in re.findall(
                                r"([A-Z][a-z]?)\s+(\d+)", head.group(3)
                            )
                        ],
                        "loc_theta": values[0],
                        "loc_phi": values[1],
                        "h1_theta": values[2],
                        "h1_phi": values[3],
                        "h1_dev": values[4],
                        "h2_theta": values[5],
                        "h2_phi": values[6],
                        "h2_dev": values[7],
                    }
                )
            if entries:
                break
        self.data["nbo_bending"] = entries

    def _charge_rows(self, start, spin_column=False):
        """Rows of an "atom symbol value [spin]" table beginning at *start*."""
        rows = []
        for line in self.lines[start:]:
            match = _ATOM_CHARGE.match(line)
            if not match:
                if rows:
                    break
                if line.strip() and not re.match(r"^\s*\d+(\s+\d+)*\s*$", line):
                    break
                continue
            row = {
                "atom_idx": int(match.group(1)) - 1,
                "atom_sym": match.group(2),
                "charge": float(match.group(3)),
            }
            if spin_column and match.group(4) is not None:
                row["spin"] = float(match.group(4))
            rows.append(row)
        return rows

    def _esp_scheme(self):
        """Display name of the electrostatic-potential charges this job requested."""
        route = (_route_text(self.lines[:400]) or "").upper()
        if "CHELPG" in route:
            return "CHELPG"
        if re.search(r"\bMK\b", route):
            return "Merz-Kollman"
        return "ESP"

    def parse_charges(self):
        """Every charge scheme the job printed; the last print of each wins."""
        self.data["charges"] = {}
        charges = self.data["charges"]
        esp_name = self._esp_scheme()
        npa_read = False
        for i, line in enumerate(self.lines):
            text = line.strip()
            if text in (
                "Mulliken charges:",
                "Mulliken charges and spin densities:",
                "Mulliken atomic charges:",
            ):
                rows = self._charge_rows(i + 1, spin_column=True)
                if rows:
                    charges["Mulliken"] = rows
            elif text == "ESP charges:":
                rows = self._charge_rows(i + 1)
                if rows:
                    charges[esp_name] = rows
            elif text in ("APT charges:", "APT atomic charges:"):
                rows = self._charge_rows(i + 1)
                if rows:
                    charges["APT"] = rows
            elif text.startswith("Hirshfeld charges, spin densities, dipoles, and CM5"):
                self._read_hirshfeld(i + 2)
            elif text.startswith("SCF Done:"):
                npa_read = False
            elif text.startswith("Summary of Natural Population Analysis") and not npa_read:
                # An open-shell job repeats this table per spin afterwards, with
                # spin populations rather than charges: only the first is charges.
                npa_read = self._read_natural_population(i)

    def _read_hirshfeld(self, start):
        """Hirshfeld and CM5 columns of the one table that prints them both."""
        hirshfeld, cm5 = [], []
        for line in self.lines[start:]:
            parts = line.split()
            if len(parts) != 8 or not parts[0].isdigit():
                break
            try:
                base = {"atom_idx": int(parts[0]) - 1, "atom_sym": parts[1]}
                hirshfeld.append(dict(base, charge=float(parts[2]), spin=float(parts[3])))
                cm5.append(dict(base, charge=float(parts[7])))
            except ValueError as exc:
                logging.warning("Hirshfeld: bad row %r: %s", line, exc)
                return
        if hirshfeld:
            self.data["charges"]["Hirshfeld"] = hirshfeld
            self.data["charges"]["CM5"] = cm5

    def _read_natural_population(self, start):
        """NPA charges from the table at *start*; True when any were read."""
        rows = []
        for line in self.lines[start + 1 : start + 20000]:
            if "Total *" in line:
                break
            parts = line.split()
            if len(parts) >= 7 and parts[1].isdigit():
                try:
                    rows.append(
                        {
                            "atom_idx": int(parts[1]) - 1,
                            "atom_sym": parts[0],
                            "charge": float(parts[2]),
                            "core": float(parts[3]),
                            "valence": float(parts[4]),
                            "rydberg": float(parts[5]),
                            "total": float(parts[6]),
                        }
                    )
                except ValueError as exc:
                    logging.warning("NPA: bad row %r: %s", line, exc)
                    return False
        if rows:
            self.data["charges"]["NBO"] = rows
        return bool(rows)
