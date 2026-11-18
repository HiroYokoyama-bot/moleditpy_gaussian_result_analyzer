# Changelog

## 1.2.0 — 2026-11-18

- New command line, with MoleditPy 5.3.0: `moleditpy -r gaussian_result_analyzer parse FILE` reads a log or checkpoint without opening a window and prints a summary (version, termination, formula, SCF energy, <S**2>, HOMO/LUMO, dipole, steps, optimisation convergence, frequencies, Gibbs energy, excited states, NMR, charges, NBO). `--json` prints the parsed data, `--section NAME` limits it to some keys, `-o PATH` writes to a file. A checkpoint is paired with its output, as in the window.

## 1.1.1 — 2026-11-18

- Read NBO tables of molecules past 100 atoms. From atom 100 on, NBO prints the number fused to the element (`H100`, `Cl101`), and those rows were skipped: NBO charges, the NPA breakdown, configurations, NAO occupancies, bending rows and the atom indices and hybrids of NBOs on those atoms were all incomplete.
- Read the NPA breakdown, natural electron configurations, NAO occupancies, bond-bending rows and NBO summary to their last row; fixed windows cut them off at 54 atoms, 78 atoms, 199 rows and 5000 lines.
- The Bond Analysis Lewis / non-Lewis line no longer fails on an output that has the Lewis total but not every non-Lewis line.

## 1.1.0 — 2026-11-15

- Read NBO 7 output. A Gaussian 16 linked with NBO 7 prints all-caps headers, E(2) rows without the "/" separator and RY (not RY*) labels; every NBO section previously came back empty.
- Bond Analysis now shows NBO hybrids; the parser never produced them before.
- New NBO data: Lewis / non-Lewis totals, the per-atom NPA breakdown, natural electron configurations, NAO occupancies and NHO bond-bending deviations.
- A log with several NBO summaries now resolves to the last one instead of one concatenated list with repeating indices.

## 1.0.1 — 2026-11-14

- Declare the Qt dependency as `PySide6-Essentials`, the package MoleditPy itself ships, instead of the full `PySide6`.

## 1.0.0 — 2026-11-14

- First stable release. Orbitals read from a log are now checked against an exact overlap matrix in the test suite; no change in behaviour since 0.3.1.

## 0.3.1 — 2026-11-12

- Open `Freq=HPModes` logs; their high-precision frequency table stopped the load.
- Counterpoise and ONIOM jobs report the whole system's charge and multiplicity, not the last fragment's.
- Read Gaussian 09's `Mulliken atomic charges:` and `APT atomic charges:` tables.
- Read a route echoed over two lines: CHELPG charges are no longer labelled "ESP", and `FOpt`/`POpt` jobs show their optimization steps.
- A job whose second step is still running, or was killed, no longer reports "Terminated normally".
- Open the log of a job above 50 atoms that ends inside a force table.

## 0.3.0 — 2026-11-11

- Read rigid `Scan` jobs: one point per energy, with the scanned value.
- Read IRC paths from the reverse end to the forward end, with the TS at zero; the final structure is the last point.
- Point the forces of every optimization step the right way above 50 atoms, not only of a single-point job.
- Opt Freq no longer shows the Freq step's check as one more optimization cycle.
- Counterpoise and anharmonic jobs report the job's own energy, structure and orbitals, not those of a monomer or a displaced geometry; the counterpoise-corrected energy and BSSE are listed.
- Read the CASSCF energy.
- A Link1 step that starts without "Proceeding to internal job step" is recognised.

## 0.2.0 — 2026-11-11

- Read every orbital energy of a job that prints no MO table. Only a few were read, and with wrong values.
- Read thermochemistry from `#T` logs, and the total energy of a double hybrid such as B2PLYP.
- Do not report the `<S**2>` of a stability analysis as spin contamination of a singlet.
- Point the forces of a molecule above 50 atoms the right way; Gaussian prints no input orientation there.
- Remove the **Load .fchk** button: the one **Open** button opens logs and formatted checkpoints, and a checkpoint beside a log is picked up by itself.
- The file box uses the same blue as the ORCA Result Analyzer.
- The description no longer says under development.

## 0.1.3 — 2026-11-10

- On MoleditPy 5.2.0 and later, pass another program's file on to the next plugin by returning `False` instead of raising, so the host no longer logs a traceback for it. Older hosts keep the exception.

## 0.1.2 — 2026-11-08

- Show every optimization cycle by default when opening the optimization tab.
- Keep the scan-point view as the default for relaxed scans with multiple points.
- Add regressions using real Gaussian logs, formatted checkpoints and real Qt windows.
