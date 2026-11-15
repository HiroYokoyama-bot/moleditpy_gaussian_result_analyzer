# Release protocol

When the user requests a release of the Gaussian Result Analyzer plugin, complete
the following steps in order:

1. Add regression tests for the fix and bump the patch version in
   `gaussian_result_analyzer/__init__.py` (`PLUGIN_VERSION`).
2. Run the full test suite with coverage using
   `python run_tests.py -p pytest_cov --cov=gaussian_result_analyzer --cov-report=term-missing --cov-report=xml`.
   Codecov is enabled for this repository; retain the full coverage upload in CI.
3. Push a feature branch and open a pull request to `main`.
4. Wait for all PR CI and Codecov checks to pass. Resolve failures before merging.
5. Merge the PR to `main`, then wait for main CI on the merged commit to pass.
6. Only after main CI is green, push the matching `vX.Y.Z` tag on that merged
   commit to trigger the automatic release workflow.
7. Verify the release workflow and release artifact, and report the PR,
   coverage result, tag, and release URL to the user.

Do not bypass the PR or tag a commit before its main CI passes.

# Shared implementation

The ORCA Result Analyzer (`HiroYokoyama/moleditpy_orca_result_analyzer_plugin`)
is the main implementation for everything except the parser. The non-parser
modules (the dialogs and analysis panels, e.g. `bond_analysis.py` and
`charge_analysis.py`) are copies of the ORCA ones: make shared changes in the
ORCA repository first, then copy them here. Keep only the differences that
are specific to this plugin — the package and dialog class names,
`__init__.py`, the `gui.py` wiring, and the parser modules, which are
specific to each host program's output format.
