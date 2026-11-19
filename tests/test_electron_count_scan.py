"""The electron counts are read once per file, not once per orbital."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _parser_loader import load_standalone_parser  # noqa: E402

_SAMPLES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_outputs")
_parser_mod = load_standalone_parser()


def _read(name):
    with open(os.path.join(_SAMPLES, name), encoding="utf-8") as handle:
        return handle.read()


class _CountingParser(_parser_mod.GaussianParser):
    scans = 0

    def _scan_electron_counts(self):
        self.scans += 1
        return super()._scan_electron_counts()


def _occupied(parser, spin):
    return sum(1 for mo in parser.data["mos"] if mo["spin"] == spin and mo["occ"] > 0)


class TestElectronCountScan(unittest.TestCase):
    def test_a_file_with_many_orbitals_is_scanned_once(self):
        parser = _CountingParser()
        parser.load_from_memory(_read("c33cl_nbo7.out"), "c33cl_nbo7.out")
        self.assertGreater(len(parser.data["mos"]), 100)
        self.assertEqual(parser.scans, 1)

    def test_a_reused_parser_reads_the_next_files_counts(self):
        parser = _CountingParser()
        parser.load_from_memory(_read("c33cl_nbo7.out"), "c33cl_nbo7.out")
        parser.load_from_memory(_read("o2_fk_uhf.out"), "o2_fk_uhf.out")
        self.assertEqual(parser.scans, 2)
        # Triplet O2: two more alpha than beta electrons.
        self.assertEqual(_occupied(parser, "alpha") - _occupied(parser, "beta"), 2)


if __name__ == "__main__":
    unittest.main()
