"""Tests for fitspeek, using small synthetic FITS files."""

import csv
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fitspeek.core import compass_vectors, make, stretch  # noqa: E402


def write_image(path, data, **keywords):
    hdr = fits.Header()
    for key, value in keywords.items():
        hdr[key] = value
    fits.writeto(path, data, hdr, overwrite=True)


def read_index(out_dir):
    with open(Path(out_dir) / "index.csv", newline="") as fh:
        return list(csv.DictReader(fh))


class TestFitsPeek(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.indir = Path(self._tmp.name) / "in"
        self.indir.mkdir()
        self.outdir = Path(self._tmp.name) / "out"
        rng = np.random.default_rng(7)
        self.field = rng.normal(100, 5, (40, 50)) + 1000 * np.exp(
            -((np.arange(50) - 25) ** 2 + (np.arange(40)[:, None] - 20) ** 2) / 8.0
        )

    def test_sheets_and_index(self):
        # Three plain images -> one 2x2 sheet, index lists every input.
        names = []
        for i in range(3):
            name = f"img{i}.fits"
            write_image(self.indir / name, self.field, OBJECT=f"T{i}", EXPTIME=10.0)
            names.append(name)
        summary = make([str(self.indir)], str(self.outdir), cols=2)
        self.assertEqual(summary["rendered"], 3)
        self.assertEqual(summary["skipped"], 0)
        self.assertEqual(summary["sheets"], 1)
        sheet = self.outdir / "sheet_001.png"
        self.assertTrue(sheet.exists())
        self.assertGreater(sheet.stat().st_size, 5000)
        rows = read_index(self.outdir)
        self.assertEqual(len(rows), 3)
        self.assertEqual(
            sorted(Path(r["source_file"]).name for r in rows), sorted(names)
        )
        for row in rows:
            self.assertEqual(row["status"], "rendered")
            self.assertEqual(row["sheet"], "sheet_001.png")

    def test_table_only_skipped_batch_continues(self):
        write_image(self.indir / "good.fits", self.field, OBJECT="OK")
        table = fits.BinTableHDU.from_columns(
            [fits.Column(name="X", format="D", array=[1.0, 2.0])]
        )
        fits.HDUList([fits.PrimaryHDU(), table]).writeto(
            self.indir / "table_only.fits", overwrite=True
        )
        summary = make([str(self.indir)], str(self.outdir), cols=2)
        self.assertEqual(summary["rendered"], 1)
        self.assertEqual(summary["skipped"], 1)
        rows = {Path(r["source_file"]).name: r for r in read_index(self.outdir)}
        self.assertEqual(rows["table_only.fits"]["status"], "skipped")
        self.assertIn("no image HDU", rows["table_only.fits"]["reason"])
        self.assertEqual(rows["good.fits"]["status"], "rendered")
        self.assertTrue((self.outdir / "sheet_001.png").exists())

    def test_cube_renders(self):
        cube = np.stack([self.field * (1 + 0.1 * i) for i in range(4)])
        write_image(self.indir / "cube.fits", cube, OBJECT="CUBE")
        summary = make([str(self.indir)], str(self.outdir), cols=1)
        self.assertEqual(summary["rendered"], 1)
        self.assertEqual(summary["skipped"], 0)
        sheet = self.outdir / "sheet_001.png"
        self.assertTrue(sheet.exists())
        self.assertGreater(sheet.stat().st_size, 5000)

    def test_all_nan_skipped(self):
        write_image(self.indir / "good.fits", self.field, OBJECT="OK")
        write_image(
            self.indir / "nan.fits", np.full((20, 20), np.nan), OBJECT="BAD"
        )
        summary = make([str(self.indir)], str(self.outdir), cols=2)
        self.assertEqual(summary["rendered"], 1)
        self.assertEqual(summary["skipped"], 1)
        rows = {Path(r["source_file"]).name: r for r in read_index(self.outdir)}
        self.assertEqual(rows["nan.fits"]["status"], "skipped")
        self.assertIn("all-NaN", rows["nan.fits"]["reason"])

    def test_wcs_image_renders(self):
        hdr = {
            "OBJECT": "WCS",
            "CTYPE1": "RA---TAN",
            "CTYPE2": "DEC--TAN",
            "CRPIX1": 25.0,
            "CRPIX2": 20.0,
            "CRVAL1": 83.822,
            "CRVAL2": -5.391,
            "CD1_1": -0.00028,
            "CD1_2": 0.0,
            "CD2_1": 0.0,
            "CD2_2": 0.00028,
        }
        write_image(self.indir / "wcs.fits", self.field, **hdr)
        summary = make([str(self.indir)], str(self.outdir), cols=1)
        self.assertEqual(summary["rendered"], 1)
        self.assertTrue((self.outdir / "sheet_001.png").exists())

    def test_compass_vectors_from_cd_matrix(self):
        hdr = fits.Header()
        hdr["CTYPE1"] = "RA---TAN"
        hdr["CTYPE2"] = "DEC--TAN"
        hdr["CRPIX1"] = 25.0
        hdr["CRPIX2"] = 20.0
        hdr["CRVAL1"] = 83.822
        hdr["CRVAL2"] = -5.391
        hdr["CD1_1"] = -0.00028
        hdr["CD1_2"] = 0.0
        hdr["CD2_1"] = 0.0
        hdr["CD2_2"] = 0.00028
        vectors = compass_vectors(hdr, self.field.shape)
        self.assertIsNotNone(vectors)
        north, east = vectors
        # CD here is axis-aligned: north is +y on the image, east is -x.
        self.assertAlmostEqual(north[1], 1.0, places=5)
        self.assertAlmostEqual(east[0], -1.0, places=5)
        # A header with no CD/PC matrix gets no compass.
        self.assertIsNone(compass_vectors(fits.Header(), self.field.shape))

    def test_corrupt_file_skipped_no_crash(self):
        write_image(self.indir / "good.fits", self.field, OBJECT="OK")
        (self.indir / "corrupt.fits").write_bytes(b"this is not a FITS file" * 10)
        summary = make([str(self.indir)], str(self.outdir), cols=2)
        self.assertEqual(summary["rendered"], 1)
        self.assertEqual(summary["skipped"], 1)
        rows = {Path(r["source_file"]).name: r for r in read_index(self.outdir)}
        self.assertEqual(rows["corrupt.fits"]["status"], "skipped")
        self.assertTrue(rows["corrupt.fits"]["reason"])

    def test_truncated_file_skipped_batch_continues(self):
        # Header block intact but data missing: astropy raises only when
        # the data is touched. That must still be a skip, not a crash.
        write_image(self.indir / "good.fits", self.field, OBJECT="OK")
        write_image(self.indir / "aaa_trunc.fits", self.field, OBJECT="TRUNC")
        trunc = self.indir / "aaa_trunc.fits"
        trunc.write_bytes(trunc.read_bytes()[:2880])
        summary = make([str(self.indir)], str(self.outdir), cols=2)
        self.assertEqual(summary["rendered"], 1)
        self.assertEqual(summary["skipped"], 1)
        rows = {Path(r["source_file"]).name: r for r in read_index(self.outdir)}
        self.assertEqual(rows["aaa_trunc.fits"]["status"], "skipped")
        self.assertIn("unreadable", rows["aaa_trunc.fits"]["reason"])
        self.assertEqual(rows["good.fits"]["status"], "rendered")

    def test_zero_size_primary_valid_extension_renders(self):
        hdul = fits.HDUList(
            [fits.PrimaryHDU(data=np.empty((0, 0))), fits.ImageHDU(data=self.field)]
        )
        hdul.writeto(self.indir / "mef.fits", overwrite=True)
        summary = make([str(self.indir)], str(self.outdir), cols=1)
        self.assertEqual(summary["rendered"], 1)
        self.assertEqual(summary["skipped"], 0)

    def test_all_nan_primary_valid_extension_renders(self):
        hdul = fits.HDUList(
            [
                fits.PrimaryHDU(data=np.full((10, 10), np.nan)),
                fits.ImageHDU(data=self.field),
            ]
        )
        hdul.writeto(self.indir / "mef.fits", overwrite=True)
        summary = make([str(self.indir)], str(self.outdir), cols=1)
        self.assertEqual(summary["rendered"], 1)
        self.assertEqual(summary["skipped"], 0)

    def test_all_inf_reason_is_not_all_nan(self):
        write_image(self.indir / "inf.fits", np.full((10, 10), np.inf))
        summary = make([str(self.indir)], str(self.outdir), cols=1)
        self.assertEqual(summary["skipped"], 1)
        rows = read_index(self.outdir)
        self.assertIn("no finite pixels", rows[0]["reason"])

    def test_constant_image_renders_mid_grey_not_black(self):
        import matplotlib.image as mpimg

        write_image(self.indir / "flat.fits", np.full((40, 40), 42.0))
        summary = make([str(self.indir)], str(self.outdir), cols=1)
        self.assertEqual(summary["rendered"], 1)
        img = mpimg.imread(self.outdir / "sheet_001.png")
        centre = img[140:180, 140:180, 0]
        self.assertGreater(float(centre.mean()), 0.3)
        self.assertLess(float(centre.mean()), 0.7)

    def test_stretch_extreme_dynamic_range(self):
        out = stretch(np.array([[0.0, 1e308], [-1e308, 0.0]]))
        self.assertTrue(np.isfinite(out).all())
        self.assertGreater(float(out.max()), 0.9)
        self.assertLess(float(out.min()), 0.1)

    def test_make_rejects_bad_cols(self):
        write_image(self.indir / "a.fits", self.field)
        with self.assertRaises(ValueError):
            make([str(self.indir)], str(self.outdir), cols=0)

    def test_empty_directory_no_crash(self):
        summary = make([str(self.indir)], str(self.outdir), cols=3)
        self.assertEqual(summary["rendered"], 0)
        self.assertEqual(summary["skipped"], 0)
        self.assertEqual(summary["sheets"], 0)
        self.assertEqual(read_index(self.outdir), [])


if __name__ == "__main__":
    unittest.main()
