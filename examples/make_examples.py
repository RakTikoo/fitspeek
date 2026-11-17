"""Generate a small set of synthetic FITS files for trying fitspeek.

Creates, in examples/sample/ (next to this script):

- m42_wcs.fits      2D image with a full WCS (CD matrix) header
- plain_field.fits  2D image with no WCS
- cube.fits         3D data cube (fitspeek median-combines the planes)
- table_only.fits   a binary table and no image, to show the skip path

Run with:  python examples/make_examples.py
Then try:  fitspeek make examples/sample --out sheets --cols 2
"""

from pathlib import Path

import numpy as np
from astropy.io import fits

OUT = Path(__file__).parent / "sample"


def star_field(rng, shape=(60, 80), n_stars=25, background=100.0):
    data = rng.normal(background, 5.0, shape)
    yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
    for _ in range(n_stars):
        cx, cy = rng.uniform(3, shape[1] - 3), rng.uniform(3, shape[0] - 3)
        amp = rng.uniform(500, 3000)
        data += amp * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 1.2**2))
    return data


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)

    # 1. Image with WCS (CD matrix), stamped header keywords.
    data = star_field(rng)
    hdr = fits.Header()
    hdr["OBJECT"] = "M42"
    hdr["DATE-OBS"] = "2026-11-16T03:12:00"
    hdr["FILTER"] = "Ha"
    hdr["EXPTIME"] = 120.0
    hdr["CTYPE1"] = "RA---TAN"
    hdr["CTYPE2"] = "DEC--TAN"
    hdr["CRPIX1"] = 40.0
    hdr["CRPIX2"] = 30.0
    hdr["CRVAL1"] = 83.822
    hdr["CRVAL2"] = -5.391
    cd = 0.00028  # ~1 arcsec/pixel
    hdr["CD1_1"] = -cd
    hdr["CD1_2"] = 0.0
    hdr["CD2_1"] = 0.0
    hdr["CD2_2"] = cd
    fits.writeto(OUT / "m42_wcs.fits", data, hdr, overwrite=True)

    # 2. Plain image, no WCS.
    data = star_field(rng)
    hdr = fits.Header()
    hdr["OBJECT"] = "FLAT"
    hdr["DATE-OBS"] = "2026-11-16T02:00:00"
    hdr["EXPTIME"] = 5.0
    fits.writeto(OUT / "plain_field.fits", data, hdr, overwrite=True)

    # 3. A small cube: 5 planes of the same field, drifting brightness.
    planes = np.stack([star_field(rng) * (1.0 + 0.1 * i) for i in range(5)])
    hdr = fits.Header()
    hdr["OBJECT"] = "NGC7000"
    hdr["DATE-OBS"] = "2026-11-16T04:45:00"
    hdr["FILTER"] = "OIII"
    hdr["EXPTIME"] = 300.0
    fits.writeto(OUT / "cube.fits", planes, hdr, overwrite=True)

    # 4. Table-only file: no image HDU at all.
    table = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="SOURCE", format="20A", array=["a", "b"]),
            fits.Column(name="FLUX", format="D", array=[1.5, 2.5]),
        ]
    )
    table.header["OBJECT"] = "CATALOG"
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(
        OUT / "table_only.fits", overwrite=True
    )

    for path in sorted(OUT.glob("*.fits")):
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
