# fitspeek — independent adversarial review (2026-11-17)

Reviewer did not build this tool. All checks below were actually run;
anything not run is said so explicitly.

## Test battery

- Original suite: **8 tests, all passing** (`python run_tests.py` with the
  shared astro-tools venv).
- Final suite after fixes and added regression tests: **15 tests, all
  passing**.

## Independent oracle verification (no fitspeek code used to build ground truth)

FITS files were crafted directly with Astropy from known arrays/headers,
the tool was run on them, and the PNG/index outputs were inspected with
PIL/NumPy:

- **Sheet geometry / tiling:** 5 known images at cols=1/2/3 produced
  5/2/1 sheets of exactly 320×320 / 640×640 / 960×960 px, and `index.csv`
  sheet+tile assignments matched the expected layout in every case. PASS
- **Orientation:** a bright block placed in each of the four data corners
  (40×40 zero image, block = 10000) appeared in the matching displayed
  corner of the PNG in all 4 cases, confirming `origin='lower'` handling
  (data row 0 at the bottom). PASS
- **Scaling:** a noise-free left-black / right-white (0 / 1000) image
  rendered left ≈ 0, right ≈ 255 in the PNG. PASS
- **Cube handling:** a 3-plane cube whose planes each contain one disjoint
  bright block median-combines to a constant image; the rendered tile was
  perfectly flat (std 0), as median-combining requires — a mean would have
  kept visible blocks. PASS
- **Compass:** the same image with a CD-matrix WCS showed 130 yellow
  (compass) pixels in the PNG; the identical image without WCS showed 0.
  PASS
- **Header stamp:** identical filename+data rendered with and without
  OBJECT/DATE-OBS/FILTER/EXPTIME — the stamped version carried ~4× the
  title-area ink (307 vs 73 dark pixels), i.e. the extra stamp line is
  really rendered. Exact stamp *text* was not OCR-verified (no OCR tool
  available in this environment); presence only.

## Adversarial attempts (selection, all executed)

Handled correctly with no fix needed: constant and zero images, 1×N / N×1 /
single-pixel images, partial NaN, mixed NaN/inf/finite, int/uint images,
4D and 5D cubes, cube with one all-NaN plane, 1D-only file, empty file,
text file with a `.fits` name, mid-data byte corruption, `.fits.gz`,
uppercase `.FITS`, no WCS, CDELT/CROTA-only WCS (no compass, as
documented), PC-matrix WCS, degenerate CD matrix, EXPTIME=0, string
EXPTIME, very long OBJECT, empty directory, nonexistent input path,
duplicate inputs (rendered twice, two index rows), table-then-image HDU
order, 1D primary + 2D extension.

## Findings and fixes

### HIGH — truncated/corrupt-data FITS crashed the whole batch (FIXED)
A file whose header block is intact but whose data is missing/truncated
raises `ValueError: cannot reshape array of size 0 into shape (20,30)`
only when `hdu.data` is touched — outside the `fits.open` try/except.
Reproduced: a batch of one truncated + one good file crashed with no
sheet and no `index.csv` at all, contradicting the README's "Corrupt
files are skipped the same way — one bad file never stops the batch."
Fix: per-HDU data access in `load_image` is now wrapped; a data-access
failure becomes `SkipImage("unreadable file (…)")`. Regression test added.

### MEDIUM — unusable first HDU hid a valid later image HDU (FIXED)
`load_image` *raised* on a zero-size or all-NaN first image HDU instead
of continuing the HDU search. A zero-size primary + valid ImageHDU, and
an all-NaN primary + valid ImageHDU, were both wrongly skipped.
Fix: per-HDU skip reasons are now remembered and the search continues;
the remembered reason is reported only if no HDU renders. Two regression
tests added. (Single-HDU zero-size / all-NaN files keep their original
skip reasons, so existing behaviour/tests are unchanged.)

### MEDIUM — constant images rendered black, not flat grey (FIXED)
`stretch()` returns 0.5 for a constant image, but `render_sheets` called
`imshow` without `vmin`/`vmax`, so Matplotlib re-normalised each
thumbnail and mapped the constant 0.5 to black — a flat/bias frame was
indistinguishable from an all-zero frame. Verified by oracle: centre of a
constant-42 tile measured 0 before the fix, 128 (mid-grey) after.
Fix: `imshow(..., vmin=0.0, vmax=1.0)`. Regression test added.

### LOW — extreme dynamic range collapsed to all-black (FIXED)
`stretch()` on data spanning ±1e308: `vmax - vmin` overflows to inf and
every pixel scaled to 0 (plus RuntimeWarnings). Fix: when the span is
not finite, data and limits are halved before scaling (ordering
preserved); arithmetic is also wrapped in `errstate`. Post-fix the same
array stretches to exactly {0, 0.5, 1}. Regression test added.

### LOW — `make()` accepted invalid `cols` (FIXED)
The CLI rejects `--cols < 1`, but the library entry point crashed with
`ZeroDivisionError` (cols=0) or a Matplotlib figure-size error (cols<0).
Fix: `make()` raises `ValueError("cols must be a positive integer")`.
Regression test added. Note: there is still no *upper* bound — cols=20
was observed to be impractically slow/large; documented in the README
rather than capped, as a cap would be a new design decision.

### LOW — all-inf image misreported as "all-NaN image" (FIXED)
An all-+inf (or all-−inf) image was skipped with reason "all-NaN image".
Fix: the reason is now "all-NaN image" only when every pixel is NaN,
otherwise "no finite pixels". Regression test added.

## README assessment

Mostly honest; the "Honest limits" section is accurate as far as it goes
(non-recursive folders, exact stamp keywords, CD/PC-only compass,
triage-not-measurement, in-memory rendering all verified). Overstatements
found and corrected in the README during this review:

- Intro claimed "a compass on every frame that has WCS" — a
  CDELT/CROTA-only WCS gets no compass (the limits section itself said
  so). Intro wording corrected.
- "`index.csv` lists every input file" — non-FITS files inside a folder
  argument are silently not collected and never appear in the index.
  Wording corrected and the behaviour spelled out.
- Added the missing limit: "`--cols` has no upper bound; very large
  values produce very large, slow sheets."
- Pre-fix, "Corrupt files are skipped the same way — one bad file never
  stops the batch" was false for truncated-data files (HIGH finding
  above); it is true after the fix.
- Not verified: the `pip install git+https://github.com/RakTikoo/fitspeek`
  line — this review was local-only by instruction, so the remote repo's
  existence/state was not checked. The local-tree install (below) works.

## Fresh-install + README demo

- A full fresh-venv `pip install <local tree>` (new venv, dependencies
  resolved from PyPI: astropy 8.0.2, numpy 2.5.5, matplotlib 3.11.2)
  succeeded and the `fitspeek` entry point ran. Caveat: the venv itself
  had to live on the home filesystem — this machine's `/tmp` is a 512 MB
  tmpfs that was ~93% full with other tasks' venvs, and a from-scratch
  dependency install in `/tmp` failed with "No space left on device"
  (environment limit, not a fitspeek defect; a `--no-deps` install of
  fitspeek itself into a `/tmp` venv did succeed).
- README demo reproduced **exactly**, both before and after the fixes:
  stdout lines and `index.csv` match the README byte-for-byte, and
  `sheet_001.png` is 640×640 as a 2×2 sheet at dpi 100 implies.

## Not a finding

- Duplicate inputs (same file twice, or file + its folder) render twice
  with two index rows — reasonable behaviour for explicit arguments.
- Byte-flipped data blocks still "render" (Astropy reads the garbage as
  floats); the values are wrong but nothing crashes — inherent to FITS
  without checksums, out of scope for a triage tool.
