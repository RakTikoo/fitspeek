"""Core logic for fitspeek.

Reads FITS files, renders one thumbnail per image, tiles thumbnails into
contact-sheet PNGs, and writes an index.csv. A file that cannot be rendered
is skipped with a recorded reason; one bad file never stops the batch.
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from astropy.io import fits  # noqa: E402
from astropy.visualization import PercentileInterval, ZScaleInterval  # noqa: E402
from astropy.wcs import WCS  # noqa: E402

FITS_SUFFIXES = (".fits", ".fit", ".fts", ".fits.gz", ".fit.gz", ".fts.gz")


class SkipImage(Exception):
    """Raised when a file cannot be rendered; str(exc) is the reason."""


def is_fits_name(name: str) -> bool:
    return name.lower().endswith(FITS_SUFFIXES)


def collect_inputs(inputs) -> list[Path]:
    """Expand file/directory arguments into a sorted list of FITS paths."""
    found: list[Path] = []
    for raw in inputs:
        path = Path(raw)
        if path.is_dir():
            for child in sorted(path.iterdir()):
                if child.is_file() and is_fits_name(child.name):
                    found.append(child)
        elif path.is_file():
            found.append(path)
        else:
            print(f"note: {path} does not exist, ignoring")
    return found


def load_image(path: Path):
    """Return (data2d, header) for the first image HDU that holds data.

    Cubes (3D or more) are median-combined down to a single plane.
    Raises SkipImage with a plain reason when there is nothing to show.
    """
    try:
        hdul = fits.open(path, memmap=False)
    except Exception as exc:  # corrupt or not a FITS file at all
        raise SkipImage(f"unreadable file ({type(exc).__name__})") from exc
    pending_reason: str | None = None
    with hdul:
        for hdu in hdul:
            try:
                data = hdu.data
            except Exception as exc:  # header parses but data is corrupt/truncated
                raise SkipImage(f"unreadable file ({type(exc).__name__})") from exc
            if data is None or not isinstance(data, np.ndarray):
                continue
            if data.ndim < 2:
                continue
            header = hdu.header
            if data.size == 0:
                # Not usable, but a later HDU may still hold a real image.
                if pending_reason is None:
                    pending_reason = "zero-size image"
                continue
            if data.ndim > 2:
                reduce_axes = tuple(range(data.ndim - 2))
                with np.errstate(all="ignore"):
                    import warnings

                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        data = np.nanmedian(data, axis=reduce_axes)
            if data.size == 0:
                if pending_reason is None:
                    pending_reason = "zero-size image"
                continue
            if not np.isfinite(data).any():
                # Same idea: an all-NaN / all-inf plane must not hide a
                # usable image in a later HDU.
                if pending_reason is None:
                    pending_reason = (
                        "all-NaN image"
                        if bool(np.isnan(data).all())
                        else "no finite pixels"
                    )
                continue
            return np.asarray(data, dtype=float), header
    raise SkipImage(pending_reason or "no image HDU with data")


def stretch(data: np.ndarray) -> np.ndarray:
    """Normalise image data to 0..1 using a ZScale stretch.

    Falls back to a percentile stretch when ZScale cannot find limits
    (e.g. nearly constant images), and to a flat grey when even that fails.
    """
    try:
        vmin, vmax = ZScaleInterval().get_limits(data)
    except Exception:
        vmin = vmax = None
    if vmin is None or vmax is None or not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        try:
            vmin, vmax = PercentileInterval(99.5).get_limits(data)
        except Exception:
            vmin, vmax = None, None
    if vmin is None or vmax is None or vmax <= vmin:
        return np.full(data.shape, 0.5, dtype=float)
    with np.errstate(over="ignore", invalid="ignore"):
        span_is_finite = bool(np.isfinite(vmax - vmin))
    if not span_is_finite:
        # Extreme dynamic range (e.g. +/-1e308): the span overflows to inf
        # and every pixel would collapse to 0. Halving the data and the
        # limits preserves the ordering and keeps the arithmetic finite.
        data = np.asarray(data, dtype=float) * 0.5
        vmin, vmax = vmin * 0.5, vmax * 0.5
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        scaled = (data - vmin) / (vmax - vmin)
    return np.clip(np.nan_to_num(scaled, nan=0.0), 0.0, 1.0)


def stamp_line(header) -> str:
    """One-line stamp from OBJECT / DATE-OBS / FILTER / EXPTIME, when present."""
    parts = []
    obj = header.get("OBJECT")
    if obj:
        parts.append(str(obj).strip())
    dateobs = header.get("DATE-OBS")
    if dateobs:
        parts.append(str(dateobs).strip())
    filt = header.get("FILTER")
    if filt:
        parts.append(f"filter {str(filt).strip()}")
    exptime = header.get("EXPTIME")
    if exptime is not None:
        try:
            parts.append(f"{float(exptime):g}s")
        except (TypeError, ValueError):
            parts.append(f"{exptime}s")
    return " · ".join(parts)


def compass_vectors(header, shape):
    """Return (north_xy, east_xy) unit vectors in pixel directions, or None.

    Only attempted when the header carries a CD or PC matrix. Pixel
    directions are in image coordinates (x right, y up as displayed with
    origin='lower').
    """
    has_matrix = any(k in header for k in ("CD1_1", "CD1_2", "CD2_1", "CD2_2")) or any(
        k in header for k in ("PC1_1", "PC1_2", "PC2_1", "PC2_2")
    )
    if not has_matrix:
        return None
    try:
        import astropy.units as u

        wcs = WCS(header)
        phys_types = [t for t in wcs.world_axis_physical_types if t]
        if sum(1 for t in phys_types if t.startswith("pos.eq")) < 2:
            return None
        ny, nx = shape
        cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
        centre = wcs.pixel_to_world(cx, cy)
        north_sky = centre.directional_offset_by(0 * u.deg, 1 * u.arcmin)
        east_sky = centre.directional_offset_by(90 * u.deg, 1 * u.arcmin)
        nx_px, ny_px = wcs.world_to_pixel(north_sky)
        ex_px, ey_px = wcs.world_to_pixel(east_sky)
        north = np.array([float(nx_px) - cx, float(ny_px) - cy])
        east = np.array([float(ex_px) - cx, float(ey_px) - cy])
        nn = np.linalg.norm(north)
        ne = np.linalg.norm(east)
        if not (np.isfinite(nn) and np.isfinite(ne)) or nn == 0 or ne == 0:
            return None
        return north / nn, east / ne
    except Exception:
        return None


def _draw_compass(ax, north, east) -> None:
    """Draw N/E arrows in the lower-right corner of a thumbnail axes."""
    cx, cy, length = 0.86, 0.14, 0.09
    for vec, label in ((north, "N"), (east, "E")):
        dx, dy = float(vec[0]) * length, float(vec[1]) * length
        ax.annotate(
            "",
            xy=(cx + dx, cy + dy),
            xytext=(cx, cy),
            xycoords="axes fraction",
            arrowprops=dict(arrowstyle="->", color="yellow", lw=1.2),
        )
        ax.text(
            cx + dx * 1.45,
            cy + dy * 1.45,
            label,
            color="yellow",
            fontsize=7,
            ha="center",
            va="center",
            transform=ax.transAxes,
        )


def render_sheets(items, out_dir: Path, cols: int) -> list[tuple[str, int]]:
    """Tile thumbnails into cols x cols sheets.

    items: list of dicts with keys data, header, name.
    Returns a list parallel to items of (sheet_filename, tile_number_1_based).
    """
    per_sheet = cols * cols
    placements: list[tuple[str, int]] = []
    n_sheets = math.ceil(len(items) / per_sheet) if items else 0
    for sheet_idx in range(n_sheets):
        chunk = items[sheet_idx * per_sheet : (sheet_idx + 1) * per_sheet]
        fig, axes = plt.subplots(cols, cols, figsize=(cols * 3.2, cols * 3.2))
        axes = np.atleast_1d(axes).ravel()
        for tile_idx, item in enumerate(chunk):
            ax = axes[tile_idx]
            # vmin/vmax pin the display to the stretch output itself;
            # without them imshow re-normalises per thumbnail, which maps
            # a constant image (stretch output 0.5 everywhere) to black.
            ax.imshow(
                stretch(item["data"]), cmap="gray", origin="lower", vmin=0.0, vmax=1.0
            )
            ax.set_xticks([])
            ax.set_yticks([])
            title = item["name"]
            stamp = stamp_line(item["header"])
            if stamp:
                title = f"{title}\n{stamp}"
            ax.set_title(title, fontsize=7)
            vectors = compass_vectors(item["header"], item["data"].shape)
            if vectors is not None:
                _draw_compass(ax, vectors[0], vectors[1])
            placements.append((f"sheet_{sheet_idx + 1:03d}.png", tile_idx + 1))
        for ax in axes[len(chunk) :]:
            ax.axis("off")
        fig.tight_layout()
        fig.savefig(out_dir / f"sheet_{sheet_idx + 1:03d}.png", dpi=100)
        plt.close(fig)
    return placements


def make(inputs, out_dir, cols: int) -> dict:
    """Run the whole job. Returns a summary dict; also prints the summary."""
    if not isinstance(cols, int) or isinstance(cols, bool) or cols < 1:
        raise ValueError("cols must be a positive integer")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = collect_inputs(inputs)

    items = []
    rows = []
    skipped_reasons: Counter = Counter()
    rendered_paths = []
    for path in paths:
        try:
            data, header = load_image(path)
        except SkipImage as exc:
            reason = str(exc)
            skipped_reasons[reason] += 1
            print(f"skipped: {path.name} — {reason}")
            rows.append(
                {
                    "source_file": str(path),
                    "sheet": "",
                    "tile": "",
                    "status": "skipped",
                    "reason": reason,
                }
            )
            continue
        items.append({"data": data, "header": header, "name": path.name})
        rendered_paths.append(path)

    placements = render_sheets(items, out, cols) if items else []
    for path, (sheet, tile) in zip(rendered_paths, placements):
        rows.append(
            {
                "source_file": str(path),
                "sheet": sheet,
                "tile": tile,
                "status": "rendered",
                "reason": "",
            }
        )

    # Keep index.csv in input order: rendered rows were appended after
    # skipped ones, so re-sort by the original collection order.
    order = {str(p): i for i, p in enumerate(paths)}
    rows.sort(key=lambda r: order.get(r["source_file"], 0))

    index_path = out / "index.csv"
    with open(index_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["source_file", "sheet", "tile", "status", "reason"]
        )
        writer.writeheader()
        writer.writerows(rows)

    n_sheets = len({s for s, _ in placements})
    summary = {
        "rendered": len(items),
        "skipped": sum(skipped_reasons.values()),
        "skipped_reasons": dict(skipped_reasons),
        "sheets": n_sheets,
        "out_dir": str(out),
    }
    print(f"Rendered: {summary['rendered']} images")
    if summary["skipped"]:
        detail = ", ".join(
            f"{count} {reason}" for reason, count in sorted(skipped_reasons.items())
        )
        print(f"Skipped: {summary['skipped']} ({detail})")
    else:
        print("Skipped: 0")
    print(f"Sheets: {n_sheets} written to {out}/ (plus index.csv)")
    return summary
