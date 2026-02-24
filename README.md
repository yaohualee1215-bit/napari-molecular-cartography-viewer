# napari-molecular-cartography-viewer

A napari dock widget to **visualize exported molecular coordinate tables (CSV)** as transcript points.
It supports **multi-gene overlay**, value thresholding, and an optional gray background layer showing all transcripts.

> This plugin is intended for **exported/parsed coordinate tables** (CSV).
> It does **not** require or reverse-engineer any proprietary file formats.

## Compatibility

- **Python:** 3.10–3.13 (tested on 3.11)
- **napari:** 0.6.x (tested on 0.6.6)

## CSV format requirements

The plugin auto-detects **4 required fields**. Column names may vary, but must match one of the accepted names below.

### Required columns

| Field | Meaning | Accepted column names |
|------|---------|------------------------|
| `x`  | X coordinate (pixel) | `x`, `X` |
| `y`  | Y coordinate (pixel) | `y`, `Y` |
| `gene` | Gene / target name | `gene`, `Gene`, `target`, `ID` |
| `val` | Numeric value (intensity/score/confidence) | `val`, `value`, `intensity`, `score`, `confidence`, `qc`, `V`, `v` |

Notes:
- napari points are rendered in **(y, x)** order internally.
- `val` must be numeric; non-numeric/NaN rows are automatically dropped.
- The CSV must contain a header row.

### Minimal example

```csv
x,y,gene,val
120.5,88.2,GLYMA_01G000100,12.3
121.0,88.9,GLYMA_01G000100,8.1
500.2,410.7,NOD26,30.0
```

## Usage

1. Start napari.
2. Open the viewer: **Plugins → Molecular Cartography Viewer**
3. Click **Choose CSV…** and select your exported/parsed `.csv`.
4. Search genes by substring, multi-select candidates, click **Add →**
5. Click **Update display** to render selected genes.

### New in v0.1.1

- Added **Apply UI Color** button for active gene layers.
- Workflow for manual color changes in napari:
  1. Select a `GENE: ...` layer in the **Layers panel** (bottom-left).
  2. Change **face color** in **layer controls** (top-left).
  3. Click **Apply UI Color** in the plugin panel to apply the color to existing points.
- Added transcript/point count display for the active/selected layer (shown in the lower-right status area).
- Improved color update behavior for large transcript point layers.
- Improved background transcript layer handling to reduce unnecessary heavy loading in common use.

## Options

- **Show gray background layer**: toggles `All_transcripts`
- **Value threshold**: keep points with `val >= threshold` (unless “Ignore threshold” is checked)
- **Ignore threshold**: show all points for selected genes
- **Scale size by value**: point size mapped to `val`
- **Opacity / Base size**: affects per-gene layers
- All point layers are forced to render with **no borders** for clean visualization.

## Installation

### From PyPI

```bash
pip install napari-molecular-cartography-viewer
```

If you need a full napari install in a fresh environment:

```bash
pip install "napari[all]" napari-molecular-cartography-viewer
```

## Notes

### Reading LZW-compressed TIFF images (optional)
If your background image is a TIFF with **LZW compression**, install `imagecodecs`:

```bash
conda install -c conda-forge imagecodecs
```

## License

BSD-3-Clause. See `LICENSE`.
