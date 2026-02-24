from __future__ import annotations

from contextlib import suppress
from dataclasses import dataclass
from pathlib import Path

import napari
import numpy as np
import pandas as pd
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QDoubleSpinBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)


# ----------------------------
# Helpers
# ----------------------------
def _force_no_border(layer) -> None:
    """Force points layer border/edge to be invisible (robust across napari versions)."""
    with suppress(Exception):
        layer.border_color = [0, 0, 0, 0]
    for attr in ("border_width", "edge_width"):
        if hasattr(layer, attr):
            with suppress(Exception):
                setattr(layer, attr, 0)


def _pick_col(df: pd.DataFrame, candidates) -> str | None:
    for n in candidates:
        if n in df.columns:
            return n
    return None


@dataclass
class MCColumns:
    x: str
    y: str
    gene: str
    val: str


# ----------------------------
# Main Widget
# ----------------------------
class MolecularCartographyViewer(QWidget):
    """Dock widget: load a coordinate table (CSV) and visualize selected genes as points."""

    ALL_NAME = "All_transcripts"

    def __init__(self, viewer: napari.viewer.Viewer):
        super().__init__()
        self.viewer = viewer

        self.df: pd.DataFrame | None = None
        self.cols: MCColumns | None = None
        self.genes_all: pd.Series | None = None
        self.unique_genes: np.ndarray = np.array([], dtype=str)
        self.val_min: float = 0.0
        self.val_max: float = 1.0

        self._gene_cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        self._gene_layers: dict[str, object] = {}
        self._color_cycle = [
            "red",
            "orange",
            "yellow",
            "green",
            "cyan",
            "blue",
            "magenta",
            "white",
        ]

        self._build_ui()

    # ---------- UI ----------
    def _build_ui(self):
        root = QVBoxLayout(self)

        title = QLabel("Molecular Cartography Viewer (Exported CSV)")
        title.setStyleSheet("font-weight: bold;")
        root.addWidget(title)

        # File row
        file_row = QHBoxLayout()
        self.btn_choose = QPushButton("Choose CSV…")
        self.lbl_file = QLabel("No file loaded")
        self.lbl_file.setTextInteractionFlags(Qt.TextSelectableByMouse)
        file_row.addWidget(self.btn_choose)
        file_row.addWidget(self.lbl_file, 1)
        root.addLayout(file_row)

        self.btn_choose.clicked.connect(self.choose_csv)

        # Search + lists
        root.addWidget(QLabel("Search genes (substring):"))
        self.search = QLineEdit()
        self.search.setPlaceholderText("e.g., GLYMA_06G... or BAC45")
        root.addWidget(self.search)

        lists_box = QHBoxLayout()

        left_box = QVBoxLayout()
        left_box.addWidget(QLabel("Candidates (multi-select):"))
        self.candidates = QListWidget()
        self.candidates.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.candidates.setMinimumHeight(320)
        left_box.addWidget(self.candidates)

        right_box = QVBoxLayout()
        right_box.addWidget(QLabel("Selected (multi-select):"))
        self.selected = QListWidget()
        self.selected.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.selected.setMinimumHeight(320)
        right_box.addWidget(self.selected)

        lists_box.addLayout(left_box, 1)
        lists_box.addLayout(right_box, 1)
        root.addLayout(lists_box)

        btn_row = QHBoxLayout()
        self.btn_add = QPushButton("Add →")
        self.btn_remove = QPushButton("← Remove")
        self.btn_clear = QPushButton("Clear")
        btn_row.addWidget(self.btn_add)
        btn_row.addWidget(self.btn_remove)
        btn_row.addWidget(self.btn_clear)
        root.addLayout(btn_row)

        # Display parameters
        params = QGroupBox("Display parameters")
        p = QVBoxLayout(params)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Base size"))
        self.base_size = QSpinBox()
        self.base_size.setRange(1, 200)
        self.base_size.setValue(30)
        row1.addWidget(self.base_size)

        row1.addWidget(QLabel("Opacity"))
        self.opacity = QDoubleSpinBox()
        self.opacity.setRange(0.05, 1.0)
        self.opacity.setSingleStep(0.05)
        self.opacity.setValue(1.0)
        row1.addWidget(self.opacity)
        p.addLayout(row1)

        row2 = QHBoxLayout()
        self.lbl_thr = QLabel("Value threshold")
        row2.addWidget(self.lbl_thr)
        self.val_thr = QDoubleSpinBox()
        self.val_thr.setRange(0.0, 1.0)
        self.val_thr.setSingleStep(1.0)
        self.val_thr.setValue(0.0)
        row2.addWidget(self.val_thr)
        p.addLayout(row2)

        self.size_by_val = QCheckBox("Scale size by value")
        self.size_by_val.setChecked(True)
        p.addWidget(self.size_by_val)

        self.show_all_val = QCheckBox("Ignore threshold (show all values)")
        self.show_all_val.setChecked(False)
        p.addWidget(self.show_all_val)

        self.show_all_layer = QCheckBox("Show gray background layer")
        self.show_all_layer.setChecked(False)
        p.addWidget(self.show_all_layer)

        row3 = QHBoxLayout()
        row3.addWidget(QLabel("Background opacity"))
        self.all_opacity = QDoubleSpinBox()
        self.all_opacity.setRange(0.0, 1.0)
        self.all_opacity.setSingleStep(0.02)
        self.all_opacity.setValue(0.05)
        row3.addWidget(self.all_opacity)
        p.addLayout(row3)

        root.addWidget(params)

        self.btn_update = QPushButton("Update display")
        self.btn_update.setStyleSheet("font-weight: bold;")
        self.btn_apply_ui_color = QPushButton("Apply UI Color")
        root.addWidget(self.btn_update)
        root.addWidget(self.btn_apply_ui_color)

        self.status = QLabel("Load a CSV to begin.")
        self.status.setStyleSheet("font-family: Menlo, monospace;")
        root.addWidget(self.status)

        # Bind events
        self.search.textChanged.connect(self.refresh_candidates)
        self.btn_add.clicked.connect(self.add_selected)
        self.btn_remove.clicked.connect(self.remove_selected)
        self.btn_clear.clicked.connect(self.clear_selected)
        self.btn_update.clicked.connect(self.update_display)
        self.btn_apply_ui_color.clicked.connect(
            self.apply_ui_color_to_active_layer
        )

        # Disable controls until CSV loaded
        self._set_enabled(False)

    def _set_enabled(self, enabled: bool):
        for w in [
            self.search,
            self.candidates,
            self.selected,
            self.btn_add,
            self.btn_remove,
            self.btn_clear,
            self.base_size,
            self.opacity,
            self.val_thr,
            self.size_by_val,
            self.show_all_val,
            self.show_all_layer,
            self.all_opacity,
            self.btn_update,
        ]:
            w.setEnabled(enabled)

    # ---------- CSV loading ----------
    def choose_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Molecular Cartography exported table (CSV)",
            str(Path.home()),
            "CSV files (*.csv);;All files (*.*)",
        )
        if not path:
            return
        self.load_csv(path)

    def load_csv(self, csv_path: str):
        csv_path = str(csv_path)
        df = pd.read_csv(csv_path)

        xcol = _pick_col(df, ["x", "X"])
        ycol = _pick_col(df, ["y", "Y"])
        gcol = _pick_col(df, ["gene", "Gene", "target", "ID"])
        vcol = _pick_col(
            df,
            [
                "val",
                "value",
                "intensity",
                "score",
                "confidence",
                "qc",
                "V",
                "v",
            ],
        )

        missing = [
            k
            for k, v in {
                "x": xcol,
                "y": ycol,
                "gene": gcol,
                "val": vcol,
            }.items()
            if v is None
        ]
        if missing:
            raise ValueError(
                f"CSV missing required columns {missing}. Found columns: {df.columns.tolist()}"
            )

        # Clean numeric values
        vals = pd.to_numeric(df[vcol], errors="coerce")
        df = df.loc[~vals.isna()].copy()
        df[vcol] = df[vcol].astype(float)

        self.df = df
        self.cols = MCColumns(x=xcol, y=ycol, gene=gcol, val=vcol)

        self.genes_all = df[gcol].astype(str)
        self.unique_genes = np.sort(self.genes_all.unique())

        self.val_min = float(df[vcol].min())
        self.val_max = float(df[vcol].max())

        self._gene_cache.clear()
        self._clear_gene_layers()

        self.lbl_file.setText(Path(csv_path).name)
        self.status.setText(
            f"Loaded: rows={len(df)} genes={len(self.unique_genes)} val_range={self.val_min:.3f}–{self.val_max:.3f}"
        )

        # Set threshold widget range
        self.lbl_thr.setText(
            f"Value threshold (range {self.val_min:.2f}–{self.val_max:.2f})"
        )
        self.val_thr.setRange(self.val_min, self.val_max)
        self.val_thr.setValue(self.val_min)

        # Background layer (All_transcripts): lazy-create only when requested
        if self.show_all_layer.isChecked():
            coords_all = df[[ycol, xcol]].to_numpy(np.float32)
            if self.ALL_NAME not in self.viewer.layers:
                all_layer = self.viewer.add_points(
                    coords_all,
                    name=self.ALL_NAME,
                    size=1.5,
                    opacity=float(self.all_opacity.value()),
                    face_color="lightgray",
                )
                _force_no_border(all_layer)
            else:
                lyr = self.viewer.layers[self.ALL_NAME]
                lyr.data = coords_all
                lyr.visible = True
                lyr.opacity = float(self.all_opacity.value())
                _force_no_border(lyr)
        else:
            # Remove stale background layer to save memory / avoid visual residue
            if self.ALL_NAME in self.viewer.layers:
                self.viewer.layers.remove(self.viewer.layers[self.ALL_NAME])

        # startup cleanup: remove All_transcripts when checkbox is off
        if (not self.show_all_layer.isChecked()) and (
            self.ALL_NAME in self.viewer.layers
        ):
            self.viewer.layers.remove(self.viewer.layers[self.ALL_NAME])

        # Init candidates list
        self._all_candidates = list(self.unique_genes)
        self.refresh_candidates("")
        self._set_enabled(True)

    # ---------- gene caching / layers ----------
    def _cache_gene(self, g: str):
        g = str(g)
        if g in self._gene_cache:
            return
        assert (
            self.df is not None
            and self.cols is not None
            and self.genes_all is not None
        )
        m = self.genes_all == g
        coords = self.df.loc[m, [self.cols.y, self.cols.x]].to_numpy(
            np.float32
        )
        vals = self.df.loc[m, self.cols.val].to_numpy(np.float32)
        self._gene_cache[g] = (coords, vals)

    def _clear_gene_layers(self):
        for _g, lyr in list(self._gene_layers.items()):
            if lyr.name in self.viewer.layers:
                self.viewer.layers.remove(lyr)
        self._gene_layers.clear()

    def _ensure_layer(
        self, g: str, color: str, base_size: float, opacity: float
    ):
        lname = f"GENE: {g}"
        if (
            g in self._gene_layers
            and self._gene_layers[g].name in self.viewer.layers
        ):
            lyr = self._gene_layers[g]
            # lyr.face_color = color
            lyr.opacity = float(opacity)
            _force_no_border(lyr)
            return lyr

        lyr = self.viewer.add_points(
            np.zeros((0, 2), np.float32),
            name=lname,
            size=float(base_size),
            opacity=float(opacity),
            face_color=color,
        )
        _force_no_border(lyr)
        self._gene_layers[g] = lyr
        return lyr

    # ---------- list ops ----------
    def refresh_candidates(self, text: str):
        if not hasattr(self, "_all_candidates"):
            return
        text = (text or "").strip().lower()
        self.candidates.clear()
        if text == "":
            show = self._all_candidates
        else:
            show = [g for g in self._all_candidates if text in g.lower()]
        for g in show:
            self.candidates.addItem(QListWidgetItem(g))

    def add_selected(self):
        existing = set(self.get_selected_genes())
        for item in self.candidates.selectedItems():
            g = item.text()
            if g not in existing:
                self.selected.addItem(QListWidgetItem(g))
                existing.add(g)

    def remove_selected(self):
        for item in self.selected.selectedItems():
            row = self.selected.row(item)
            self.selected.takeItem(row)

    def clear_selected(self):
        self.selected.clear()

    def get_selected_genes(self):
        return [
            self.selected.item(i).text() for i in range(self.selected.count())
        ]

    def apply_ui_color_to_active_layer(self):
        """Apply napari UI current face color to all existing points in the active points layer."""
        lyr = self.viewer.layers.selection.active
        if lyr is None:
            self.status.setText("No active layer selected.")
            return
        if not hasattr(lyr, "face_color") or not hasattr(lyr, "data"):
            self.status.setText("Active layer is not a points layer.")
            return

        # 可选：只允许基因层，避免误改背景层（如需允许全部点层可删掉这段）
        if hasattr(lyr, "name") and lyr.name == self.ALL_NAME:
            self.status.setText(
                "All_transcripts selected. Choose a gene layer instead."
            )
            return

        try:
            c = getattr(lyr, "current_face_color", None)
            if c is None:
                fc = getattr(lyr, "face_color", None)
                if fc is not None and len(fc) > 0:
                    c = fc[0]
            if c is None:
                self.status.setText("No UI face color found on active layer.")
                return

            # 内存友好：只有颜色真的变化才广播到所有点
            need_apply = True
            fc = getattr(lyr, "face_color", None)
            if fc is not None and len(fc) > 0:
                try:
                    need_apply = not np.allclose(fc[0], c)
                except (TypeError, ValueError):
                    need_apply = True

            if need_apply:
                lyr.face_color = c
                lyr.refresh()

            npts = len(lyr.data) if hasattr(lyr, "data") else 0
            self.status.setText(
                f"Applied UI color to {lyr.name}; n_points={npts}"
            )
        except (AttributeError, TypeError, ValueError, RuntimeError) as e:
            self.status.setText(f"Apply UI Color failed: {e}")

    # ---------- render ----------
    def update_display(self):
        if self.df is None or self.cols is None:
            return
        sel = self.get_selected_genes()
        if len(sel) == 0:
            self.status.setText("No genes selected.")
            return

        # background layer visibility / lazy-create
        if self.show_all_layer.isChecked():
            if self.ALL_NAME in self.viewer.layers:
                bg = self.viewer.layers[self.ALL_NAME]
                bg.visible = True
                bg.opacity = float(self.all_opacity.value())
                _force_no_border(bg)
            else:
                if self.df is not None and self.cols is not None:
                    coords_all = self.df[[self.cols.y, self.cols.x]].to_numpy(
                        np.float32
                    )
                    bg = self.viewer.add_points(
                        coords_all,
                        name=self.ALL_NAME,
                        size=1.5,
                        opacity=float(self.all_opacity.value()),
                        face_color="lightgray",
                    )
                    _force_no_border(bg)
        else:
            if self.ALL_NAME in self.viewer.layers:
                self.viewer.layers.remove(self.viewer.layers[self.ALL_NAME])

        base_size = float(self.base_size.value())
        opacity = float(self.opacity.value())
        thr = float(self.val_thr.value())
        use_size_by_val = self.size_by_val.isChecked()
        ignore_thr = self.show_all_val.isChecked()

        # self._clear_gene_layers()

        for i, g in enumerate(sel):
            self._cache_gene(g)
            coords, vals = self._gene_cache[g]

            if ignore_thr:
                keep = np.ones_like(vals, dtype=bool)
                vals2 = vals
            else:
                keep = vals >= thr
                vals2 = vals[keep]

            coords2 = coords[keep]
            color = self._color_cycle[i % len(self._color_cycle)]

            lyr = self._ensure_layer(
                g, color=color, base_size=base_size, opacity=opacity
            )
            lyr.data = coords2

            # Apply the currently selected face color (from napari UI) to all existing points.
            # In napari Points (direct mode), changing face color in the UI may only update
            # current_face_color, not the existing face_color array.
            try:
                c = getattr(lyr, "current_face_color", None)
                if c is None:
                    fc = getattr(lyr, "face_color", None)
                    if fc is not None and len(fc) > 0:
                        c = fc[0]
                if c is not None:
                    lyr.face_color = c
            except (AttributeError, TypeError, ValueError):
                pass

            _force_no_border(lyr)

            if use_size_by_val and len(vals2) > 0:
                norm = (vals2 - self.val_min) / (
                    self.val_max - self.val_min + 1e-6
                )
                lyr.size = (0.5 + norm) * base_size
            else:
                lyr.size = base_size

        self.status.setText(
            f"Updated: genes={len(sel)}; threshold={'ALL' if ignore_thr else thr}"
        )


def make_molecular_cartography_viewer_widget():
    """Entry point for napari to create the dock widget."""
    viewer = napari.current_viewer()
    return MolecularCartographyViewer(viewer)
