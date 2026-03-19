# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
# for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.

import os
from xml.etree import ElementTree as ET
from xml.dom.minidom import parseString
from datetime import datetime

from pyqtgraph import FileDialog
from pyqtgraph.Qt.QtWidgets import (
    QDialog, QVBoxLayout, QGridLayout,
    QDialogButtonBox, QTableWidget, QTableWidgetItem,
    QComboBox, QLineEdit, QLabel, QHeaderView,
)
from pyqtgraph.Qt.QtCore import Qt

from .boxes import msgbox
from .panelparser import _CHANNEL_INDEX_TO_COLOR as _CH_COLOR

_COLOR_TO_CHANNEL = {v: k for k, v in _CH_COLOR.items()}

# ── Valid CODIS 3.2 locus names (XSD LocusNameType) ──────────────────────────
# Both TH01/THO1 and TPOX/TP0X are valid per the XSD; keep whatever the panel
# uses as long as it is in this set.
CODIS_LOCI = frozenset([
    "AMEL", "Amelogenin", "CSF1PO", "D13S317", "D16S539",
    "D18S51", "D19S433", "D21S11", "D2S1338", "D3S1358",
    "D5S818", "D7S820", "D8S1179", "FGA", "Penta D", "Penta E",
    "TH01", "THO1", "TP0X", "TPOX", "vWA",
])

# Map common panel name variants not already in CODIS_LOCI to a valid name.
# The map key is the uppercased panel marker name.
_ALIAS_MAP = {
    "AMELOGENIN": "Amelogenin",
    "VWA":        "vWA",
    "PENTA D":    "Penta D",
    "PENTA E":    "Penta E",
}

SPECIMEN_CATEGORIES = [
    "Forensic, Unknown",
    "Suspect, Known",
    "Victim, Known",
    "Convicted Offender",
    "Arrestee",
    "Unidentified Person",
    "Missing Person",
    "Deceased",
    "Elimination, Known",
    "Proficiency",
    "Population",
    "Biological Mother",
    "Biological Father",
    "Biological Sibling",
    "Biological Child",
    "Alleged Mother",
    "Alleged Father",
    "Maternal Relative",
    "Paternal Relative",
    "Deduced Victim Known",
    "Deduced Suspect",
    "Forensic Mixture",
    "CO Duplicate",
    "Staff",
    "Juvenile",
    "Volunteer",
    "Spouse",
    "Legal",
    "Other",
]


def to_codis_locus(name):
    # Return a valid CODIS locus name for *name*, or None if not recognised.
    if name in CODIS_LOCI:
        return name
    upper = name.upper()
    mapped = _ALIAS_MAP.get(upper)
    if mapped:
        return mapped
    # Case-insensitive fallback
    for locus in CODIS_LOCI:
        if locus.upper() == upper:
            return locus
    return None


def extract_loci(state):
    """Return {codis_locus_name: [allele_str, ...]} for one FileState.

    Requires panel data and sized peaks; returns {} otherwise.
    OL (off-ladder) peaks within a locus range are exported with their numeric
    size as the allele value.  ILS peaks and blank labels are skipped.
    """
    if not state.panel_data or len(state.peaksizes) == 0:
        return {}
    panel_name = state.panel_combo.currentText()
    panel = state.panel_data.get(panel_name, {})
    if not panel:
        return {}

    result = {}
    for marker, info in panel.items():
        codis_name = to_codis_locus(marker)
        if codis_name is None:
            continue
        ch_idx = _COLOR_TO_CHANNEL.get(info["dye"].lower())
        if ch_idx is None:
            continue
        alleles = []
        for ch, sz, allele in zip(state.peakchannels, state.peaksizes,
                                  state.peakalleles):
            if int(ch) != ch_idx or (info["min_size"] <= float(sz) <= info["max_size"]):
                continue
            if not allele or allele == "ILS":
                continue
            alleles.append(str(sz) if allele == "OL" else str(allele))
        if alleles:
            result[codis_name] = alleles
    return result


def _pretty_xml(root):
    """Return indented XML string (Python 3.8-compatible via minidom)."""
    raw = ET.tostring(root, encoding="unicode")
    dom = parseString(raw.encode("utf-8"))
    lines = [ln for ln in dom.toprettyxml(indent="  ").splitlines() if ln.strip()]
    return "\n".join(lines)


def build_codis_xml(rows, dest_ori, source_lab, submit_user,
                    submit_dt, batch_id, kit):
    """Build a CODIS 3.2 CMF XML string.

    rows — list of dicts:
        specimen_id  str
        category     str  (one of SPECIMEN_CATEGORIES)
        comment      str  (may be empty)
        loci         dict  {locus_name: [allele_str, ...]}
    """
    NS = "urn:CODISImportFile-schema"
    XSI = "http://www.w3.org/2001/XMLSchema-instance"
    ET.register_namespace("", NS)
    ET.register_namespace("xsi", XSI)

    root = ET.Element(
        f"{{{NS}}}CODISImportFile",
        {f"{{{XSI}}}schemaLocation":
         f"{NS} http://www.ncbi.nlm.nih.gov/projects/SNP/osiris/"
         "CODIS-32.Appendix-B.xsd"},
    )
    ET.SubElement(root, f"{{{NS}}}HEADERVERSION").text = "3.2"
    ET.SubElement(root, f"{{{NS}}}MESSAGETYPE").text = "Import"
    ET.SubElement(root, f"{{{NS}}}DESTINATIONORI").text = dest_ori
    ET.SubElement(root, f"{{{NS}}}SOURCELAB").text = source_lab
    ET.SubElement(root, f"{{{NS}}}SUBMITBYUSERID").text = submit_user
    ET.SubElement(root, f"{{{NS}}}SUBMITDATETIME").text = submit_dt
    if batch_id:
        ET.SubElement(root, f"{{{NS}}}BATCHID").text = batch_id
    if kit:
        ET.SubElement(root, f"{{{NS}}}KIT").text = kit

    for row in rows:
        s_el = ET.SubElement(root, f"{{{NS}}}SPECIMEN")
        ET.SubElement(s_el, f"{{{NS}}}SPECIMENID").text = row["specimen_id"][:24]
        ET.SubElement(s_el, f"{{{NS}}}SPECIMENCATEGORY").text = row["category"]
        if row.get("comment"):
            ET.SubElement(s_el, f"{{{NS}}}SPECIMENCOMMENT").text = \
                row["comment"][:255]
        for locus_name, alleles in row["loci"].items():
            loc_el = ET.SubElement(s_el, f"{{{NS}}}LOCUS")
            ET.SubElement(loc_el, f"{{{NS}}}LOCUSNAME").text = locus_name
            ET.SubElement(loc_el, f"{{{NS}}}READINGBY").text = submit_user
            ET.SubElement(loc_el, f"{{{NS}}}READINGDATETIME").text = submit_dt
            for av in alleles:
                al_el = ET.SubElement(loc_el, f"{{{NS}}}ALLELE")
                ET.SubElement(al_el, f"{{{NS}}}ALLELEVALUE").text = str(av)[:10]

    return _pretty_xml(root)


class CODISExportDialog(QDialog):
    def __init__(self, file_states, tab_names, iface, parent=None):
        super().__init__(parent)
        self._states = file_states
        self._tab_names = tab_names
        self._msg = iface
        self.setWindowTitle(iface["codisdlgtitle"])
        self.setMinimumWidth(700)
        self._build_ui()

    # ── UI construction ──────────────────────────────────────────────────────

    def _build_ui(self):
        root_layout = QVBoxLayout(self)
        root_layout.setSpacing(8)

        # Header / session fields
        grid = QGridLayout()
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(6)

        self._dest_ori = QLineEdit(); self._dest_ori.setMaxLength(10)
        self._source_lab = QLineEdit(); self._source_lab.setMaxLength(10)
        self._submit_user = QLineEdit(); self._submit_user.setMaxLength(20)
        self._batch_id = QLineEdit(); self._batch_id.setMaxLength(32)
        self._kit = QLineEdit(); self._kit.setMaxLength(32)
        self._dt = QLineEdit()
        self._dt.setText(datetime.now().strftime("%Y-%m-%dT%H:%M:%S"))

        fields = [
            (0, 0, self._msg["codisdestori"],   self._dest_ori),
            (0, 2, self._msg["codissourcelab"], self._source_lab),
            (1, 0, self._msg["codisanalyst"],   self._submit_user),
            (1, 2, self._msg["codisbatch"],     self._batch_id),
            (2, 0, self._msg["codiskit"],       self._kit),
            (2, 2, self._msg["codisdt"],        self._dt),
        ]
        for row, col, label, widget in fields:
            grid.addWidget(QLabel(label), row, col)
            grid.addWidget(widget,        row, col + 1)

        root_layout.addLayout(grid)
        root_layout.addWidget(QLabel(self._msg["codisspecimens"]))

        # Specimen table
        self._table = QTableWidget(len(self._states), 4)
        self._table.setHorizontalHeaderLabels([
            "✓",
            self._msg["codisfile"],
            self._msg["codisspecimenid"],
            self._msg["codiscategory"],
        ])
        hdr = self._table.horizontalHeader()
        try:
            hdr.setSectionResizeMode(0, QHeaderView.ResizeToContents)
            hdr.setSectionResizeMode(1, QHeaderView.Stretch)
            hdr.setSectionResizeMode(2, QHeaderView.Stretch)
            hdr.setSectionResizeMode(3, QHeaderView.Stretch)
        except AttributeError:
            hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
            hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
            hdr.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
            hdr.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        self._table.verticalHeader().setVisible(False)

        try:
            _checked = Qt.Checked
            _checkable = Qt.ItemIsUserCheckable | Qt.ItemIsEnabled
            _no_edit = Qt.ItemIsEditable
        except AttributeError:
            _checked = Qt.CheckState.Checked
            _checkable = Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled
            _no_edit = Qt.ItemFlag.ItemIsEditable

        self._cat_combos = []
        for i, (state, name) in enumerate(zip(self._states, self._tab_names)):
            chk_item = QTableWidgetItem()
            chk_item.setCheckState(_checked)
            chk_item.setFlags(_checkable)
            self._table.setItem(i, 0, chk_item)

            name_item = QTableWidgetItem(name)
            name_item.setFlags(name_item.flags() & ~_no_edit)
            self._table.setItem(i, 1, name_item)

            spec_id = os.path.splitext(name)[0][:24]
            self._table.setItem(i, 2, QTableWidgetItem(spec_id))

            cat_combo = QComboBox()
            cat_combo.addItems(SPECIMEN_CATEGORIES)
            self._table.setCellWidget(i, 3, cat_combo)
            self._cat_combos.append(cat_combo)

        root_layout.addWidget(self._table)

        # Warnings for tabs missing panel or sized peaks
        warn_parts = []
        if any(not s.panel_data for s in self._states):
            warn_parts.append(self._msg["codisnopanel"])
        if any(len(s.peaksizes) == 0 for s in self._states):
            warn_parts.append(self._msg["codisnosize"])
        if warn_parts:
            warn_lbl = QLabel("\n".join(warn_parts))
            warn_lbl.setStyleSheet("color: orange;")
            warn_lbl.setWordWrap(True)
            root_layout.addWidget(warn_lbl)

        # Dialog buttons
        try:
            btns = QDialogButtonBox(
                QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            export_btn = btns.button(QDialogButtonBox.Ok)
        except AttributeError:
            btns = QDialogButtonBox(
                QDialogButtonBox.StandardButton.Ok |
                QDialogButtonBox.StandardButton.Cancel)
            export_btn = btns.button(QDialogButtonBox.StandardButton.Ok)
        export_btn.setText(self._msg["codisexport"])
        btns.accepted.connect(self._do_export)
        btns.rejected.connect(self.reject)
        root_layout.addWidget(btns)

    # ── Export logic ─────────────────────────────────────────────────────────

    def _do_export(self):
        dest_ori = self._dest_ori.text().strip()
        source_lab = self._source_lab.text().strip()
        submit_user = self._submit_user.text().strip()
        submit_dt = self._dt.text().strip()

        if not all([dest_ori, source_lab, submit_user, submit_dt]):
            msgbox("", self._msg["codisvalidation"], 1)
            return

        try:
            _checked = Qt.Checked
        except AttributeError:
            _checked = Qt.CheckState.Checked

        rows = []
        for i, (state, name) in enumerate(zip(self._states, self._tab_names)):
            if self._table.item(i, 0).checkState() != _checked:
                continue
            spec_id = (self._table.item(i, 2).text() or "").strip()
            if not spec_id:
                msgbox("", self._msg["codisemptyid"], 1)
                return
            rows.append({
                "specimen_id": spec_id,
                "category":    self._cat_combos[i].currentText(),
                "comment":     "",
                "loci":        extract_loci(state),
            })

        if not rows:
            msgbox("", self._msg["codisnorows"], 1)
            return

        xml_str = build_codis_xml(
            rows, dest_ori, source_lab, submit_user, submit_dt,
            self._batch_id.text().strip(),
            self._kit.text().strip(),
        )

        fname, _ = FileDialog.getSaveFileName(
            self, self._msg["codissave"], "", "CODIS XML (*.xml)")
        if not fname:
            return
        if not fname.lower().endswith(".xml"):
            fname += ".xml"
        with open(fname, "w", encoding="UTF-8") as f:
            f.write(xml_str)
        self.accept()
