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

from os.path import splitext
from xml.etree import ElementTree as ET
from xml.dom.minidom import parseString
from datetime import datetime
from pyqtgraph import FileDialog
from pyqtgraph.Qt.QtWidgets import (
    QDialog, QVBoxLayout, QGridLayout,
    QDialogButtonBox, QTableWidget, QTableWidgetItem,
    QComboBox, QLineEdit, QLabel, QHeaderView, QWidget,
)
from pyqtgraph.Qt.QtCore import Qt

from .boxes import msgbox

# ── Export formats ───────────────────────────────────────────────────────────
# The standard CMF 3.2/3.3 family shares the urn:CODISImportFile-schema layout
# (the only difference is HEADERVERSION and the set of valid loci); the Rapid
# Import CMF is a distinct urn:CODISRapidImportFile-schema layout described in
# the CODIS Rapid Import CMF Interface Specification (R17).
FORMAT_CMF_32 = "CMF 3.2"
FORMAT_CMF_33 = "CMF 3.3"
FORMAT_RAPID = "Rapid Import CMF"
EXPORT_FORMATS = (FORMAT_CMF_32, FORMAT_CMF_33, FORMAT_RAPID)

# ── Valid CMF 3.2 locus names (CODIS-32 Appendix-B LocusNameType) ────────────
# Both TH01/THO1 and TPOX/TP0X are valid per the 3.2 XSD; keep whatever the
# panel uses as long as it is in this set.
CODIS_LOCI_32 = frozenset([
    "AMEL", "Amelogenin", "CSF1PO", "D13S317", "D16S539",
    "D18S51", "D19S433", "D21S11", "D2S1338", "D3S1358",
    "D5S818", "D7S820", "D8S1179", "FGA", "Penta D", "Penta E",
    "TH01", "THO1", "TP0X", "TPOX", "vWA",
])

# ── Valid CMF 3.3 / Rapid Import locus names (R17 Appendix B LocusNameType,
#    a.k.a. STR Import CMF 3.3 loci; Appendix C) ─────────────────────────────
_AUTOSOMAL_33 = frozenset([
    "Amelogenin", "CSF1PO", "D10S1248", "D12S391", "D13S317", "D16S539",
    "D18S51", "D19S433", "D1S1656", "D21S11", "D22S1045", "D2S1338",
    "D2S441", "D3S1358", "D5S818", "D6S1043", "D7S820", "D8S1179", "FGA",
    "Penta D", "Penta E", "SE33", "TH01", "TPOX", "vWA",
])
_YSTR_33 = frozenset([
    "DYF387S1", "DYS19", "DYS385", "DYS389 I", "DYS389 II", "DYS390",
    "DYS391", "DYS392", "DYS393", "DYS437", "DYS438", "DYS439", "DYS448",
    "DYS449", "DYS456", "DYS458", "DYS460", "DYS481", "DYS518", "DYS533",
    "DYS549", "DYS570", "DYS576", "DYS627", "DYS635", "DYS643", "YGATAH4",
    "Yindel",
])
CODIS_LOCI_33 = _AUTOSOMAL_33 | _YSTR_33

# Superset of all recognised loci — used for permissive aliasing/lookups.
CODIS_LOCI = CODIS_LOCI_32 | CODIS_LOCI_33

# Map common panel marker-name variants to a canonical CODIS name. The map key
# is the uppercased panel marker name; the mapped value is only used when it is
# valid for the targeted CMF version.
_ALIAS_MAP = {
    "AMELOGENIN": "Amelogenin",
    "AMEL":       "Amelogenin",
    "VWA":        "vWA",
    "PENTA D":    "Penta D",
    "PENTAD":     "Penta D",
    "PENTA E":    "Penta E",
    "PENTAE":     "Penta E",
    "THO1":       "TH01",
    "TP0X":       "TPOX",
    "DYS389I":    "DYS389 I",
    "DYS389 1":   "DYS389 I",
    "DYS389II":   "DYS389 II",
    "DYS389 2":   "DYS389 II",
}

# ── Valid specimen categories ────────────────────────────────────────────────
# Standard CMF (CODIS-32 Appendix-B SpecimenCategoryType).
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

# Rapid Import CMF (R17 Appendix D / SpecimenCategoryType enumeration).
RAPID_SPECIMEN_CATEGORIES = [
    "Arrestee",
    "Convicted Offender",
    "Detainee",
    "Juvenile",
    "Legal",
]


def valid_loci_for(export_format):
    # Return the frozenset of valid locus names for an export format.
    return CODIS_LOCI_32 if export_format == FORMAT_CMF_32 else CODIS_LOCI_33


def to_codis_locus(name, valid=CODIS_LOCI):
    # Return a valid CODIS locus name for *name* within *valid*, or None.
    if name in valid:
        return name
    upper = name.upper()
    mapped = _ALIAS_MAP.get(upper)
    if mapped and mapped in valid:
        return mapped
    # Case-insensitive fallback against the target set.
    for locus in valid:
        if locus.upper() == upper:
            return locus
    return None


def _allele_sort_key(value):
    # SWGDAM-ish ordering: numeric repeats ascending, then everything else.
    try:
        return (0, float(value), "")
    except (TypeError, ValueError):
        return (1, 0.0, str(value))


def extract_loci(state, valid_loci=CODIS_LOCI):
    # Return {codis_locus_name: [allele_str, ...]} for one FileState.
    #
    # Peak labels use the "MARKER:ALLELE" format produced by allele binning;
    # ILS, OL and blank labels are skipped. Marker names are mapped to a valid
    # CODIS locus for the requested CMF version (markers that do not map are
    # skipped). Duplicate (homozygous) alleles are collapsed, alleles are sorted
    # by nomenclature, and each locus is capped at the CODIS maximum of 8.
    if len(state.peaksizes) == 0 or not state.peakalleles:
        return {}

    heights = state.peakheights
    by_locus = {}
    for i, label in enumerate(state.peakalleles):
        if not label or label in ("OL", "ILS") or ":" not in label:
            continue
        marker, allele = label.split(":", 1)
        if not allele:
            continue
        codis_name = to_codis_locus(marker, valid_loci)
        if codis_name is None:
            continue
        height = float(heights[i]) if i < len(heights) else 0.0
        by_locus.setdefault(codis_name, {})
        # Keep the tallest occurrence of each distinct allele value.
        if allele not in by_locus[codis_name] or height > by_locus[codis_name][allele]:
            by_locus[codis_name][allele] = height

    result = {}
    for locus, allele_heights in by_locus.items():
        # Cap at 8 alleles per locus (CODIS limit): keep the tallest, then
        # restore nomenclature order for the output.
        tallest = sorted(allele_heights.items(),
                         key=lambda kv: kv[1], reverse=True)[:8]
        alleles = sorted((a for a, _ in tallest), key=_allele_sort_key)
        result[locus] = alleles
    return result


def _pretty_xml(root):
    # Return indented XML string (Python 3.8-compatible via minidom).
    raw = ET.tostring(root, encoding="unicode")
    dom = parseString(raw.encode("utf-8"))
    lines = [ln for ln in dom.toprettyxml(indent="  ").splitlines() if ln.strip()]
    return "\n".join(lines)


def build_codis_xml(rows, dest_ori, source_lab, submit_user,
                    submit_dt, batch_id, kit, version="3.2"):
    # Build a standard CODIS Import CMF XML string (urn:CODISImportFile-schema).
    #
    # version — "3.2" or "3.3"; only HEADERVERSION and the set of valid loci
    #           (validated by the caller) differ between the two.
    # rows — list of dicts:
    #     specimen_id  str
    #     category     str  (one of SPECIMEN_CATEGORIES)
    #     comment      str  (may be empty)
    #     loci         dict  {locus_name: [allele_str, ...]}
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
    ET.SubElement(root, f"{{{NS}}}HEADERVERSION").text = version
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


def build_rapid_cmf_xml(rows, dest_ori, source_ori, msg_creator, msg_dt,
                        msg_id, instrument, batch_id, kit):
    # Build a CODIS Rapid Import CMF XML string
    # (urn:CODISRapidImportFile-schema, R17 Appendix B).
    #
    # instrument — dict: id (required), manufacturer, model, software_version.
    # rows — list of dicts (standard keys plus the Rapid-required fields):
    #     specimen_id          str   (required)
    #     category             str   (one of RAPID_SPECIMEN_CATEGORIES)
    #     unique_event         str   (required, UNIQUEEVENTID)
    #     fingerprint_dt       str   (required, FINGERPRINTDATE)
    #     arrest_offense       str   (required, ARRESTOFFENSECATEGORY)
    #     comment              str   (optional, SPECIMENCOMMENT)
    #     sid / ucn / arrest_dt / booking_id / arresting_id  (optional)
    #     loci                 dict  {locus_name: [allele_str, ...]}
    NS = "urn:CODISRapidImportFile-schema"
    ET.register_namespace("", NS)

    def sub(parent, tag, text):
        el = ET.SubElement(parent, f"{{{NS}}}{tag}")
        el.text = text
        return el

    root = ET.Element(f"{{{NS}}}CODISRapidImportFile")

    header = ET.SubElement(root, f"{{{NS}}}HEADER")
    sub(header, "MESSAGEVERSION", "1.0")
    sub(header, "MESSAGETYPE", "Rapid Import")
    sub(header, "MESSAGEID", str(msg_id))
    sub(header, "MESSAGEDATETIME", msg_dt)
    sub(header, "MSGCREATORUSERID", msg_creator[:20])
    sub(header, "DESTINATIONORI", dest_ori[:10])
    sub(header, "SOURCEORI", source_ori[:10])

    device = ET.SubElement(root, f"{{{NS}}}DEVICE")
    sub(device, "INSTRUMENTID", instrument["id"][:32])
    if instrument.get("manufacturer"):
        sub(device, "MANUFACTURER", instrument["manufacturer"][:32])
    if instrument.get("model"):
        sub(device, "MODEL", instrument["model"][:32])
    if instrument.get("software_version"):
        sub(device, "SOFTWAREVERSION", instrument["software_version"][:32])

    for row in rows:
        s_el = ET.SubElement(root, f"{{{NS}}}SPECIMEN")
        sub(s_el, "SPECIMENID", row["specimen_id"][:24])
        sub(s_el, "SPECIMENCATEGORY", row["category"])
        if row.get("sid"):
            sub(s_el, "SID", row["sid"][:32])
        if row.get("ucn"):
            sub(s_el, "FBI_NUMBER_UCN", row["ucn"][:9])
        sub(s_el, "UNIQUEEVENTID", row["unique_event"][:32])
        if row.get("booking_id"):
            sub(s_el, "BOOKINGCUSTOMID", row["booking_id"][:32])
        if row.get("arresting_id"):
            sub(s_el, "ARRESTINGCUSTOMID", row["arresting_id"][:32])
        if row.get("arrest_dt"):
            sub(s_el, "ARRESTDATE", row["arrest_dt"])
        sub(s_el, "FINGERPRINTDATE", row["fingerprint_dt"])
        sub(s_el, "ARRESTOFFENSECATEGORY", row["arrest_offense"][:300])
        if row.get("comment"):
            sub(s_el, "SPECIMENCOMMENT", row["comment"][:512])
        for locus_name, alleles in row["loci"].items():
            loc_el = ET.SubElement(s_el, f"{{{NS}}}LOCUS")
            sub(loc_el, "LOCUSNAME", locus_name)
            if kit:
                sub(loc_el, "KIT", kit[:32])
            if batch_id:
                sub(loc_el, "BATCHID", batch_id[:32])
            for av in alleles:
                al_el = ET.SubElement(loc_el, f"{{{NS}}}ALLELE")
                sub(al_el, "ALLELEVALUE", str(av)[:10])

    return _pretty_xml(root)


class CODISExportDialog(QDialog):
    def __init__(self, file_states, tab_names, iface, parent=None):
        super().__init__(parent)
        self._states = file_states
        self._tab_names = tab_names
        self._msg = iface
        self.setWindowTitle(self._t("codiscmftitle", "CODIS CMF Export"))
        self.setMinimumWidth(760)
        self._build_ui()
        self._on_format_changed()

    # ── Localisation helper ────────────────────────────────────────────────
    def _t(self, key, default):
        # New CMF-3.3/Rapid strings fall back to English when a translation is
        # absent for the active language.
        return self._msg.get(key, default)

    # ── Qt enum helpers (PyQt5 vs PyQt6/PySide6) ───────────────────────────
    @staticmethod
    def _checked_state():
        try:
            return Qt.Checked
        except AttributeError:
            return Qt.CheckState.Checked

    # ── UI construction ────────────────────────────────────────────────────
    def _build_ui(self):
        root_layout = QVBoxLayout(self)
        root_layout.setSpacing(8)

        # Format selector
        fmt_grid = QGridLayout()
        fmt_grid.setHorizontalSpacing(8)
        self._format = QComboBox()
        self._format.addItems(EXPORT_FORMATS)
        self._format.currentIndexChanged.connect(self._on_format_changed)
        fmt_grid.addWidget(QLabel(self._t("codisformat", "Format:")), 0, 0)
        fmt_grid.addWidget(self._format, 0, 1)
        fmt_grid.setColumnStretch(2, 1)
        root_layout.addLayout(fmt_grid)

        # Header / session fields (shared by all formats)
        grid = QGridLayout()
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(6)

        self._dest_ori = QLineEdit()
        self._dest_ori.setMaxLength(10)
        self._source_lab = QLineEdit()
        self._source_lab.setMaxLength(10)
        self._submit_user = QLineEdit()
        self._submit_user.setMaxLength(20)
        self._batch_id = QLineEdit()
        self._batch_id.setMaxLength(32)
        self._kit = QLineEdit()
        self._kit.setMaxLength(32)
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

        # Rapid-only header block (Message ID + device fields)
        self._rapid_box = QWidget()
        rgrid = QGridLayout(self._rapid_box)
        rgrid.setContentsMargins(0, 0, 0, 0)
        rgrid.setHorizontalSpacing(8)
        rgrid.setVerticalSpacing(6)
        self._msg_id = QLineEdit("1")
        self._instr_id = QLineEdit()
        self._instr_id.setMaxLength(32)
        self._manuf = QLineEdit()
        self._manuf.setMaxLength(32)
        self._model = QLineEdit()
        self._model.setMaxLength(32)
        self._softver = QLineEdit()
        self._softver.setMaxLength(32)
        rapid_fields = [
            (0, 0, self._t("codismsgid", "Message ID:"),         self._msg_id),
            (0, 2, self._t("codisinstrid", "Instrument ID:"),    self._instr_id),
            (1, 0, self._t("codismanuf", "Manufacturer:"),       self._manuf),
            (1, 2, self._t("codismodel", "Model:"),              self._model),
            (2, 0, self._t("codissoftver", "Software version:"), self._softver),
        ]
        for row, col, label, widget in rapid_fields:
            rgrid.addWidget(QLabel(label), row, col)
            rgrid.addWidget(widget,        row, col + 1)
        rgrid.setColumnStretch(1, 1)
        rgrid.setColumnStretch(3, 1)
        root_layout.addWidget(self._rapid_box)

        root_layout.addWidget(QLabel(self._msg["codisspecimens"]))

        # Specimen table — columns 4–6 are Rapid-only.
        self._COL_UNIQUE, self._COL_FPDT, self._COL_OFFENSE = 4, 5, 6
        self._table = QTableWidget(len(self._states), 7)
        self._table.setHorizontalHeaderLabels([
            "✓",
            self._msg["codisfile"],
            self._msg["codisspecimenid"],
            self._msg["codiscategory"],
            self._t("codisuniqueevent", "Unique Event ID"),
            self._t("codisfingerprintdt", "Fingerprint date/time"),
            self._t("codisarrestoffense", "Arrest offense"),
        ])
        hdr = self._table.horizontalHeader()
        try:
            _to_contents = QHeaderView.ResizeToContents
            _stretch = QHeaderView.Stretch
        except AttributeError:
            _to_contents = QHeaderView.ResizeMode.ResizeToContents
            _stretch = QHeaderView.ResizeMode.Stretch
        hdr.setSectionResizeMode(0, _to_contents)
        for c in range(1, 7):
            hdr.setSectionResizeMode(c, _stretch)
        self._table.verticalHeader().setVisible(False)

        try:
            _checked = Qt.Checked
            _checkable = Qt.ItemIsUserCheckable | Qt.ItemIsEnabled
            _no_edit = Qt.ItemIsEditable
        except AttributeError:
            _checked = Qt.CheckState.Checked
            _checkable = Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled
            _no_edit = Qt.ItemFlag.ItemIsEditable

        dt_default = self._dt.text()
        self._cat_combos = []
        for i, (state, name) in enumerate(zip(self._states, self._tab_names)):
            chk_item = QTableWidgetItem()
            chk_item.setCheckState(_checked)
            chk_item.setFlags(_checkable)
            self._table.setItem(i, 0, chk_item)

            name_item = QTableWidgetItem(name)
            name_item.setFlags(name_item.flags() & ~_no_edit)
            self._table.setItem(i, 1, name_item)

            spec_id = splitext(name)[0][:24]
            self._table.setItem(i, 2, QTableWidgetItem(spec_id))

            cat_combo = QComboBox()
            cat_combo.addItems(SPECIMEN_CATEGORIES)
            self._table.setCellWidget(i, 3, cat_combo)
            self._cat_combos.append(cat_combo)

            # Rapid-only cells, prefilled with sensible defaults.
            self._table.setItem(i, self._COL_UNIQUE, QTableWidgetItem(spec_id))
            self._table.setItem(i, self._COL_FPDT, QTableWidgetItem(dt_default))
            self._table.setItem(i, self._COL_OFFENSE, QTableWidgetItem(""))

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

    # ── Format switching ───────────────────────────────────────────────────
    def _on_format_changed(self, *_):
        is_rapid = self._format.currentText() == FORMAT_RAPID
        self._rapid_box.setVisible(is_rapid)
        for col in (self._COL_UNIQUE, self._COL_FPDT, self._COL_OFFENSE):
            self._table.setColumnHidden(col, not is_rapid)

        # Specimen-category list differs between standard and Rapid CMF.
        categories = RAPID_SPECIMEN_CATEGORIES if is_rapid else SPECIMEN_CATEGORIES
        for combo in self._cat_combos:
            previous = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(categories)
            if previous in categories:
                combo.setCurrentText(previous)
            combo.blockSignals(False)

    # ── Export logic ─────────────────────────────────────────────────────────
    def _selected_rows(self):
        # Yield (table_index, FileState) for every checked, non-empty-loci row.
        checked = self._checked_state()
        for i, state in enumerate(self._states):
            if self._table.item(i, 0).checkState() == checked:
                yield i, state

    def _cell_text(self, row, col):
        item = self._table.item(row, col)
        return (item.text() if item else "").strip()

    def _do_export(self):
        export_format = self._format.currentText()
        if export_format == FORMAT_RAPID:
            xml_str = self._build_rapid()
        else:
            xml_str = self._build_standard(export_format)
        if xml_str is None:
            return  # a validation error was already reported

        fname, _ = FileDialog.getSaveFileName(
            self, self._msg["codissave"], "", "CODIS XML (*.xml)")
        if not fname:
            return
        if not fname.lower().endswith(".xml"):
            fname += ".xml"
        with open(fname, "w", encoding="UTF-8") as f:
            f.write(xml_str)
        self.accept()

    def _build_standard(self, export_format):
        dest_ori = self._dest_ori.text().strip()
        source_lab = self._source_lab.text().strip()
        submit_user = self._submit_user.text().strip()
        submit_dt = self._dt.text().strip()
        if not all([dest_ori, source_lab, submit_user, submit_dt]):
            msgbox("", self._msg["codisvalidation"], 1)
            return None

        valid_loci = valid_loci_for(export_format)
        rows = []
        for i, state in self._selected_rows():
            spec_id = self._cell_text(i, 2)
            if not spec_id:
                msgbox("", self._msg["codisemptyid"], 1)
                return None
            loci = extract_loci(state, valid_loci)
            if not loci:
                continue  # a specimen with no loci is not a valid CMF specimen
            rows.append({
                "specimen_id": spec_id,
                "category": self._cat_combos[i].currentText(),
                "comment": "",
                "loci": loci,
            })

        if not rows:
            msgbox("", self._msg["codisnorows"], 1)
            return None

        version = "3.3" if export_format == FORMAT_CMF_33 else "3.2"
        return build_codis_xml(
            rows, dest_ori, source_lab, submit_user, submit_dt,
            self._batch_id.text().strip(), self._kit.text().strip(),
            version=version)

    def _build_rapid(self):
        dest_ori = self._dest_ori.text().strip()
        source_ori = self._source_lab.text().strip()
        msg_creator = self._submit_user.text().strip()
        msg_dt = self._dt.text().strip()
        msg_id = self._msg_id.text().strip()
        instr_id = self._instr_id.text().strip()

        rapid_msg = self._t(
            "codisrapidvalidation",
            "Rapid Import requires Destination ORI, Source ORI, Message "
            "Creator User ID, Message date/time, Message ID and Instrument "
            "ID, plus per specimen: Unique Event ID, Fingerprint date/time "
            "and Arrest offense.")
        if not all([dest_ori, source_ori, msg_creator, msg_dt, instr_id]) \
                or not msg_id.isdigit() or int(msg_id) < 1:
            msgbox("", rapid_msg, 1)
            return None

        instrument = {
            "id": instr_id,
            "manufacturer": self._manuf.text().strip(),
            "model": self._model.text().strip(),
            "software_version": self._softver.text().strip(),
        }

        rows = []
        for i, state in self._selected_rows():
            spec_id = self._cell_text(i, 2)
            if not spec_id:
                msgbox("", self._msg["codisemptyid"], 1)
                return None
            unique_event = self._cell_text(i, self._COL_UNIQUE)
            fingerprint_dt = self._cell_text(i, self._COL_FPDT)
            arrest_offense = self._cell_text(i, self._COL_OFFENSE)
            if not all([unique_event, fingerprint_dt, arrest_offense]):
                msgbox("", rapid_msg, 1)
                return None
            loci = extract_loci(state, CODIS_LOCI_33)
            if not loci:
                continue
            rows.append({
                "specimen_id": spec_id,
                "category": self._cat_combos[i].currentText(),
                "unique_event": unique_event,
                "fingerprint_dt": fingerprint_dt,
                "arrest_offense": arrest_offense,
                "comment": "",
                "loci": loci,
            })

        if not rows:
            msgbox("", self._msg["codisnorows"], 1)
            return None

        return build_rapid_cmf_xml(
            rows, dest_ori, source_ori, msg_creator, msg_dt, int(msg_id),
            instrument, self._batch_id.text().strip(),
            self._kit.text().strip())
