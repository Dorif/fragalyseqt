# FFT Low-Pass Denoising — Implementation Plan

## Source algorithm

`Fragman/R/transfft.R` — brick-wall low-pass FFT filter:

```r
transfft <- function(sn, top=0.3){
  sn.fft = fft(sn)
  qq <- length(sn) * top
  sn.fft[qq:(length(sn)-qq)] = 0 + 0i
  sn.ifft = fft(sn.fft, inverse=TRUE) / length(sn.fft)
  return(Re(sn.ifft))
}
```

Keeps the lowest `top` fraction of the frequency spectrum on each side;
zeros the rest. With `top=0.3` and `Rate1=15 Hz`: cutoff = 4.5 Hz.

## Python translation

Use `rfft`/`irfft` instead of full `fft` — real-valued input only needs
the positive-frequency half, ~2× faster, no manual conjugate-mirror handling.

Fragman equivalence: `cutoff_hz = top * scan_rate_hz`

```python
def _fft_lowpass(signal, cutoff_hz, scan_rate_hz):
    n = len(signal)
    spectrum = rfft(signal)
    freqs = rfftfreq(n, d=1.0 / scan_rate_hz)
    spectrum[freqs > cutoff_hz] = 0.0
    return irfft(spectrum, n=n)
```

Default cutoff: **3.0 Hz** — more conservative than Fragman's 4.5 Hz at 15 Hz
scan rate; preserves all genuine peaks (seconds wide) while removing detector
shot noise concentrated above ~3 Hz.

**Order rationale:** jBCD first → FFT second is the safe default because the
brick-wall FFT filter produces Gibbs ringing at sharp transitions; jBCD seeing
ringing artifacts would distort its baseline estimate. However, FFT first →
jBCD second may perform better on signals with a very noisy baseline that
impedes jBCD convergence. Both orders are exposed for testing.

---

## Implementation steps

### Step 1 — Add numpy.fft imports in `fragalyseqt.py`

```python
# extend the existing numpy import line:
from numpy import around, multiply, array, concatenate, transpose, where
from numpy.fft import rfft, irfft, rfftfreq
```

Three specific functions only — no new library.

### Step 2 — Add module-level helpers after `_refine_peak_positions`

```python
_FFT_CUTOFF_HZ = 3.0


def _fft_lowpass(signal, cutoff_hz, scan_rate_hz):
    n = len(signal)
    spectrum = rfft(signal)
    freqs = rfftfreq(n, d=1.0 / scan_rate_hz)
    spectrum[freqs > cutoff_hz] = 0.0
    return irfft(spectrum, n=n)


def _get_scan_rate(abif_raw):
    return abif_raw.get("Rate1", 15)
```

### Step 3 — Replace the BCD checkbox with a ComboBox in `FileState`

`FileState.do_BCD` (bool) → `FileState.denoise_mode` (str).

```python
# in FileState.__init__:
self.denoise_mode = "off"
self.bcd = None          # widget reference; now holds ComboBox, not QCheckBox
```

### Step 4 — Replace QCheckBox with ComboBox in `_create_tab_content`

`ComboBox` is already imported from pyqtgraph.

```python
bcd_label = QLabel(ifacemsg["bcd"])
bcd_label.setStyleSheet(''' font-size: 10pt; ''')
controls_layout.addWidget(bcd_label, 8, 0)
bcd = ComboBox()
bcd.addItems(["Off", "jBCD", "FFT", "jBCD \u2192 FFT", "FFT \u2192 jBCD"])
bcd.currentIndexChanged.connect(self.setbcd)
controls_layout.addWidget(bcd, 8, 1)
```

`QLabel` — check if already imported; add to Qt imports if not.

### Step 5 — Update `setbcd()`

```python
def setbcd(self):
    s = self._state
    if s is None:
        return
    mode_map = {0: "off", 1: "jbcd", 2: "fft", 3: "jbcd_fft", 4: "fft_jbcd"}
    s.denoise_mode = mode_map[self.sender().currentIndex()]
    self.reanalyse()
```

### Step 6 — Update `findpeaks()` denoising dispatch

Replace both `if s.do_BCD:` blocks with a single `_denoise()` helper.
`_sr` is computed once and captured by closure — no repeated dict lookups.

```python
_sr = _get_scan_rate(s.abif_raw)
half_win = (s.winwidth - 1) // 2

def _denoise(sig):
    mode = s.denoise_mode
    if mode == "jbcd":
        _, p = jbcd(sig, half_window=half_win)
        return list(p['signal'])
    if mode == "fft":
        return list(_fft_lowpass(sig, _FFT_CUTOFF_HZ, _sr))
    if mode == "jbcd_fft":
        _, p = jbcd(sig, half_window=half_win)
        return list(_fft_lowpass(p['signal'], _FFT_CUTOFF_HZ, _sr))
    if mode == "fft_jbcd":
        sig = _fft_lowpass(sig, _FFT_CUTOFF_HZ, _sr)
        _, p = jbcd(sig, half_window=half_win)
        return list(p['signal'])
    return list(sig)   # "off"
```

ILS channel:
```python
if s.denoise_mode != "off":
    ils_data = _denoise(ils_data)
```

Sample channels (replaces `_bcd_channel` closure):
```python
def _denoise_channel(chnum):
    return _denoise(s.abif_raw[s.udatac[chnum]])

if s.denoise_mode != "off":
    with ThreadPoolExecutor() as executor:
        s.ch = list(executor.map(_denoise_channel, s.dyerange))
else:
    for chnum in s.dyerange:
        s.ch.append(list(s.abif_raw[s.udatac[chnum]]))
```

### Step 7 — Localization in `localize.py`

The existing `"bcd"` key becomes a ComboBox section label. Update its value
in all 6 languages from the current checkbox label to e.g. "Denoising" /
"Шумоподавление" / etc.

---

## What does NOT change

- `reanalyse()`, `replot()`, `export_csv()` — untouched
- jBCD parameters — untouched
- Peak detection logic — untouched
- `_FFT_CUTOFF_HZ` is a module-level constant; can become a SpinBox later

---

## Rate1 field notes

| Instrument               | Rate1 (scans/sec) | Nyquist  | Cutoff at 3 Hz |
|--------------------------|-------------------|----------|----------------|
| ABI 310 / 3100 / 3130    | 15                | 7.5 Hz   | keeps 40%      |
| ABI 3500 / SeqStudio     | 15–20             | ~10 Hz   | keeps 30–40%   |
| RapidHIT ID              | 15                | 7.5 Hz   | keeps 40%      |
| CEQ 8000 (when supported)| unknown (est. 17) | ~8.5 Hz  | keeps ~35%     |

`abif_raw.get("Rate1", 15)` is safe for files where the tag is absent.

---

## Files changed

| File | Change |
|------|--------|
| `src/fragalyseqt/fragalyseqt.py` | +1 import line, +2 helpers, +1 constant, checkbox→ComboBox, `setbcd` updated, `findpeaks` dispatch rewritten |
| `src/fragalyseqt/localize.py` | `"bcd"` label text updated in all 6 languages |
