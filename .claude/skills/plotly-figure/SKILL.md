---
name: plotly-figure
description: Build publication-quality interactive Plotly.js charts with figure controls (line width, marker, fonts, axis ranges, aspect ratio, dark/light/pub themes, PNG/SVG export). Use when the user needs to create or configure scientific plots in a web UI.
---

# Plotly.js Publication-Quality Figure Skill

When building interactive scientific plots with Plotly.js in a web UI, follow these patterns and best practices.

## Architecture: Reactive Settings Pattern

All figure controls should trigger a single `onSettingChange()` handler that re-renders the entire chart with current settings. This avoids partial update bugs.

```
User changes control → onSettingChange() → getSettings() → renderResults(data)
                                                ↓
                                    buildTraces() + buildLayout()
                                                ↓
                                        Plotly.react(div, traces, layout, config)
```

### 1. Settings Object

Collect all figure settings from UI controls into a single object:

```javascript
function getSettings() {
  const ni = id => { const v = parseFloat($(id)?.value); return isNaN(v) ? null : v; };
  return {
    lineWidth:   ni('inp_lw')       ?? 2.5,   // px
    markerSize:  ni('inp_ms')       ?? 8,     // px
    markerEvery: ni('inp_mi')       ?? 20,    // show 1 marker per N data points
    frameWidth:  ni('inp_fw')       ?? 1.5,   // axis border width
    fsTitle:     ni('inp_fs_title') ?? 14,    // pt
    fsAxis:      ni('inp_fs_axis')  ?? 12,    // pt
    fsTick:      ni('inp_fs_tick')  ?? 11,    // pt
    fsLeg:       ni('inp_fs_leg')   ?? 11,    // pt
    fontWeight:  parseInt($('sel_fw')?.value ?? '400'),
    exportWidth: ni('inp_er')       ?? 1600,  // px for export
    pubMode:     $('chk_pub')?.checked ?? false,
    aspectRatio: $('sel_ar')?.value ?? '4/3',
    // axis range overrides (null = auto)
    xMin: ni('inp_xmin'), xMax: ni('inp_xmax'),
    yMin: ni('inp_ymin'), yMax: ni('inp_ymax'),
  };
}
```

### 2. Trace Building

Key patterns for scientific traces:

```javascript
function buildTraces(defs, dataObj, xArr, cfg, isLog) {
  const { lineWidth, markerSize, markerEvery, pubMode } = cfg;
  // Sparse markers: show ~1 per markerEvery points, but legend marker always shows
  const nShow = Math.max(2, Math.ceil(xArr.length / markerEvery));

  return defs.map(def => ({
    x: xArr, y: dataObj[def.key],
    name: def.name, type: 'scatter', mode: 'lines+markers',
    line: { color: def.color, width: lineWidth, dash: def.dash ?? 'solid' },
    marker: {
      symbol: def.sym, size: markerSize,
      color: def.fillColor,
      line: { width: def.isOpen ? Math.max(1.5, lineWidth * 0.7) : 0, color: def.color },
      maxdisplayed: nShow,  // limits markers in chart; legend marker unaffected
    },
    hovertemplate: isLog
      ? `x = %{x:.2f}<br><b>${def.name}</b> = %{y:.3e}<extra></extra>`
      : `x = %{x:.2f}<br><b>${def.name}</b> = %{y:.2f}%<extra></extra>`,
  }));
}
```

**Important:**
- Use `maxdisplayed` (not `nth` filtering) — it keeps the legend marker visible while reducing chart clutter.
- Open markers need `marker.line.width > 0` and `marker.color = background color`.
- Use `hovertemplate` with `<extra></extra>` to suppress the trace name box.

### 3. Layout Building (Theme-Aware)

```javascript
function buildLayout(title, yTitle, isLog, cfg, yRange) {
  const { fsTitle, fsAxis, fsTick, fsLeg, frameWidth, pubMode, xMin, xMax } = cfg;
  // Theme-dependent colors
  const paperBg = pubMode ? 'white' : 'rgba(0,0,0,0)';  // transparent for dark mode
  const plotBg  = pubMode ? 'white' : 'rgba(0,0,0,0)';
  const textC   = pubMode ? '#111' : (isDark ? '#b0b4cc' : '#3a3e58');
  const gridC   = pubMode ? '#ccc' : (isDark ? '#2e3250' : '#d8daea');

  const axBase = {
    gridcolor: gridC, color: textC,
    linecolor: textC, mirror: true, linewidth: frameWidth,  // frame around plot
    tickfont: { size: fsTick },
  };

  const yaxis = {
    ...axBase,
    title: { text: yTitle, font: { size: fsAxis } },
    type: isLog ? 'log' : 'linear',
    exponentformat: 'e',  // scientific notation: 1e+6
    ...(isLog ? { dtick: 1 } : {}),  // one tick per decade on log scale
  };
  // Log scale range must be in log10 units
  if (yRange && isLog) {
    yaxis.range = yRange.map(v => v !== null ? Math.log10(v) : undefined);
  }

  return {
    paper_bgcolor: paperBg, plot_bgcolor: plotBg,
    title: { text: title, font: { size: fsTitle, color: textC } },
    xaxis: { ...axBase, title: { text: 'X Label', font: { size: fsAxis } } },
    yaxis,
    legend: {
      orientation: 'v', x: 1.02, y: 1, xanchor: 'left', yanchor: 'top',
      bgcolor: 'rgba(255,255,255,0.85)', borderwidth: 1,
      font: { size: fsLeg },
    },
    margin: { t: 46, b: 52, l: 72, r: 110 },  // right margin for legend
  };
}
```

**Important:**
- `mirror: true` + `linewidth: frameWidth` creates a box frame around the plot.
- Log axis `range` must be in `Math.log10()` units, not raw values.
- Legend outside plot: `x: 1.02, xanchor: 'left'` with enough right margin.
- Use `rgba(0,0,0,0)` for transparent background in dark mode (not `'transparent'`).

### 4. Font Weight Override (CSS Injection)

Plotly.js doesn't natively support font-weight per layout. Use CSS injection:

```javascript
function applyFontWeight(weight) {
  document.getElementById('plot_fw_style').textContent =
    `.js-plotly-plot text { font-weight: ${weight} !important; }`;
}
```

Requires a `<style id="plot_fw_style"></style>` in `<head>`.

### 5. Aspect Ratio

```javascript
function applyAspectRatio(ratio) {
  const css = ratio === 'sqrt(2)' ? String(Math.SQRT2) : ratio;
  document.getElementById('chart').style.aspectRatio = css;
}
```

Common ratios: `4/3` (default), `3/2`, `16/9`, `sqrt(2)` (for A4 paper).

### 6. Export / Save Figure

```javascript
function downloadChart(divId, filename, format = 'png') {
  const cfg = getSettings();
  const w = document.getElementById(divId).offsetWidth || 800;
  const scale = cfg.exportWidth / w;
  Plotly.downloadImage(divId, { format, scale, filename });
}
```

**Tips:**
- `scale` ensures exported image has the desired pixel width regardless of screen size.
- Support both PNG (raster, for slides) and SVG (vector, for papers).
- Default export width ~1600px gives good quality for presentations.

### 7. Per-Trace Controls (Color, Symbol, Visibility)

Allow users to override each trace individually:

```html
<div class="trace-ctrl-row">
  <input type="checkbox" id="vis_0" checked onchange="onSettingChange()">
  <input type="color" id="col_0" value="#4f8ef7" oninput="onSettingChange()">
  <select id="sym_0" onchange="onSettingChange()">
    <option>circle</option><option>square</option><option>diamond</option>
    <option>circle-open</option><option>square-open</option>
  </select>
  <span>Trace Name</span>
</div>
```

Set `visible: true` or `'legendonly'` based on checkbox state.

### 8. Publication Mode Checklist

When `pubMode` is enabled:
- White background (not transparent)
- Black text and axis colors
- Gray grid lines (#ccc)
- No theme-dependent styling
- Result: clean figure suitable for journal submission

### 9. Plotly Config

```javascript
const PLOTLY_CONFIG = { responsive: true, displayModeBar: false };
```

- `responsive: true` — chart resizes with container.
- `displayModeBar: false` — hides Plotly toolbar for clean UI. Users export via custom buttons instead.

### 10. Common Gotchas

1. **Use `Plotly.react()` not `Plotly.newPlot()`** for updates — it's more efficient and preserves zoom state.
2. **Call `Plotly.purge()` before switching datasets** to avoid ghost traces.
3. **Log axis range in log10 units** — `range: [4, 8]` means 10⁴ to 10⁸.
4. **`maxdisplayed` vs filtering data** — always prefer `maxdisplayed` for sparse markers; it keeps hover working on all points.
5. **Font weight requires CSS injection** — Plotly ignores `font.weight` in layout config.
6. **Aspect ratio via CSS `aspect-ratio`** on the chart container div, not via Plotly layout.
