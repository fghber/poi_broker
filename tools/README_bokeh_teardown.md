# Bokeh Teardown & Benchmark — README

This document describes the client-side Bokeh teardown logic used by the POI broker, the runtime toggles exposed for testing, and how to run the headless teardown benchmark.

Purpose
- Explain the teardown modes implemented in `poi_broker/templates/main.html` and how the automatic escalation works.
- Show the runtime toggles you can change from the browser console.
- Describe how to run the Playwright benchmark in `tools/bokeh_teardown_benchmark.py`.

Key files
- poi_broker/templates/main.html — main teardown and escalation logic (client JS).
- tools/bokeh_teardown_benchmark.py — headless Playwright benchmark harness.

Runtime notes (browser console)
- `window.BOKEH_TEARDOWN_DIAGNOSTICS` (boolean): When true, runtime stats are collected to `window._bokehTeardownStats`.

- The client uses a single aggressive teardown mode and performs a page-wide `globalWipe()` on modal close (250ms after `resetModalContent()`).
  

Data collected (when diagnostics enabled)
- `window._bokehTeardownStats` is an ordered array of stat objects. Typical fields:
  - `container`: selector cleaned (e.g. '#locus-plot')
  - `mode`: mode used for that cleanup
  - `preIndexCount`: number of keys in `Bokeh.index` before cleanup
  - `postIndexCount`: number of keys after cleanup
  - `docsCleared`, `modelsRemoved`, `viewsRemoved`, `indexDeleted`, `orphansRemoved` — counters from the per-container sweep
  - For `global`-wipe entries the stat object contains `mode: 'global-wipe'` and fields such as `docsCleared`, `viewsRemoved`, and `indexCleared` emitted by `globalWipe()`.
  - `globalEscalated` (for escalation entries): true if escalation ran
  - `reason` (for escalation entries): e.g. 'external-views', 'no-entries', 'no-escalation-needed'
  - `timestamp` (ms since epoch)

How the teardown works
1. `resetModalContent()` calls per-container `teardownBokehInContainer()` for the modal's plot containers (e.g. `#locus-plot`, `#feature-plot`).
2. On modal `hidden.bs.modal`, after a 250ms timeout, `globalWipe()` runs unconditionally; it attempts to clear `Bokeh.documents` and deletes remaining keys from `Bokeh.index`, calling `view.remove()` and `document.remove_root()` where available.
3. If your page hosts persistent, non-modal Bokeh plots, guard `globalWipe()` or revert to a conditional escalation approach; `globalWipe()` is destructive and intended for modal-only pages.

Safety notes
- `global` is destructive: it will clear any Bokeh documents/views on the page, including ones outside the modal. The helper refuses to escalate if it detects other on-page Bokeh views.
 - The new default is aggressive teardown with unconditional global wipe on modal close. If your page contains other persistent Bokeh plots, consider guarding `globalWipe()` or reverting to conditional logic.

Quick console recipe (development)

1. Enable diagnostics and inspect stats:

```js
window.BOKEH_TEARDOWN_DIAGNOSTICS = true;
// inspect collected stats (pretty-print)
console.log(window._bokehTeardownStats);
copy(JSON.stringify(window._bokehTeardownStats, null, 2));
```

2. Force a manual global wipe from console (use only for debugging):

```js
teardownBokehInContainer('#locus-plot');
globalWipe();
```

Benchmark (headless Playwright)

Prerequisites (local dev machine):

```bash
python -m pip install -r requirements-dev.txt
# or if you install playwright directly:
python -m pip install playwright
python -m playwright install chromium
```

Start the app (Windows `cmd` example):

```cmd
set SECRET_KEY=devsecret
python -m flask --app wsgi:app run --no-debug --no-reload
```

Run the benchmark (example):

```bash
python tools/bokeh_teardown_benchmark.py --url http://127.0.0.1:5000 --iterations 50 --wait 30
```

Options (benchmark script)
- `--url`: Base URL of the running app.
- `--iterations`: Number of open/close iterations.
- `--wait`: Wait timeout for server readiness (seconds).

Interpreting results
- The benchmark collects per-iteration diagnostics. Look at `preIndexCount` and `postIndexCount` across iterations. If `postIndexCount` grows linearly you still have a leak.
- Prefer `postIndexCount === 0` after `global` runs when the page contains only modal plots.
- Use Chrome DevTools heap snapshots (baseline, mid-run, end-run) and search for retained Bokeh model instances or typed arrays from `ColumnDataSource`.

Notes & next steps
- Long-term: replace inline `components()` HTML (server-side) with `bokeh.embed.json_item` / `Bokeh.embed.embed_item()` to avoid repeated inline function compilation and closure retention. That change greatly simplifies client-side teardown.
- If the benchmark still shows leaks with `global`, capture heap snapshots and share them for deeper analysis.

Contact
- If you need the benchmark run or help interpreting heap snapshots, ping the maintainer on the repo.
