# Bokeh Teardown & Benchmark — README

This document describes the client-side Bokeh teardown logic used by the POI broker, the runtime toggles exposed for testing, and how to run the headless teardown benchmark.

Purpose
- Explain the teardown modes implemented in `poi_broker/templates/main.html` and how the automatic escalation works.
- Show the runtime toggles you can change from the browser console.
- Describe how to run the Playwright benchmark in `tools/bokeh_teardown_benchmark.py`.

Key files
- poi_broker/templates/main.html — main teardown and escalation logic (client JS).
- tools/bokeh_teardown_benchmark.py — headless Playwright benchmark harness.

Runtime toggles (browser console)
- `window.BOKEH_TEARDOWN_MODE` (string): 'safe' | 'aggressive' | 'global'
  - 'safe': minimal cleanup (prefer remove_root for shared docs)
  - 'aggressive': per-document clears when a document is local to the modal + orphan sweep + delete index entries
  - 'global': last-resort full wipe of `Bokeh.documents` and empty `Bokeh.index`
  - Default: 'aggressive'

- `window.BOKEH_TEARDOWN_DIAGNOSTICS` (boolean): When true, runtime stats are collected to `window._bokehTeardownStats`.

- The client now uses a single aggressive teardown mode and performs a page-wide `globalWipe()` on modal close (250ms after `resetModalContent()`).
  - The previous toggles `window.BOKEH_TEARDOWN_AUTO_GLOBAL` and `window.BOKEH_TEARDOWN_AUTO_GLOBAL_FALLBACK` have been removed.

Data collected (when diagnostics enabled)
- `window._bokehTeardownStats` is an ordered array of stat objects. Typical fields:
  - `container`: selector cleaned (e.g. '#locus-plot')
  - `mode`: mode used for that cleanup
  - `preIndexCount`: number of keys in `Bokeh.index` before cleanup
  - `postIndexCount`: number of keys after cleanup
  - `docsCleared`, `modelsRemoved`, `viewsRemoved`, `indexDeleted`, `orphansRemoved` — counters from the per-container sweep
  - `globalDocsWiped`, `globalIndexCleared` — counters emitted when `global` wipe ran
  - `globalEscalated` (for escalation entries): true if escalation ran
  - `reason` (for escalation entries): e.g. 'external-views', 'no-entries', 'no-escalation-needed'
  - `timestamp` (ms since epoch)

How the conditional escalation works
1. `resetModalContent()` calls `teardownBokehInContainer('#locus-plot')` and `teardownBokehInContainer('#feature-plot')` to attempt local cleanup.
2. On modal `hidden.bs.modal`, after a 250ms timeout, `performGlobalWipeIfNeeded('#objectIdModal')` runs.
3. `performGlobalWipeIfNeeded` checks:
   - If there are no `Bokeh.index` entries, skip escalation.
   - If any `Bokeh.index` view is attached to the document outside the modal container, skip escalation to avoid clearing globally-used plots.
   - If the most recent teardown diagnostic (`window._bokehTeardownStats.slice(-1)[0]`) shows `postIndexCount > 0` (i.e. the aggressive cleanup left entries), escalation proceeds.
  - The new behavior always runs `globalWipe()` on modal close; if you rely on persistent non-modal Bokeh plots, you must change this behavior.
4. Escalation runs the same 'global' branch used by `teardownBokehInContainer` (clear `Bokeh.documents`, iterate and clear `Bokeh.index`, call `view.remove()` where available), then appends a diagnostic stat.

Safety notes
- `global` is destructive: it will clear any Bokeh documents/views on the page, including ones outside the modal. The helper refuses to escalate if it detects other on-page Bokeh views.
 - The new default is aggressive teardown with unconditional global wipe on modal close. If your page contains other persistent Bokeh plots, consider guarding `globalWipe()` or reverting to conditional logic.

Quick console recipe (development)

1. Enable diagnostics and inspect stats:

```js
window.BOKEH_TEARDOWN_DIAGNOSTICS = true;
// optionally set mode
window.BOKEH_TEARDOWN_MODE = 'aggressive'; // or 'global' for manual testing
// inspect collected stats (pretty-print)
console.log(window._bokehTeardownStats);
copy(JSON.stringify(window._bokehTeardownStats, null, 2));
```

2. Force a manual global wipe from console (use only for debugging):

```js
window.BOKEH_TEARDOWN_MODE = 'global';
teardownBokehInContainer('#locus-plot');
performGlobalWipeIfNeeded('#objectIdModal');
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
python tools/bokeh_teardown_benchmark.py --url http://127.0.0.1:5000 --iterations 50 --mode global --wait 30
```

Options (benchmark script)
- `--url`: Base URL of the running app.
- `--iterations`: Number of open/close iterations.
- `--mode`: `safe|aggressive|global` — sets `window.BOKEH_TEARDOWN_MODE` on the page before starting iterations.
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
