**Findings (Prioritized)**

### Critical

- [x] Stored XSS risk in watchlist rendering on profile page  

   - Evidence: profile.html, profile.html, profile.html, visual_query.py, visual_query.py  
   - Why this matters: user-controlled fields like watchlist name and SQL snippet are inserted into innerHTML without escaping. A crafted value can execute script in the user’s session.  
   - Refactor suggestion:
     1. Escape all interpolated content before rendering (same pattern as existing escapeHtml in this file).
     2. Prefer textContent and element creation over innerHTML for user data.
     3. Keep innerHTML only for static markup wrappers.

- [x] XSS surface from API data inserted as raw HTML in modal panels  
   - Evidence: main.html, main.html, main.html, main.html  
   - Why this matters: JSON fields from backend data sources are directly templated into HTML strings; if any field contains markup/script payloads, client-side execution is possible.  
   - Refactor suggestion:
     1. Escape all dynamic cell values before concatenation.
     2. Build table rows with DOM APIs and set textContent.
     3. If rich HTML must be supported, sanitize with a strict allowlist sanitizer.

- [x] The only remaining caveat is the crossmatches footer path, which still builds a small HTML string with values like ztf_object_id and locus_id and inserts it via .html(...).

### Medium

- [x] N+1 query pattern on main list when loading classification  
   - Evidence: app.py, main.html, relationship definition at models.py  
   - Why this matters: page query loads Ztf rows, then template access to alert.classification can trigger one query per row.  
   - Refactor suggestion:
     1. Use eager loading with joinedload/selectinload for Ztf.classification.
     2. Or select Classification.prob_class in the base query and render from flat row data.
     3. Validate with SQL logging to confirm query count drops.

- [x] Incomplete transaction error handling around commits

   - Evidence: auth.py, auth.py, auth.py, auth.py, auth.py, filter_bookmarks.py, filter_bookmarks.py  
   - Why this matters: on IntegrityError/DB errors, routes can 500 without structured response and without local rollback handling in some paths.  
   - Refactor suggestion:
     1. Wrap write operations in try/except around commit.
     2. Always rollback on exception.
     3. Return controlled error payloads/messages (avoid exposing raw exception text).

- [x] Error details leaked to clients  
   - Evidence: app.py, app.py  
   - Why this matters: raw exception messages can disclose internals (schema, file paths, query context).  
   - Refactor suggestion:
     1. Log full exception server-side.
     2. Return generic user-safe error messages and consistent JSON error schema for APIs.

- [x] Frontend state fragility in feature query handling  
   - Evidence: main.html, main.html  
   - Why this matters: JSON.parse is called in complete callback regardless of HTTP status/content type; failures can throw and break modal state.  
   - Refactor suggestion:
     1. Use success/error callbacks (or fetch with response.ok checks).
     2. Guard JSON parsing and display a fallback UI message.
     3. Keep previous successful state if refresh fails.

### Low

- [x] 10. Repetitive DOM queries in hot UI paths  
   - Evidence: main.html, main.html, profile.html, profile.html  
   - Why this matters: repeated selector lookups and full-list scans are small individually but add overhead in frequent interactions.  
   - Refactor suggestion:
     1. Cache frequently used selectors (CSRF token, checkbox collections, key container nodes).
     2. Prefer incremental state updates over full recount/requery when possible.

> partially resolved

main.html still has some repeated selector work in a few hot spots, including repeated lookups for the modal labels, CSRF meta tag, and the multiselect widget, plus a full $('.sortable') scan on each sort interaction.
So this looks like a "partially stale" finding: the worst of it has been reduced, but there are still a few minor DOM-query inefficiencies left.
---

**Security Spot Check Summary**

- SQL injection: no obvious direct SQL injection found in reviewed paths; ORM usage is generally safe and raw SQL in classification uses bound parameterization (classification.py, classification.py).  
- CSRF: generally good coverage for forms and AJAX headers (base.html, login.html, profile.html); logout GET remains a weak point.  
- XSS: main concern area, especially profile watchlist and modal HTML composition paths.

---

**Open Questions / Assumptions**

2. Are DB indexes guaranteed by migration/bootstrap scripts in deployed environments (not just comments)? Performance assumptions depend on that.  
3. Do you want API errors standardized to JSON everywhere, including legacy endpoints returning HTML/text?

---

**Testing Gaps**

- Existing tests cover route smoke/security basics, but I did not see tests specifically for output-escaping/XSS in watchlist/profile rendering or modal table rendering from API payloads.  
- I also did not see query-count regression tests for N+1 hotspots (main list classification and favorite group counts).