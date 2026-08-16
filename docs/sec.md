# Security Headers for Flask Applications

This Flask implementation uses the `@app.after_request` decorator to inject security-focused HTTP headers into every response, ensuring a **"secure by default"** posture. This approach hardens the web application against common vulnerabilities such as **Cross-Site Scripting (XSS)**, **Clickjacking**, and **Protocol Downgrade attacks**.

---

## Security Headers Breakdown

### 1. HTTP Strict Transport Security (HSTS)

**Header:**  
`'Strict-Transport-Security': 'max-age=31536000; includeSubDomains; preload'`

- **Purpose:** Forces browsers to use HTTPS for all connections to your site for one year.
- **Benefit:** Prevents **Man-in-the-Middle (MITM)** attacks by blocking unencrypted HTTP connections. The `preload` directive allows inclusion in browser preload lists for built-in HTTPS enforcement.

---

### 2. MIME Sniffing Prevention

**Header:**  
`'X-Content-Type-Options': 'nosniff'`

- **Purpose:** Instructs browsers to strictly respect the `Content-Type` header.
- **Benefit:** Mitigates **MIME-sniffing attacks**, where browsers might execute files (e.g., `.txt` or images) as scripts.

---

### 3. Clickjacking Protection

**Header:**  
`'X-Frame-Options': 'DENY'`

- **Purpose:** Prevents the site from being embedded in `<iframe>`, `<frame>`, or `<object>` elements.
- **Benefit:** Eliminates **Clickjacking** risks, where attackers overlay your site to trick users into unintended actions.

---

### 4. Referrer Policy

**Header:**  
`'Referrer-Policy': 'strict-origin-when-cross-origin'`

- **Purpose:** Controls the information sent in the `Referer` header during navigation.
- **Benefit:** Protects user privacy by sending:
  - Full URL for same-origin requests.
  - Only the domain for cross-origin HTTPS requests.
  - No referrer for HTTPS-to-HTTP transitions.

---

### 5. Permissions Policy

**Header:**  
`'Permissions-Policy': 'geolocation=(), microphone=(), camera=()'`

- **Purpose:** Disables browser access to specific hardware/APIs for your site.
- **Benefit:** Reduces the **attack surface** by blocking access to camera, microphone, or location, even if malicious scripts are injected.

---

## Content Security Policy (CSP)

The CSP header acts as a gatekeeper, specifying trusted sources for content types. Below is the configuration and its purpose:


| Directive                | Configuration Summary                                                                    | Benefit                                                                                    |
| ------------------------ | ---------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| `default-src 'self'`     | Only allows content from your domain by default.                                         | Blocks unauthorized data loading.                                                          |
| `script-src`             | Allows your domain, jQuery/jsDelivr CDNs, and `unsafe-inline`/`eval`. Bokeh JS is same-origin. | Enables complex libraries (e.g., Bokeh) while blocking unknown external scripts.           |
| `style-src`              | Allows your domain, Google Fonts, and jQuery styles.                                     | Ensures UI loads from trusted design repositories.                                         |
| `connect-src`            | Allows connections to scientific/astronomical data origins (CDS, IRSA, ESAC, etc.).      | **Critical for Data Apps:** Permits real-time data fetching from trusted global databases. |
| `img-src`                | Allows images from your domain, data URIs, and any HTTPS source.                         | Flexibility for dynamic images and embedded icons.                                         |
| `frame-ancestors 'none'` | Modern alternative to `X-Frame-Options`.                                                 | Secondary defense against Clickjacking.                                                    |


---

## Summary of Benefits

- **XSS Mitigation:** CSP restricts script execution to trusted sources, making XSS attacks significantly harder.
- **Protocol Security:** HSTS ensures users remain on encrypted connections.
- **Data Integrity:** Whitelisting domains in `connect-src` ensures interactions only with trusted data providers.
- **Compliance:** Aligns with security best practices (e.g., OWASP Top 10) and audit requirements.

> **Note:** The use of `'unsafe-inline'` and `'unsafe-eval'` in `script-src` is a necessary compromise for compatibility with older libraries or CDNs. Monitor these closely to minimize XSS risks.

Plot endpoints (`/query_classification`, `/query_lightcurve_data`, `/query_featureplot_data`) return JSON `{div, script}` rather than concatenated HTML. Empty plot data is HTTP 200 with a static warning `div` and empty `script` via `bokeh_warning_payload()` — **not** 404/`{error}` (that paints a red modal error). The catalog modal inserts the Bokeh `div` and evaluates only the JSON `script` field. Classification never echoes `alertId`. Observing plots return JSON `{image, moonHtml}` (moon up), `{image, moonMessage}` (moon down), or `{message}` (not visible). The client styles `message`/`moonMessage` via `textContent` (not innerHTML) and sets `img.src` only when the value is a `data:image/png;base64,` URL. Last-selected observatory is written only by authenticated `POST /api/last-observatory` with `X-CSRFToken`, not by `GET /query_observing_plot`. Export status errors use `textContent`.

Signup (`POST /signup`) uses one generic flash and a `/login` redirect for a new account, a duplicate email, and a uniqueness race, so the response is not an account-existence oracle. Forgot-password (`POST /forgot-password`) uses one generic flash for an existing email, a missing email, and a mail/commit failure. Invalid email/password/name stay distinct validation errors; a signup mail send failure stays an operational error. Verification tokens expire after 24 hours (`email_verification_token_expires`); reset tokens expire after 1 hour. Set `PUBLIC_BASE_URL` (e.g. `https://poibroker.example.edu`) so verification/reset emails do not take their host from `Host` / `X-Forwarded-Host`. `RATELIMIT_STORAGE_URI` defaults to `memory://` (per-process); multi-worker Gunicorn must set a shared backend URI. HTML CSRF failures redirect to `request.path`, not `Referer`. Catalog modal query strings encode `locusId` / `alertId` with `encodeURIComponent`. Bokeh JS is served same-origin from the installed package (`GET /bokeh.min.js`); QueryBuilder JS/CSS are local static snapshots with no CDN fallback. Do not attach SRI to a jsDelivr minify URL for those files — jsDelivr may re-minify and the hash will not stay valid.

jQuery is **3.7.1** (CVE-2020-11022 / CVE-2020-11023 in `.html()` / `.append()` were fixed in 3.5.0). jQuery 4 is blocked by Bootstrap 4. Frontend versions and the Bootstrap 4 / Popper v1+v2 split are in `docs/spec.md` §1.5.