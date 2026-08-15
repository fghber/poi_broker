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
| `script-src`             | Allows your domain, specific CDNs (jQuery, Bokeh, jsDelivr), and `unsafe-inline`/`eval`. | Enables complex libraries (e.g., Bokeh) while blocking unknown external scripts.           |
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

jQuery is **3.7.1** (CVE-2020-11022 / CVE-2020-11023 in `.html()` / `.append()` were fixed in 3.5.0). jQuery 4 is blocked by Bootstrap 4. Frontend versions and the Bootstrap 4 / Popper v1+v2 split are in `docs/spec.md` §1.5.