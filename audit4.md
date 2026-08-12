**Critical Issues**

| Severity | Finding | Evidence | Risk |
|---|---|---|---|
| Critical | Client-side XSS sink in modal population via html() with dynamic values from data-* attributes | main.html, main.html, main.html, main.html, main.html, main.html | If attacker-controlled content reaches alert fields, DOM-based XSS can execute script in authenticated user sessions. |

| High | Password-reset and email-verification tokens are stored in plaintext in DB | Token columns: models.py, models.py; tokens generated/used directly: auth.py, auth.py, auth.py | DB read exposure immediately enables account takeover for active tokens. |
| High | ~~Internal exception details can leak to clients on some endpoints~~ | | Information disclosure aids recon and exploit chaining. |
| Medium | N+1 API/data-loading pattern in profile favorites UI | Per-group fetch loop: profile.html, profile.html, profile.html | For G groups, page triggers about 1+G requests and DB queries, creating avoidable latency/load. |
| Medium | Synchronous blocking work in request path (SMTP, heavy plotting/astronomy computations) | SMTP in auth flow: auth.py, auth.py, auth.py; CPU-heavy plotting in request: observing_tool.py, features.py, plotting_service.py | Slower tail latency, worker starvation, and poor scalability under concurrent load. |
| Medium | Rate-limit backend defaults to in-memory store | settings.py | In multi-process/multi-instance deployment, limits are inconsistent and easier to bypass. |
| Medium | Dependency vulnerability scan could not complete due TLS interception/cert-chain failure | Audit run failed resolving PyPI advisories (SSL cert verify failure). Manifests use wildcard pins: requirements.txt, requirements-dev.txt | Known CVEs cannot be confirmed from this environment; wildcard ranges weaken reproducibility and patch governance. |

**Flask Best Practices Assessment**

- Good:
1. Factory pattern is correctly used in __init__.py.
2. Blueprints are modular and registered centrally in app.py.
3. CSRF and login protections are present for state-changing routes, with explicit CSRF test coverage in test_security_regressions.py.

- Needs improvement:
1. Error handling is inconsistent; some routes return raw exception strings.
2. No unified JSON/HTML error boundary strategy for 400/401/403/404/429/500.
3. Security headers are not globally enforced.

**Refactoring Suggestions**

[x] Harden authentication/reset flows (token hashing + signed tokens + stronger password hashing)
    ```
        import hashlib
        import secrets
        from itsdangerous import URLSafeTimedSerializer, BadSignature, SignatureExpired
        from werkzeug.security import generate_password_hash

        def hash_token(token: str) -> str:
            return hashlib.sha256(token.encode('utf-8')).hexdigest()

        # generate token for email link
        raw_token = secrets.token_urlsafe(32)
        user.reset_token = hash_token(raw_token)

        # verify token from URL
        incoming = request.view_args['token']
        user = User.query.filter_by(reset_token=hash_token(incoming)).first()

        # use stronger KDF explicitly
        user.password = generate_password_hash(new_password, method='scrypt')
    ```

5. Remove N+1 profile loading by returning grouped favorites in one response

    @favorites_bp.route('/favorites-grouped', methods=['GET'])
    @login_required
    def favorites_grouped():
        rows = (
            db.session.query(Favorite.id, Favorite.locus_id, Favorite.group_id)
            .filter(Favorite.user_id == current_user.id)
            .all()
        )
        grouped = {}
        for fid, locus_id, group_id in rows:
            key = 'ungrouped' if group_id is None else str(group_id)
            grouped.setdefault(key, []).append({'id': fid, 'locusId': locus_id})
        return jsonify({'favoritesByGroup': grouped})

Then load once in profile.html and render tabs from cached payload.

6. Offload blocking SMTP/plot jobs

    # enqueue in request thread
    email_queue.enqueue('poi_broker.jobs.send_email_job', to_email, subject, body)

    # route returns immediately
    return jsonify({'status': 'accepted'}), 202

For plot-heavy endpoints, precompute/cache by locus_id and selected feature set, or move to background job + polling.


8. Dependency governance

- Replace wildcard major pins with bounded ranges or lock file.
- Run pip-audit in CI in a trusted cert environment.
- Fail CI on known high/critical advisories.

    # Example CI step
    python -m pip install -r requirements-dev.txt
    python -m pip install pip-audit
    pip-audit -r requirements.txt
    pip-audit -r requirements-dev.txt

If you want, I can implement the top 3 remediations directly now: XSS sink fixes in main.html, global security headers in __init__.py, and uniform 500 error responses in the affected routes.