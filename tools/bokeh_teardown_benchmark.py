#!/usr/bin/env python3
"""
Headless Bokeh teardown benchmark using Playwright.

Usage:
    python tools/bokeh_teardown_benchmark.py --url http://127.0.0.1:5000 --iterations 50

Notes:
- Requires `playwright` to be installed and browsers downloaded: `python -m pip install playwright` and
    `python -m playwright install chromium`
- The script enables `window.BOKEH_TEARDOWN_DIAGNOSTICS = true` and clears `window._bokehTeardownStats`
    on the page. The client now forces a single aggressive teardown mode and performs a page-global
    `globalWipe()` on modal close; the CLI no longer controls teardown mode.
"""
import json
import time
import sys
import argparse
from urllib.parse import urlparse
import socket

try:
    from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeout
except Exception as e:
    print("Playwright is not installed. Install with: python -m pip install playwright", file=sys.stderr)
    raise


def wait_for_server(url, timeout=30, interval=0.5):
    parsed = urlparse(url)
    host = parsed.hostname or '127.0.0.1'
    port = parsed.port or (443 if parsed.scheme == 'https' else 80)
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with socket.create_connection((host, port), timeout=2):
                return True
        except Exception:
            time.sleep(interval)
    return False


def run(url="http://127.0.0.1:5000", iterations=50, headless=True, wait_timeout=30):
    results = []
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=headless)
        context = browser.new_context()
        page = context.new_page()

        # Collect browser console lines (helpful for debugging)
        console_lines = []
        def on_console(msg):
            try:
                console_lines.append(msg.text())
            except Exception:
                pass
        page.on("console", on_console)

        print(f"Waiting for server {url} (timeout {wait_timeout}s)")
        if not wait_for_server(url, timeout=wait_timeout):
            print(f"Timed out waiting for server at {url}", file=sys.stderr)
            browser.close()
            return 3

        print(f"Navigating to {url}")
        try:
            page.goto(url, wait_until='networkidle')
        except Exception as e:
            print(f"Page.goto failed: {e}", file=sys.stderr)
            browser.close()
            return 2

        # Ensure modal opener exists
        try:
            page.wait_for_selector('a[data-target="#objectIdModal"]', timeout=5000)
        except PlaywrightTimeout:
            print('No opener link found on the page; aborting.', file=sys.stderr)
            browser.close()
            return 2

        # Turn on diagnostics and set teardown mode on the page
        page.evaluate("() => { window.BOKEH_TEARDOWN_DIAGNOSTICS = true; window._bokehTeardownStats = []; }")

        opener_selector = 'a[data-target="#objectIdModal"]'
        close_selector = '#objectIdModal button.close'
        tab_selectors = ['#overview-tab', '#thumbnails-tab', '#features-tab', '#observing-tab',  '#classification-tab']  # Add more tabs as needed

        for i in range(iterations):
            print(f"Iteration {i+1}/{iterations}")
            # Record index before
            pre = page.evaluate("() => (window.Bokeh ? Object.keys(Bokeh.index||{}).length : 0)")

            # Click opener and wait for modal + nav links to be visible
            try:
                time.sleep(0.5)  # small delay to ensure page is ready for next interaction
                page.click(opener_selector)
                page.wait_for_selector('#objectIdModal', state='visible', timeout=5000)
                page.wait_for_selector('#objectIdModal .nav-link', timeout=5000)
            except Exception as e:
                print('Failed to open modal or find modal nav links:', e, file=sys.stderr)
                break

            # Click through tabs and wait for each to become active
            for tab_selector in tab_selectors:
                try:
                    time.sleep(0.5)
                    page.click(tab_selector)
                except Exception as e:
                    print(f'Failed to click tab {tab_selector}:', e, file=sys.stderr)
                    break

                # Wait for the tab to have the 'active' class
                try:
                    page.wait_for_selector(f"{tab_selector}.active", timeout=10000)
                except PlaywrightTimeout:
                    print(f'Timeout waiting for tab {tab_selector} to become active', file=sys.stderr)

            # Close modal
            try:
                page.click(close_selector)
            except Exception as e:
                print('Failed to click close button:', e, file=sys.stderr)

            # Wait for teardown to run (either registry empty or container empty)
            try:
                page.wait_for_function("() => (window.Bokeh && Object.keys(Bokeh.index||{}).length === 0) || (document.getElementById('locus-plot') && document.getElementById('locus-plot').children.length === 0)", timeout=10000)
            except PlaywrightTimeout:
                # Not fatal; collect what we have
                pass

            post = page.evaluate("() => (window.Bokeh ? Object.keys(Bokeh.index||{}).length : 0)")
            stats = page.evaluate("() => (window._bokehTeardownStats || []).slice(-1)[0] || null")
            results.append({ 
                'iteration': i+1, 
                'preIndex': pre, 
                'postIndex': post, 
                'stat': stats 
            })

            # Small sleep to give the browser some time
            time.sleep(0.25)

        # Collect any console messages and the full stats array
        all_stats = page.evaluate("() => (window._bokehTeardownStats || [])")
        browser.close()

    output = {
        'url': url,
        'iterations': iterations,
        'per_iteration': results,
        'all_stats': all_stats,
        'console': console_lines
    }

    print(json.dumps(output, indent=2, default=str))
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--url', default='http://127.0.0.1:5000')
    parser.add_argument('--iterations', type=int, default=50)
    parser.add_argument('--headless', action='store_true', default=True)
    # The script no longer accepts a teardown mode; the client forces aggressive teardown.
    parser.add_argument('--wait', type=int, default=30, help='Seconds to wait for target server to accept connections')
    args = parser.parse_args()
    rc = run(url=args.url, iterations=args.iterations, headless=args.headless, wait_timeout=args.wait)
    sys.exit(rc)
