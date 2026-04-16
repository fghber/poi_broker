import json
import re

from poi_broker.constants.features import FEATURE_COLUMNS
from werkzeug.security import generate_password_hash


def _extract_csrf_token(page_html):
    """Extract CSRF token from rendered HTML form."""
    match = re.search(r'name="csrf_token"\s+value="([^"]+)"', page_html)
    return match.group(1) if match else None


def test_settings_requires_login(client):
    r = client.get('/settings', follow_redirects=False)
    assert r.status_code in (302, 303)
    assert '/login' in r.headers.get('Location', '')


def test_settings_page_loads(auth_client):
    r = auth_client.get('/settings')
    assert r.status_code == 200
    assert b'Feature Plot Defaults' in r.data
    assert b'default_feature_plot_columns' in r.data
    assert b'csrf_token' in r.data


def test_save_and_apply_default_feature_columns(auth_client):
    # Retrieve settings page to get CSRF token
    r = auth_client.get('/settings')
    assert r.status_code == 200
    page_html = r.get_data(as_text=True)
    csrf_token = _extract_csrf_token(page_html)
    assert csrf_token is not None, 'CSRF token not found in form'

    # Save a custom default feature selection
    feature_names = [
        'feature_amplitude_magn_r',
        'feature_beyond_1_std_magn_r',
    ]
    r = auth_client.post(
        '/settings',
        data={
            'default_feature_plot_columns': feature_names,
            'csrf_token': csrf_token,
        },
        follow_redirects=True,
    )
    assert r.status_code == 200

    # The settings page should reflect the saved values.
    r = auth_client.get('/settings')
    assert r.status_code == 200
    page_text = r.get_data(as_text=True)
    assert 'feature_amplitude_magn_r" selected' in page_text
    assert 'feature_beyond_1_std_magn_r" selected' in page_text

    # The main page should include the saved list in the rendered JavaScript context.
    r = auth_client.get('/')
    assert r.status_code == 200
    page_text = r.get_data(as_text=True)
    assert json.dumps(feature_names) in page_text


def test_save_truncates_to_ten_feature_columns(auth_client):
    r = auth_client.get('/settings')
    assert r.status_code == 200
    csrf_token = _extract_csrf_token(r.get_data(as_text=True))
    assert csrf_token is not None

    keys = list(FEATURE_COLUMNS.keys())[:11]
    assert len(keys) == 11

    r = auth_client.post(
        '/settings',
        data={
            'default_feature_plot_columns': keys,
            'csrf_token': csrf_token,
        },
        follow_redirects=True,
    )
    assert r.status_code == 200

    r = auth_client.get('/settings')
    page_text = r.get_data(as_text=True)
    for k in keys[:10]:
        assert f'{k}" selected' in page_text
    assert f'{keys[10]}" selected' not in page_text
