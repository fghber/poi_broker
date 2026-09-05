from poi_broker import db
from poi_broker.models import Classification


def test_classification_plot_requires_alert_id(client):
    response = client.get('/query_classification')

    assert response.status_code == 400
    assert response.is_json
    assert response.get_json()['error'] == 'Missing alertId'


def test_classification_plot_rejects_oversized_alert_id(client):
    response = client.get('/query_classification', query_string={'alertId': 'x' * 129})

    assert response.status_code == 400
    assert response.is_json
    assert response.get_json()['error'] == 'alertId is too long'


def test_classification_plot_does_not_reflect_alert_id(client):
    xss = '"><script>alert(1)</script>'
    response = client.get('/query_classification', query_string={'alertId': xss})

    assert response.status_code == 200
    assert response.is_json
    payload = response.get_json()
    assert payload['script'] == ''
    assert 'No classification data found' in payload['div']
    assert 'alert-warning' in payload['div']
    body = response.get_data(as_text=True)
    assert '<script>' not in body
    assert xss not in body
    assert 'alert_id=' not in body


def test_classification_plot_returns_no_data_for_unknown_alert_id(client):
    response = client.get('/query_classification', query_string={'alertId': 'missing-alert'})

    assert response.status_code == 200
    assert response.is_json
    payload = response.get_json()
    assert payload['script'] == ''
    assert payload['div'].startswith('<div')
    assert 'alert-warning' in payload['div']
    assert 'No classification data found' in payload['div']
    assert 'missing-alert' not in payload['div']


def test_classification_plot_renders_bokeh_components_for_existing_record(client, app):
    alert_id = 'alert-123'
    classification_row = Classification(
        alert_id=alert_id,
        p_cvnova=None,
        p_e=0.1,
        p_lpv=None,
        p_puls=None,
        p_periodic_other=None,
        p_quas=0.2,
        p_sn=0.87,
        p_yso=0.4,
    )

    with app.app_context():
        db.session.add(classification_row)
        db.session.commit()

    response = client.get('/query_classification', query_string={'alertId': alert_id})

    assert response.status_code == 200
    assert response.is_json
    payload = response.get_json()
    assert payload['div'].startswith('<div')
    assert '<script' not in payload['script']
    assert 'Bokeh' in payload['script']
    assert 'Classified as SN(0.87)' in payload['script'] or 'Classified as SN(0.87)' in payload['div']
