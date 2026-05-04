from poi_broker import db
from poi_broker.classification import classification_plot
from poi_broker.models import Classification


def test_classification_plot_requires_alert_id(app):
    with app.test_request_context('/query_classification'):
        response = classification_plot()

    assert isinstance(response, tuple)
    assert "Missing alertId" in response[0]
    assert response[1] == ""


def test_classification_plot_returns_no_data_for_unknown_alert_id(app):
    alert_id = "missing-alert"

    with app.test_request_context(f'/query_classification?alertId={alert_id}'):
        response = classification_plot()

    assert isinstance(response, tuple)
    assert f"No classification data found for alert_id={alert_id}" in response[0]
    assert response[1] == ""


def test_classification_plot_renders_bokeh_components_for_existing_record(app):
    alert_id = "alert-123"
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

    with app.test_request_context(f'/query_classification?alertId={alert_id}'):
        response = classification_plot()

    assert isinstance(response, str)
    assert response.startswith('<div')
    assert '<script' in response
    assert 'Classified as SN(0.87)' in response
