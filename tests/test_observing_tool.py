from astropy import units as u
from astropy.coordinates import EarthLocation


def test_query_observing_plot_requires_all_parameters(client):
    response = client.get('/query_observing_plot')

    assert response.status_code == 400
    assert 'Missing required query parameters' in response.get_data(as_text=True)


def test_query_observing_plot_rejects_invalid_parameter_format(client):
    response = client.get(
        '/query_observing_plot',
        query_string={
            'obs_loc': 'Any Observatory',
            'obs_date': '2024-99-99',
            'obs_tz': 'option_utc',
            'ra': 'not-a-number',
            'dec': 'not-a-number',
        },
    )

    assert response.status_code == 400
    assert 'Invalid parameter format' in response.get_data(as_text=True)


def test_query_observing_plot_returns_400_for_unknown_observatory(client, monkeypatch):
    import poi_broker.observing_tool as observing_tool

    def fake_of_site(_site_name):
        raise Exception('site unknown')

    monkeypatch.setattr(observing_tool.EarthLocation, 'of_site', fake_of_site)

    response = client.get(
        '/query_observing_plot',
        query_string={
            'obs_loc': 'Unknown Observatory',
            'obs_date': '2025-01-01',
            'obs_tz': 'option_utc',
            'ra': '101.28715533',
            'dec': '16.71611586',
        },
    )

    assert response.status_code == 400
    assert 'Unknown observatory location' in response.get_data(as_text=True)


def test_query_observing_plot_generates_image_and_moon_panel(client, monkeypatch):
    import poi_broker.observing_tool as observing_tool

    location = EarthLocation(lat=19.8261 * u.deg, lon=-155.4700 * u.deg, height=4145 * u.m)

    monkeypatch.setattr(observing_tool.EarthLocation, 'of_site', lambda site_name: location)
    monkeypatch.setattr(observing_tool.TimezoneFinder, 'timezone_at', lambda self, lng, lat: 'UTC')
    from datetime import timezone
    monkeypatch.setattr(observing_tool, 'ZoneInfo', lambda tz_name: timezone.utc)

    response = client.get(
        '/query_observing_plot',
        query_string={
            'obs_loc': 'Keck Observatory',
            'obs_date': '2025-01-01',
            'obs_tz': 'option_utc',
            'ra': '101.28715533',
            'dec': '16.71611586',
        },
    )

    text = response.get_data(as_text=True)
    assert response.status_code == 200
    assert '<img src="data:image/png;base64,' in text
    assert '<div class="col-md-7">' in text
    assert '<div class="col-md-5">' in text
