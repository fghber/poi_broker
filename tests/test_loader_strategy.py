"""O1/O2: classification lazy='raise' and plot column projection."""

import pytest
from sqlalchemy.dialects import sqlite as sqlite_dialect
from sqlalchemy.exc import InvalidRequestError

from poi_broker import db
from poi_broker.models import Classification, Ztf
from poi_broker.routes.lightcurve import (
    build_lightcurve_csv_query,
    build_lightcurve_plot_query,
)
from poi_broker.services.feature_service import (
    build_feature_plot_query,
    query_feature_plot_data,
)

SELECTED_FEATURE = 'feature_amplitude_magn_r'


def _seed_locus_rows() -> str:
    locus_id = 'L-plot-1'
    db.session.add_all(
        [
            Ztf(
                alert_id='ztf_candidate:2000000000000000001',
                date_alert_mjd=60010.0,
                locus_id=locus_id,
                ztf_object_id='Z-plot-1',
                ant_passband='g',
                ant_mag_corrected=18.5,
                feature_amplitude_magn_r=0.12,
            ),
            Ztf(
                alert_id='ztf_candidate:2000000000000000002',
                date_alert_mjd=60011.0,
                locus_id=locus_id,
                ztf_object_id='Z-plot-1',
                ant_passband='R',
                ant_mag_corrected=18.7,
                feature_amplitude_magn_r=0.15,
            ),
            Classification(
                alert_id='ztf_candidate:2000000000000000001',
                prob_class='sn',
                p_sn=0.9,
            ),
        ]
    )
    db.session.commit()
    return locus_id


def test_ztf_classification_lazy_raise(app):
    with app.app_context():
        _seed_locus_rows()
        db.session.expunge_all()
        ztf = db.session.query(Ztf).first()
        assert ztf is not None
        with pytest.raises(InvalidRequestError):
            _ = ztf.classification


def test_plot_routes_return_data_for_seeded_locus(client, app):
    with app.app_context():
        locus_id = _seed_locus_rows()

    plot_response = client.get('/query_lightcurve_data', query_string={'locusId': locus_id})
    assert plot_response.status_code == 200
    assert plot_response.is_json
    plot_payload = plot_response.get_json()
    assert 'No lightcurve data' not in plot_payload.get('div', '')

    csv_response = client.get('/locus_plot_csv', query_string={'locusId': locus_id})
    assert csv_response.status_code == 200
    assert 'text/csv' in csv_response.content_type
    csv_text = csv_response.get_data(as_text=True)
    assert csv_text.startswith('locus_id,date_alert_mjd,ant_mag_corrected\n')
    assert f'{locus_id},60010.0,18.5' in csv_text
    assert f'{locus_id},60011.0,18.7' in csv_text

    feature_response = client.get(
        '/query_featureplot_data',
        query_string={'locusId': locus_id, 'features': SELECTED_FEATURE},
    )
    assert feature_response.status_code == 200
    assert feature_response.is_json
    feature_payload = feature_response.get_json()
    feature_div = feature_payload.get('div', '')
    feature_script = feature_payload.get('script', '')
    assert 'no feature data' not in feature_div.lower()
    assert 'No features selected' not in feature_div
    assert '<script' not in feature_script
    assert 'Bokeh' in feature_script


def test_lightcurve_query_select_list_is_projected(app):
    with app.app_context():
        plot_query = build_lightcurve_plot_query('L-plot-1')
        plot_sql = str(plot_query.statement.compile(dialect=sqlite_dialect.dialect()))
        plot_names = [desc['name'] for desc in plot_query.column_descriptions]
        assert plot_names == ['date_alert_mjd', 'ant_mag_corrected', 'ant_passband']
        assert 'feature_amplitude_magn_r' not in plot_sql
        assert 'prob_class' not in plot_sql

        csv_query = build_lightcurve_csv_query('L-plot-1')
        csv_sql = str(csv_query.statement.compile(dialect=sqlite_dialect.dialect()))
        csv_names = [desc['name'] for desc in csv_query.column_descriptions]
        assert csv_names == ['locus_id', 'date_alert_mjd', 'ant_mag_corrected']
        assert 'ant_passband' not in csv_sql
        assert 'prob_class' not in csv_sql


def test_feature_plot_query_projects_columns_not_full_entity(app):
    with app.app_context():
        locus_id = _seed_locus_rows()
        db.session.expunge_all()

        selected = [SELECTED_FEATURE]
        query, used = build_feature_plot_query(locus_id, selected)
        compiled = str(query.statement.compile(dialect=sqlite_dialect.dialect()))
        names = [desc['name'] for desc in query.column_descriptions]
        assert used == selected
        assert names == [
            'date_alert_mjd',
            'ant_mag_corrected',
            SELECTED_FEATURE,
        ]
        assert 'ant_passband' not in compiled
        assert 'prob_class' not in compiled

        data, used_from_query = query_feature_plot_data(locus_id, selected)
        assert used_from_query == selected
        assert len(data) == 2
        for row in data:
            assert not isinstance(row, Ztf)
        identities = list(db.session.identity_map.values())
        assert not any(isinstance(obj, Ztf) for obj in identities)
