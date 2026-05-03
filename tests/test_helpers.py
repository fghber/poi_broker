import json

from poi_broker import db
from poi_broker.helpers import (
    extract_dates,
    extract_float_filter,
    extract_int_filter,
    extract_numbers,
    extract_mjd_filter,
    object_as_dict,
    result_to_dict,
    safe_serialize,
    serialize_fallback,
)
from poi_broker.models import Ztf


def test_extract_numbers_returns_none_when_no_numeric_text():
    assert extract_numbers('no numbers here') is None


def test_extract_numbers_preserves_single_value_with_comparator():
    assert extract_numbers('>42') == ['>42']
    assert extract_numbers('<-1.5') == ['<-1.5']


def test_extract_numbers_returns_two_values_without_comparators():
    assert extract_numbers('>1.2 and <3.4') == ['1.2', '3.4']
    assert extract_numbers('5 10') == ['5', '10']


def test_extract_dates_parses_yyyymmdd_and_iso_formats():
    assert extract_dates('20250115') == ['2025-01-15']
    assert extract_dates('2025-01-15T12:30:00') == ['2025-01-15T12:30:00']
    assert extract_dates('<2025-01-15 00:00:00') == ['<2025-01-15 00:00:00']


def test_extract_dates_falls_back_when_yyyymmdd_is_invalid_or_prefixed():
    assert extract_dates('>20250115') == ['>20250115']
    assert extract_dates('20251301') == ['20251301']
    assert extract_dates('not-a-date') == []


def test_extract_float_filter_creates_equality_filter_for_exact_value(app):
    with app.app_context():
        query = db.session.query(Ztf)
        filtered = extract_float_filter(['1.23'], Ztf.date_alert_mjd, query)
        dialect = query.session.get_bind().dialect
        compiled = str(filtered.statement.compile(dialect=dialect, compile_kwargs={'literal_binds': True}))

    assert 'date_alert_mjd' in compiled
    assert '= 1.23' in compiled or '== 1.23' in compiled


def test_extract_int_filter_handles_two_values_and_orders_inputs(app):
    with app.app_context():
        query = db.session.query(Ztf)
        filtered = extract_int_filter(['5', '1'], Ztf.date_alert_mjd, query)
        dialect = query.session.get_bind().dialect
        compiled = str(filtered.statement.compile(dialect=dialect, compile_kwargs={'literal_binds': True}))

    assert 'date_alert_mjd >= 1' in compiled
    assert 'date_alert_mjd <= 5' in compiled


def test_extract_mjd_filter_converts_iso_date_to_mjd_range(app):
    with app.app_context():
        query = db.session.query(Ztf)
        filtered = extract_mjd_filter(['2025-01-15'], Ztf.date_alert_mjd, query)
        dialect = query.session.get_bind().dialect
        compiled = str(filtered.statement.compile(dialect=dialect, compile_kwargs={'literal_binds': True}))

    assert 'date_alert_mjd' in compiled
    assert '>=' in compiled and '<=' in compiled


def test_safe_serialize_falls_back_for_non_json_objects():
    payload = {'data': b'test-bytes'}
    result = safe_serialize(payload)

    assert isinstance(result, str)
    assert json.loads(result)['data'] == 'test-bytes'


def test_serialize_fallback_handles_nested_collections():
    obj = {
        'bytes': b'foo',
        'nested': [b'bar', {'inner': b'baz'}],
    }

    serialized = serialize_fallback(obj)
    assert serialized['bytes'] == 'foo'
    assert serialized['nested'][0] == 'bar'
    assert serialized['nested'][1]['inner'] == 'baz'


def test_object_as_dict_and_result_to_dict_with_model_instance(app):
    sample = Ztf(alert_id='test-1', date_alert_mjd=59000.5)

    as_dict = object_as_dict(sample)
    assert as_dict['alert_id'] == 'test-1'
    assert as_dict['date_alert_mjd'] == 59000.5

    results = result_to_dict([sample])
    assert isinstance(results, list)
    assert results[0]['alert_id'] == 'test-1'
    assert results[0]['date_alert_mjd'] == 59000.5
