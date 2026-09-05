"""Tests for hybrid catalog pagination and cheap / on-demand counts."""

from astropy.time import Time

from poi_broker import db
from poi_broker.models import Classification, Ztf
from poi_broker.services.catalog_list import (
    PAGE_SIZE,
    cached_catalog_total,
    clear_catalog_total_cache,
    count_matches,
    cursor_from_row,
    fetch_keyset,
    fetch_page,
    last_page_from_total,
    parse_keyset_request,
    resolve_catalog_total,
    should_use_keyset,
)
from poi_broker.services.catalog_query import build_catalog_query, catalog_href, catalog_query_string


def _add_alert(*, alert_id: str, mjd: float, locus_id: str, object_id: str, passband: str = 'g'):
    row = Ztf(
        alert_id=alert_id,
        date_alert_mjd=mjd,
        locus_id=locus_id,
        ztf_object_id=object_id,
        ant_passband=passband,
        locus_ra=10.0,
        locus_dec=20.0,
        ant_mag_corrected=18.0,
        num_alerts=1,
        num_mag_values=1,
    )
    db.session.add(row)
    return row


def test_fetch_page_has_next_from_extra_row(app):
    with app.app_context():
        for index in range(5):
            _add_alert(
                alert_id=f'ztf_candidate:{index:019d}',
                mjd=60000.0 + index,
                locus_id=f'L{index}',
                object_id=f'Z{index}',
            )
        db.session.commit()
        query = db.session.query(Ztf).order_by(Ztf.date_alert_mjd.asc())
        page1, has_next = fetch_page(query, 1, page_size=2)
        assert len(page1) == 2
        assert has_next is True
        page3, has_next = fetch_page(query, 3, page_size=2)
        assert len(page3) == 1
        assert has_next is False


def test_fetch_keyset_no_overlap_with_tied_mjd(app):
    with app.app_context():
        _add_alert(alert_id='ztf_candidate:1000000000000000001', mjd=60000.0, locus_id='L-a', object_id='Za')
        _add_alert(alert_id='ztf_candidate:1000000000000000002', mjd=60000.0, locus_id='L-b', object_id='Zb')
        _add_alert(alert_id='ztf_candidate:1000000000000000003', mjd=60000.0, locus_id='L-c', object_id='Zc')
        db.session.commit()
        query = db.session.query(Ztf)
        first, has_next, has_prev = fetch_keyset(
            query, cursor=None, direction='first', sort_desc=True, page_size=2
        )
        assert has_prev is False
        assert has_next is True
        assert len(first) == 2
        nxt, has_next, has_prev = fetch_keyset(
            query,
            cursor=cursor_from_row(first[-1]),
            direction='before',
            sort_desc=True,
            page_size=2,
        )
        assert has_prev is True
        assert has_next is False
        first_ids = {row.alert_id for row in first}
        next_ids = {row.alert_id for row in nxt}
        assert first_ids.isdisjoint(next_ids)
        assert first_ids | next_ids == {
            'ztf_candidate:1000000000000000001',
            'ztf_candidate:1000000000000000002',
            'ztf_candidate:1000000000000000003',
        }


def test_cached_catalog_total_ttl(app, monkeypatch):
    with app.app_context():
        clear_catalog_total_cache()
        _add_alert(alert_id='ztf_candidate:2000000000000000001', mjd=60001.0, locus_id='Lc', object_id='Zc')
        db.session.commit()
        clock = {'now': 1000.0}
        monkeypatch.setattr('poi_broker.services.catalog_list.time.monotonic', lambda: clock['now'])
        assert cached_catalog_total(ttl_seconds=30) == 1
        _add_alert(alert_id='ztf_candidate:2000000000000000002', mjd=60002.0, locus_id='Ld', object_id='Zd')
        db.session.commit()
        assert cached_catalog_total(ttl_seconds=30) == 1
        clock['now'] = 1031.0
        assert cached_catalog_total(ttl_seconds=30) == 2
        clear_catalog_total_cache()


def test_resolve_catalog_total_defers_unselective_filter(app):
    with app.app_context():
        build = build_catalog_query({'ant_passband': 'g'})
        total, deferred = resolve_catalog_total(
            build,
            items=[object()],
            has_next=True,
            page=1,
            use_keyset=True,
            keyset_direction='first',
        )
        assert total is None
        assert deferred is True


def test_resolve_catalog_total_uses_cheap_equality(app):
    with app.app_context():
        _add_alert(
            alert_id='ztf_candidate:3000000000000000001',
            mjd=60010.0,
            locus_id='L-cheap',
            object_id='Zcheap',
        )
        db.session.commit()
        build = build_catalog_query({'locus_id': 'L-cheap'})
        assert build.has_cheap_equality is True
        total, deferred = resolve_catalog_total(
            build,
            items=[object(), object()],
            has_next=True,
            page=1,
            use_keyset=True,
            keyset_direction='first',
        )
        assert deferred is False
        assert total == 1


def test_should_use_keyset_falls_back_to_offset_for_page_jump():
    assert should_use_keyset(True, {'page': '2'}) is False
    assert should_use_keyset(True, {}) is True
    assert should_use_keyset(False, {}) is False
    assert should_use_keyset(True, {'before_mjd': '1', 'before_alert_id': 'a', 'before_locus_id': 'l'}) is True


def test_parse_keyset_request_last_wins():
    direction, cursor = parse_keyset_request({'last': '1'})
    assert direction == 'last'
    assert cursor is None


def test_last_page_from_total():
    assert last_page_from_total(0) == 1
    assert last_page_from_total(100) == 1
    assert last_page_from_total(101) == 2


def test_catalog_href_and_query_string_strip_pagination():
    qs = catalog_query_string({'ant_passband': 'g', 'page': '2', 'before_mjd': '1'})
    assert qs == 'ant_passband=g'
    assert catalog_href(qs, {'page': '3'}) == '/?ant_passband=g&page=3'
    assert catalog_href('', None) == '/'


def test_unfiltered_list_shows_exact_short_total(client, app):
    with app.app_context():
        clear_catalog_total_cache()
        _add_alert(alert_id='ztf_candidate:4000000000000000001', mjd=60020.0, locus_id='L1', object_id='Z1')
        db.session.commit()
    response = client.get('/')
    assert response.status_code == 200
    assert b'id="catalog-rows-value"' in response.data
    assert b'>1<' in response.data
    assert b'id="btn-catalog-count"' not in response.data


def test_cheap_locus_filter_shows_exact_count(client, app):
    with app.app_context():
        _add_alert(alert_id='ztf_candidate:5000000000000000001', mjd=60021.0, locus_id='L-only', object_id='Z1')
        _add_alert(alert_id='ztf_candidate:5000000000000000002', mjd=60022.0, locus_id='L-other', object_id='Z2')
        db.session.commit()
    response = client.get('/?locus_id=L-only')
    assert response.status_code == 200
    assert b'id="btn-catalog-count"' not in response.data
    assert b'>1<' in response.data


def test_unselective_filter_defers_count_when_page_is_full(client, app, monkeypatch):
    monkeypatch.setattr('poi_broker.app.PAGE_SIZE', 2)
    monkeypatch.setattr('poi_broker.services.catalog_list.PAGE_SIZE', 2)
    with app.app_context():
        for index in range(3):
            _add_alert(
                alert_id=f'ztf_candidate:{6000000000000000000 + index}',
                mjd=60100.0 + index,
                locus_id=f'Lp{index}',
                object_id=f'Zp{index}',
                passband='g',
            )
        db.session.commit()
    response = client.get('/?ant_passband=g')
    assert response.status_code == 200
    assert b'2+' in response.data
    assert b'id="btn-catalog-count"' in response.data


def test_catalog_count_api_returns_exact_match_count(client, app):
    with app.app_context():
        _add_alert(
            alert_id='ztf_candidate:7000000000000000001',
            mjd=60200.0,
            locus_id='Lc1',
            object_id='Zc1',
            passband='g',
        )
        _add_alert(
            alert_id='ztf_candidate:7000000000000000002',
            mjd=60201.0,
            locus_id='Lc2',
            object_id='Zc2',
            passband='R',
        )
        db.session.commit()
    response = client.get('/api/catalog-count?ant_passband=g&sort__date=desc&page=2')
    assert response.status_code == 200
    assert response.get_json() == {'count': 1}


def test_catalog_count_includes_classification_only_when_needed(client, app):
    with app.app_context():
        _add_alert(
            alert_id='ztf_candidate:8000000000000000001',
            mjd=60300.0,
            locus_id='Lsn',
            object_id='Zsn',
        )
        db.session.add(Classification(alert_id='ztf_candidate:8000000000000000001', prob_class='sn'))
        db.session.commit()
        build = build_catalog_query({'prob_class': 'sn'}, include_sort=False)
        assert count_matches(build.count_query) == 1


def test_keyset_next_link_uses_cursor_not_page(client, app, monkeypatch):
    monkeypatch.setattr('poi_broker.app.PAGE_SIZE', 2)
    monkeypatch.setattr('poi_broker.services.catalog_list.PAGE_SIZE', 2)
    with app.app_context():
        for index in range(3):
            _add_alert(
                alert_id=f'ztf_candidate:{9000000000000000000 + index}',
                mjd=Time(60400.0 + index, format='mjd').mjd,
                locus_id=f'Lk{index}',
                object_id=f'Zk{index}',
            )
        db.session.commit()
    response = client.get('/')
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert 'before_mjd=' in body
    assert 'before_alert_id=' in body
    assert 'page=2' not in body


def test_sort_locus_id_is_applied_and_disables_keyset(app):
    with app.app_context():
        _add_alert(alert_id='ztf_candidate:1300000000000000003', mjd=60700.0, locus_id='L-c', object_id='Zc')
        _add_alert(alert_id='ztf_candidate:1300000000000000001', mjd=60700.0, locus_id='L-a', object_id='Za')
        _add_alert(alert_id='ztf_candidate:1300000000000000002', mjd=60700.0, locus_id='L-b', object_id='Zb')
        db.session.commit()

        asc_args = {'sort__locus_id': 'asc'}
        asc_build = build_catalog_query(asc_args)
        assert asc_build.is_date_only_sort is False
        assert should_use_keyset(asc_build.is_date_only_sort, asc_args) is False
        asc_items, _ = fetch_page(asc_build.list_query, 1, page_size=10)
        assert [row.locus_id for row in asc_items] == ['L-a', 'L-b', 'L-c']

        desc_args = {'sort__locus_id': 'desc'}
        desc_build = build_catalog_query(desc_args)
        assert desc_build.is_date_only_sort is False
        desc_items, _ = fetch_page(desc_build.list_query, 1, page_size=10)
        assert [row.locus_id for row in desc_items] == ['L-c', 'L-b', 'L-a']


def test_offset_page_past_end_is_404(client, app):
    with app.app_context():
        _add_alert(alert_id='ztf_candidate:1100000000000000001', mjd=60500.0, locus_id='Lx', object_id='Zx')
        db.session.commit()
    response = client.get('/?sort__locus_ra=asc&page=9')
    assert response.status_code == 404


def test_unselective_list_does_not_call_count_matches(client, app, monkeypatch):
    called = []

    def _fail_count(_query):
        called.append('count')
        return 0

    monkeypatch.setattr('poi_broker.app.count_matches', _fail_count)
    monkeypatch.setattr('poi_broker.services.catalog_list.count_matches', _fail_count)
    with app.app_context():
        for index in range(3):
            _add_alert(
                alert_id=f'ztf_candidate:{1200000000000000000 + index}',
                mjd=60600.0 + index,
                locus_id=f'Ln{index}',
                object_id=f'Zn{index}',
                passband='g',
            )
        db.session.commit()
    monkeypatch.setattr('poi_broker.app.PAGE_SIZE', 2)
    monkeypatch.setattr('poi_broker.services.catalog_list.PAGE_SIZE', 2)
    response = client.get('/?ant_passband=g')
    assert response.status_code == 200
    assert called == []


def test_catalog_alert_id_prefix_escapes_like_metacharacters(app):
    with app.app_context():
        build = build_catalog_query({'alert_id': 'ztf_'})
        compiled = build.list_query.statement.compile()
        sql = str(compiled)
        assert 'LIKE' in sql.upper()
        assert 'ESCAPE' in sql.upper()
        assert any(v == 'ztf\\_%' for v in compiled.params.values())
