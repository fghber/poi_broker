import pytest

import poi_broker.models as models_module
from poi_broker import db
from poi_broker.models import Classification, Ztf
from poi_broker.querybuilder_translator import Filter, TableNotFoundError


def test_filter_module_init_scans_models(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter(models_module, base_query)

        assert 'featuretable' in filter_obj.models
        assert 'classification' in filter_obj.models


def test_querybuilder_with_empty_rules_returns_query(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf}, base_query)

        filtered_query = filter_obj.querybuilder({'rules': []})
        assert filtered_query.whereclause is None


def test_querybuilder_unsupported_operator_raises(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf}, base_query)
        rules = {'rules': [{'field': 'featuretable.alert_id', 'operator': 'not_a_real_op', 'value': 'x'}]}

        with pytest.raises(NotImplementedError):
            filter_obj.querybuilder(rules)


def test_querybuilder_invalid_field_format_raises(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf}, base_query)
        rules = {'rules': [{'field': 'alert_id', 'operator': 'equal', 'value': 'x'}]}

        with pytest.raises(ValueError, match='Invalid field format'):
            filter_obj.querybuilder(rules)


def test_querybuilder_unknown_table_name_raises(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf}, base_query)
        rules = {'rules': [{'field': 'unknown.alert_id', 'operator': 'equal', 'value': 'x'}]}

        with pytest.raises(TableNotFoundError):
            filter_obj.querybuilder(rules)


def test_querybuilder_unknown_column_name_raises(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf}, base_query)
        rules = {'rules': [{'field': 'featuretable.no_column', 'operator': 'equal', 'value': 'x'}]}

        with pytest.raises(ValueError, match='Unknown field'):
            filter_obj.querybuilder(rules)


def test_querybuilder_between_invalid_payload_raises(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf}, base_query)
        rules = {'rules': [{'field': 'featuretable.date_alert_mjd', 'operator': 'between', 'value': [1]}]}

        with pytest.raises(ValueError, match='Operator "between" requires a two-item array value'):
            filter_obj.querybuilder(rules)


def test_querybuilder_in_invalid_payload_raises(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf}, base_query)
        rules = {'rules': [{'field': 'featuretable.alert_id', 'operator': 'in', 'value': []}]}

        with pytest.raises(ValueError, match='Operator "in" requires a non-empty array value'):
            filter_obj.querybuilder(rules)


def test_querybuilder_missing_value_raises(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf}, base_query)
        rules = {'rules': [{'field': 'featuretable.alert_id', 'operator': 'equal', 'value': None}]}

        with pytest.raises(ValueError, match='requires a value'):
            filter_obj.querybuilder(rules)


def test_querybuilder_unary_operator_is_null(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf}, base_query)
        rules = {'rules': [{'field': 'featuretable.ztf_object_id', 'operator': 'is_null'}]}

        filtered_query = filter_obj.querybuilder(rules)
        sql = str(filtered_query.whereclause.compile(compile_kwargs={'literal_binds': True}))

        assert 'IS NULL' in sql.upper()


def test_querybuilder_add_entity_and_nested_rules(app):
    with app.app_context():
        base_query = db.session.query(Ztf)
        filter_obj = Filter({'featuretable': Ztf, 'classification': Classification}, base_query)
        rules = {
            'condition': 'AND',
            'rules': [
                {'field': 'featuretable.locus_ra', 'operator': 'greater', 'value': 1},
                {
                    'condition': 'OR',
                    'rules': [
                        {'field': 'featuretable.locus_dec', 'operator': 'less', 'value': 2},
                        {'field': 'classification.p_cvnova', 'operator': 'greater_or_equal', 'value': 0.5},
                    ],
                },
            ],
        }

        filtered_query = filter_obj.querybuilder(rules)
        sql = str(filtered_query.whereclause.compile(compile_kwargs={'literal_binds': True}))

        assert 'featuretable.locus_ra > 1' in sql
        assert 'featuretable.locus_dec < 2' in sql
        assert 'classification.p_cvnova >= 0.5' in sql
        assert ' OR ' in sql.upper()
