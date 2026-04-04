import types
from inspect import signature

from sqlalchemy import and_, or_
from sqlalchemy import inspect as sa_inspect
from sqlalchemy.exc import NoInspectionAvailable
from sqlalchemy.orm import DeclarativeMeta


class TableNotFoundError(Exception):
    pass


OPERATORS = {
    'equal': lambda f, a: f.__eq__(a),
    'not_equal': lambda f, a: f.__ne__(a),
    'less': lambda f, a: f.__lt__(a),
    'greater': lambda f, a: f.__gt__(a),
    'less_or_equal': lambda f, a: f.__le__(a),
    'greater_or_equal': lambda f, a: f.__ge__(a),
    'in': lambda f, a: f.in_(a),
    'not_in': lambda f, a: f.notin_(a),
    'ends_with': lambda f, a: f.like('%' + a),
    'begins_with': lambda f, a: f.like(a + '%'),
    'contains': lambda f, a: f.like('%' + a + '%'),
    'not_contains': lambda f, a: f.notlike('%' + a + '%'),
    'not_begins_with': lambda f, a: f.notlike(a + '%'),
    'not_ends_with': lambda f, a: f.notlike('%' + a),
    'is_empty': lambda f: f.__eq__(''),
    'is_not_empty': lambda f: f.__ne__(''),
    'is_null': lambda f: f.is_(None),
    'is_not_null': lambda f: f.is_not(None),
    'between': lambda f, a: f.between(a[0], a[1]),
}


class Filter(object):
    """Translate jQuery QueryBuilder rules to SQLAlchemy ORM filters."""

    def __init__(self, models, query, operators=None):
        if isinstance(models, types.ModuleType):
            model_dict = {}
            for attr in models.__dict__.values():
                if isinstance(attr, DeclarativeMeta):
                    try:
                        inspected = sa_inspect(attr)
                        table_name = getattr(inspected.persist_selectable, 'name', None)
                        if table_name:
                            model_dict[table_name] = attr
                    except NoInspectionAvailable:
                        pass
            self.models = model_dict
        else:
            self.models = dict(models)

        self.query = query
        self.operators = operators if operators else OPERATORS

    def querybuilder(self, rules):
        query, cond_list = self._make_query(self.query, rules)
        if not cond_list:
            return query

        condition = rules.get('condition', 'AND').upper()
        operator = or_ if condition == 'OR' else and_
        return query.filter(operator(*cond_list))

    def _make_query(self, query, rules):
        cond_list = []
        for cond in rules.get('rules', []):
            if 'condition' not in cond:
                operator_name = cond.get('operator')
                if operator_name not in self.operators:
                    raise NotImplementedError(f'Unsupported operator: {operator_name}')

                field_name = cond.get('field', '')
                parts = field_name.split('.', 1)
                if len(parts) != 2:
                    raise ValueError(f'Invalid field format: {field_name}')

                table_name, column_name = parts
                try:
                    model = self.models[table_name]
                except KeyError:
                    raise TableNotFoundError(table_name)

                for table in query.column_descriptions:
                    if table.get('entity') == model:
                        break
                else:
                    query = query.add_entity(model)

                try:
                    field = getattr(model, column_name)
                except AttributeError as exc:
                    raise ValueError(f'Unknown field: {field_name}') from exc

                function = self.operators[operator_name]
                arity = len(signature(function).parameters)
                if arity == 1:
                    cond_list.append(function(field))
                elif arity == 2:
                    value = cond.get('value')

                    # Validate known multi-value operators so malformed payloads
                    # become client errors (ValueError) instead of server errors.
                    if operator_name == 'between':
                        if not isinstance(value, (list, tuple)) or len(value) != 2:
                            raise ValueError('Operator "between" requires a two-item array value')
                    elif operator_name in ('in', 'not_in'):
                        if not isinstance(value, (list, tuple)) or len(value) == 0:
                            raise ValueError(f'Operator "{operator_name}" requires a non-empty array value')
                    else:
                        if value is None:
                            raise ValueError(f'Operator "{operator_name}" requires a value')

                    try:
                        cond_list.append(function(field, value))
                    except (TypeError, ValueError, IndexError) as exc:
                        raise ValueError(f'Invalid value for operator "{operator_name}"') from exc
                else:
                    raise NotImplementedError(f'Unsupported operator arity: {arity}')
            else:
                query, cond_subrule = self._make_query(query, cond)
                nested_condition = cond.get('condition', 'AND').upper()
                nested_operator = or_ if nested_condition == 'OR' else and_
                if cond_subrule:
                    cond_list.append(nested_operator(*cond_subrule))

        return query, cond_list
