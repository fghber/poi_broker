# Flask Astronomy App

Python 3.12, Flask, SQLAlchemy, pytest.

Architecture:
- App factory pattern
- Blueprints (`routes/`)
- Service layer (`services/`)

Key files in `poi_broker/`:
- `models.py` → data models
- `routes/` → endpoints
- `services/` → business logic
- `tests/` → pytest tests
- `wsgi.py` → entrypoint

Conventions:
- Type hints
- Descriptive names
- `lowercase_with_underscores`

Commands:
- Environment: `source /c/dev/python/poi_broker/.venv/Scripts/activate`
- Run: `python -m flask --app wsgi:app run --debug`
- Test: `python -m pytest -q`

For new features: models → routes → services → tests.