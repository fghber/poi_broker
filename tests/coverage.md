# tests coverage report generated with `python -m pytest --cov=poi_broker --cov-report=term`

|Name                                       |Stmts  | Miss | Cover|
|-------------------------------------------|-------|------|------|
|poi_broker\__init__.py                     |   89  |    6 |   93%|
|poi_broker\app.py                          |  269  |   92 |   66%|
|poi_broker\auth.py                         |  227  |   94 |   59%|
|poi_broker\classification.py               |   55  |    1 |   98%|
|poi_broker\constants\__init__.py           |    0  |    0 |  100%|
|poi_broker\constants\features.py           |    5  |    0 |  100%|
|poi_broker\helpers.py                      |   25  |    1 |   96%|
|poi_broker\models.py                       |  404  |   19 |   95%|
|poi_broker\observing_tool.py               |  226  |   39 |   83%|
|poi_broker\querybuilder_translator.py      |   82  |    5 |   94%|
|poi_broker\routes\__init__.py              |    7  |    0 |  100%|
|poi_broker\routes\favorites.py             |   66  |    2 |   97%|
|poi_broker\routes\features.py              |   41  |   11 |   73%|
|poi_broker\routes\filter_bookmarks.py      |  105  |   27 |   74%|
|poi_broker\routes\lightcurve.py            |   40  |   10 |   75%|
|poi_broker\routes\user_observatories.py    |   96  |   15 |   84%|
|poi_broker\routes\visual_query.py          |  109  |   35 |   68%|
|poi_broker\services\__init__.py            |    0  |    0 |  100%|
|poi_broker\services\email_service.py       |   56  |    0 |  100%|
|poi_broker\services\favorites_service.py   |  125  |   15 |   88%|
|poi_broker\services\feature_service.py     |   39  |   12 |   69%|
|poi_broker\services\filter_service.py      |   54  |    1 |   98%|
|poi_broker\services\input_parser.py        |   73  |    3 |   96%|
|poi_broker\services\plotting_service.py    |   58  |    2 |   97%|
|poi_broker\services\query_service.py       |   27  |    3 |   89%|
|poi_broker\services\search_service.py      |   30  |    0 |  100%|
|poi_broker\settings.py                     |   50  |    7 |   86%|
|poi_broker\user_settings.py                |   98  |   17 |   83%|
|-------------------------------------------|-------|------|------|
|TOTAL                                      | 2456  |  417 |   83%|


# Fast parallel execution (recommended for development)
`python -m pytest -n auto -q --tb=no`

# Skip coverage for even faster runs
`python -m pytest -n auto -q --tb=no --no-cov -m "not slow"`

# Run only unit tests (if you add markers)
`python -m pytest -n auto -q -m "not slow"`