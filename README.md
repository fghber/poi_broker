# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/fghber/poi_broker/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                       |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| poi\_broker/\_\_init\_\_.py                |       79 |        5 |       10 |        3 |     89% |64, 89-\>92, 136-139 |
| poi\_broker/app.py                         |      225 |      126 |      102 |       19 |     36% |28-30, 41-\>176, 45-49, 52-56, 59-66, 69, 74-75, 78, 81-85, 88-92, 95-99, 102-104, 108-112, 119-122, 127-130, 135-138, 143-146, 151-154, 211-216, 228-239, 244-282, 291, 305-307, 312-319, 323-324, 330-337 |
| poi\_broker/auth.py                        |      277 |      205 |       68 |        3 |     22% |19, 22-31, 34-45, 65-66, 71-72, 83-84, 92, 108-109, 114-183, 188-210, 216-227, 231, 236-266, 270-276, 281-314, 320, 330-372, 386-436 |
| poi\_broker/classification.py              |       55 |       45 |       10 |        0 |     15% |    15-142 |
| poi\_broker/constants/\_\_init\_\_.py      |        0 |        0 |        0 |        0 |    100% |           |
| poi\_broker/constants/features.py          |        5 |        0 |        0 |        0 |    100% |           |
| poi\_broker/helpers.py                     |       82 |       61 |       26 |        0 |     19% |12, 17, 24-31, 35-55, 58-60, 63-64, 67-69, 73-96, 103-105, 109-116 |
| poi\_broker/models.py                      |      389 |       18 |        2 |        0 |     95% |304, 349, 352, 355, 359-367, 382, 398, 415, 429, 443 |
| poi\_broker/observing\_tool.py             |      157 |      134 |       28 |        0 |     12% |40-205, 215-272 |
| poi\_broker/querybuilder\_translator.py    |       82 |       31 |       38 |       11 |     57% |42-52, 62, 74, 79, 84-85, 88-\>87, 91, 95-96, 114-128 |
| poi\_broker/routes/\_\_init\_\_.py         |        6 |        0 |        0 |        0 |    100% |           |
| poi\_broker/routes/favorites.py            |       66 |        2 |       14 |        2 |     95% |    40, 83 |
| poi\_broker/routes/features.py             |       41 |       11 |       10 |        3 |     69% |32-35, 46, 54-56, 66-68 |
| poi\_broker/routes/filter\_bookmarks.py    |      105 |       27 |       34 |       10 |     73% |50, 57, 60, 62, 67, 82-85, 91, 118, 124, 133, 144-150, 165-171 |
| poi\_broker/routes/lightcurve.py           |       40 |       10 |        8 |        3 |     69% |22, 39-41, 51, 65-66, 73-75 |
| poi\_broker/routes/visual\_query.py        |      109 |       35 |       20 |        8 |     65% |32, 36-\>39, 43-45, 55, 59-\>62, 64-69, 79, 83, 85, 90-93, 111-119, 133-135, 149-154 |
| poi\_broker/services/\_\_init\_\_.py       |        0 |        0 |        0 |        0 |    100% |           |
| poi\_broker/services/favorites\_service.py |      125 |       46 |       42 |       13 |     62% |27, 30, 49, 54-\>57, 74, 77, 88-92, 95-\>99, 101-107, 122, 140-146, 157, 186-188, 202, 205, 220-226, 240, 253-259 |
| poi\_broker/services/feature\_service.py   |       39 |       12 |        6 |        2 |     64% |21, 45-48, 67-71, 82-84 |
| poi\_broker/services/plotting\_service.py  |       58 |       45 |       32 |        2 |     17% |45-67, 84-128 |
| poi\_broker/services/query\_service.py     |       27 |        3 |        6 |        3 |     82% |26, 29, 42 |
| poi\_broker/settings.py                    |       50 |        7 |       12 |        5 |     81% |21, 42, 54, 58, 62-63, 65 |
| poi\_broker/user\_settings.py              |       62 |       15 |       16 |        5 |     74% |18, 22-\>24, 37-38, 60, 63-64, 67-\>69, 74-82 |
| **TOTAL**                                  | **2079** |  **838** |  **484** |   **92** | **54%** |           |


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/fghber/poi_broker/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/fghber/poi_broker/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/fghber/poi_broker/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/fghber/poi_broker/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2Ffghber%2Fpoi_broker%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/fghber/poi_broker/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.