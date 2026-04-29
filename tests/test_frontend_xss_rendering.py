from pathlib import Path


TEMPLATES_DIR = Path(__file__).resolve().parents[1] / "poi_broker" / "templates"


def _read_template(name: str) -> str:
    return (TEMPLATES_DIR / name).read_text(encoding="utf-8")


def test_profile_watchlist_renderer_escapes_user_fields():
    profile_html = _read_template("profile.html")

    # Regression guard: watchlist fields must be escaped before insertion.
    assert "${escapeHtml(item.name)}" in profile_html
    assert "${escapeHtml(item.sql_where || '')}" in profile_html


def test_profile_toast_renderer_avoids_html_injection_for_messages():
    profile_html = _read_template("profile.html")

    # Regression guard: toast messages must not be interpolated into HTML strings.
    assert "toastContainer.insertAdjacentHTML" not in profile_html
    assert "showToast.template" in profile_html
    assert "body.textContent = message == null ? '' : String(message);" in profile_html


def test_main_modal_table_renderers_bind_dynamic_values_with_textcontent():
    main_html = _read_template("main.html")

    # Regression guard: central cell binder must use textContent for all dynamic values.
    assert "function setCells(root, values)" in main_html
    assert "cell.textContent = value == null ? '' : String(value);" in main_html

    # Renderers should flow through the safe helper path.
    assert "appendTemplateRows(table, 'features-table-row-template'" in main_html
    assert "appendTemplateRows(tbody, 'crossmatches-table-row-template'" in main_html


def test_main_auth_required_toast_uses_inline_template_clone():
    main_html = _read_template("main.html")

    # Regression guard: auth-required toast should avoid string HTML injection.
    assert "showAuthRequiredToast.template" in main_html
    assert "toast = showAuthRequiredToast.template.content.firstElementChild.cloneNode(true);" in main_html
    assert "toast.innerHTML =" not in main_html


def test_main_modal_table_renderers_no_legacy_html_string_rendering():
    main_html = _read_template("main.html")

    # Regression guard: ensure legacy unsafe string-template renderers are gone.
    assert "json2featurestable" not in main_html
    assert "json2crossmatchestable" not in main_html

    # Guard the specific endpoints now using safe renderers.
    assert "renderFeaturesTable('#features_table_output', response);" in main_html
    assert "renderCrossmatchesTable('#crossmatches_table_output', response);" in main_html


def test_main_feature_query_uses_success_error_instead_of_complete_parse():
    main_html = _read_template("main.html")
    function_block = main_html.split("function query_features_table(alert_id) {", 1)[1].split("function query_featureplot", 1)[0]

    assert "function query_features_table(alert_id)" in main_html
    assert "success: function(response)" in function_block
    assert "error: function(xhr, status, error)" in function_block
    assert "var json_obj = JSON.parse(r.responseText);" not in function_block
    assert "complete: function(r){" not in function_block
