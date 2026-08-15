from pathlib import Path


TEMPLATES_DIR = Path(__file__).resolve().parents[1] / "poi_broker" / "templates"


def _read_template(name: str) -> str:
    return (TEMPLATES_DIR / name).read_text(encoding="utf-8")


def test_profile_watchlist_renderer_escapes_user_fields():
    profile_html = _read_template("profile.html")

    # Regression guard: watchlist fields must be escaped before insertion.
    assert "${escapeHtml(item.name)}" in profile_html
    assert "${escapeHtml(item.sql_where || '')}" in profile_html


def test_profile_group_tabs_escape_labels():
    profile_html = _read_template("profile.html")

    # Regression guard: favorite group names must be escaped before tab insertion.
    assert "${escapeHtml(label)}" in profile_html


def test_profile_favorites_list_escapes_locus_ids():
    profile_html = _read_template("profile.html")

    # Regression guard: locusId must be escaped in HTML and encoded in hrefs.
    assert "${encodeURIComponent(fav.locusId)}" in profile_html
    assert "${escapeHtml(fav.locusId)}" in profile_html
    assert 'data-locus="${escapeHtml(fav.locusId)}"' in profile_html


def test_profile_filter_bookmarks_escape_path_display():
    profile_html = _read_template("profile.html")

    # Regression guard: bookmark path must not re-enter markup after decodeURIComponent.
    assert 'href="${escapeHtml(item.path)}"' in profile_html
    assert "${escapeHtml(decodeURIComponent(item.path || '/'))}" in profile_html


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


def test_pagination_spinner_uses_delegated_click():
    main_html = _read_template("main.html")
    assert "$(document).on('click', 'a[href*=\"page=\"]'" in main_html
    assert 'document.querySelectorAll(\'a[href*="page="]\')' not in main_html


def test_sort_restore_does_not_leak_id_global():
    main_html = _read_template("main.html")
    assert "$id = $('#'+key.substring(6))" not in main_html
    assert "var $th = $('#' + key.substring(6));" in main_html
    restore = main_html.split("adhere to sort order provided in URL params", 1)[1]
    restore = restore.split("</script>", 1)[0]
    assert "hasClass('asc')" not in restore
    assert "hasClass('desc')" not in restore


def test_favorite_click_uses_inflight_lock():
    main_html = _read_template("main.html")
    assert "if (btn.data('favXhr')) return;" in main_html
    assert "btn.removeData('favXhr')" in main_html
    assert "btn.data('favXhr', xhr);" in main_html


def test_observing_retrieve_is_not_inline_onclick():
    main_html = _read_template("main.html")
    assert 'onclick="query_observing_plot()"' not in main_html
    assert 'data-role="observing-retrieve"' in main_html
    assert "$objectIdModal.on('click', '[data-role=\"observing-retrieve\"]'" in main_html


def test_profile_delete_buttons_use_closest():
    profile_html = _read_template("profile.html")
    assert "clickEl.closest('.remove-fav-btn')" in profile_html
    assert "clickEl.closest('.delete-group-btn')" in profile_html
    assert "clickEl.closest('.remove-watchlist-btn')" in profile_html
    assert "clickEl.closest('.remove-filter-bookmark-btn')" in profile_html
    assert "e.target.classList.contains('remove-watchlist-btn')" not in profile_html
    assert "e.target.classList.contains('delete-group-btn')" not in profile_html


def test_query_builder_invalidates_stale_preview_and_count():
    qb_js = _read_template("_query_builder_js.html")
    assert "function abortPreviewRequests()" in qb_js
    assert "if (seq !== previewSeq) return;" in qb_js
    assert "Query is valid, but row count failed." in qb_js
    assert "timeout: 15000" in qb_js


def test_query_builder_destroy_splits_tooltip_and_plugin_try():
    qb_js = _read_template("_query_builder_js.html")
    block = qb_js.split("function destroyBuilderPlugin()", 1)[1].split("function resetState", 1)[0]
    assert block.count("try {") >= 2
    assert "tooltip('dispose')" in block
    assert "queryBuilder('destroy')" in block
    dispose_try, _, rest = block.partition("tooltip('dispose')")
    assert "queryBuilder('destroy')" not in dispose_try
    assert "queryBuilder('destroy')" in rest
    profile_html = _read_template("profile.html")
    assert "$('body > .tooltip').remove()" in profile_html
    assert "$('.tooltip.show, .tooltip').remove()" not in profile_html


def test_main_query_urls_encode_locus_and_alert_ids():
    main_html = _read_template("main.html")

    assert "'/locus_plot_csv?locusId=' + encodeURIComponent(locus_id)" in main_html
    assert "`/api/favorite?locusId=${encodeURIComponent(locus_id)}`" in main_html
    assert "`/query_features?alert_id=${encodeURIComponent(alert_id)}`" in main_html
    assert "`/query_crossmatches?locusId=${encodeURIComponent(locusId)}`" in main_html
    assert "cdn.bokeh.org" not in main_html
    assert "url_for('main.bokeh_js'" in main_html


def test_query_builder_assets_are_local_only():
    qb_js = _read_template("_query_builder_js.html")
    qb_css = _read_template("_query_builder_assets.html")

    assert "cdn.jsdelivr.net" not in qb_js
    assert "onerror=" not in qb_js
    assert "cdn.jsdelivr.net" not in qb_css
    assert "onerror=" not in qb_css


def test_classification_ajax_has_fail_handler():
    main_html = _read_template("main.html")
    block = main_html.split("function query_classification(alertId, sessionId) {", 1)[1]
    block = block.split("function renderAlert", 1)[0]
    assert "dataType: 'json'" in block
    assert "encodeURIComponent(alertId)" in block
    assert "insertBokehFromPayload" in block
    assert "tmp.innerHTML = data" not in block
    assert "new Function(rawJs)" not in block
    assert "plotRequestFail(xhr, statusText, sessionId, '#classification_output'" in block
    assert "Error loading classification data." in block
    assert "timeout: 15000" in block


def test_plot_ajax_uses_json_bokeh_payload():
    main_html = _read_template("main.html")
    assert "function insertBokehFromPayload(container, payload)" in main_html
    assert "function plotRequestFail(xhr, statusText, sessionId, selector, fallback)" in main_html
    assert "insertBokehFromPayload(document.getElementById('locus-plot'), payload)" in main_html
    assert "insertBokehFromPayload(document.getElementById('feature-plot'), payload)" in main_html
    featureplot_fail = (
        main_html.split("function query_featureplot", 1)[1]
        .split("const queryObservingPlot", 1)[0]
        .split(".fail(function(xhr, statusText) {", 1)[1]
    )
    before_hide = featureplot_fail.split("loadingSpinnerController.hide()", 1)[0]
    assert "statusText !== 'abort'" in before_hide or "statusText === 'abort'" in before_hide
    assert "plotRequestFail(xhr, statusText, sessionId, '#feature-plot'" in featureplot_fail
    assert "insertBokehFromPayload(document.getElementById('classification_output'), payload)" in main_html
    assert "$(\"#observing_output\").html(data)" not in main_html
    assert "function renderObservingOutput(payload)" in main_html
    block = main_html.split("function renderObservingOutput(payload) {", 1)[1]
    block = block.split("function persistLastObservatory", 1)[0]
    assert "renderAlert('#observing_output', 'warning', payload.message)" in block
    assert "moonAlert.textContent = payload.moonMessage" in block
    assert "alert alert-info" in block
    assert "moonCol.innerHTML = payload.moonHtml" in block
    assert "innerHTML = payload.message" not in block
    assert "innerHTML = payload.moonMessage" not in block
    assert "data:image/png;base64," in main_html
    assert "fetch('/api/last-observatory'" in main_html
    assert "'X-CSRFToken': csrfToken" in main_html


def test_export_error_uses_textcontent():
    export_html = _read_template("export.html")
    start_fn = export_html.split("function startExport(rules)", 1)[1].split("// Auto-refresh", 1)[0]
    error_branch = start_fn.split("} else {", 1)[1].split(".catch(", 1)[0]
    assert "res.data.error" in error_branch
    assert "innerHTML" not in error_branch
    assert "(res.data.error || 'Export failed') + '</span>'" not in export_html
    assert "createElement('i')" in error_branch
    assert "glyphicon-exclamation-sign" in error_branch
    assert "createTextNode" in error_branch or "textContent" in error_branch


def test_favorite_post_has_timeout_and_logs_non_auth_errors():
    main_html = _read_template("main.html")
    block = main_html.split("url: '/api/favorite'", 1)[1].split("$('#btn_save_filter_bookmark')", 1)[0]
    assert "timeout: 10000" in block
    assert "Failed to save favorite" in block


def test_dead_lightcurve_helpers_removed():
    main_html = _read_template("main.html")
    assert "function UrlExists" not in main_html
    assert "function CsvExists" not in main_html
    assert "function generate_lightcurveJS" not in main_html


def test_export_count_failure_does_not_start_export():
    export_html = _read_template("export.html")
    assert "Could not check export size. Please retry." in export_html
    assert export_html.count("startExport(rules);") == 1
    assert "pagehide" in export_html
    start_fn = export_html.split("function startExport(rules)", 1)[1].split("// Auto-refresh", 1)[0]
    assert start_fn.rstrip().endswith("}")
    assert not start_fn.rstrip().endswith("});")


def test_site_footer_block_is_not_nested_in_body():
    site_html = _read_template("site.html")
    without_content = site_html.replace("{% block content %}{% endblock %}", "")
    body_pos = without_content.find("{% block body %}")
    foot_pos = without_content.find("{% block foot %}")
    body_end = without_content.find("{% endblock %}", body_pos)
    assert body_pos != -1 and foot_pos != -1
    assert body_end < foot_pos


def test_profile_group_create_and_move_check_response_ok():
    profile_html = _read_template("profile.html")
    assert "if (!r.ok) throw new Error('Failed to create group');" in profile_html
    assert "if (!r.ok) throw new Error('Failed to move favorite');" in profile_html
    assert "confirmMoveBtn.disabled = true;" in profile_html
    assert "confirmMoveBtn.disabled = false;" in profile_html


def test_bootstrap_datepicker_assets_are_gone():
    root = Path(__file__).resolve().parents[1]
    leftover = list((root / "poi_broker" / "static").glob("**/bootstrap-datepicker*"))
    assert leftover == [], leftover
    templates_dir = root / "poi_broker" / "templates"
    for path in templates_dir.glob("*.html"):
        assert "bootstrap-datepicker" not in path.read_text(encoding="utf-8")
