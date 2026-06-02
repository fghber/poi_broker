(function (global) {
  'use strict';

  var DEFAULT_SPINNER_TEMPLATE = [
    '<div id="loading-spinner" class="d-none">',
    '  <div class="spinner-overlay">',
    '    <div class="spinner-border text-primary" role="status">',
    '      <span class="sr-only">Loading...</span>',
    '    </div>',
    '    <div class="spinner-text mt-2">Loading...</div>',
    '  </div>',
    '</div>'
  ].join('');

  function createLoadingSpinnerController(options) {
    var config = options || {};
    var spinnerElementId = config.spinnerElementId || 'loading-spinner';
    var spinnerContainerSelector = config.spinnerContainerSelector || null;
    var hiddenClass = config.hiddenClass || 'd-none';
    var defaultDelayMs = typeof config.defaultDelayMs === 'number' ? config.defaultDelayMs : 1500;
    var safetyTimeoutMs = typeof config.safetyTimeoutMs === 'number' ? config.safetyTimeoutMs : 30000;
    var spinnerTemplateHtml = config.spinnerTemplateHtml || DEFAULT_SPINNER_TEMPLATE;

    var spinnerTimeout = null;
    var spinnerSafetyTimeout = null;
    var spinnerVisible = false;

    function getSpinnerElement() {
      return document.getElementById(spinnerElementId);
    }

    function getSpinnerContainer() {
      if (spinnerContainerSelector) {
        return document.querySelector(spinnerContainerSelector);
      }
      return document.body;
    }

    function ensureSpinnerElement() {
      var existing = getSpinnerElement();
      if (existing) {
        return existing;
      }

      var container = getSpinnerContainer();
      if (!container) {
        return null;
      }

      var template = document.createElement('template');
      template.innerHTML = spinnerTemplateHtml.trim();
      var spinnerElement = template.content.firstElementChild;
      if (!spinnerElement) {
        return null;
      }

      if (!spinnerElement.id) {
        spinnerElement.id = spinnerElementId;
      }
      if (!spinnerElement.classList.contains(hiddenClass)) {
        spinnerElement.classList.add(hiddenClass);
      }

      container.appendChild(spinnerElement);
      return spinnerElement;
    }

    function show(delayMs) {
      if (spinnerVisible) return;

      if (spinnerTimeout) {
        clearTimeout(spinnerTimeout);
      }
      if (spinnerSafetyTimeout) {
        clearTimeout(spinnerSafetyTimeout);
      }

      var delay = typeof delayMs === 'number' ? delayMs : defaultDelayMs;
      spinnerTimeout = setTimeout(function () {
        var spinnerElement = ensureSpinnerElement();
        if (!spinnerElement) return;

        spinnerElement.classList.remove(hiddenClass);
        spinnerVisible = true;

        spinnerSafetyTimeout = setTimeout(function () {
          hide();
        }, safetyTimeoutMs);
      }, delay);
    }

    function hide() {
      if (spinnerTimeout) {
        clearTimeout(spinnerTimeout);
        spinnerTimeout = null;
      }
      if (spinnerSafetyTimeout) {
        clearTimeout(spinnerSafetyTimeout);
        spinnerSafetyTimeout = null;
      }

      if (!spinnerVisible) return;

      var spinnerElement = getSpinnerElement();
      if (spinnerElement) {
        spinnerElement.classList.add(hiddenClass);
      }
      spinnerVisible = false;
    }

    return {
      show: show,
      hide: hide
    };
  }

  global.createLoadingSpinnerController = createLoadingSpinnerController;
})(window);