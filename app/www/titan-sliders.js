/* TITAN — bipolar (centre-fill) weight sliders
 *
 * Shiny uses ionRangeSlider for sliderInput(). The filled bar (.irs-bar) always
 * runs from the left edge to the handle. For ±1 dims we want it to fill from
 * the centre outward (green = positive, red = negative).
 *
 * Strategy:
 *  1. Wrap each bipolar slider in a .titan-bipolar-slider div (done in R).
 *  2. On init, inject a custom fill <span> inside .irs-line.
 *  3. On shiny:inputchanged, update the fill based on the new value.
 *  4. Use MutationObserver to catch sliders added after page load (tab switch).
 */

(function () {
  'use strict';

  // Map inputId → container element (populated as sliders are initialised)
  var registry = {};

  // ── Position the centre fill for a given value (-1..1) ─────────────────────
  function setFill(container, value) {
    var fill = container._titanFill;
    if (!fill) return;

    var v = parseFloat(value) || 0;

    if (Math.abs(v) < 0.005) {
      // Exactly at zero — hide fill, just show centre marker
      fill.style.left  = '50%';
      fill.style.width = '0px';
      fill.className   = 'titan-bipolar-fill';
      return;
    }

    var halfPct = Math.min(Math.abs(v), 1) * 50; // 0..50

    if (v > 0) {
      fill.style.left  = '50%';
      fill.style.width = halfPct + '%';
      fill.className   = 'titan-bipolar-fill pos';
    } else {
      fill.style.left  = (50 - halfPct) + '%';
      fill.style.width = halfPct + '%';
      fill.className   = 'titan-bipolar-fill neg';
    }
  }

  // ── Initialise one .titan-bipolar-slider container ──────────────────────────
  function initContainer(container) {
    if (container._titanInit) return;

    var irsLine = container.querySelector('.irs-line');
    if (!irsLine) return; // ionRangeSlider not yet rendered — will retry via MO

    container._titanInit = true;

    // Centre tick mark
    var marker = document.createElement('span');
    marker.className = 'titan-bipolar-center';
    irsLine.appendChild(marker);

    // Colour fill
    var fill = document.createElement('span');
    fill.className = 'titan-bipolar-fill';
    irsLine.appendChild(fill);
    container._titanFill = fill;

    // Register by input id (the hidden <input> ionRangeSlider controls)
    var input = container.querySelector('input.js-range-slider');
    if (input && input.id) {
      registry[input.id] = container;
    }

    // Set initial fill (default value = 0)
    setFill(container, 0);
  }

  // ── Scan DOM for uninitialised bipolar containers ───────────────────────────
  function scanAll() {
    document.querySelectorAll('.titan-bipolar-slider').forEach(initContainer);
  }

  // ── React to Shiny input changes (user drag or server updateSliderInput) ────
  $(document).on('shiny:inputchanged', function (e) {
    var container = registry[e.name];
    if (container) setFill(container, e.value);
  });

  // ── MutationObserver: catch sliders rendered after page load ─────────────────
  var mo = new MutationObserver(scanAll);

  $(document).on('shiny:connected', function () {
    mo.observe(document.body, { childList: true, subtree: true });
    scanAll();
  });

})();

/* ── Bootstrap popovers for peptide hits in the sequence viewer ──────────────
 * Uses a MutationObserver so popovers auto-init whenever the sequence card
 * re-renders (different ORF selected).  The data-tt-init sentinel prevents
 * double-initialisation of elements that are already set up.
 */
(function () {
  'use strict';

  function initPepPopovers() {
    document.querySelectorAll('[data-bs-toggle="popover"]:not([data-tt-init])').forEach(function (el) {
      el.setAttribute('data-tt-init', '1');
      if (typeof bootstrap !== 'undefined' && bootstrap.Popover) {
        new bootstrap.Popover(el, {
          container  : 'body',
          html       : true,
          sanitize   : false,
          trigger    : 'hover focus',
          delay      : { show: 180, hide: 120 },
          customClass: 'pep-popover'
        });
      }
    });
  }

  var popMO = new MutationObserver(initPepPopovers);

  $(document).on('shiny:connected', function () {
    popMO.observe(document.body, { childList: true, subtree: true });
    initPepPopovers();
  });

})();
