// Checkbox-dropdown ("tickbox") column filters for the small-variant and SV tables.
//
// DT gives every column a free-text search box (filter = "top"), which only helps a reader
// who already knows what the column contains — nobody remembers the ~30 VEP consequence
// terms. For the categorical columns this replaces that box with a dropdown listing the
// values actually present in this sample, with per-value row counts, so the filter doubles
// as a summary of the column.
//
// R publishes the values (window.SNV_FACETS / window.SV_FACETS, from js_facet_defs() in
// R/utils.R) and the column-name -> index maps (window.SNV_COLS / window.SV_COLS) that
// every by-name lookup here goes through. Four things about the DOM this works in, each
// learned from DT/DataTables internals rather than assumed:
//
//   1. DT's own `div.form-group` in the filter cell is KEPT and the <input> inside it only
//      hidden. DT binds the handler that un-clips `overflow:hidden` on
//      .dataTables_scrollHead to that div, so triggering its "show"/"hide" events is what
//      lets an open menu escape the scroll head under scrollX. The widget is appended
//      inside that div so DT's own $td.children('div').first()/.last() lookups still hold.
//   2. No `id` attributes anywhere in this markup, and no handlers delegated from
//      `document` onto widget classes. DataTables rebuilds a zero-height sizing clone of
//      the whole header in the scroll body on every draw, refilling each cell by innerHTML
//      after stripping ids — an id here would reappear duplicated once per draw. Hence
//      aria-expanded rather than aria-controls, handlers bound to the nodes we create, and
//      every query rooted at api.table().header(), never at `document`.
//   3. The R-side column index is NOT the DOM position: DataTables detaches the filter
//      cells of visible:false columns. columns(':visible').indexes() translates it, and the
//      header label above the chosen cell is compared to the column name — a reorder that
//      broke the translation would otherwise attach a working dropdown to the wrong column.
//   4. Ticked state lives here, in a closure, not in the DOM and not in DT's per-column
//      search slots. That is what makes it survive a gene-panel redraw, Scroller's node
//      recycling and the per-draw header clone.
(function () {
  "use strict";

  // consequence is ~30 terms; impact is 4. Below this a find box is just clutter.
  var SEARCH_THRESHOLD = 12;
  // One draw over 167k rows runs two ext.search predicates per row, so coalesce rapid ticks.
  var DRAW_DEBOUNCE_MS = 120;
  var LABEL_MAX = 16;

  // Which tables get facets, keyed on the globals they publish from initComplete — the same
  // identity test the gene-panel predicate in per_sample.qmd uses, which excludes the QC,
  // ASCAT and phasing DT tables without any extra condition.
  var TABLES = [
    { key: "snv",
      elem:      function () { return window.snvTableElem; },
      cols:      function () { return window.SNV_COLS; },
      facets:    function () { return window.SNV_FACETS; },
      requested: function () { return window.SNV_FACET_COLS; } },
    { key: "sv",
      elem:      function () { return window.svTableElem; },
      cols:      function () { return window.SV_COLS; },
      facets:    function () { return window.SV_FACETS; },
      requested: function () { return window.SV_FACET_COLS; } }
  ];

  // key -> { colName -> Set of ticked tokens }. `null` in a Set is the "(none)" bucket.
  // A column absent from a table's object is unconstrained.
  var STATE    = { snv: {}, sv: {} };
  var ATTACHED = {};
  var WIDGETS  = [];      // every built widget, for the "clear all" control
  var PENDING  = {};      // per-table draw debounce timers
  var OPEN     = null;    // the one open menu

  function tableFor(node) {
    for (var i = 0; i < TABLES.length; i++) if (TABLES[i].elem() === node) return TABLES[i];
    return null;
  }

  function esc(s) {
    return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;")
                    .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }

  // ---- matching ------------------------------------------------------------
  // Mirrors how the cell was built, which is why the separator comes from R rather than
  // being guessed here: `consequence` is VEP's "&"-joined terms rewritten to commas
  // (R/parse_smallvariants.R) and `callers` is paste(sort(unique(caller)), collapse = ",").
  // NA and "" both collapse to the single null bucket the dropdown labels "(none)".
  function tokens(v, sep) {
    if (v === null || v === undefined) return [null];
    var s = String(v).trim();
    if (s === "") return [null];
    if (!sep) return [s];
    var out = s.split(sep).map(function (x) { return x.trim(); })
                          .filter(function (x) { return x.length > 0; });
    return out.length ? out : [null];
  }

  // One predicate for both tables, keyed on the table node exactly as the gene-panel
  // predicate is. OR within a column, AND across columns. Both predicates are ext.search
  // entries, so DataTables ANDs them: a tick composes with the panel filter and with every
  // remaining per-column text box for free.
  //
  // Chosen over column().search() with a regex because the token rule is the load-bearing
  // part: as a regex it would be re-encoded in a second language with an escaping
  // obligation (svclass values carry spaces, callers carries "clairs-to"), and an
  // unescaped metacharacter yields a filter that matches MORE rows and says nothing.
  $.fn.dataTable.ext.search.push(function (settings, searchData, index, rowData) {
    var t = tableFor(settings.nTable);
    if (!t) return true;                                  // another DT table in the report
    var sel = STATE[t.key], names = Object.keys(sel);
    if (names.length === 0) return true;
    var cols = t.cols() || {}, defs = t.facets() || {};
    // rowData is the original (typed) row, searchData the rendered strings. Always the
    // former: the SV table renders `locus`/`size` orthogonally, and mixing the two sources
    // per column is the trap variant_key() documents on the R side.
    var row = rowData || searchData;
    for (var i = 0; i < names.length; i++) {
      var chosen = sel[names[i]];
      if (!chosen || chosen.size === 0) continue;
      var idx = cols[names[i]];
      // Defensive only — widgets are built solely for columns that resolved. A vanished
      // column must not empty the table; the loud failure is the attach-time warning.
      if (idx === undefined) continue;
      var def  = defs[names[i]];
      var toks = tokens(row[idx], def ? def.sep : null);
      var hit  = false;
      for (var k = 0; k < toks.length; k++) {
        if (chosen.has(toks[k])) { hit = true; break; }
      }
      if (!hit) return false;                             // AND across columns
    }
    return true;
  });

  function scheduleDraw(w) {
    clearTimeout(PENDING[w.key]);
    PENDING[w.key] = setTimeout(function () {
      // resetPaging default: a filter change should return to the top, which is also what
      // Scroller does with the scroll position.
      w.api.draw();
    }, DRAW_DEBOUNCE_MS);
  }

  // ---- the open menu -------------------------------------------------------
  function closeOpen() {
    if (!OPEN) return;
    var w = OPEN;
    w.$menu.prop("hidden", true).removeClass("facet__menu--right");
    w.$btn.attr("aria-expanded", "false");
    w.$wrap.trigger("hide");                 // DT re-clips the scroll head
    w.$wrapper.removeClass("facet-open");    // and the report's own wrapper overflow
    OPEN = null;
  }

  function openMenu(w) {
    if (OPEN === w) { closeOpen(); return; }
    closeOpen();
    w.$wrap.trigger("show");
    w.$wrapper.addClass("facet-open");
    w.$menu.prop("hidden", false);
    w.$btn.attr("aria-expanded", "true");
    OPEN = w;
    // Flip to the right edge when the menu would run off the viewport. Cheap, and the
    // faceted columns sit far enough right in the SV table for it to matter.
    var box = w.$menu[0].getBoundingClientRect();
    if (box.right > window.innerWidth - 8) w.$menu.addClass("facet__menu--right");
    if (w.$find.length) w.$find.trigger("focus");
  }

  // ---- widget -------------------------------------------------------------
  function markup(name, def) {
    var withFind = def.values.length > SEARCH_THRESHOLD;
    var opts = def.values.map(function (v) {
      var tok = v[0], n = v[1];
      var isNone = (tok === null);
      var label = isNone ? "(none)" : tok;
      return '<label class="facet__opt' + (isNone ? " facet__opt--none" : "") + '"' +
             ' title="' + esc(label) + '">' +
             '<input type="checkbox"' +
             (isNone ? ' data-none="1"' : ' data-tok="' + esc(tok) + '"') + '>' +
             '<span class="facet__v">' + esc(label) + '</span>' +
             '<span class="facet__n">' + Number(n).toLocaleString() + '</span>' +
             '</label>';
    }).join("");

    return '<div class="facet">' +
             '<button type="button" class="facet__btn" aria-expanded="false"' +
             ' aria-haspopup="true" title="Filter ' + esc(name) + ' by value">' +
               '<span class="facet__label">All</span>' +
               '<span class="facet__caret" aria-hidden="true">&#9662;</span>' +
             '</button>' +
             '<div class="facet__menu" role="group" aria-label="' + esc(name) +
             ' values" hidden>' +
               '<div class="facet__tools">' +
                 (withFind ? '<input type="text" class="facet__find"' +
                             ' placeholder="Find value…" autocomplete="off">' : "") +
                 '<button type="button" class="facet__all" title="Tick all shown">All' +
                 '</button>' +
                 '<button type="button" class="facet__none" title="Untick all shown">None' +
                 '</button>' +
               '</div>' +
               '<div class="facet__list">' + opts + '</div>' +
               '<p class="facet__foot">counts over all rows</p>' +
             '</div>' +
           '</div>';
  }

  function tokenOf($cb) {
    return $cb.data("none") ? null : String($cb.attr("data-tok"));
  }

  // A column with nothing ticked is deliberately absent from STATE rather than holding an
  // empty Set, so the predicate's `names.length === 0` fast path stays genuinely free for
  // the untouched table — it runs per row, 167k times per draw.
  var NO_SELECTION = new Set();

  function readSel(w) {
    return STATE[w.key][w.name] || NO_SELECTION;
  }

  function ensureSel(w) {
    var sel = STATE[w.key][w.name];
    if (!sel) { sel = new Set(); STATE[w.key][w.name] = sel; }
    return sel;
  }

  function pruneSel(w) {
    var sel = STATE[w.key][w.name];
    if (sel && sel.size === 0) delete STATE[w.key][w.name];
  }

  function trunc(s) {
    return s.length > LABEL_MAX ? s.slice(0, LABEL_MAX - 1) + "…" : s;
  }

  function refreshLabel(w) {
    var sel = readSel(w), n = sel.size, text, title = "";
    if (n === 0) {
      text = "All";
    } else if (n === 1) {
      var only = sel.values().next().value;
      var lab  = (only === null) ? "(none)" : only;
      text  = trunc(lab);
      title = lab;
    } else {
      text  = n + " selected";
      title = Array.from(sel).map(function (v) { return v === null ? "(none)" : v; })
                             .join(", ");
    }
    w.$label.text(text);
    w.$btn.attr("title", title || ("Filter " + w.name + " by value"));
    w.$container.toggleClass("facet--active", n > 0);
  }

  function anyActive() {
    return TABLES.some(function (t) {
      return Object.keys(STATE[t.key]).some(function (c) {
        return STATE[t.key][c] && STATE[t.key][c].size > 0;
      });
    });
  }

  // The report's own control bar owns the reset, because with scrollY = "400px" a ticked
  // filter can be scrolled out of sight and "the table is empty and I don't know why" is
  // the predictable question. Absent in an older template, hence the guard.
  function syncClearButton() {
    var el = document.getElementById("facet-clear");
    if (el) el.hidden = !anyActive();
  }

  function clearAll() {
    WIDGETS.forEach(function (w) {
      delete STATE[w.key][w.name];
      w.$menu.find('input[type="checkbox"]').prop("checked", false);
      refreshLabel(w);
    });
    syncClearButton();
    TABLES.forEach(function (t) {
      var node = t.elem();
      if (node && $.fn.dataTable.isDataTable(node)) $(node).DataTable().draw();
    });
  }

  function build(t, api, $td, name, def) {
    var $wrap  = $td.children("div").first();      // DT's div.form-group.has-feedback
    var $input = $wrap.children("input");

    // Hidden, not removed: DT holds a reference to this input, and the scroll-head
    // un-clipping handler is bound to $wrap itself. We never write to the input, so this
    // column's DataTables search stays "" and cannot collide with the facet filter.
    $input.hide().attr("tabindex", "-1").attr("aria-hidden", "true");
    $wrap.children("span.glyphicon").hide();

    var $container = $(markup(name, def));
    $wrap.append($container);                      // inside $wrap, see note 1 at the top

    // Both of these clip: report.scss sets overflow:hidden on .dataTables_wrapper AND on
    // the htmlwidget's own .datatables div (for the card's border radius), and DT's
    // show/hide handler only un-clips the scroll head inside them.
    var $cont = $(api.table().container());

    var w = {
      key: t.key, name: name, api: api,
      $wrap: $wrap, $wrapper: $cont.add($cont.closest(".datatables")),
      $container: $container,
      $btn:   $container.find(".facet__btn"),
      $label: $container.find(".facet__label"),
      $menu:  $container.find(".facet__menu"),
      $find:  $container.find(".facet__find"),
      $list:  $container.find(".facet__list")
    };
    WIDGETS.push(w);

    w.$btn.on("click", function (e) {
      e.stopPropagation();
      openMenu(w);
    });

    // Delegated within our own node — never from `document`, which the per-draw header
    // clone would also match.
    w.$menu.on("change", 'input[type="checkbox"]', function () {
      var sel = ensureSel(w), tok = tokenOf($(this));
      if (this.checked) sel.add(tok); else sel.delete(tok);
      pruneSel(w);
      refreshLabel(w);
      syncClearButton();
      scheduleDraw(w);
    });

    w.$menu.find(".facet__all, .facet__none").on("click", function () {
      var on  = $(this).hasClass("facet__all");
      var sel = ensureSel(w);
      // Only what is currently shown, which is what the button titles say — once the find
      // box has narrowed the list, acting on the hidden values too would be a surprise.
      w.$list.children(".facet__opt").not(".facet__opt--hidden")
        .find('input[type="checkbox"]').each(function () {
          this.checked = on;
          var tok = tokenOf($(this));
          if (on) sel.add(tok); else sel.delete(tok);
        });
      pruneSel(w);
      refreshLabel(w);
      syncClearButton();
      scheduleDraw(w);
    });

    if (w.$find.length) {
      w.$find.on("input", function () {
        var q = this.value.trim().toLowerCase();
        w.$list.children(".facet__opt").each(function () {
          var txt = $(this).find(".facet__v").text().toLowerCase();
          $(this).toggleClass("facet__opt--hidden", q.length > 0 && txt.indexOf(q) < 0);
        });
      });
    }

    w.$container.on("keydown", function (e) {
      if (e.key === "Escape" || e.keyCode === 27) {
        closeOpen();
        w.$btn.trigger("focus");
      }
    });

    // A click anywhere inside the menu must not reach the outside-click handler below.
    w.$menu.on("click", function (e) { e.stopPropagation(); });

    refreshLabel(w);
    return w;
  }

  // ---- attach -------------------------------------------------------------
  function attach(t) {
    if (ATTACHED[t.key]) return;
    var node = t.elem();
    if (!node || !$.fn.dataTable.isDataTable(node)) return;
    var defs = t.facets(), cols = t.cols();
    if (!defs || !cols) return;

    var api    = $(node).DataTable();
    var $head  = $(api.table().header());          // the live thead, wherever DT moved it
    var $rows  = $head.children("tr");
    if ($rows.length < 2) return;                  // no filter row: nothing to replace
    ATTACHED[t.key] = true;

    var $label  = $rows.first().children("th,td");
    var $filter = $rows.last().children("td");
    var visIdx  = api.columns(":visible").indexes().toArray();

    if ($filter.length !== visIdx.length) {
      console.warn("facet_filter: " + t.key + " filter row has " + $filter.length +
                   " cells for " + visIdx.length +
                   " visible columns — leaving the text filters alone.");
      return;
    }

    // R drops a column it cannot offer a dropdown for — fewer than two distinct values, or
    // gone from the frame entirely — and knitr swallows its message(), so this is where that
    // decision becomes visible. console.info rather than warn: for `callers` on the VEP text
    // path and the SV `caller` column it is the expected outcome, and warnings here are
    // reserved for something actually being wrong.
    var requested = t.requested() || [];
    var missing = requested.filter(function (n) { return !(n in defs); });
    if (missing.length) {
      console.info("facet_filter: no value list for " + t.key + " column(s) " +
                   missing.join(", ") + " (fewer than " +
                   "2 distinct values, or not in the table) — they keep their text filter.");
    }

    Object.keys(defs).forEach(function (name) {
      var idx = cols[name];
      if (idx === undefined) {
        console.warn("facet_filter: '" + name + "' is not in " + t.key.toUpperCase() +
                     "_COLS — was the column renamed?");
        return;
      }
      var pos = visIdx.indexOf(idx);
      if (pos < 0) {
        console.warn("facet_filter: '" + name + "' is a hidden column.");
        return;
      }
      // Decisive cross-check on the translation: the label above the cell we picked must be
      // the column we think it is. Attaching an svtype dropdown to the locus column would
      // otherwise look like a filter that works and filters the wrong thing.
      var got = $.trim($label.eq(pos).text());
      if (got !== name) {
        console.warn("facet_filter: header mismatch for '" + name + "' at visible column " +
                     pos + " (found '" + got + "') — skipping.");
        return;
      }
      var $td = $filter.eq(pos);
      if ($td.attr("data-type") !== "character") {
        // A factor/logical column gets DT's own selectize control; two widgets in one cell
        // would fight over the same input.
        console.warn("facet_filter: '" + name + "' has data-type '" +
                     $td.attr("data-type") + "' — leaving DT's own filter in place.");
        return;
      }
      build(t, api, $td, name, defs[name]);
    });

    syncClearButton();
  }

  // DT wires its own filter-row handlers after $().DataTable() returns, i.e. after init.dt
  // has fired, and one of those is the handler that un-clips the scroll head for an open
  // menu. Defer a tick so it exists before we touch the cell. The load sweep is the
  // idempotent fallback for a table that initialised before this script was parsed.
  $(document).on("init.dt", function (e, settings) {
    var t = tableFor(settings.nTable);
    if (t) setTimeout(function () { attach(t); }, 0);
  });

  $(window).on("load", function () {
    TABLES.forEach(function (t) { if (t.elem()) attach(t); });
    var clear = document.getElementById("facet-clear");
    if (clear) clear.addEventListener("click", clearAll);
    syncClearButton();
  });

  // The one outside-click handler. Bound to document on purpose — it tests containment
  // against the open menu, not against a widget class, so the header's sizing clone is
  // irrelevant to it.
  $(document).on("mousedown", function (e) {
    if (OPEN && !OPEN.$container[0].contains(e.target)) closeOpen();
  });

  // Which values are ticked, for inspection from the browser console and for
  // tests/js/test-facet-predicate.js — the ticked state is otherwise unreachable without a
  // DOM. Same reason window.svPanelState is exposed in per_sample.qmd: one shared object
  // beats re-deriving the state from the rendered header.
  window.facetFilterState = STATE;
})();
