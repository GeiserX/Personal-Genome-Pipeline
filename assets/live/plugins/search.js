// search.js: find a gene symbol or rsID across one run's results (a plugin, bin/pgp_live_search.py).
// When state.json says the run's <sample>_search_index.json exists, it is fetched once from the
// run folder (same origin, nothing else; again only when the report is written anew). Typing
// filters that index here in the browser: no request per keystroke, so the query is in no
// request and no server log. A hit lists each section that names the term, with a link to it in
// the report and, for a mention with a locus, to the reads in the genome viewer. The query rule
// is bin/results_index.py's query(): an rsID, a gene symbol (a prefix also lists the terms that
// start with it), or a position, which the index cannot answer (scripts/pgp-lookup.sh can).
// tests/test_results_index.py runs match() under node against the Python rule.
(function (root) {
  'use strict';
  var SUGGESTIONS = 20;

  function match(index, q) {
    q = String(q || '').trim();
    var terms = (index && index.terms) || {};
    if (/^rs\d+$/i.test(q)) {
      var t = q.toLowerCase();
      return { kind: 'rsid', term: t, items: terms[t] || [], suggestions: [] };
    }
    var site = /^(?:chr)?([0-9]{1,2}|X|Y|M|MT)[:_-]([0-9][0-9,]*)$/i.exec(q);
    if (site) {
      // Position 0 is no position (genome_lookup.parse_query refuses it too).
      var kind = parseInt(site[2].replace(/,/g, ''), 10) < 1 ? 'error' : 'site';
      return { kind: kind, term: null, items: [], suggestions: [] };
    }
    if (/^[A-Za-z0-9][A-Za-z0-9._-]*$/.test(q)) {
      var g = q.toUpperCase();
      var sugg = Object.keys(terms).sort().filter(function (k) { return k !== g && k.indexOf(g) === 0; });
      return { kind: 'gene', term: g, items: terms[g] || [], suggestions: sugg.slice(0, SUGGESTIONS) };
    }
    return { kind: 'error', term: null, items: [], suggestions: [] };
  }

  if (typeof module !== 'undefined' && module.exports) { module.exports = { match: match }; }
  if (!root || !root.PGP || typeof document === 'undefined') { return; }

  var index = null, shownKey = null, input = null, out = null, status = null;

  function el(tag, attrs, text) {
    var e = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    if (text !== undefined && text !== null) { e.textContent = String(text); }
    return e;
  }

  function panel(box) {
    box.appendChild(el('h2', {}, 'Search the results'));
    var form = el('form', { role: 'search' });
    form.addEventListener('submit', function (e) { e.preventDefault(); show(); });
    input = el('input', { type: 'search', 'aria-label': 'Gene symbol or rsID', placeholder: 'BRCA2 or rs429358',
                          autocomplete: 'off', spellcheck: 'false' });
    input.disabled = true;
    input.addEventListener('input', show);
    form.appendChild(input);
    box.appendChild(form);
    status = el('p', { 'class': 'muted' }, 'The search opens when the report step has written its index.');
    box.appendChild(status);
    out = el('div', { 'class': 'search-results', 'aria-live': 'polite' });
    box.appendChild(out);
  }

  function item(it) {
    var li = el('li', { 'data-section': it.section });
    var head = el('p', {});
    head.appendChild(el('a', { href: 'report.html#' + (it.anchor || ''), 'class': 'report' }, it.title));
    head.appendChild(document.createTextNode(' (step ' + it.step + ')' + (it.state === 'stale' ? ', stale: from an earlier run' : '')
      + ': ' + it.count + (it.count === 1 ? ' mention' : ' mentions')));
    li.appendChild(head);
    var ul = el('ul', {});
    (it.mentions || []).forEach(function (m) {
      var x = el('li', {}, m.text);
      if (m.locus) {
        x.appendChild(document.createTextNode(' · '));
        x.appendChild(el('a', { href: 'viewer/#' + m.locus, 'class': 'viewer' }, 'reads at ' + m.locus));
      }
      ul.appendChild(x);
    });
    if (it.count > (it.mentions || []).length) {
      ul.appendChild(el('li', { 'class': 'muted' }, '... and ' + (it.count - it.mentions.length) + ' more in ' + it.source));
    }
    li.appendChild(ul);
    return li;
  }

  function show() {
    if (!out) { return; }
    out.textContent = '';
    var q = input.value;
    if (!index || !q.trim()) { document.documentElement.removeAttribute('data-pgp-search'); return; }
    var r = match(index, q);
    document.documentElement.setAttribute('data-pgp-search', String(r.items.length));
    if (r.kind === 'site') {
      out.appendChild(el('p', {}, 'The search finds gene symbols and rsIDs. A position is a question about the genotype: '
        + 'run scripts/pgp-lookup.sh on this computer.'));
      return;
    }
    if (r.kind === 'error') {
      out.appendChild(el('p', {}, 'Type a gene symbol (CFTR) or an rsID (rs123).'));
      return;
    }
    if (r.items.length) {
      out.appendChild(el('h3', {}, r.term));
      var ul = el('ul', { 'class': 'search-hits' });
      r.items.forEach(function (it) { ul.appendChild(item(it)); });
      out.appendChild(ul);
    } else {
      out.appendChild(el('p', {}, 'None of these results mentions ' + r.term + '.'));
    }
    if (r.suggestions.length) {
      var p = el('p', { 'class': 'suggestions' }, 'Starting with ' + r.term + ': ');
      r.suggestions.forEach(function (t, i) {
        if (i) { p.appendChild(document.createTextNode(', ')); }
        var b = el('button', { type: 'button' }, t);
        b.addEventListener('click', function () { input.value = t; show(); });
        p.appendChild(b);
      });
      out.appendChild(p);
    }
  }

  function onState(state) {
    var c = state && state.plugins && state.plugins.search;
    if (!out || !c) { return; }
    if (!c.file) {
      // The file went away (step 24 deletes it before it writes the report again).
      if (shownKey !== null) { shownKey = null; index = null; input.disabled = true; show(); }
      return;
    }
    var key = c.file + '@' + c.mtime;
    if (key === shownKey) { return; }
    shownKey = key;
    fetch(encodeURIComponent(c.file), { credentials: 'same-origin', cache: 'no-store' })
      // A refused or missing file (an expired session, a file mid-rewrite) clears the key, so the next state retries.
      .then(function (r) {
        if (key !== shownKey) { return null; }
        if (!r.ok) { shownKey = null; return null; }
        return r.json();
      })
      .then(function (data) {
        if (key !== shownKey || !data || !data.terms) { return; }
        index = data;
        input.disabled = false;
        status.textContent = Object.keys(data.terms).length + ' genes and rsIDs in these results. Typing searches them here, '
          + 'in this browser.';
        show();
      })
      .catch(function () { if (key === shownKey) { shownKey = null; } });
  }

  root.PGP.register({ name: 'search', panel: panel, onState: onState });
}(typeof window !== 'undefined' ? window : null));
