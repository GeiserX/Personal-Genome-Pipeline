// charts.js: the report's charts on the live page (a plugin, bin/pgp_live_charts.py).
// When state.json says the run's <sample>_charts.json exists, it is fetched from the run
// folder (same origin, nothing else) and each chart's SVG is shown with its caveat. A click on
// a mark opens a panel with the mark's detail, the chart's caveat, a link to the reads in the
// genome viewer (viewer/#<locus>) and a link to the section of the report. No inline script,
// no style attribute, no eval: the SVG is parsed as XML and holds no script.
(function () {
  'use strict';
  if (!window.PGP) { return; }
  var list = null, side = null, shownKey = null;

  function el(tag, attrs, text) {
    var e = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    if (text !== undefined && text !== null) { e.textContent = String(text); }
    return e;
  }

  function panel(box) {
    box.appendChild(el('h2', {}, 'Charts'));
    list = el('div', { 'class': 'charts' });
    list.appendChild(el('p', { 'class': 'muted' }, 'The charts appear when the report step has written them.'));
    side = el('div', { 'class': 'chart-detail', 'aria-live': 'polite' });
    side.hidden = true;
    box.appendChild(side);
    box.appendChild(list);
  }

  function openMark(sp, m) {
    side.textContent = '';
    side.hidden = false;
    side.setAttribute('data-mark', m.id);
    side.appendChild(el('h3', {}, m.label));
    if (m.detail) { side.appendChild(el('p', {}, m.detail)); }
    side.appendChild(el('p', { 'class': 'caveat' }, m.caveat || sp.caveat));
    var links = el('p', {});
    if (m.locus) { links.appendChild(el('a', { href: 'viewer/#' + m.locus, 'class': 'viewer' }, 'See the reads at ' + m.locus)); links.appendChild(document.createTextNode(' · ')); }
    links.appendChild(el('a', { href: 'report.html#' + (m.anchor || ''), 'class': 'report' }, 'Open it in the report'));
    side.appendChild(links);
    var close = el('button', { type: 'button' }, 'Close');
    close.addEventListener('click', function () { side.hidden = true; side.textContent = ''; });
    side.appendChild(close);
  }

  function chart(sp) {
    var fig = el('figure', { 'class': 'chart', 'data-chart': sp.id });
    fig.appendChild(el('h3', {}, sp.title));
    fig.appendChild(el('p', { 'class': 'caveat muted' }, sp.caveat));
    var doc = new DOMParser().parseFromString(sp.svg, 'image/svg+xml');
    if (doc.documentElement && doc.documentElement.nodeName === 'svg') {
      fig.appendChild(document.importNode(doc.documentElement, true));
    }
    var byId = {};
    (sp.marks || []).forEach(function (m) { byId[m.id] = m; });
    Array.prototype.forEach.call(fig.querySelectorAll('a[data-mark]'), function (a) {
      a.addEventListener('click', function (e) {
        e.preventDefault();
        var m = byId[a.getAttribute('data-mark')];
        if (m) { openMark(sp, m); }
      });
    });
    return fig;
  }

  function draw(data) {
    if (!data || !data.charts) { return; }
    list.textContent = '';
    data.charts.forEach(function (sp) { list.appendChild(chart(sp)); });
    document.documentElement.setAttribute('data-pgp-charts', String(data.charts.length));
  }

  function onState(state) {
    var c = state && state.plugins && state.plugins.charts;
    if (!c || !c.file || !list) { return; }
    var key = c.file + '@' + c.mtime;
    if (key === shownKey) { return; }
    shownKey = key;
    fetch(encodeURIComponent(c.file), { credentials: 'same-origin', cache: 'no-store' })
      .then(function (r) { return r.ok ? r.json() : null; })
      .then(draw)
      .catch(function () { shownKey = null; });
  }

  window.PGP.register({ name: 'charts', panel: panel, onState: onState });
}());
