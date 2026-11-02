// recheck.js: the database re-check on the live page (a plugin, bin/pgp_live_recheck.py).
// state.json's plugins.recheck holds a banner line and, per step, its badges ("newer ClinVar
// available", "redone 2026-10-05"). The panel shows the banner and one line per marked step; on
// the step map a marked step gets a small dot with the badges as its tooltip. No inline script,
// no style attribute, no eval; nothing is fetched beyond state.json.
(function () {
  'use strict';
  if (!window.PGP) { return; }
  var NS = 'http://www.w3.org/2000/svg';
  var banner = null, list = null, shownKey = null;

  function el(tag, attrs, text) {
    var e = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    if (text !== undefined && text !== null) { e.textContent = String(text); }
    return e;
  }

  function panel(box) {
    box.appendChild(el('h2', {}, 'Database re-check'));
    banner = el('p', { 'class': 'muted' }, 'No re-check yet: bin/pgp_recheck.py writes one (docs/recheck.md).');
    list = el('ul', { 'class': 'plain recheck-steps' });
    box.appendChild(banner);
    box.appendChild(list);
  }

  function mark(step, badges) {
    var g = document.querySelector('#map g.box[data-step="' + step + '"]');
    if (!g) { return; }
    var old = g.querySelector('.recheck-mark');
    if (old) { g.removeChild(old); }
    if (!badges || !badges.length) { return; }
    var r = g.querySelector('rect');
    var x = Number(r.getAttribute('x')) + Number(r.getAttribute('width')) - 7;
    var y = Number(r.getAttribute('y')) + 7;
    var c = document.createElementNS(NS, 'circle');
    c.setAttribute('class', 'recheck-mark');
    c.setAttribute('cx', String(x));
    c.setAttribute('cy', String(y));
    c.setAttribute('r', '4');
    c.setAttribute('fill', '#e0b000');
    var t = document.createElementNS(NS, 'title');
    t.textContent = badges.join('; ');
    c.appendChild(t);
    g.appendChild(c);
  }

  function onState(state) {
    var r = state && state.plugins && state.plugins.recheck;
    if (!banner) { return; }
    var steps = (r && r.steps) || {};
    // The map can be drawn anew: mark every box on each state, not only when the check changed.
    Array.prototype.forEach.call(document.querySelectorAll('#map g.box'), function (g) {
      mark(g.getAttribute('data-step'), steps[g.getAttribute('data-step')]);
    });
    var key = r && r.file ? r.file + '@' + r.mtime : null;
    if (key === shownKey) { return; }
    shownKey = key;
    list.textContent = '';
    if (!key) {
      banner.textContent = 'No re-check yet: bin/pgp_recheck.py writes one (docs/recheck.md).';
      document.documentElement.removeAttribute('data-pgp-recheck');
      return;
    }
    banner.textContent = r.banner;
    Object.keys(steps).sort().forEach(function (n) { list.appendChild(el('li', {}, 'Step ' + n + ': ' + steps[n].join('; '))); });
    document.documentElement.setAttribute('data-pgp-recheck', String(Object.keys(steps).length));
  }

  window.PGP.register({ name: 'recheck', panel: panel, onState: onState });
}());
