// family.js: the family card of the live page (bin/pgp_live_family.py fills plugins.family).
// On a family's page: the pedigree, each child's chromosome map, the relationship check, the
// shared segments and the couple carrier screen, as steps 38 and 39 wrote them. On a sample's
// page: a link to its family. Text only through textContent; images from this server only.
(function () {
  'use strict';
  var box = null, shownKey = null;

  function el(tag, attrs, text) {
    var e = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    if (text !== undefined && text !== null) { e.textContent = String(text); }
    return e;
  }
  function link(url, text) { return el('a', { href: url }, text); }
  function para(text, cls) { return el('p', cls ? { 'class': cls } : {}, text); }
  function table(head, rows) {
    var t = el('table', { 'class': 'family-table' }), tr = el('tr');
    head.forEach(function (h) { tr.appendChild(el('th', {}, h)); });
    t.appendChild(tr);
    rows.forEach(function (r) {
      var row = el('tr');
      r.forEach(function (c) { row.appendChild(el('td', {}, c === null || c === undefined ? '' : c)); });
      t.appendChild(row);
    });
    return t;
  }

  function family(f) {
    box.appendChild(el('h2', {}, 'Family'));
    if (f.error) { box.appendChild(para(f.error)); return; }
    var img = el('img', { src: f.pedigree, alt: 'Pedigree of ' + f.name, 'class': 'pedigree' });
    img.style.maxWidth = '100%';
    box.appendChild(img);

    box.appendChild(el('h3', {}, 'Chromosome maps'));
    if (!f.children.length) { box.appendChild(para('The pedigree names no child.', 'muted')); }
    f.children.forEach(function (c) {
      var d = el('details', { 'data-child': c.child });
      d.appendChild(el('summary', {}, c.child + ': father ' + (c.father || 'not named') + ', mother ' + (c.mother || 'not named')));
      if (c.segments) {
        var m = el('img', { src: c.segments, alt: 'Chromosome map of ' + c.child });
        m.style.maxWidth = '100%';
        d.appendChild(m);
      } else {
        d.appendChild(para('No chromosome map yet: step 38 draws one for a child with a sequenced parent.', 'muted'));
      }
      box.appendChild(d);
    });

    box.appendChild(el('h3', {}, 'Relationship check'));
    if (f.relationships.length) {
      box.appendChild(table(['Pair', 'Declared', 'Verdict', 'Reason'], f.relationships.map(function (r) {
        return [r.a + ' / ' + r.b, r.declared, r.verdict, r.reason];
      })));
    } else { box.appendChild(para('Not run yet (step 38).', 'muted')); }

    box.appendChild(el('h3', {}, 'Shared segments'));
    if (f.ibd_pairs.length) {
      box.appendChild(table(['Pair', 'Declared', 'Shared cM', 'Degree'], f.ibd_pairs.map(function (r) {
        return [r.a + ' / ' + r.b, r.declared, r.shared_cm, r.degree];
      })));
    } else { box.appendChild(para('Not run yet (step 38 with the panel of step 26b).', 'muted')); }

    box.appendChild(el('h3', {}, 'Couple carrier screen'));
    if (!f.couples.length) { box.appendChild(para('The pedigree names no partners.', 'muted')); }
    f.couples.forEach(function (c) {
      var p = el('p', { 'data-couple': c.a + ' ' + c.b });
      if (!c.screened) {
        p.textContent = c.a + ' and ' + c.b + ': not screened yet (step 39).';
      } else {
        // A count of zero is not an all-clear: the screen sees only what it lists as assessed.
        p.textContent = c.a + ' and ' + c.b + ': ' + c.shared + ' gene(s) where both carry a finding, ' + c.x_linked +
          ' X-linked. Read the report for what the screen does not see. ';
        if (c.report) { p.appendChild(link(c.report, 'Report')); }
      }
      box.appendChild(p);
    });

    if (f.plan.length) {
      box.appendChild(el('h3', {}, 'Plan'));
      box.appendChild(table(['Analysis', 'Decision', 'Reason'], f.plan.map(function (r) {
        return [r.analysis, r.decision, r.reason];
      })));
    }
    if (f.reports.length) {
      var ul = el('ul');
      f.reports.forEach(function (r) { var li = el('li'); li.appendChild(link(r.url, r.name)); ul.appendChild(li); });
      box.appendChild(el('h3', {}, 'Reports'));
      box.appendChild(ul);
    }
  }

  function member(f) {
    box.appendChild(el('h2', {}, 'Family'));
    f.families.forEach(function (x) {
      var p = el('p', {}, 'In the family ');
      p.appendChild(link(x.url, x.name));
      box.appendChild(p);
    });
  }

  function onState(state) {
    if (!box) { return; }
    var f = state && state.plugins ? state.plugins.family : null;
    var key = JSON.stringify(f || null);
    if (key === shownKey) { return; }  // state.json comes after every message; redraw only on change
    shownKey = key;
    while (box.firstChild) { box.removeChild(box.firstChild); }
    box.hidden = !f;
    if (!f) { return; }
    if (f.kind === 'family') { family(f); } else if (f.kind === 'member') { member(f); }
  }

  window.PGP.register({ name: 'family', panel: function (b) { box = b; box.hidden = true; }, onState: onState });
}());
