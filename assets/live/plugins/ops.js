// ops.js: the "Run health" card of the live page, drawn from state.json's plugins.ops
// (bin/pgp_live_ops.py): the alerts, the hours left as a range, and the memory and disk headroom.
// No outside request, no inline script, text only (never innerHTML).
(function () {
  'use strict';
  var box = null, last = null;

  function el(tag, cls, text) {
    var e = document.createElement(tag);
    if (cls) { e.className = cls; }
    if (text !== undefined && text !== null) { e.textContent = String(text); }
    return e;
  }
  function hours(s) {
    if (s === null || s === undefined) { return null; }
    var h = s / 3600;
    return h < 1 ? Math.round(s / 60) + ' min' : (h < 10 ? h.toFixed(1) : Math.round(h)) + ' h';
  }
  function pct(r) { return r === null || r === undefined ? 'unknown' : Math.round(r * 100) + '%'; }
  function size(b) {
    if (b === null || b === undefined) { return 'unknown'; }
    var u = ['B', 'KB', 'MB', 'GB', 'TB'], i = 0;
    while (b >= 1024 && i < u.length - 1) { b /= 1024; i++; }
    return (i ? b.toFixed(1) : Math.round(b)) + ' ' + u[i];
  }

  function draw() {
    if (!box) { return; }
    while (box.firstChild) { box.removeChild(box.firstChild); }
    box.appendChild(el('h2', null, 'Run health'));
    var o = last && last.plugins && last.plugins.ops;
    if (!o) { box.appendChild(el('p', 'muted', 'No numbers yet.')); return; }
    if (o.error) { box.appendChild(el('p', 'tone-attention', o.error)); return; }

    var ul = el('ul', 'found ops-alerts');
    (o.alerts || []).forEach(function (a) {
      var li = el('li', a.level === 'info' ? 'tone-info' : 'tone-attention', a.text);
      li.setAttribute('data-kind', a.kind);
      li.setAttribute('data-level', a.level);
      ul.appendChild(li);
    });
    if (ul.firstChild) { box.appendChild(ul); } else { box.appendChild(el('p', 'muted', 'No alert.')); }

    var h = o.hours_left;
    if (h) {
      var lo = hours(h.low_s), hi = hours(h.high_s);
      var t = hi === null ? 'at least ' + lo + ' left; the high end is unknown while a step runs late'
        : (lo === hi ? lo : lo + ' to ' + hi) + ' left';
      var p = el('p', 'headline ops-left', t);
      p.title = h.basis || '';
      box.appendChild(p);
      if (h.unknown_steps && h.unknown_steps.length) {
        box.appendChild(el('p', 'muted', 'Not counted (no usual time): step ' + h.unknown_steps.join(', ')));
      }
    }

    (o.disk || []).forEach(function (d) {
      var s = 'Disk ' + d.fs + ': ' + pct(d.used_ratio) + ' used, ' + size(d.free_bytes) + ' free';
      if (d.inodes_used_ratio !== null && d.inodes_used_ratio !== undefined) { s += ', inodes ' + pct(d.inodes_used_ratio); }
      box.appendChild(el('p', d.warn ? 'tone-attention ops-disk' : 'muted ops-disk', s));
    });
    var m = o.memory || {};
    var ms = 'Memory: ' + size(m.host_bytes) + ' on this machine';
    if (m.max_memory_bytes !== null && m.max_memory_bytes !== undefined) { ms += ', max_memory ' + size(m.max_memory_bytes); }
    if (m.retries) { ms += ', ' + m.retries + ' retried after a memory or time kill'; }
    box.appendChild(el('p', m.max_memory_above_ram || m.retries ? 'tone-attention ops-mem' : 'muted ops-mem', ms));
    (m.tasks || []).forEach(function (t) {
      box.appendChild(el('p', 'muted', (t.step ? 'Step ' + t.step + ', ' : '') + t.process + ': peak_rss ' +
        size(t.peak_rss_bytes) + ' of ' + size(t.request_bytes) + ' asked (' + pct(t.ratio) + ')'));
    });
  }

  window.PGP.register({
    name: 'ops',
    panel: function (el_) { box = el_; draw(); },
    onState: function (state) { last = state; draw(); }
  });
}());
