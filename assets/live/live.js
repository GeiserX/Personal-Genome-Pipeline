// live.js: the live page of one run (assets/live/index.html, served by bin/pgp_live.py).
// One data source feeds the page: in live mode state.json plus the run's event stream; with
// <html data-mode="replay"> a frames.json of recorded states played on a timer. No outside
// request, no inline script, no eval. Plugins call window.PGP.register({name, onState, onEvent, panel}).
(function () {
  'use strict';
  var root = document.documentElement;
  var blocked = [];
  document.addEventListener('securitypolicyviolation', function (e) {
    blocked.push(e.violatedDirective + ' ' + (e.blockedURI || ''));
    root.setAttribute('data-pgp-blocked', blocked.join(' | '));
  });

  // One table maps each run_state state to its class and the words the legend shows.
  var STATES = {
    'waiting': ['s-waiting', 'waiting'],
    'running': ['s-running', 'running'],
    'done': ['s-done', 'done'],
    'failed': ['s-failed', 'failed'],
    'skipped': ['s-skipped', 'skipped'],
    'not finished': ['s-notfinished', 'not finished'],
    'not planned': ['s-notplanned', 'not planned']
  };
  var SVGNS = 'http://www.w3.org/2000/svg';
  var BOX_W = 170, BOX_H = 34, GAP_X = 36, GAP_Y = 10, PAD = 8;
  var plugins = [];
  var last = null, selected = null, layoutKey = null, boxes = {};

  function el(tag, attrs, text) {
    var e = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    if (text !== undefined && text !== null) { e.textContent = String(text); }
    return e;
  }
  function svg(tag, attrs) {
    var e = document.createElementNS(SVGNS, tag);
    Object.keys(attrs || {}).forEach(function (k) { e.setAttribute(k, attrs[k]); });
    return e;
  }
  function stateOf(s) { return STATES[s] || ['s-other', s]; }

  function shown(steps) {
    var v = steps.filter(function (s) { return s.state !== 'not planned'; });
    return v.length ? v : steps;
  }

  // Columns by depth: a step sits one column right of the deepest step it reads from.
  function layout(steps) {
    var by = {}, depth = {};
    steps.forEach(function (s) { by[s.step] = s; });
    function d(s, seen) {
      if (depth[s.step] !== undefined) { return depth[s.step]; }
      if (seen[s.step]) { return 0; }
      seen[s.step] = true;
      var m = -1;
      (s.upstream || []).forEach(function (u) { if (by[u]) { m = Math.max(m, d(by[u], seen)); } });
      depth[s.step] = m + 1;
      return depth[s.step];
    }
    steps.forEach(function (s) { d(s, {}); });
    var rows = {}, pos = {};
    steps.forEach(function (s) {
      var c = depth[s.step];
      rows[c] = (rows[c] || 0) + 1;
      pos[s.step] = { x: PAD + c * (BOX_W + GAP_X), y: PAD + (rows[c] - 1) * (BOX_H + GAP_Y) };
    });
    var cols = Math.max.apply(null, Object.keys(rows).map(Number)) + 1;
    var maxRows = Math.max.apply(null, Object.keys(rows).map(function (k) { return rows[k]; }));
    return { pos: pos, w: 2 * PAD + cols * BOX_W + (cols - 1) * GAP_X, h: 2 * PAD + maxRows * BOX_H + (maxRows - 1) * GAP_Y };
  }

  function drawMap(steps) {
    var key = steps.map(function (s) { return s.step; }).join(',');
    var map = document.getElementById('map');
    if (key !== layoutKey) {
      layoutKey = key;
      boxes = {};
      var L = layout(steps);
      var s = svg('svg', { viewBox: '0 0 ' + L.w + ' ' + L.h, width: L.w, height: L.h, role: 'img', 'aria-label': 'Steps of the run' });
      var edges = svg('g', { 'class': 'edges' });
      s.appendChild(edges);
      steps.forEach(function (st) {
        (st.upstream || []).forEach(function (u) {
          var a = L.pos[u], b = L.pos[st.step];
          if (!a) { return; }
          edges.appendChild(svg('line', { x1: a.x + BOX_W, y1: a.y + BOX_H / 2, x2: b.x, y2: b.y + BOX_H / 2 }));
        });
      });
      steps.forEach(function (st) {
        var p = L.pos[st.step];
        var g = svg('g', { 'class': 'box', 'data-step': st.step, tabindex: '0', role: 'button' });
        g.appendChild(svg('rect', { x: p.x, y: p.y, width: BOX_W, height: BOX_H, rx: 5 }));
        var t = svg('text', { x: p.x + 8, y: p.y + BOX_H / 2 + 4 });
        var label = st.step + ' ' + st.label;
        t.textContent = label.length > 24 ? label.slice(0, 23) + '…' : label;
        g.appendChild(t);
        var title = svg('title', {});
        g.appendChild(title);
        g.addEventListener('click', function () { selected = st.step; detail(); });
        g.addEventListener('keydown', function (e) { if (e.key === 'Enter' || e.key === ' ') { selected = st.step; detail(); } });
        s.appendChild(g);
        boxes[st.step] = g;
      });
      map.textContent = '';
      map.appendChild(s);
    }
    steps.forEach(function (st) { paint(st.step, st.state, st.label); });
  }

  function paint(step, state, label) {
    var g = boxes[step];
    if (!g) { return; }
    var cls = stateOf(state);
    g.setAttribute('class', 'box ' + cls[0] + (step === selected ? ' selected' : ''));
    g.setAttribute('data-state', state);
    var title = g.querySelector('title');
    if (title && label) { title.textContent = step + ' ' + label + ': ' + cls[1]; }
  }

  function legend() {
    var l = document.getElementById('legend');
    if (l.childNodes.length) { return; }
    Object.keys(STATES).forEach(function (k) {
      var i = el('span', { 'class': 'key ' + STATES[k][0] });
      i.appendChild(el('span', { 'class': 'swatch' }));
      i.appendChild(document.createTextNode(STATES[k][1]));
      l.appendChild(i);
    });
  }

  function detail() {
    var box = document.getElementById('detail');
    if (!last || !selected) { return; }
    var st = last.steps.filter(function (s) { return s.step === selected; })[0];
    if (!st) { return; }
    box.textContent = '';
    box.appendChild(el('h2', {}, 'Step ' + st.step + ': ' + st.label));
    var dl = el('dl', {});
    function row(k, v) { if (v === null || v === undefined || v === '') { return; } dl.appendChild(el('dt', {}, k)); dl.appendChild(el('dd', {}, v)); }
    row('State', stateOf(st.state)[1] + (st.reason ? ' (' + st.reason + ')' : ''));
    row('Note', st.note);
    row('Started', st.started_utc);
    row('Ended', st.ended_utc);
    row('Usual hours', st.hours);
    box.appendChild(dl);
    var procs = Object.keys(st.processes || {});
    if (procs.length) {
      var t = el('table', { 'class': 'tasks' });
      var h = el('tr', {});
      ['Process', 'waiting', 'running', 'done', 'cached', 'failed'].forEach(function (c) { h.appendChild(el('th', {}, c)); });
      t.appendChild(h);
      procs.forEach(function (p) {
        var r = el('tr', {});
        r.appendChild(el('td', {}, p));
        ['waiting', 'running', 'done', 'cached', 'failed'].forEach(function (c) { r.appendChild(el('td', {}, st.processes[p][c])); });
        t.appendChild(r);
      });
      box.appendChild(t);
    }
    (st.sections || []).forEach(function (k) {
      var hl = (last.headlines || {})[k];
      if (hl) { box.appendChild(el('p', { 'class': 'headline tone-' + hl.tone }, hl.text)); }
      if (root.getAttribute('data-mode') === 'live') {
        box.appendChild(el('a', { href: 'report.html#' + encodeURIComponent(k) }, 'Open in the report: ' + k));
      }
    });
    drawMap(shown(last.steps));
  }

  function found(state) {
    var ul = document.getElementById('found');
    var hl = state.headlines || {};
    var keys = [];
    state.steps.forEach(function (s) { (s.sections || []).forEach(function (k) { if (hl[k]) { keys.push(k); } }); });
    ul.textContent = '';
    if (!keys.length) { ul.appendChild(el('li', { 'class': 'muted' }, 'Nothing finished yet.')); return; }
    keys.forEach(function (k) {
      var li = el('li', { 'class': 'tone-' + hl[k].tone, 'data-section': k });
      li.appendChild(el('span', { 'class': 'step' }, hl[k].step));
      li.appendChild(document.createTextNode(' ' + hl[k].text));
      li.addEventListener('click', function () { selected = hl[k].step; detail(); });
      ul.appendChild(li);
    });
  }

  function show(state) {
    last = state;
    window.PGP.state = state;
    legend();
    var badge = document.getElementById('run-state');
    badge.textContent = state.run.state;
    badge.setAttribute('class', 'badge run-' + state.run.state.replace(/\s+/g, '-'));
    drawMap(shown(state.steps));
    found(state);
    detail();
    plugins.forEach(function (p) { if (p.onState) { try { p.onState(state); } catch (e) { console.error(p.name, e); } } });
    root.setAttribute('data-pgp-status', 'ready');
  }

  // --- live: state.json and the event stream ---
  var timer = null;
  function refresh() {
    fetch('state.json', { credentials: 'same-origin', cache: 'no-store' }).then(function (r) {
      if (r.status === 401) { window.location.href = '/login?next=' + encodeURIComponent(window.location.pathname); return null; }
      return r.ok ? r.json() : null;
    }).then(function (s) { if (s) { show(s); } }).catch(function () { root.setAttribute('data-pgp-status', 'offline'); });
  }
  function later() { if (timer) { clearTimeout(timer); } timer = setTimeout(refresh, 800); }
  function live() {
    refresh();
    var es = new EventSource('events');
    es.onmessage = function (e) {
      var m;
      try { m = JSON.parse(e.data); } catch (x) { return; }
      // Paint at once what a new message says, then let state.json settle it. A message the
      // last state.json already counted (the stream replays from the start) is not painted.
      var fresh = !last || Number(e.lastEventId) > (last.events || 0);
      if (fresh && m.step && (m.event === 'process_started' || m.event === 'process_submitted' || m.state === 'started')) { paint(m.step, 'running'); }
      plugins.forEach(function (p) { if (p.onEvent) { try { p.onEvent(m); } catch (x) { console.error(p.name, x); } } });
      later();
    };
  }

  // --- replay: frames.json on a timer ---
  function replay() {
    var bar = document.getElementById('replay');
    bar.hidden = false;
    var report = document.getElementById('report-link');
    if (report) { report.hidden = true; }
    fetch('frames.json').then(function (r) { return r.json(); }).then(function (data) {
      var frames = data.frames || [];
      var slider = document.getElementById('replay-slider'), play = document.getElementById('replay-play');
      var at = 0, tick = null;
      slider.max = String(Math.max(frames.length - 1, 0));
      function go(i) {
        at = Math.min(Math.max(i, 0), frames.length - 1);
        slider.value = String(at);
        document.getElementById('replay-time').textContent = frames[at].t || '';
        show(frames[at].state);
        root.setAttribute('data-pgp-frame', String(at));
      }
      function stop() { if (tick) { clearInterval(tick); tick = null; } play.textContent = 'Play'; }
      play.addEventListener('click', function () {
        if (tick) { stop(); return; }
        if (at >= frames.length - 1) { go(0); }
        play.textContent = 'Pause';
        tick = setInterval(function () { if (at >= frames.length - 1) { stop(); } else { go(at + 1); } }, Number(data.interval_ms) || 1000);
      });
      slider.addEventListener('input', function () { stop(); go(Number(slider.value)); });
      if (frames.length) { go(0); }
      if (data.autoplay) { play.click(); }
    });
  }

  window.PGP = {
    state: null,
    register: function (p) {
      plugins.push(p);
      if (p.panel) {
        var box = el('section', { 'class': 'card plugin', 'data-plugin': p.name || '' });
        document.getElementById('plugin-panels').appendChild(box);
        try { p.panel(box); } catch (e) { console.error(p.name, e); }
      }
      if (last && p.onState) { p.onState(last); }
    }
  };

  function startPage() {
    if (root.getAttribute('data-mode') === 'replay') { replay(); } else { live(); }
  }
  if (document.readyState === 'loading') { document.addEventListener('DOMContentLoaded', startPage); } else { startPage(); }
}());
