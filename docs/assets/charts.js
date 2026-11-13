/* Phase Equilibrium Calculator: small SVG line-chart library shared by all pages.
   PEC.drawChart(container, cfg) draws responsive line charts with a hover
   read-out; PEC.legend and PEC.dataTable build the matching legend and table. */
(function () {
  "use strict";
  var NS = "http://www.w3.org/2000/svg";
  var PEC = (window.PEC = window.PEC || {});

  PEC.fmt = function (v, d) {
    if (v == null || !isFinite(v)) return "–";
    if (Math.abs(v) < 0.5 * Math.pow(10, -d)) v = 0;
    return Number(v).toLocaleString("en-US", { minimumFractionDigits: d, maximumFractionDigits: d });
  };
  /* significant-figure formatting for read-outs */
  PEC.sig = function (v, n) {
    if (v == null || !isFinite(v)) return "–";
    if (v === 0) return "0";
    var e = Math.floor(Math.log10(Math.abs(v)));
    var d = Math.max(0, Math.min(8, (n || 4) - 1 - e));
    return PEC.fmt(v, d);
  };
  PEC.sci = function (v, d) {
    if (v === 0) return "0";
    var e = Math.floor(Math.log10(Math.abs(v))), m = v / Math.pow(10, e);
    return PEC.fmt(m, d) + " × 10<sup>" + String(e).replace("-", "−") + "</sup>";
  };

  function el(tag, attrs, parent) {
    var e = document.createElementNS(NS, tag);
    for (var a in attrs) e.setAttribute(a, attrs[a]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function niceTicks(lo, hi, n) {
    var span = hi - lo;
    if (!(span > 0)) { span = Math.abs(lo) || 1; lo -= span / 2; hi += span / 2; span = hi - lo; }
    var step = Math.pow(10, Math.floor(Math.log10(span / n)));
    var err = span / n / step;
    if (err >= 7.5) step *= 10; else if (err >= 3.5) step *= 5; else if (err >= 1.5) step *= 2;
    var t = [], v = Math.ceil(lo / step - 1e-9) * step;
    for (; v <= hi + step * 1e-9; v += step) t.push(Math.round(v / step) * step);
    return { ticks: t, step: step };
  }
  function decimals(step) { return Math.max(0, -Math.floor(Math.log10(step) + 1e-9)); }

  /* a domain that contains the data with a little headroom, on round numbers */
  PEC.niceDomain = function (vals, opts) {
    opts = opts || {};
    var lo = Infinity, hi = -Infinity;
    vals.forEach(function (v) { if (v != null && isFinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; } });
    if (!isFinite(lo)) return [0, 1];
    if (opts.include != null) { lo = Math.min(lo, opts.include); hi = Math.max(hi, opts.include); }
    var span = hi - lo;
    if (span < 1e-9 * Math.max(1, Math.abs(hi))) { span = Math.abs(hi) * 0.1 || 1; }
    var pad = (opts.pad == null ? 0.06 : opts.pad) * span;
    var a = lo - pad, b = hi + pad;
    if (opts.floor != null) a = Math.max(a, opts.floor);
    if (lo >= 0 && a < 0) a = 0;
    var t = niceTicks(a, b, 10);        // round outward to a fine step: little wasted space
    return [Math.floor(a / t.step + 1e-9) * t.step, Math.ceil(b / t.step - 1e-9) * t.step];
  };

  PEC.legend = function (target, items) {
    var ul = typeof target === "string" ? document.getElementById(target) : target;
    ul.innerHTML = items.map(function (it) {
      var cls = "key" + (it.area ? " area" : "") + (it.dot ? " dot" : "");
      return "<li><span class='" + cls + "' style='background:var(" + it.color + ")'></span>" + it.name + "</li>";
    }).join("");
  };

  PEC.dataTable = function (target, head, rows) {
    var host = typeof target === "string" ? document.getElementById(target) : target;
    var h = "<table><thead><tr>" + head.map(function (c, i) {
      return "<th" + (i ? " class='num'" : "") + ">" + c + "</th>"; }).join("") + "</tr></thead><tbody>";
    rows.forEach(function (r) {
      h += "<tr>" + r.map(function (c, i) { return "<td" + (i ? " class='num'" : "") + ">" + c + "</td>"; }).join("") + "</tr>";
    });
    host.innerHTML = h + "</tbody></table>";
  };

  var tip;
  function tooltip() {
    if (!tip) {
      tip = document.getElementById("tooltip");
      if (!tip) { tip = document.createElement("div"); tip.id = "tooltip"; tip.className = "tooltip"; tip.setAttribute("role", "presentation"); document.body.appendChild(tip); }
    }
    return tip;
  }

  /* cfg: {
       series: [{ name, color, pts: [[x, y, payload?] | null], width, opacity, label: {i, text, dx, dy, anchor, minW} }],
       areas:  [{ color, poly: [[x, y]], opacity }],
       xDomain, yDomain (optional: computed from the data), xTitle, yTitle, xDec, yDec,
       ratio, height, hlines: [{y, text}], vlines: [{x, text}], diagonal: true,
       regions: [{x, y, text, anchor, minW}], points: [{x, y, color, text, dx, dy, anchor, minW}],
       tip: function (series, point) -> html,  aria }                                   */
  PEC.drawChart = function (container, cfg) {
    container.__surface = null;                 // leave 3D mode, if it was on
    var T = tooltip();
    function render() {
      var W = Math.max(280, container.clientWidth);
      var H = cfg.height || Math.round(Math.min(cfg.maxHeight || 400, Math.max(240, W * (cfg.ratio || 0.62))));
      var m = { t: 12, r: 18, b: 44, l: 56 };
      var iw = W - m.l - m.r, ih = H - m.t - m.b;
      container.innerHTML = "";
      var svg = el("svg", { viewBox: "0 0 " + W + " " + H, width: W, height: H, role: "img", "aria-label": cfg.aria || "" }, container);

      var xs = [], ys = [];
      cfg.series.forEach(function (s) { s.pts.forEach(function (p) { if (p) { xs.push(p[0]); ys.push(p[1]); } }); });
      var xD = cfg.xDomain || PEC.niceDomain(xs, { pad: 0.02 });
      var yD = cfg.yDomain || PEC.niceDomain(ys, cfg.yOpts);
      var x0 = xD[0], x1 = xD[1], y0 = yD[0], y1 = yD[1];
      function X(v) { return m.l + (v - x0) / (x1 - x0) * iw; }
      function Y(v) { return m.t + ih - (v - y0) / (y1 - y0) * ih; }

      var g = el("g", { "class": "tick" }, svg);
      var xt = niceTicks(x0, x1, Math.max(3, Math.floor(iw / 80))), yt = niceTicks(y0, y1, Math.max(3, Math.floor(ih / 48)));
      var xd = cfg.xDec != null ? cfg.xDec : decimals(xt.step), yd = cfg.yDec != null ? cfg.yDec : decimals(yt.step);
      yt.ticks.forEach(function (v) {
        el("line", { x1: m.l, x2: m.l + iw, y1: Y(v), y2: Y(v), "class": "gridline" }, g);
        el("text", { x: m.l - 8, y: Y(v) + 4, "text-anchor": "end" }, g).textContent = PEC.fmt(v, yd);
      });
      xt.ticks.forEach(function (v) {
        el("line", { x1: X(v), x2: X(v), y1: m.t + ih, y2: m.t + ih + 4, "class": "baseline" }, g);
        el("text", { x: X(v), y: m.t + ih + 18, "text-anchor": "middle" }, g).textContent = PEC.fmt(v, xd);
      });
      el("line", { x1: m.l, x2: m.l + iw, y1: m.t + ih, y2: m.t + ih, "class": "baseline" }, svg);
      el("line", { x1: m.l, x2: m.l, y1: m.t, y2: m.t + ih, "class": "baseline" }, svg);
      el("text", { x: m.l + iw / 2, y: H - 6, "text-anchor": "middle", "class": "axis-title" }, svg).textContent = cfg.xTitle || "";
      el("text", { x: 14, y: m.t + ih / 2, "text-anchor": "middle", "class": "axis-title",
        transform: "rotate(-90 14 " + (m.t + ih / 2) + ")" }, svg).textContent = cfg.yTitle || "";

      var clipId = "c" + Math.random().toString(36).slice(2);
      var cp = el("clipPath", { id: clipId }, svg);
      el("rect", { x: m.l, y: m.t - 4, width: iw, height: ih + 5 }, cp);
      var plot = el("g", { "clip-path": "url(#" + clipId + ")" }, svg);

      (cfg.areas || []).forEach(function (a) {
        if (a.poly.length < 3) return;
        el("path", { d: "M" + a.poly.map(function (p) { return X(p[0]).toFixed(1) + "," + Y(p[1]).toFixed(1); }).join("L") + "Z",
          fill: "var(" + a.color + ")", "fill-opacity": a.opacity || 0.1, stroke: "none" }, plot);
      });
      if (cfg.diagonal) el("line", { x1: X(x0), y1: Y(x0), x2: X(x1), y2: Y(x1), "class": "annot-line" }, plot);
      (cfg.hlines || []).forEach(function (h) {
        if (h.y < y0 || h.y > y1) return;
        el("line", { x1: m.l, x2: m.l + iw, y1: Y(h.y), y2: Y(h.y), "class": "annot-line" }, svg);
        if (h.text) el("text", { x: m.l + iw - 4, y: Y(h.y) - 6, "text-anchor": "end", "class": "annot-text" }, svg).textContent = h.text;
      });
      (cfg.vlines || []).forEach(function (v) {
        if (v.x < x0 || v.x > x1) return;
        el("line", { x1: X(v.x), x2: X(v.x), y1: m.t, y2: m.t + ih, "class": "annot-line" }, svg);
        if (v.text) el("text", { x: X(v.x) + 5, y: m.t + 12, "class": "annot-text" }, svg).textContent = v.text;
      });
      cfg.series.forEach(function (s) {
        var segs = [[]];
        s.pts.forEach(function (p) { if (p == null) segs.push([]); else segs[segs.length - 1].push(p); });
        segs.forEach(function (sg) {
          if (sg.length < 2) return;
          el("path", { d: "M" + sg.map(function (p) { return X(p[0]).toFixed(1) + "," + Y(p[1]).toFixed(1); }).join("L"),
            "class": "series-line", stroke: "var(" + s.color + ")", "stroke-width": s.width || 2, opacity: s.opacity || 1,
            "stroke-dasharray": s.dash || "none" }, plot);
        });
      });
      (cfg.regions || []).forEach(function (r) {
        if (r.minW && iw < r.minW) return;
        el("text", { x: X(r.x), y: Y(r.y), "text-anchor": r.anchor || "middle", "class": "region-label" }, svg).textContent = r.text;
      });
      cfg.series.forEach(function (s) {
        if (!s.label || (s.label.minW && iw < s.label.minW)) return;
        var p = s.pts[s.label.i < 0 ? s.pts.length + s.label.i : s.label.i];
        if (!p) return;
        el("text", { x: X(p[0]) + (s.label.dx || 0), y: Y(p[1]) + (s.label.dy || 0),
          "text-anchor": s.label.anchor || "start", "class": "direct-label" }, svg).textContent = s.label.text;
      });
      (cfg.points || []).forEach(function (p) {
        if (p.x < x0 || p.x > x1 || p.y < y0 || p.y > y1) return;
        el("circle", { cx: X(p.x), cy: Y(p.y), r: 4.5, fill: "var(" + (p.color || "--ink-2") + ")", "class": "hover-dot" }, svg);
        if (p.text && !(p.minW && iw < p.minW)) {
          el("text", { x: X(p.x) + (p.dx || 0), y: Y(p.y) + (p.dy || 0), "text-anchor": p.anchor || "start", "class": "annot-text" }, svg).textContent = p.text;
        }
      });

      /* hover layer: nearest point of any series */
      var cross = el("line", { "class": "crosshair", y1: m.t, y2: m.t + ih, opacity: 0 }, svg);
      var dot = el("circle", { r: 5, "class": "hover-dot", opacity: 0 }, svg);
      var hit = el("rect", { x: m.l, y: m.t, width: iw, height: ih, fill: "transparent" }, svg);
      function hide() { T.classList.remove("on"); dot.setAttribute("opacity", 0); cross.setAttribute("opacity", 0); }
      function move(ev) {
        var b = svg.getBoundingClientRect(), sx = W / b.width;
        var px = (ev.clientX - b.left) * sx, py = (ev.clientY - b.top) * sx;
        var best = null, bd = Infinity;
        cfg.series.forEach(function (s) {
          if (s.noHover) return;
          s.pts.forEach(function (p) {
            if (!p) return;
            var dx = X(p[0]) - px, dy = Y(p[1]) - py, d = dx * dx + dy * dy;
            if (d < bd) { bd = d; best = { s: s, p: p }; }
          });
        });
        if (!best || bd > 70 * 70) { hide(); return; }
        var cx = X(best.p[0]), cy = Y(best.p[1]);
        dot.setAttribute("cx", cx); dot.setAttribute("cy", cy);
        dot.setAttribute("fill", "var(" + best.s.color + ")"); dot.setAttribute("opacity", 1);
        cross.setAttribute("x1", cx); cross.setAttribute("x2", cx); cross.setAttribute("opacity", 1);
        var body = cfg.tip ? cfg.tip(best.s, best.p)
          : (cfg.fx ? cfg.fx(best.p[0]) : PEC.sig(best.p[0], 4)) + "<br>" + (cfg.fy ? cfg.fy(best.p[1]) : PEC.sig(best.p[1], 4));
        T.innerHTML = "<span class='sw' style='background:var(" + best.s.color + ")'></span>" + best.s.name + "<br>" + body;
        var host = container.getBoundingClientRect();
        var left = window.scrollX + host.left + cx / sx, half = T.offsetWidth / 2;
        left = Math.max(window.scrollX + half + 4, Math.min(left, window.scrollX + document.documentElement.clientWidth - half - 4));
        T.style.left = left + "px";
        T.style.top = (window.scrollY + host.top + cy / sx) + "px";
        T.classList.add("on");
      }
      hit.addEventListener("pointermove", move);
      hit.addEventListener("pointerdown", move);
      hit.addEventListener("pointerleave", hide);
    }
    render();
    if (container.__ro) container.__ro.disconnect();
    var last = container.clientWidth;
    container.__ro = new ResizeObserver(function () {
      if (Math.abs(container.clientWidth - last) > 4) { last = container.clientWidth; render(); }
    });
    container.__ro.observe(container);
  };

  /* 3D surface, drawn in SVG with the painter's algorithm and rotated by dragging.
     cfg = { surfaces: [{ name, color, grid: [[ [x, y, z] | null ]] }]   (grid rows: layers in y)
             lines:    [{ color, width, pts: [[x, y, z] | null] }]       (e.g. the current slice)
             x, y, z:  { title, domain, log }                            (z is vertical)
             aria }                                                                     */
  PEC.drawSurface = function (container, cfg) {
    var view = container.__view || (container.__view = { yaw: -0.85, pitch: 0.38 });
    var css = getComputedStyle(document.documentElement);
    function rgb(name) {
      var h = css.getPropertyValue(name).trim();
      if (h.length === 4) h = "#" + h[1] + h[1] + h[2] + h[2] + h[3] + h[3];
      var n = parseInt(h.slice(1), 16);
      return isFinite(n) ? [(n >> 16) & 255, (n >> 8) & 255, n & 255] : [120, 120, 120];
    }
    function lg(v) { return Math.log(v) / Math.LN10; }
    function nrm(v, ax) {
      var d = ax.domain;
      if (ax.log) return 2 * (lg(v) - lg(d[0])) / (lg(d[1]) - lg(d[0])) - 1;
      return 2 * (v - d[0]) / (d[1] - d[0]) - 1;
    }
    function N(p) { return p ? [nrm(p[0], cfg.x), nrm(p[1], cfg.y), nrm(p[2], cfg.z)] : null; }
    var surf = cfg.surfaces.map(function (s) {
      return { s: s, base: rgb(s.color), g: s.grid.map(function (row) { return row.map(N); }) };
    });
    var lines = (cfg.lines || []).map(function (l) { return { l: l, pts: l.pts.map(N) }; });
    var light = (function (v) { var m = Math.hypot(v[0], v[1], v[2]); return v.map(function (c) { return c / m; }); })([-0.5, -0.7, 0.9]);

    function ticksFor(ax) {
      if (!ax.log) {
        var t = niceTicks(ax.domain[0], ax.domain[1], 5);
        return t.ticks.filter(function (v) { return v >= ax.domain[0] - 1e-12 && v <= ax.domain[1] + 1e-12; })
          .map(function (v) { return [v, PEC.fmt(v, decimals(t.step))]; });
      }
      var out = [];
      for (var e = Math.floor(lg(ax.domain[0])); e <= Math.ceil(lg(ax.domain[1])); e++)
        [1, 2, 5].forEach(function (m) {
          var v = m * Math.pow(10, e);
          if (v >= ax.domain[0] * 0.999 && v <= ax.domain[1] * 1.001) out.push([v, PEC.sig(v, 2)]);
        });
      return out;
    }

    function render() {
      var W = Math.max(280, container.clientWidth);
      var H = Math.round(Math.min(600, Math.max(340, W * 0.72)));
      var cy = Math.cos(view.yaw), sy = Math.sin(view.yaw), cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
      var sc = Math.min(W * 0.245, H * 0.3), ox = W / 2, oy = H / 2 - sc * 0.04;
      function pr(q) {                      // -> [screen x, screen y, depth (larger: farther)]
        var a = q[0] * cy - q[1] * sy, d = q[0] * sy + q[1] * cy;
        return [ox + sc * a, oy - sc * (q[2] * cp + d * sp), d * cp - q[2] * sp];
      }
      container.innerHTML = "";
      /* surfaces go on a canvas with a depth buffer (correct hiding even where the
         two sides of a two-phase lens nearly touch); text goes on an SVG above it */
      var dpr = view.fast ? 1 : Math.min(2, window.devicePixelRatio || 1);   // coarser while dragging
      var cw = Math.round(W * dpr), ch = Math.round(H * dpr);
      var cv = document.createElement("canvas");
      cv.width = cw; cv.height = ch;
      cv.style.cssText = "position:absolute;left:0;top:0;width:" + W + "px;height:" + H + "px";
      container.appendChild(cv);
      var svg = el("svg", { viewBox: "0 0 " + W + " " + H, width: W, height: H, role: "img", "aria-label": cfg.aria || "" }, container);
      svg.style.cssText = "position:relative;cursor:grab";
      var ctx = cv.getContext("2d"), img = ctx.createImageData(cw, ch), px = img.data;
      var zb = new Float32Array(cw * ch).fill(Infinity), ZB = zb, PX = px;   // ZB, PX: the buffers being drawn into
      function S(q) { var p = pr(q); return [p[0] * dpr, p[1] * dpr, p[2]]; }

      function tri(a, b, c, col) {                 // flat-shaded, depth-tested triangle
        var x0 = Math.max(0, Math.floor(Math.min(a[0], b[0], c[0]))), x1 = Math.min(cw - 1, Math.ceil(Math.max(a[0], b[0], c[0])));
        var y0 = Math.max(0, Math.floor(Math.min(a[1], b[1], c[1]))), y1 = Math.min(ch - 1, Math.ceil(Math.max(a[1], b[1], c[1])));
        var den = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]);
        if (Math.abs(den) < 1e-9) return;
        for (var y = y0; y <= y1; y++) {
          var yc = y + 0.5;
          for (var x = x0; x <= x1; x++) {
            var xc = x + 0.5;
            var w0 = ((b[1] - c[1]) * (xc - c[0]) + (c[0] - b[0]) * (yc - c[1])) / den;
            var w1 = ((c[1] - a[1]) * (xc - c[0]) + (a[0] - c[0]) * (yc - c[1])) / den;
            var w2 = 1 - w0 - w1;
            if (w0 < -1e-6 || w1 < -1e-6 || w2 < -1e-6) continue;
            var z = w0 * a[2] + w1 * b[2] + w2 * c[2], k = y * cw + x;
            if (z < ZB[k]) { ZB[k] = z; k *= 4; PX[k] = col[0]; PX[k + 1] = col[1]; PX[k + 2] = col[2]; PX[k + 3] = 255; }
          }
        }
      }
      function seg(a, b, col, width, alpha, tol, hidden) {  // hidden parts: alpha 'hidden' (or not drawn)
        var n = Math.max(1, Math.ceil(Math.hypot(b[0] - a[0], b[1] - a[1]) / 0.7)), r = width * dpr / 2;
        for (var s = 0; s <= n; s++) {
          var t = s / n, x = a[0] + (b[0] - a[0]) * t, y = a[1] + (b[1] - a[1]) * t, z = a[2] + (b[2] - a[2]) * t;
          var cx0 = Math.round(x), cy0 = Math.round(y);
          if (cx0 < 0 || cy0 < 0 || cx0 >= cw || cy0 >= ch) continue;
          var al = z > ZB[cy0 * cw + cx0] + tol ? (hidden || 0) : alpha;
          if (!al) continue;
          for (var yy = Math.floor(y - r); yy <= Math.ceil(y + r); yy++)
            for (var xx = Math.floor(x - r); xx <= Math.ceil(x + r); xx++) {
              if (xx < 0 || yy < 0 || xx >= cw || yy >= ch || (xx + 0.5 - x) * (xx + 0.5 - x) + (yy + 0.5 - y) * (yy + 0.5 - y) > r * r) continue;
              var k = (yy * cw + xx) * 4;
              if (PX[k + 3] === 0) { PX[k] = col[0]; PX[k + 1] = col[1]; PX[k + 2] = col[2]; PX[k + 3] = Math.round(255 * al); }
              else {                                      // blend over what is below
                PX[k] += (col[0] - PX[k]) * al; PX[k + 1] += (col[1] - PX[k + 1]) * al; PX[k + 2] += (col[2] - PX[k + 2]) * al;
                PX[k + 3] = Math.max(PX[k + 3], Math.round(255 * al));
              }
            }
        }
      }

      /* each surface into its own depth buffer, with its mesh lines */
      var bufs = [];
      surf.forEach(function (Sf) {
        ZB = new Float32Array(cw * ch).fill(Infinity); PX = new Uint8ClampedArray(cw * ch * 4);
        var g = Sf.g;
        for (var i = 0; i + 1 < g.length; i++)
          for (var j = 0; j + 1 < g[i].length; j++) {
            var q = [g[i][j], g[i][j + 1], g[i + 1][j + 1], g[i + 1][j]];
            if (!q[0] || !q[1] || !q[2] || !q[3]) continue;
            var u = [q[2][0] - q[0][0], q[2][1] - q[0][1], q[2][2] - q[0][2]],
                v = [q[3][0] - q[1][0], q[3][1] - q[1][1], q[3][2] - q[1][2]];
            var n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
            var m = Math.hypot(n[0], n[1], n[2]) || 1;
            var sh = 0.62 + 0.38 * Math.abs((n[0] * light[0] + n[1] * light[1] + n[2] * light[2]) / m);
            var col = Sf.base.map(function (c) { return Math.round(c * sh); });
            var P = q.map(S);
            tri(P[0], P[1], P[2], col); tri(P[0], P[2], P[3], col);
          }
        var dark = Sf.base.map(function (c) { return Math.round(c * 0.45); });
        for (var i = 0; i < g.length; i++)
          for (var j = 0; j < g[i].length; j++) {
            if (g[i][j] && j + 1 < g[i].length && g[i][j + 1] && (i % 1 === 0)) seg(S(g[i][j]), S(g[i][j + 1]), dark, 0.6, 0.35, 0.012);
            if (g[i][j] && i + 1 < g.length && g[i + 1][j] && j % 3 === 0) seg(S(g[i][j]), S(g[i + 1][j]), dark, 0.6, 0.35, 0.012);
          }
        bufs.push({ zb: ZB, px: PX });
      });
      /* composite: at each pixel the surfaces' front fragments, far to near, the nearer
         ones see-through, so a surface behind another still shows */
      var NA = 0.62, FA = 0.94, nb = bufs.length, fr = [];
      for (var q = 0, Np = cw * ch; q < Np; q++) {
        var nf = 0;
        for (var b = 0; b < nb; b++) if (bufs[b].zb[q] < Infinity) fr[nf++] = b;
        if (!nf) continue;
        if (nf === 2 && bufs[fr[0]].zb[q] < bufs[fr[1]].zb[q]) { var t0 = fr[0]; fr[0] = fr[1]; fr[1] = t0; }
        else if (nf > 2) fr.length = nf, fr.sort(function (u, v) { return bufs[v].zb[q] - bufs[u].zb[q]; });
        var r = 0, g2 = 0, bl = 0, A = 0, o = q * 4;
        for (var f = 0; f < nf; f++) {
          // a surface lying almost on the one behind it is drawn solid (no muddy mix)
          var gap = f > 0 ? bufs[fr[f - 1]].zb[q] - bufs[fr[f]].zb[q] : 0, tg = Math.min(1, gap / 0.08);
          var B = bufs[fr[f]].px, al = (f === 0 ? FA : 1 - (1 - NA) * tg) * B[o + 3] / 255;
          r = B[o] * al + r * (1 - al); g2 = B[o + 1] * al + g2 * (1 - al); bl = B[o + 2] * al + bl * (1 - al);
          A = al + A * (1 - al);
          if (bufs[fr[f]].zb[q] < zb[q]) zb[q] = bufs[fr[f]].zb[q];
        }
        px[o] = r / A; px[o + 1] = g2 / A; px[o + 2] = bl / A; px[o + 3] = Math.round(255 * A);
      }
      ZB = zb; PX = px;
      /* see-through planes (e.g. phase regions on one diagram), hidden behind nearer surfaces */
      (cfg.planes || []).forEach(function (pl) {
        var col = rgb(pl.color), al = pl.alpha || 0.3, mask = new Uint8Array(cw * ch);
        function cover(a, b, c) {
          var x0 = Math.max(0, Math.floor(Math.min(a[0], b[0], c[0]))), x1 = Math.min(cw - 1, Math.ceil(Math.max(a[0], b[0], c[0])));
          var y0 = Math.max(0, Math.floor(Math.min(a[1], b[1], c[1]))), y1 = Math.min(ch - 1, Math.ceil(Math.max(a[1], b[1], c[1])));
          var den = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]);
          if (Math.abs(den) < 1e-9) return;
          for (var y = y0; y <= y1; y++)
            for (var x = x0; x <= x1; x++) {
              var xc = x + 0.5, yc = y + 0.5;
              var w0 = ((b[1] - c[1]) * (xc - c[0]) + (c[0] - b[0]) * (yc - c[1])) / den;
              var w1 = ((c[1] - a[1]) * (xc - c[0]) + (a[0] - c[0]) * (yc - c[1])) / den;
              var w2 = 1 - w0 - w1;
              if (w0 < -1e-6 || w1 < -1e-6 || w2 < -1e-6) continue;
              var k = y * cw + x;
              if (w0 * a[2] + w1 * b[2] + w2 * c[2] <= zb[k] + 0.01) mask[k] = 1;
            }
        }
        pl.quads.forEach(function (q) {
          var P = q.map(function (p) { return S(N(p)); });
          cover(P[0], P[1], P[2]); cover(P[0], P[2], P[3]);
        });
        for (var k = 0, n = cw * ch; k < n; k++) {
          if (!mask[k]) continue;
          var o = k * 4;
          if (px[o + 3] === 0) { px[o] = col[0]; px[o + 1] = col[1]; px[o + 2] = col[2]; px[o + 3] = Math.round(255 * al); }
          else { px[o] += (col[0] - px[o]) * al; px[o + 1] += (col[1] - px[o + 1]) * al; px[o + 2] += (col[2] - px[o + 2]) * al; }
        }
      });
      /* the box and the floor grid, then the lines on the surface */
      var C = [], k;
      for (k = 0; k < 8; k++) C.push([k & 1 ? 1 : -1, k & 2 ? 1 : -1, k & 4 ? 1 : -1]);
      var E = [[0, 1], [2, 3], [4, 5], [6, 7], [0, 2], [1, 3], [4, 6], [5, 7], [0, 4], [1, 5], [2, 6], [3, 7]];
      var cAxis = rgb("--axis"), cHair = rgb("--hair");
      [cfg.x, cfg.y].forEach(function (ax, which) {
        ticksFor(ax).forEach(function (t) {
          var u = nrm(t[0], ax);
          seg(S(which ? [-1, u, -1] : [u, -1, -1]), S(which ? [1, u, -1] : [u, 1, -1]), cHair, 1, 1, 0.01);
        });
      });
      E.forEach(function (e) { seg(S(C[e[0]]), S(C[e[1]]), cAxis, 1, 1, 0.01); });
      /* the current diagram and the azeotrope line, on the surface */
      lines.forEach(function (Ln) {
        var col = rgb(Ln.l.color);
        for (var i = 0; i + 1 < Ln.pts.length; i++)
          if (Ln.pts[i] && Ln.pts[i + 1]) seg(S(Ln.pts[i]), S(Ln.pts[i + 1]), col, Ln.l.width || 2.5, 1, 0.03, 0.45);
      });
      ctx.putImageData(img, 0, 0);

      /* labels on the edges nearest the viewer */
      var ctr = pr([0, 0, 0]);
      function label(q, text, cls, push) {
        var p = pr(q), dx = p[0] - ctr[0], dy = p[1] - ctr[1], m = Math.hypot(dx, dy) || 1;
        var x = p[0] + dx / m * push, y = p[1] + dy / m * push;
        var t = el("text", { x: x, y: y + 4, "text-anchor": Math.abs(dx / m) < 0.3 ? "middle" : dx > 0 ? "start" : "end", class: cls }, svg);
        t.textContent = text;
        var bb = t.getBBox ? t.getBBox() : null;          // keep it inside the picture
        if (bb && bb.width) {
          if (bb.x < 2) t.setAttribute("x", x + 2 - bb.x);
          else if (bb.x + bb.width > W - 2) t.setAttribute("x", x - (bb.x + bb.width - W + 2));
        }
      }
      var yEdge = pr([0, -1, -1])[2] < pr([0, 1, -1])[2] ? -1 : 1;
      var xEdge = pr([-1, 0, -1])[2] < pr([1, 0, -1])[2] ? -1 : 1;
      ticksFor(cfg.x).forEach(function (t) { label([nrm(t[0], cfg.x), yEdge, -1], t[1], "tick-3d", 12); });
      ticksFor(cfg.y).forEach(function (t) { label([xEdge, nrm(t[0], cfg.y), -1], t[1], "tick-3d", 12); });
      var zc = [[-1, -1], [1, -1], [-1, 1], [1, 1]].reduce(function (best, c) {
        return pr([c[0], c[1], 0])[0] < pr([best[0], best[1], 0])[0] ? c : best; });
      ticksFor(cfg.z).forEach(function (t) {          // the lowest one would sit on the floor labels
        if (nrm(t[0], cfg.z) > -0.97) label([zc[0], zc[1], nrm(t[0], cfg.z)], t[1], "tick-3d", 10); });
      label([0, yEdge, -1], cfg.x.title, "axis-title", 40);
      label([xEdge, 0, -1], cfg.y.title, "axis-title", 46);
      var zt = pr([zc[0], zc[1], 1.16]);
      var tz = el("text", { x: zt[0], y: zt[1], "text-anchor": "middle", class: "axis-title" }, svg); tz.textContent = cfg.z.title;
    }

    /* drag to rotate, double-click to reset; listeners live on the container */
    if (!container.__drag) {
      container.__drag = true;
      var start = null, raf = 0;
      container.addEventListener("pointerdown", function (e) {
        if (!container.__surface) return;
        start = { x: e.clientX, y: e.clientY, yaw: container.__view.yaw, pitch: container.__view.pitch };
        container.setPointerCapture(e.pointerId);
      });
      container.addEventListener("pointermove", function (e) {
        if (!start || !container.__surface) return;
        var v = container.__view;
        v.fast = true;
        v.yaw = start.yaw + (e.clientX - start.x) * 0.01;
        v.pitch = Math.max(-0.15, Math.min(1.45, start.pitch + (e.clientY - start.y) * 0.01));
        if (!raf) raf = requestAnimationFrame(function () { raf = 0; if (container.__surface) container.__surface(); });
      });
      ["pointerup", "pointercancel"].forEach(function (t) {
        container.addEventListener(t, function () {
          start = null;
          if (container.__view.fast) { container.__view.fast = false; if (container.__surface) container.__surface(); }
        });
      });
      container.addEventListener("dblclick", function () {
        if (!container.__surface) return;
        container.__view.yaw = -0.85; container.__view.pitch = 0.38; container.__surface();
      });
    }
    container.__surface = render;
    render();
    if (container.__ro) container.__ro.disconnect();
    var last = container.clientWidth;
    container.__ro = new ResizeObserver(function () {
      if (Math.abs(container.clientWidth - last) > 4) { last = container.clientWidth; render(); }
    });
    container.__ro.observe(container);
  };

  PEC.zip = function (a, b) {
    return a.map(function (v, i) { return v == null || b[i] == null ? null : [v, b[i]]; });
  };
  PEC.interp = function (xs, ys, x) {
    for (var i = 0; i < xs.length - 1; i++) {
      if (xs[i] == null || xs[i + 1] == null) continue;
      if ((xs[i] - x) * (xs[i + 1] - x) <= 0 && xs[i] !== xs[i + 1]) {
        return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i]);
      }
    }
    return null;
  };

  /* ---------- page furniture: math, contents highlighting, data-bound numbers ---------- */
  PEC.typeset = function (root) {
    if (!window.renderMathInElement) return setTimeout(function () { PEC.typeset(root); }, 60);
    renderMathInElement(root || document.body, {
      delimiters: [{ left: "\\[", right: "\\]", display: true }, { left: "\\(", right: "\\)", display: false }],
      throwOnError: false
    });
  };
  PEC.bindNumbers = function (data) {
    function get(path) {
      return path.split(".").reduce(function (o, key) {
        if (o == null) return undefined;
        if (Array.isArray(o) && /^-?\d+$/.test(key)) { var i = +key; return o[i < 0 ? o.length + i : i]; }
        return o[key];
      }, data);
    }
    document.querySelectorAll("[data-v]").forEach(function (node) {
      var v = get(node.getAttribute("data-v"));
      if (v == null) return;
      if (node.hasAttribute("data-e")) node.innerHTML = PEC.sci(v, +node.getAttribute("data-e"));
      else node.textContent = PEC.fmt(v, +(node.getAttribute("data-d") || 3));
    });
  };
  document.addEventListener("DOMContentLoaded", function () {
    var links = Array.prototype.slice.call(document.querySelectorAll("nav.toc a"));
    if (links.length) {
      var targets = links.map(function (a) { return document.querySelector(a.getAttribute("href")); });
      var onScroll = function () {
        var cur = 0;
        targets.forEach(function (t, i) { if (t && t.getBoundingClientRect().top < 140) cur = i; });
        links.forEach(function (a, i) { a.setAttribute("aria-current", i === cur ? "true" : "false"); });
      };
      window.addEventListener("scroll", onScroll, { passive: true });
      onScroll();
    }
    if (document.querySelector("[data-math]")) PEC.typeset();
  });
})();
