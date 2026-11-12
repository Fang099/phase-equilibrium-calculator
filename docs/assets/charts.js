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
            "class": "series-line", stroke: "var(" + s.color + ")", "stroke-width": s.width || 2, opacity: s.opacity || 1 }, plot);
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
