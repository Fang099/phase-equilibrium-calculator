/* Calculator page. All thermodynamics is done by the project's Python modules
   (web_api.py), running in a Web Worker through Pyodide. This file only
   collects the choices, sends requests and draws what comes back. */
(function () {
  "use strict";
  var PEC = window.PEC, $ = function (id) { return document.getElementById(id); };
  var query = new URLSearchParams(location.search);

  /* ------------------------------------------------------------ engine */
  /* Two ways to reach the same Python code (web_api.handle):
     - normally, a Web Worker runs it in the browser with Pyodide;
     - when the page is served by scripts/dev_server.py on this computer, requests
       go to that server instead, so the calculator also works with no internet.
     ?engine=pyodide or ?engine=http forces one of them. */
  function Engine(onStatus, onFatal) {
    var self = this, forced = query.get("engine");
    var onLocalhost = /^(localhost|127\.0\.0\.1|\[::1\])$/.test(location.hostname);
    this.local = false;
    function viaHttp(msg) {
      return fetch("/api", { method: "POST", body: JSON.stringify(msg) }).then(function (r) {
        if (!r.ok) throw new Error("no local engine");
        return r.json();
      });
    }
    function startWorker() {
      var worker = new Worker("assets/worker.js"), pending = {}, seq = 0;
      worker.onmessage = function (e) {
        var d = e.data;
        if (d.type === "status") onStatus(d.text);
        else if (d.type === "fatal") onFatal(d.error);
        else if (d.type === "result" && pending[d.id]) { pending[d.id](d.result); delete pending[d.id]; }
      };
      worker.onerror = function (e) { onFatal(e.message || "the worker stopped"); };
      return function (msg) {
        return new Promise(function (res) { var id = ++seq; pending[id] = res; worker.postMessage({ id: id, msg: msg }); });
      };
    }
    var impl;
    if (forced === "http" || (onLocalhost && forced !== "pyodide")) {
      impl = viaHttp({ op: "ping" }).then(function (r) {
        if (!r || !r.ok) throw new Error("no local engine");
        self.local = true;
        return viaHttp;
      }).catch(function () { return startWorker(); });
    } else {
      impl = Promise.resolve(startWorker());
    }
    this.call = function (msg) { return impl.then(function (f) { return f(msg); }); };
  }

  /* ------------------------------------------------------------ definitions */
  var EOS = ["PR", "SRK", "vdW"];
  var MODEL_SHORT = { PR: "Peng–Robinson", SRK: "SRK", vdW: "van der Waals", ideal: "Ideal", margules: "Margules", vanlaar: "van Laar", nrtl: "NRTL" };
  var X_OPTS = [["xy", "Composition, liquid x and vapor y"], ["x1", "Liquid composition x"], ["y1", "Vapor composition y"], ["T", "Temperature"], ["P", "Pressure"]];
  var Y_OPTS = [["P", "Pressure"], ["T", "Temperature"], ["y1", "Vapor composition y"], ["K", "K-values, y/x"],
    ["alpha", "Relative volatility"], ["gamma", "Activity coefficients γ"], ["lngamma", "ln γ"], ["f", "Fugacities"],
    ["phiL", "Fugacity coefficients, liquid"], ["phiV", "Fugacity coefficients, vapor"], ["GE", "Excess Gibbs energy Gᴱ/RT"],
    ["Gmix", "Gibbs energy of mixing ΔG/RT"], ["stab", "Stability of the liquid"]];
  var PRESETS = [["Pxy", "P", "xy", "T"], ["Txy", "T", "xy", "P"], ["xy", "y1", "x1", null], ["P–T", "P", "T", "x"],
    ["K-values", "K", "x1", null], ["Activity coefficients", "gamma", "x1", null], ["Fugacity", "f", "x1", "T"],
    ["Excess Gibbs energy", "GE", "x1", null], ["Relative volatility", "alpha", "x1", null]];
  var PAIRED = { K: ["K1", "K2"], gamma: ["gamma1", "gamma2"], lngamma: ["gamma1", "gamma2"], f: ["f1", "f2"], phiL: ["phiL1", "phiL2"], phiV: ["phiV1", "phiV2"] };
  var SINGLE = { alpha: "alpha", GE: "GE_RT", Gmix: "Gmix_RT", stab: "stability", y1: "y1", P: "P", T: "T" };
  var SOURCE_LABEL = { chemsep_nrtl: "ChemSep NRTL library", chemsep_pr: "ChemSep Peng–Robinson library",
    nrtl_infinite_dilution: "matched to the NRTL model at infinite dilution", none: "no database value",
    vanlaar_fit: "database, A = a + b/T fitted to the DECHEMA NRTL parameters" };

  /* the edition (set at build time in site/edition.json) limits what is offered */
  var ED = window.PEC_EDITION || {};
  function offered(kind, v) { return !ED[kind] || ED[kind].indexOf(v) >= 0; }
  var ACT = ["ideal", "margules", "vanlaar", "nrtl"].filter(function (m) { return offered("models", m); });
  EOS = EOS.filter(function (m) { return offered("models", m); });
  PRESETS = PRESETS.filter(function (p) { return offered("presets", p[0]); });
  X_OPTS = X_OPTS.filter(function (o) { return offered("x_axes", o[0]); });
  Y_OPTS = Y_OPTS.filter(function (o) { return offered("y_axes", o[0]); });
  if (!offered("models", "nrtl"))
    SOURCE_LABEL.nrtl_infinite_dilution = "infinite-dilution activity coefficients from the ChemSep interaction data";
  var DEF = ED["default"] || { a: "ethanol", b: "water", model: "nrtl" };

  var state = { a: DEF.a, b: DEF.b, model: DEF.model, compare: "", x: "xy", y: "T", hold: "P", T: 350, P: 1.01325, x1: 0.5, params: {} };
  var catalog = null, pair = null, last = null, engine, busy = false, dirty = false, timer = null;

  /* ------------------------------------------------------------ URL state */
  function readHash() {
    var h = new URLSearchParams(location.hash.slice(1));
    ["a", "b", "model", "compare", "x", "y", "hold"].forEach(function (k) { if (h.get(k) != null) state[k] = h.get(k); });
    ["T", "P", "x1"].forEach(function (k) { if (h.get(k) != null && isFinite(+h.get(k))) state[k] = +h.get(k); });
    var p = {};
    ["kij", "A12", "A21", "tdep", "dg12", "dg21", "alpha"].forEach(function (k) { if (h.get(k) != null && isFinite(+h.get(k))) p[k] = +h.get(k); });
    p.__hasT = h.get("T") != null;
    return p;
  }
  function writeHash() {
    var h = new URLSearchParams();
    ["a", "b", "model", "y", "x", "hold"].forEach(function (k) { h.set(k, state[k]); });
    h.set(state.hold === "x" ? "x1" : state.hold, state[state.hold === "x" ? "x1" : state.hold]);
    if (state.compare) h.set("compare", state.compare);
    activeParams().forEach(function (k) { if (state.params[k] !== defaults()[k]) h.set(k, state.params[k]); });
    history.replaceState(null, "", "#" + h.toString());
  }

  /* ------------------------------------------------------------ helpers */
  function comp(id) { return catalog.compounds.filter(function (c) { return c.id === id; })[0]; }
  function nameA() { return comp(state.a).name; }
  function nameB() { return comp(state.b).name; }
  function lower(s) { return s.toLowerCase(); }
  function isComp(x) { return x === "xy" || x === "x1" || x === "y1"; }
  function activeParams() {
    if (EOS.indexOf(state.model) >= 0) return ["kij"];
    if (state.model === "margules") return ["A12", "A21"];
    if (state.model === "vanlaar") return pair.vanlaar ? ["A12", "A21", "tdep"] : [];
    if (state.model === "nrtl") return ["dg12", "dg21", "alpha"];
    return [];
  }
  function defaults() {
    var n = pair.nrtl || { dg12: 0, dg21: 0, alpha: 0.3 };
    var A = (state.model === "vanlaar" && pair.vanlaar) ? pair.vanlaar : pair.A;
    return { kij: state.model === "PR" ? pair.kij.value : 0, A12: A.A12, A21: A.A21, tdep: 1, dg12: n.dg12, dg21: n.dg21, alpha: n.alpha };
  }
  function seg(host, items, current, onPick, disabled) {
    host.innerHTML = "";
    items.forEach(function (it) {
      var b = document.createElement("button");
      b.type = "button"; b.textContent = it[1];
      b.setAttribute("aria-pressed", String(it[0] === current));
      if (disabled && disabled(it[0])) b.disabled = true;
      b.addEventListener("click", function () { onPick(it[0]); });
      host.appendChild(b);
    });
  }

  /* axis choices decide what is held constant: a binary two-phase system has
     two degrees of freedom, so a diagram fixes one variable and sweeps another */
  function resolve() {
    if (isComp(state.x)) {
      if (state.x === "xy" && state.y !== "P" && state.y !== "T") state.x = "x1";
      if (state.x === "y1" && state.y === "y1") state.x = "x1";
      if (state.y === "P") state.hold = "T";
      else if (state.y === "T") state.hold = "P";
      else if (state.hold === "x") state.hold = "T";
    } else {
      if (state.y === state.x) state.y = state.x === "T" ? "P" : "T";
      state.hold = "x";
    }
  }
  function holdAllowed(h) {
    if (!isComp(state.x)) return h === "x";
    if (state.y === "P") return h === "T";
    if (state.y === "T") return h === "P";
    return h !== "x";
  }

  /* ------------------------------------------------------------ controls */
  function fillCompounds() {
    var fam = {}, order = [];
    catalog.compounds.forEach(function (c) { if (!fam[c.family]) { fam[c.family] = []; order.push(c.family); } fam[c.family].push(c); });
    var plural = { alkane: "Alkanes", aromatic: "Aromatics", alcohol: "Alcohols", water: "Water", ketone: "Ketones", ester: "Esters", ether: "Ethers", chlorinated: "Chlorinated solvents" };
    [["comp-a", "a"], ["comp-b", "b"]].forEach(function (pr) {
      var sel = $(pr[0]); sel.innerHTML = "";
      order.forEach(function (f) {
        var g = document.createElement("optgroup"); g.label = plural[f] || f;
        fam[f].forEach(function (c) { var o = document.createElement("option"); o.value = c.id; o.textContent = c.name + "  (" + c.formula + ")"; g.appendChild(o); });
        sel.appendChild(g);
      });
      sel.value = state[pr[1]];
      sel.addEventListener("change", function () {
        var other = pr[1] === "a" ? "b" : "a";
        if (sel.value === state[other]) state[other] = state[pr[1]];      // picking the same one swaps them
        state[pr[1]] = sel.value;
        $("comp-a").value = state.a; $("comp-b").value = state.b;
        loadPair({});
      });
    });
    $("swap").addEventListener("click", function () {
      var t = state.a; state.a = state.b; state.b = t;
      $("comp-a").value = state.a; $("comp-b").value = state.b;
      state.x1 = Math.round((1 - state.x1) * 1000) / 1000;
      loadPair({});
    });
  }
  function drawModels() {
    function pick(m) { state.model = m; state.params = defaults(); if (state.compare === m) state.compare = ""; drawModels(); drawParams(); drawCompare(); schedule(0); }
    seg($("models-eos"), EOS.map(function (m) { return [m, MODEL_SHORT[m]]; }), state.model, pick);
    seg($("models-act"), ACT.map(function (m) { return [m, MODEL_SHORT[m]]; }), state.model, pick,
      function (m) { return m === "vanlaar" && !pair.vanlaar; });
    var vb = [].slice.call($("models-act").children).filter(function (b) { return b.textContent === MODEL_SHORT.vanlaar; })[0];
    if (vb && vb.disabled) vb.title = "The database has no van Laar parameters for this pair";
  }
  function drawCompare() {
    var sel = $("compare");
    sel.innerHTML = "<option value=''>Nothing</option>" + catalog.models.filter(function (m) { return m.id !== state.model && offered("models", m.id) && (m.id !== "vanlaar" || pair.vanlaar); })
      .map(function (m) { return "<option value='" + m.id + "'>" + m.label + "</option>"; }).join("");
    sel.value = state.compare;
  }
  var PARAM_META = {
    kij: { label: function () { return "Binary interaction parameter k₁₂"; }, min: -0.3, max: 0.3, step: 0.005 },
    A12: { label: function () { return "A₁₂ = ln γ∞ of " + lower(nameA()) + " in " + lower(nameB()) + atRef(); }, min: -3, max: 5, step: 0.01 },
    A21: { label: function () { return "A₂₁ = ln γ∞ of " + lower(nameB()) + " in " + lower(nameA()) + atRef(); }, min: -3, max: 5, step: 0.01 },
    dg12: { label: function () { return "Δg₁₂ (J/mol)"; }, min: -6000, max: 16000, step: 10 },
    dg21: { label: function () { return "Δg₂₁ (J/mol)"; }, min: -6000, max: 16000, step: 10 },
    alpha: { label: function () { return "Non-randomness α"; }, min: 0.05, max: 0.9, step: 0.005 }
  };
  function atRef() {
    return (state.model === "vanlaar" && pair.vanlaar && state.params.tdep) ? ", at " + PEC.fmt(pair.vanlaar.T_ref, 1) + " K" : "";
  }
  function paramSource() {
    if (state.model === "PR") return pair.kij.source;
    if (state.model === "vanlaar") return pair.vanlaar ? pair.vanlaar.source : "none";
    if (state.model === "SRK" || state.model === "vdW") return "none";
    if (state.model === "nrtl") return pair.nrtl ? pair.nrtl.source : "none";
    if (state.model === "ideal") return null;
    return pair.A.source;
  }
  function drawParams() {
    var host = $("params"), keys = activeParams();
    host.innerHTML = "";
    keys.forEach(function (k) {
      if (k === "tdep") {                                   // a switch, not a slider
        var row0 = document.createElement("div"); row0.className = "param";
        row0.innerHTML = "<label class='check'><input type='checkbox' id='p-tdep'" + (state.params.tdep ? " checked" : "") +
          "> Vary with temperature</label>";
        host.appendChild(row0);
        row0.querySelector("input").addEventListener("change", function () {
          state.params.tdep = this.checked ? 1 : 0; drawParams(); schedule(0);
        });
        return;
      }
      var meta = PARAM_META[k], row = document.createElement("div");
      row.className = "param";
      row.innerHTML = "<label for='p-" + k + "'>" + meta.label() + "</label>" +
        "<input type='range' id='r-" + k + "' min='" + meta.min + "' max='" + meta.max + "' step='" + meta.step + "' aria-label='" + k + "'>" +
        "<input type='number' id='p-" + k + "' step='" + meta.step + "'>";
      host.appendChild(row);
      var r = row.querySelector("input[type=range]"), n = row.querySelector("input[type=number]");
      r.value = n.value = state.params[k];
      function set(v) { if (!isFinite(v)) return; state.params[k] = v; r.value = v; n.value = v; drawParamNote(); schedule(120); }
      r.addEventListener("input", function () { set(+r.value); });
      n.addEventListener("change", function () { set(+n.value); });
    });
    var note = document.createElement("p"); note.className = "source"; note.id = "param-note"; host.appendChild(note);
    drawParamNote();
  }
  function drawParamNote() {
    var note = $("param-note"), src = paramSource(), def = defaults();
    if (!note) return;
    if (src == null) { note.textContent = ""; return; }
    var edited = activeParams().some(function (k) { return state.params[k] !== def[k]; });
    var extra = "";
    if (state.model === "vanlaar" && pair.vanlaar) {
      var v = pair.vanlaar;
      extra = state.params.tdep
        ? " A changes with temperature as a + b/T (b₁₂ = " + PEC.fmt(v.b12, 1) + " K, b₂₁ = " + PEC.fmt(v.b21, 1) + " K), fitted from " +
          PEC.fmt(v.Tmin, 0) + " to " + PEC.fmt(v.Tmax, 0) + " K."
        : " A is held constant at these values.";
    }
    note.innerHTML = (edited ? "Edited. " : "") + "Starting values: " + (SOURCE_LABEL[src] || src) + "." + extra +
      (edited ? " <button class='linkbtn' type='button' id='reset-params'>Reset</button>" : "");
    var b = $("reset-params");
    if (b) b.addEventListener("click", function () { state.params = defaults(); drawParams(); schedule(0); });
  }
  function drawAxes() {
    resolve();
    var ys = $("axis-y"), xs = $("axis-x");
    ys.innerHTML = Y_OPTS.map(function (o) { return "<option value='" + o[0] + "'>" + o[1] + "</option>"; }).join("");
    xs.innerHTML = X_OPTS.map(function (o) {
      var bad = o[0] === "xy" && state.y !== "P" && state.y !== "T";
      return "<option value='" + o[0] + "'" + (bad ? " disabled" : "") + ">" + o[1] + "</option>";
    }).join("");
    ys.value = state.y; xs.value = state.x;
    seg($("presets"), PRESETS.map(function (p, i) { return [String(i), p[0]]; }),
      String(PRESETS.findIndex(function (p) { return p[1] === state.y && p[2] === state.x; })),
      function (i) { var p = PRESETS[+i]; state.y = p[1]; state.x = p[2]; if (p[3]) state.hold = p[3]; drawAxes(); schedule(0); });
    seg($("hold"), [["T", "Temperature"], ["P", "Pressure"], ["x", "Composition"]], state.hold,
      function (h) { state.hold = h; drawAxes(); schedule(0); }, function (h) { return !holdAllowed(h); });
    var r = $("hold-range"), n = $("hold-value"), u = $("hold-unit");
    if (state.hold === "T") {
      r.min = pair.ranges.T[0]; r.max = pair.ranges.T[1]; r.step = 1; r.value = state.T;
      n.value = state.T; n.step = 1; u.textContent = "K";
    } else if (state.hold === "P") {          // logarithmic slider
      r.min = Math.log10(pair.ranges.P[0]); r.max = Math.log10(pair.ranges.P[1]); r.step = 0.005; r.value = Math.log10(state.P);
      n.value = +state.P.toPrecision(5); n.step = 0.1; u.textContent = "bar";
    } else {
      r.min = 0.02; r.max = 0.98; r.step = 0.01; r.value = state.x1;
      n.value = state.x1; n.step = 0.01; u.textContent = "mole fraction of " + lower(nameA());
    }
  }
  function bindAxes() {
    /* choosing the same variable for both axes moves the other axis to a sensible partner */
    $("axis-y").addEventListener("change", function () {
      state.y = this.value;
      if (state.y === state.x) state.x = (state.y === "T" || state.y === "P") ? "xy" : "x1";
      drawAxes(); schedule(0);
    });
    $("axis-x").addEventListener("change", function () {
      state.x = this.value;
      if (state.x === state.y) state.y = state.x === "T" ? "P" : state.x === "P" ? "T" : "P";
      drawAxes(); schedule(0);
    });
    $("compare").addEventListener("change", function () { state.compare = this.value; schedule(0); });
    $("hold-range").addEventListener("input", function () {
      var v = +this.value;
      if (state.hold === "T") state.T = v;
      else if (state.hold === "P") state.P = +Math.pow(10, v).toPrecision(4);
      else state.x1 = v;
      $("hold-value").value = state.hold === "x" ? state.x1 : state[state.hold];
      schedule(120);
    });
    $("hold-value").addEventListener("change", function () {
      var v = +this.value; if (!isFinite(v) || v <= 0) return;
      if (state.hold === "T") state.T = v;
      else if (state.hold === "P") state.P = v;
      else state.x1 = Math.min(0.999, Math.max(0.001, v));
      drawAxes(); schedule(0);
    });
    var tabs = ["data", "input", "sources"];
    tabs.forEach(function (t) {
      $("tab-" + t).addEventListener("click", function () {
        tabs.forEach(function (o) { $("tab-" + o).setAttribute("aria-selected", String(o === t)); $("panel-" + o).hidden = o !== t; });
      });
    });
    $("download").addEventListener("click", downloadCsv);
  }

  /* ------------------------------------------------------------ requests */
  function loadPair(urlParams) {
    overlay("Loading " + lower(nameA()) + " and " + lower(nameB()), "");
    return engine.call({ op: "pair", a: state.a, b: state.b }).then(function (p) {
      pair = p;
      if (state.model === "vanlaar" && !pair.vanlaar) state.model = "PR";      // no database values for this pair
      if (state.compare === "vanlaar" && !pair.vanlaar) state.compare = "";
      if (urlParams.__hasT) state.T = Math.min(Math.max(state.T, p.ranges.T[0]), p.ranges.T[1]);
      else state.T = p.ranges.T_default;
      state.params = defaults();
      Object.keys(urlParams).forEach(function (k) { if (k.indexOf("__") !== 0) state.params[k] = urlParams[k]; });
      drawModels(); drawParams(); drawCompare(); drawAxes(); drawSources();
      schedule(0);
    });
  }
  function schedule(delay) {
    clearTimeout(timer);
    timer = setTimeout(run, delay);
  }
  function request() {
    var p = {};
    activeParams().forEach(function (k) { p[k] = state.params[k]; });
    var req = { components: [state.a, state.b], model: state.model, params: p, hold: state.hold,
      npts: state.hold === "x" ? 31 : 41, compare: state.compare || null };
    if (state.hold === "T") req.T = state.T; else if (state.hold === "P") req.P = state.P; else req.x1 = state.x1;
    return req;
  }
  function run() {
    if (busy) { dirty = true; return; }
    busy = true; dirty = false;
    writeHash();
    $("timing").textContent = "Calculating…"; $("timing").className = "busy";
    var req = request(), view = { x: state.x, y: state.y, hold: state.hold };
    engine.call({ op: "compute", request: req }).then(function (res) {
      busy = false;
      if (dirty) { run(); if (!res.ok) return; }
      show(res, req, view);
    });
  }

  /* ------------------------------------------------------------ drawing */
  function overlay(title, text, isError) {
    $("overlay").hidden = false;
    $("spinner").style.display = isError ? "none" : "";
    $("overlay-title").textContent = title; $("overlay-text").textContent = text || "";
  }
  function show(res, req, view) {
    $("input-text").textContent = res.input_text || "";
    $("report-text").textContent = res.report || "";
    if (!res.ok) {
      overlay("This case could not be calculated", res.error || "", true);
      $("timing").textContent = ""; $("route").textContent = "";
      $("facts").innerHTML = ""; $("notice").hidden = true;      // nothing from the previous case stays on screen
      return;
    }
    $("overlay").hidden = true;
    last = { res: res, req: req, view: view };
    plot(res, view);
    $("about-plot").href = "concepts.html#" + aboutAnchor(view);
    $("route").textContent = res.route.charAt(0).toUpperCase() + res.route.slice(1);
    $("timing").className = "";
    $("timing").textContent = res.columns.x1.length + " points · " + res.ms + " ms · " +
      (engine.local ? "local Python engine" : "Python running in your browser");
    drawFacts(res);
    var notes = (res.notes || []).filter(function (n) { return n.indexOf("azeotrope at") !== 0; });
    $("notice").hidden = !notes.length; $("notice").textContent = notes.join(" ");
    drawTable(res);
  }
  /* the section of the Concepts page that explains the current plot */
  function aboutAnchor(view) {
    var Y = view.y;
    if (!isComp(view.x)) return "pt";
    if (Y === "gamma" || Y === "lngamma" || Y === "GE") return "activity";
    if (Y === "f" || Y === "phiL" || Y === "phiV") return "fugacity";
    if (Y === "K" || Y === "alpha") return "kvalues";
    if (Y === "Gmix" || Y === "stab") return "two-liquids";
    return "diagrams";
  }
  function axisTitle(v) {
    var a = lower(nameA());
    return { xy: "Mole fraction of " + a + " (liquid x, vapor y)", x1: "Liquid mole fraction of " + a, y1: "Vapor mole fraction of " + a,
      T: "Temperature (K)", P: "Pressure (bar)", K: "K-value, y/x", alpha: "Relative volatility of " + a + " to " + lower(nameB()),
      gamma: "Activity coefficient", lngamma: "ln γ", f: "Fugacity (bar)", phiL: "Liquid fugacity coefficient",
      phiV: "Vapor fugacity coefficient", GE: "Excess Gibbs energy, Gᴱ/RT", Gmix: "Gibbs energy of mixing, ΔG/RT",
      stab: "Stability function, x₁x₂ ∂²(G/RT)/∂x₁²" }[v];
  }
  function seriesFor(cols, view, hold, suffix, gray) {
    var out = [], idx = cols.x1.map(function (_, i) { return i; });
    function pts(xc, yc, tf) {
      return idx.map(function (i) {
        var x = xc[i], y = yc[i];
        if (x == null || y == null) return null;
        if (tf) { if (y <= 0) return null; y = Math.log(y); }
        return [x, y, i];
      });
    }
    function add(name, color, p) { out.push({ name: name + suffix, color: gray ? "--muted" : color, pts: p, width: gray ? 1.5 : 2 }); }
    var n1 = lower(nameA()), n2 = lower(nameB()), Y = view.y, X = view.x;
    if (hold === "x" && ((X === "T" && Y === "P") || (X === "P" && Y === "T"))) {
      if (X === "T") { add("Bubble curve", "--liquid", pts(cols.T, cols.P)); add("Dew curve", "--vapor", pts(cols.T, cols.P_dew)); }
      else { add("Bubble curve", "--liquid", pts(cols.P, cols.T)); add("Dew curve", "--vapor", pts(cols.P_dew, cols.T)); }
      out.lens = true;
      return out;
    }
    if (X === "xy") {
      add("Bubble curve (liquid)", "--liquid", pts(cols.x1, cols[Y]));
      add("Dew curve (vapor)", "--vapor", pts(cols.y1, cols[Y]));
      out.lens = true;
      return out;
    }
    var xc = cols[X];
    if (PAIRED[Y]) {
      var sym = { K: "K", gamma: "γ", lngamma: "ln γ", f: "f", phiL: "φ liquid", phiV: "φ vapor" }[Y];
      add(sym + " of " + n1, "--liquid", pts(xc, cols[PAIRED[Y][0]], Y === "lngamma"));
      add(sym + " of " + n2, "--vapor", pts(xc, cols[PAIRED[Y][1]], Y === "lngamma"));
    } else {
      var nm = { alpha: "Relative volatility", GE: "Gᴱ/RT", Gmix: "ΔG/RT", stab: "Stability function (negative: the liquid splits)", y1: "Vapor against liquid composition",
        P: X === "y1" ? "Dew curve" : "Bubble curve", T: X === "y1" ? "Dew curve" : "Bubble curve" }[Y];
      add(nm, Y === "P" || Y === "T" ? (X === "y1" ? "--vapor" : "--liquid") : "--liquid", pts(xc, cols[SINGLE[Y]]));
    }
    return out;
  }
  function lensArea(series) {
    var a = series[0].pts.filter(Boolean), b = series[1].pts.filter(Boolean);
    return a.concat(b.slice().reverse());
  }
  function plot(res, view) {
    var cols = res.columns, hold = res.hold, X = view.x, Y = view.y;
    var series = seriesFor(cols, view, hold, "", false), items = [];
    var cfg = { series: [], areas: [], points: [], hlines: [], xTitle: axisTitle(X), yTitle: axisTitle(Y), maxHeight: 460, ratio: 0.6,
      aria: axisTitle(Y) + " against " + axisTitle(X) + " for " + nameA() + " and " + nameB() };
    if (res.compare) {
      var cs = seriesFor(res.compare.columns, view, hold, "", true);
      cs.forEach(function (s) { s.name += " · " + MODEL_SHORT[res.compare.model]; cfg.series.push(s); });
    }
    /* ideal-solution straight lines for the fugacity plot at constant temperature */
    if (Y === "f" && X === "x1" && hold === "T" && cols.x1[0] < 1e-6 && cols.x1[cols.x1.length - 1] > 1 - 1e-6) {
      var f1p = cols.f1[cols.f1.length - 1], f2p = cols.f2[0];
      cfg.series.push({ name: "Ideal solution", color: "--muted", width: 1, noHover: true, pts: [[0, 0], [1, f1p]] });
      cfg.series.push({ name: "Ideal solution", color: "--muted", width: 1, noHover: true, pts: [[0, f2p], [1, 0]] });
      items.push({ name: "Ideal solution, f = x · f of the pure liquid", color: "--muted" });
    }
    series.forEach(function (s) { cfg.series.push(s); items.push({ name: s.name, color: s.color }); });
    if (series.lens) {
      cfg.areas.push({ color: "--liquid", poly: lensArea(series) });
      items.push({ name: "Two-phase region", color: "--liquid", area: true });
    }
    if (res.compare) items.push({ name: MODEL_SHORT[res.compare.model] + " for comparison", color: "--muted" });
    if (isComp(X)) cfg.xDomain = [0, 1];
    if (Y === "y1" && isComp(X)) { cfg.yDomain = [0, 1]; cfg.diagonal = true; cfg.ratio = 0.7; }
    if (Y === "gamma" || Y === "K" || Y === "alpha" || Y === "phiL" || Y === "phiV") { cfg.yOpts = { include: 1 }; cfg.hlines.push({ y: 1 }); }
    if (Y === "lngamma" || Y === "GE" || Y === "f") cfg.yOpts = { include: 0 };
    if (Y === "lngamma" || Y === "GE" || Y === "stab") cfg.hlines.push({ y: 0 });
    if (Y === "stab") cfg.yOpts = { include: 0 };
    var tp = res.three_phase;
    if (tp && isComp(X) && (Y === "P" || Y === "T" || Y === "y1")) {
      var lv = Y === "y1" ? tp.y : (tp.P != null ? tp.P : tp.T);
      if (X !== "y1") cfg.series.push({ name: "Two liquids and vapor", color: "--third", width: 3, noHover: true, pts: [[tp.x_alpha, lv], [tp.x_beta, lv]] });
      if (Y !== "y1") cfg.points.push({ x: tp.y, y: lv, color: "--third" });
      items.push({ name: "Two liquids and vapor coexist", color: "--third" });
    }
    var az = res.azeotrope;
    if (az && isComp(X) && (Y === "P" || Y === "T" || Y === "y1")) {
      cfg.points.push({ x: az.x1, y: Y === "y1" ? az.x1 : (az.P != null ? az.P : az.T), color: "--ink-2" });
      items.push({ name: "Azeotrope", color: "--ink-2", dot: true });
    }
    cfg.tip = function (s, p) {
      var i = p[2], c = (s.color === "--muted" && res.compare) ? res.compare.columns : cols;
      if (i == null) return "";
      var val = PEC.sig(p[1], 5), unit = Y === "P" || Y === "f" ? " bar" : Y === "T" ? " K" : "";
      return "<b>" + val + unit + "</b><br><span class='dim'>x " + PEC.fmt(c.x1[i], 3) + " · y " + PEC.fmt(c.y1[i], 3) +
        "<br>" + PEC.fmt(c.T[i], 2) + " K · " + PEC.sig(c.P[i], 4) + " bar</span>";
    };
    PEC.legend("legend", items);
    PEC.drawChart($("chart"), cfg);
  }
  function drawFacts(res) {
    var c = res.columns, n = c.x1.length, out = [], a = lower(nameA()), b = lower(nameB());
    var full = c.x1[0] < 1e-6 && c.x1[n - 1] > 1 - 1e-6;
    function add(k, v) { out.push("<li><span class='k'>" + k + "</span><span class='v'>" + v + "</span></li>"); }
    if (res.azeotrope) {
      add("Azeotrope", "x = " + PEC.fmt(res.azeotrope.x1, 3) + " at " + (res.azeotrope.P != null ? PEC.sig(res.azeotrope.P, 4) + " bar" : PEC.fmt(res.azeotrope.T, 2) + " K"));
    } else if (res.hold !== "x" && full && !res.three_phase) add("Azeotrope", "none");
    if (res.hold === "T" && full) {
      add("Vapor pressure of " + a, PEC.sig(c.P[n - 1], 4) + " bar"); add("Vapor pressure of " + b, PEC.sig(c.P[0], 4) + " bar");
    } else if (res.hold === "P" && full) {
      add("Boiling point of " + a, PEC.fmt(c.T[n - 1], 2) + " K"); add("Boiling point of " + b, PEC.fmt(c.T[0], 2) + " K");
    }
    if (res.hold !== "x" && full && c.gamma1[0] != null && c.gamma2[n - 1] != null) {
      add("γ∞ of " + a + " in " + b, PEC.sig(c.gamma1[0], 4)); add("γ∞ of " + b + " in " + a, PEC.sig(c.gamma2[n - 1], 4));
    }
    if (res.hold === "x") add("Composition", PEC.fmt(res.fixed.x1, 3) + " " + a);
    if (res.unstable) add("Liquid unstable", "x from " + PEC.fmt(res.unstable[0], 2) + " to " + PEC.fmt(res.unstable[1], 2));
    if (res.three_phase) {
      var t3 = res.three_phase;
      add("Two liquid phases", "x = " + PEC.fmt(t3.x_alpha, 3) + " and " + PEC.fmt(t3.x_beta, 3));
      add(t3.heteroazeotrope ? "Heterogeneous azeotrope" : "Vapor with both liquids",
        "y = " + PEC.fmt(t3.y, 3) + " at " + (t3.P != null ? PEC.sig(t3.P, 4) + " bar" : PEC.fmt(t3.T, 2) + " K"));
    }
    $("facts").innerHTML = out.join("");
  }
  var TABLE_COLS = [["x1", "x₁", 4], ["y1", "y₁", 4], ["T", "T (K)", 2], ["P", "P (bar)", 4], ["P_dew", "P dew (bar)", 4], ["K1", "K₁", 4], ["K2", "K₂", 4],
    ["alpha", "α₁₂", 4], ["gamma1", "γ₁", 4], ["gamma2", "γ₂", 4], ["f1", "f₁ (bar)", 4], ["f2", "f₂ (bar)", 4],
    ["phiL1", "φ₁ liquid", 4], ["phiL2", "φ₂ liquid", 4], ["phiV1", "φ₁ vapor", 4], ["phiV2", "φ₂ vapor", 4], ["GE_RT", "Gᴱ/RT", 4], ["Gmix_RT", "ΔG/RT", 4], ["stability", "Stability", 4]];
  function drawTable(res) {
    var c = res.columns, use = TABLE_COLS.filter(function (t) { return c[t[0]]; });
    var rows = c.x1.map(function (_, i) { return use.map(function (t) { return t[0] === "T" ? PEC.fmt(c[t[0]][i], t[2]) : PEC.sig(c[t[0]][i], 5); }); });
    PEC.dataTable("data-table", use.map(function (t) { return t[1]; }), rows);
    $("data-caption").textContent = "Subscript 1 is " + lower(nameA()) + ", 2 is " + lower(nameB()) + ". " +
      (res.hold === "x" ? "Values along the bubble curve at fixed liquid composition." : "Values along the bubble and dew curves.");
  }
  function downloadCsv() {
    if (!last) return;
    var c = last.res.columns, use = TABLE_COLS.filter(function (t) { return c[t[0]]; });
    var lines = ["# " + nameA() + " (1) / " + nameB() + " (2); model " + last.req.model + "; T in K, P and f in bar",
      use.map(function (t) { return t[0]; }).join(",")];
    c.x1.forEach(function (_, i) { lines.push(use.map(function (t) { var v = c[t[0]][i]; return v == null ? "" : v; }).join(",")); });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([lines.join("\n") + "\n"], { type: "text/csv" }));
    a.download = state.a + "_" + state.b + "_" + state.model + ".csv";
    document.body.appendChild(a); a.click(); a.remove();
  }
  function drawSources() {
    var A = comp(state.a), B = comp(state.b), used = {};
    function src(k) { used[k] = true; return k; }
    var rows = [
      ["Formula", A.formula, B.formula, ""], ["CAS number", A.cas, B.cas, ""],
      ["Molar mass (g/mol)", PEC.fmt(A.M, 3), PEC.fmt(B.M, 3), src(A.sources.M)],
      ["Critical temperature (K)", PEC.fmt(A.Tc, 2), PEC.fmt(B.Tc, 2), src(A.sources.Tc)],
      ["Critical pressure (bar)", PEC.fmt(A.Pc_bar, 2), PEC.fmt(B.Pc_bar, 2), src(A.sources.Pc)],
      ["Acentric factor", PEC.fmt(A.omega, 3), PEC.fmt(B.omega, 3), src(A.sources.omega)],
      ["Antoine A, B, C", [A.antoine.A, A.antoine.B, A.antoine.C].join(", "), [B.antoine.A, B.antoine.B, B.antoine.C].join(", "), src(A.sources.antoine)],
      ["Antoine range (K)", A.antoine.Tmin + " to " + A.antoine.Tmax, B.antoine.Tmin + " to " + B.antoine.Tmax, A.sources.antoine],
      ["Normal boiling point (K)", PEC.fmt(A.Tb, 2), PEC.fmt(B.Tb, 2), src(A.sources.Tb)]];
    var h = "<table><thead><tr><th>Property</th><th>" + A.name + "</th><th>" + B.name + "</th><th>Source</th></tr></thead><tbody>" +
      rows.map(function (r) { return "<tr><td>" + r[0] + "</td><td>" + r[1] + "</td><td>" + r[2] + "</td><td>" + r[3] + "</td></tr>"; }).join("");
    var pr = [];
    var V = pair.vanlaar;
    if (V) {
      var tr = function (a, b) { return PEC.fmt(a, 4) + (b < 0 ? " − " : " + ") + PEC.fmt(Math.abs(b), 1) + "/T"; };
      pr.push(["van Laar A₁₂(T), A₂₁(T), T in K", tr(V.A12 - V.b12 / V.T_ref, V.b12) + ";  " + tr(V.A21 - V.b21 / V.T_ref, V.b21) +
        " (fitted " + PEC.fmt(V.Tmin, 0) + "–" + PEC.fmt(V.Tmax, 0) + " K, rms error in ln γ " + V.rms_lngamma + ")", src("vanlaar_fit")]);
      src("chemsep_nrtl");
    } else pr.push(["van Laar parameters", "not in the database", ""]);
    if (!offered("models", "nrtl")) {
    } else if (pair.nrtl) { pr.push(["NRTL Δg₁₂, Δg₂₁ (J/mol), α", pair.nrtl.dg12 + ", " + pair.nrtl.dg21 + ", " + pair.nrtl.alpha, src(pair.nrtl.source)]); }
    else pr.push(["NRTL parameters", "not in the database", ""]);
    if (pair.kij.source !== "none") pr.push(["Peng–Robinson k₁₂", pair.kij.value, src(pair.kij.source)]);
    else pr.push(["Peng–Robinson k₁₂", "not in the database (0 used)", ""]);
    h += pr.map(function (r) { return "<tr><td>" + r[0] + "</td><td colspan='2'>" + r[1] + "</td><td>" + r[2] + "</td></tr>"; }).join("") + "</tbody></table>";
    $("source-table").innerHTML = h;
    used.libraries = true;
    $("source-notes").innerHTML = "<p class='cite'>Antoine equation: " + A.antoine.form + ".</p>" +
      Object.keys(used).filter(function (k) { return catalog.sources[k]; })
        .map(function (k) { return "<p class='cite'><b>" + k + "</b>: " + catalog.sources[k] + "</p>"; }).join("");
  }

  /* ------------------------------------------------------------ start */
  engine = new Engine(function (text) { overlay(text, "The first visit downloads Python (about 12 MB) and can take up to half a minute. Later visits start in a few seconds."); },
    function (err) { overlay("Python could not start in this browser", "Check your connection and reload the page. (" + err + ")", true); });
  engine.call({ op: "catalog" }).then(function (cat) {
    catalog = cat;
    catalog.compounds = catalog.compounds.filter(function (c) { return offered("compounds", c.id); });
    var urlParams = readHash();
    if (!comp(state.a)) state.a = DEF.a;
    if (!comp(state.b) || state.b === state.a) state.b = state.a === DEF.b ? DEF.a : DEF.b;
    if (!MODEL_SHORT[state.model] || !offered("models", state.model)) state.model = DEF.model;
    if (state.compare && !offered("models", state.compare)) state.compare = "";
    if (!Y_OPTS.some(function (o) { return o[0] === state.y; })) { state.y = "T"; state.x = "xy"; }
    if (!X_OPTS.some(function (o) { return o[0] === state.x; })) state.x = "xy";
    fillCompounds(); bindAxes();
    urlParams.__first = true;
    return loadPair(urlParams);
  });
})();
