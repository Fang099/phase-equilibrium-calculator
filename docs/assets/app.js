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

  var state = { a: DEF.a, b: DEF.b, model: DEF.model, compare: "", x: "xy", y: "T", hold: "P", T: 350, P: 1.01325, x1: 0.5, view: "2d", shade: true, params: {},
    /* external field: kind, its parameters, the overall composition z1 at the reference
       position, and which profile is plotted */
    field: "none", g: 9.80665, H: 500, rpm: 20000, ri: 0.05, ro: 0.10, z1: 0.4, fy: "P" };
  var FIELD_PRESETS = [["P", "Pressure"], ["x1", "Composition"], ["rho", "Density"]];
  var LAYERS = 15, LAYER_PTS = 31;                 // 3D surface: diagrams stacked along the fixed variable
  var surfaceCache = { key: null, res: null };
  var catalog = null, pair = null, last = null, engine, busy = false, dirty = false, timer = null;

  /* ------------------------------------------------------------ URL state */
  function readHash() {
    var h = new URLSearchParams(location.hash.slice(1));
    ["a", "b", "model", "compare", "x", "y", "hold", "view"].forEach(function (k) { if (h.get(k) != null) state[k] = h.get(k); });
    ["T", "P", "x1"].forEach(function (k) { if (h.get(k) != null && isFinite(+h.get(k))) state[k] = +h.get(k); });
    if (h.get("shade") === "0") state.shade = false;
    ["field", "fy"].forEach(function (k) { if (h.get(k) != null) state[k] = h.get(k); });
    ["g", "H", "rpm", "ri", "ro", "z1"].forEach(function (k) { if (h.get(k) != null && isFinite(+h.get(k))) state[k] = +h.get(k); });
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
    if (state.view === "3d") h.set("view", "3d");
    if (!state.shade) h.set("shade", "0");
    if (state.field !== "none") {
      h.set("field", state.field); h.set("fy", state.fy); h.set("z1", state.z1); h.set("T", state.T); h.set("P", state.P);
      (state.field === "gravity" ? ["g", "H"] : ["rpm", "ri", "ro"]).forEach(function (k) { h.set(k, state[k]); });
    }
    activeParams().forEach(function (k) { if (state.params[k] !== defaults()[k]) h.set(k, state.params[k]); });
    history.replaceState(null, "", "#" + h.toString());
  }

  /* ------------------------------------------------------------ helpers */
  function comp(id) { return catalog.compounds.filter(function (c) { return c.id === id; })[0]; }
  function nameA() { return comp(state.a).name; }
  function nameB() { return comp(state.b).name; }
  function lower(s) { return s.toLowerCase(); }
  function isComp(x) { return x === "xy" || x === "x1" || x === "y1"; }
  /* a surface stacks Pxy or Txy-type diagrams along the fixed T or P */
  /* liquid, vapor and two-phase regions exist on diagrams of T or P */
  function canShade() {
    if (state.field !== "none") return true;
    return (isComp(state.x) && state.x !== "y1" && (state.y === "T" || state.y === "P")) ||
      (state.hold === "x" && ((state.x === "T" && state.y === "P") || (state.x === "P" && state.y === "T")));
  }
  function can3d() { return state.field === "none" && (state.hold === "T" || state.hold === "P") && isComp(state.x); }
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
      function (m) { return !hasData(m); });
    [].slice.call($("models-act").children).forEach(function (b) {
      if (b.disabled) b.title = "The database has no " + b.textContent + " parameters for this pair";
    });
  }
  function drawCompare() {
    var sel = $("compare");
    sel.innerHTML = "<option value=''>Nothing</option>" + catalog.models.filter(function (m) { return m.id !== state.model && offered("models", m.id) && hasData(m.id); })
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
  /* van Laar and NRTL need parameters from the database for the pair */
  function hasData(m) {
    if (m === "vanlaar") return !!pair.vanlaar;
    if (m === "nrtl") return !!pair.nrtl;
    return true;
  }
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
        : " A does not change with temperature.";
    }
    note.innerHTML = (edited ? "Edited. " : "") + "Starting values: " + (SOURCE_LABEL[src] || src) + "." + extra +
      (edited ? " <button class='linkbtn' type='button' id='reset-params'>Reset</button>" : "");
    var b = $("reset-params");
    if (b) b.addEventListener("click", function () { state.params = defaults(); drawParams(); schedule(0); });
  }
  function drawAxes() {
    var fm = state.field !== "none";
    document.querySelector(".axes").hidden = fm;
    $("view").hidden = fm;
    if (fm) {
      seg($("presets"), FIELD_PRESETS, state.fy, function (k) { state.fy = k; drawAxes(); if (last) redraw(); });
      var sh0 = $("shade"); sh0.checked = state.shade; sh0.disabled = false; sh0.parentNode.title = "";
      return;
    }
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
    if (!can3d()) state.view = "2d";
    seg($("view"), [["2d", "Diagram"], ["3d", "3D surface"]], state.view,
      function (v) { state.view = v; drawAxes(); schedule(0); }, function (v) { return v === "3d" && !can3d(); });
    var sh = $("shade");
    sh.checked = state.shade; sh.disabled = !canShade();
    sh.parentNode.title = sh.disabled ? "Phase regions are shaded on temperature and pressure diagrams" : "";
    var b3 = $("view").children[1];
    if (b3 && b3.disabled) b3.title = "3D needs a composition axis and a fixed temperature or pressure";
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
  /* ------------------------------------------------------------ is there equilibrium here? */
  var check = { mode: "point", mark: null };
  function phaseWord(p) { return { liquid: "Liquid only", vapor: "Vapor only", "two-phase": "Liquid and vapor in equilibrium" }[p] || p; }
  function drawCheckMode() {
    seg($("check-mode"), [["point", "A state of this mixture"], ["file", "Your own input"]], check.mode, function (m) {
      check.mode = m; $("check-point").hidden = m !== "point"; $("check-file").hidden = m !== "file"; drawCheckMode();
      if (m === "file" && !$("ck-text").value && last) $("ck-text").value = last.res.input_text || "";
    });
  }
  function fillCheckDefaults() {
    if (!$("ck-T").value) $("ck-T").value = state.T;
    if (!$("ck-P").value) $("ck-P").value = +(+state.P).toPrecision(5);
    if (!$("ck-z").value) $("ck-z").value = 0.5;
    $("ck-z-label").textContent = "Overall mole fraction of " + lower(nameA());
  }
  function runPoint() {
    var T = +$("ck-T").value, P = +$("ck-P").value, z = +$("ck-z").value, box = $("ck-result");
    if (!(T > 0) || !(P > 0) || !(z > 0 && z < 1)) { box.innerHTML = "<p class='notice'>Give T > 0, P > 0 and a mole fraction between 0 and 1.</p>"; return; }
    var p = {}; activeParams().forEach(function (k) { p[k] = state.params[k]; });
    box.innerHTML = "<p class='source'>Calculating…</p>";
    engine.call({ op: "point", request: { components: [state.a, state.b], model: state.model, params: p, T: T, P: P, z1: z } }).then(function (r) {
      if (!r.ok) { box.innerHTML = "<p class='notice'>" + (r.error || "This state could not be calculated.") + "</p>"; check.mark = null; return; }
      var a = lower(nameA()), b = r.bounds || {}, f = function (v, d) { return v == null ? "–" : PEC.fmt(v, d); };
      var s = function (v) { return v == null ? "–" : PEC.sig(v, 4); };
      var h = "<p class='check-verdict ph-" + r.phase + "'>" + phaseWord(r.phase) + "</p>";
      if (r.phase === "two-phase") {
        h += "<p>Yes: at " + f(T, 2) + " K and " + s(P) + " bar this mixture splits into a liquid with x = " + f(r.x1, 4) + " and a vapor with y = " + f(r.y1, 4) +
          " (mole fractions of " + a + "), " + f(r.beta * 100, 1) + "% of it vapor. The two compositions are the ends of the tie line through your state.</p>";
      } else {
        h += "<p>No two-phase equilibrium at this state: the mixture is all " + (r.phase === "liquid" ? "liquid (the pressure is above its bubble pressure)" : "vapor (the pressure is below its dew pressure)") + ".</p>";
      }
      h += "<p><b>Where equilibrium occurs for this composition:</b> at " + f(T, 2) + " K, between the dew pressure " + s(b.P_dew) + " bar and the bubble pressure " + s(b.P_bubble) +
        " bar; at " + s(P) + " bar, between the bubble temperature " + f(b.T_bubble, 2) + " K and the dew temperature " + f(b.T_dew, 2) + " K." +
        (b.P_dew != null && b.P_bubble != null && Math.abs(b.P_bubble - b.P_dew) < 1e-3 * b.P_bubble ? " The two coincide: this composition is (very nearly) the azeotrope." : "") + "</p>";
      if (r.K) h += "<p class='source'>K = " + r.K.map(s).join(", ") + (r.gamma ? " · γ = " + r.gamma.map(s).join(", ") : "") + " · " + r.ms + " ms</p>";
      h += "<p><button type='button' class='linkbtn' id='ck-show'>Show it on the Txy diagram at " + s(P) + " bar</button>" +
        " · <button type='button' class='linkbtn' id='ck-show-p'>on the Pxy diagram at " + f(T, 1) + " K</button></p>";
      box.innerHTML = h;
      check.mark = { T: T, P: P, z1: z, x1: r.x1, y1: r.y1, phase: r.phase };
      function go(hold) {
        state.field = "none"; state.view = "2d"; state.hold = hold; state.x = "xy"; state.y = hold === "P" ? "T" : "P";
        if (hold === "P") state.P = P; else state.T = T;
        drawField(); drawAxes(); schedule(0);
        document.querySelector(".workspace").scrollIntoView({ behavior: "smooth", block: "start" });
      }
      $("ck-show").addEventListener("click", function () { go("P"); });
      $("ck-show-p").addEventListener("click", function () { go("T"); });
      if (last) redraw();
    });
  }
  function runFile() {
    var box = $("ck-file-result"), text = $("ck-text").value;
    if (!text.trim()) { box.innerHTML = "<p class='notice'>Write an input first, or start from the current calculation.</p>"; return; }
    box.innerHTML = "<p class='source'>Checking and solving…</p>";
    engine.call({ op: "run_input", text: text }).then(function (r) {
      if (!r.ok) { box.innerHTML = "<p class='notice'>" + (r.error || "Failed") + "</p>"; return; }
      var esc = function (t) { return t.replace(/&/g, "&amp;").replace(/</g, "&lt;"); };
      box.innerHTML = "<p class='check-verdict st-" + r.status.toLowerCase() + "'>" + r.status + "</p>" +
        "<pre>" + esc(r.report) + "</pre>" + (r.results ? "<pre>" + esc(r.results) + "</pre>" : "") +
        "<p class='source'>" + r.ms + " ms</p>";
    });
  }
  function bindCheck() {
    drawCheckMode();
    $("ck-run").addEventListener("click", runPoint);
    ["ck-T", "ck-P", "ck-z"].forEach(function (id) { $(id).addEventListener("keydown", function (e) { if (e.key === "Enter") runPoint(); }); });
    $("ck-file-run").addEventListener("click", runFile);
    $("ck-file-fill").addEventListener("click", function () { if (last) $("ck-text").value = last.res.input_text || ""; });
  }

  /* the settings column: scrolls on its own beside the plot (a sheet at the bottom
     on small screens), and can be hidden to give the plot more room */
  function bindControls() {
    var shell = $("main"), KEY = "pec-controls-hidden";
    function set(hidden) {
      shell.classList.toggle("controls-hidden", hidden);
      $("toggle-controls").textContent = hidden ? "Show" : "Hide";
      $("toggle-controls").setAttribute("aria-expanded", String(!hidden));
      $("show-controls").hidden = !hidden;
      try { localStorage.setItem(KEY, hidden ? "1" : "0"); } catch (e) { /* not remembered */ }
    }
    var start = false;
    try { start = localStorage.getItem(KEY) === "1"; } catch (e) { /* default: shown */ }
    set(start);
    $("toggle-controls").addEventListener("click", function () { set(!shell.classList.contains("controls-hidden")); });
    $("show-controls").addEventListener("click", function () { set(false); });
  }
  function bindAxes() {
    bindControls();
    bindCheck();
    $("shade").addEventListener("change", function () {          // a redraw, not a new calculation
      state.shade = this.checked; writeHash();
      if (last) show(last.res, last.req, last.view, last.sres);
    });
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
      if (!hasData(state.model)) state.model = "PR";      // no database values for this pair
      if (state.compare && !hasData(state.compare)) state.compare = "";
      if (urlParams.__hasT) state.T = Math.min(Math.max(state.T, p.ranges.T[0]), p.ranges.T[1]);
      else state.T = p.ranges.T_default;
      state.params = defaults();
      Object.keys(urlParams).forEach(function (k) { if (k.indexOf("__") !== 0) state.params[k] = urlParams[k]; });
      drawModels(); drawParams(); drawCompare(); drawField(); drawAxes(); drawSources();
      check.mark = null; $("ck-result").innerHTML = ""; fillCheckDefaults();
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
    if (state.field !== "none")
      return { components: [state.a, state.b], model: state.model, params: p, hold: "field", npts: 61,
        T: state.T, P: state.P, z1: state.z1,
        field: state.field === "gravity" ? { kind: "gravity", g: state.g, height: state.H }
                                         : { kind: "centrifugal", rpm: state.rpm, r_in: state.ri, r_out: state.ro } };
    var req = { components: [state.a, state.b], model: state.model, params: p, hold: state.hold,
      npts: state.hold === "x" ? 41 : 81, compare: state.compare || null };
    if (state.hold === "T") req.T = state.T; else if (state.hold === "P") req.P = state.P; else req.x1 = state.x1;
    return req;
  }
  function run() {
    if (busy) { dirty = true; return; }
    busy = true; dirty = false;
    writeHash();
    $("timing").textContent = "Calculating…"; $("timing").className = "busy";
    var req = request(), view = { x: state.x, y: state.y, hold: state.hold };
    if (state.view === "3d" && can3d()) req.compare = null;
    engine.call({ op: "compute", request: req }).then(function (res) {
      if (!(state.view === "3d" && can3d()) || !res.ok) {
        busy = false;
        if (dirty) { run(); if (!res.ok) return; }
        show(res, req, view);
        return;
      }
      /* 3D: the layers depend on everything except the fixed value, so they are cached */
      var sreq = surfaceRequest(req), key = JSON.stringify(sreq);
      function done(sres) {
        busy = false;
        if (dirty) { run(); return; }
        show(res, req, view, sres);
      }
      if (surfaceCache.key === key) { done(surfaceCache.res); return; }
      sourceVersion.then(function (ver) {
        var kept = ver && storeGet(key, ver);
        if (kept) { kept.source = "remembered"; surfaceCache.key = key; surfaceCache.res = kept; done(kept); return; }
        precomputed(sreq).then(function (pre) {
          if (pre) { surfaceCache.key = key; surfaceCache.res = pre; done(pre); return; }
          $("timing").textContent = "Calculating " + LAYERS + " diagrams for the surface…";
          overlay("Calculating the surface", LAYERS + " diagrams of " + LAYER_PTS + " points each. This takes a few seconds, once: this browser remembers the result.", false);
          engine.call({ op: "surface", request: sreq }).then(function (sres) {
            if (sres.ok) {
              surfaceCache.key = key; surfaceCache.res = sres;
              if (ver) storePut(key, ver, sres);
            }
            done(sres);
          });
        });
      });
    });
  }

  function surfaceRequest(req) {
    var lo, hi, vals = [], k, r = {};
    Object.keys(req).forEach(function (k2) { if (k2 !== "T" && k2 !== "P" && k2 !== "compare") r[k2] = req[k2]; });
    if (state.hold === "T") { lo = pair.ranges.T[0]; hi = pair.ranges.T[1]; }
    else { lo = Math.log10(pair.ranges.P[0]); hi = Math.log10(pair.ranges.P[1]); }
    for (k = 0; k < LAYERS; k++) {
      var u = lo + (hi - lo) * k / (LAYERS - 1);
      vals.push(state.hold === "T" ? +u.toFixed(2) : +Math.pow(10, u).toPrecision(5));
    }
    r.values = vals; r.npts = LAYER_PTS;
    /* only the columns the 3D plots use, to keep the result small */
    var cols = ["x1", "y1", "T", "P"], need = { gamma: ["gamma1", "gamma2"], lngamma: ["gamma1", "gamma2"], K: ["K1", "K2"],
      f: ["f1", "f2"], phiL: ["phiL1", "phiL2"], phiV: ["phiV1", "phiV2"], alpha: ["alpha"], GE: ["GE_RT"], Gmix: ["Gmix_RT"], stab: ["stability"] };
    Y_OPTS.forEach(function (o) { (need[o[0]] || []).forEach(function (c) { if (cols.indexOf(c) < 0) cols.push(c); }); });
    r.keep = cols;
    return r;
  }

  /* Surfaces are slow to compute in the browser, so the page looks in three places
     first: this visit's last surface, surfaces this browser computed before, and the
     precomputed files in surfaces/ (scripts/precompute_surfaces.py, for the database
     parameters). All must belong to the current version of the Python code. */
  var sourceVersion = fetch("py/manifest.json", { cache: "no-cache" }).then(function (r) { return r.json(); })
    .then(function (m) { return m.source || null; }).catch(function () { return null; });
  var surfaceIndex = sourceVersion.then(function (v) {
    if (!v) return null;
    return fetch("surfaces/index.json", { cache: "no-cache" }).then(function (r) { return r.ok ? r.json() : null; })
      .then(function (ix) { return ix && ix.version === v && ix.layers === LAYERS && ix.points === LAYER_PTS ? ix : null; })
      .catch(function () { return null; });
  });
  function sameNumbers(a, b, tol) {
    var ka = Object.keys(a || {}), kb = Object.keys(b || {});
    if (ka.length !== kb.length) return false;
    return ka.every(function (k) { var x = +a[k], y = +b[k]; return isFinite(x) && isFinite(y) && Math.abs(x - y) <= tol * Math.max(1, Math.abs(x)); });
  }
  function precomputed(sreq) {
    return surfaceIndex.then(function (ix) {
      if (!ix) return null;
      var key = [sreq.components[0], sreq.components[1], sreq.model, sreq.hold].join("~");
      if (ix.keys.indexOf(key) < 0) return null;
      return fetch("surfaces/" + encodeURIComponent(key) + ".json?v=" + ix.version).then(function (r) { return r.ok ? r.json() : null; })
        .then(function (f) {
          if (!f || f.version !== ix.version || !f.result || !f.result.ok) return null;
          var q = f.request;
          var ok = q.model === sreq.model && q.hold === sreq.hold && q.npts === sreq.npts &&
            q.components.join() === sreq.components.join() && sameNumbers(q.params, sreq.params, 1e-9) &&
            sameNumbers(q.values, sreq.values, 1e-6) && sreq.keep.every(function (c) { return q.keep.indexOf(c) >= 0; });
          if (ok) f.result.source = "precomputed";
          return ok ? f.result : null;
        });
    }).catch(function () { return null; });
  }
  var STORE = "pec-surface:";
  function storeGet(key, version) {
    try {
      var v = localStorage.getItem(STORE + key);
      if (!v) return null;
      v = JSON.parse(v);
      return v.version === version ? v.res : null;
    } catch (e) { return null; }
  }
  function storePut(key, version, res) {
    var item = JSON.stringify({ version: version, res: res });
    try { localStorage.setItem(STORE + key, item); return; } catch (e) { /* full: make room */ }
    try {
      Object.keys(localStorage).filter(function (k) { return k.indexOf(STORE) === 0; })
        .forEach(function (k) { localStorage.removeItem(k); });
      localStorage.setItem(STORE + key, item);
    } catch (e) { /* storage unavailable: nothing is remembered */ }
  }

  /* ------------------------------------------------------------ drawing */
  function overlay(title, text, isError) {
    $("overlay").hidden = false;
    $("spinner").style.display = isError ? "none" : "";
    $("overlay-title").textContent = title; $("overlay-text").textContent = text || "";
  }
  function show(res, req, view, sres) {
    $("input-text").textContent = res.input_text || "";
    $("report-text").textContent = res.report || "";
    if (!res.ok) {
      overlay("This case could not be calculated", res.error || "", true);
      $("timing").textContent = ""; $("route").textContent = "";
      $("facts").innerHTML = ""; $("notice").hidden = true;      // nothing from the previous case stays on screen
      return;
    }
    $("overlay").hidden = true;
    last = { res: res, req: req, view: view, sres: sres };
    var is3d = !!sres;
    $("chart").classList.toggle("is-3d", is3d);
    $("chart-hint").hidden = !is3d;
    var isField = res.hold === "field";
    if (is3d && !sres.ok) { overlay("The surface could not be calculated", sres.error || "", true); }
    else if (is3d) plot3d(res, view, sres);
    else if (isField) plotField(res);
    else plot(res, view);
    $("about-plot").href = "concepts.html#" + (isField ? "beyond" : aboutAnchor(view));
    $("route").textContent = res.route.charAt(0).toUpperCase() + res.route.slice(1);
    $("timing").className = "";
    $("timing").textContent = (is3d && sres.ok
      ? sres.layers.length + " diagrams × " + LAYER_PTS + " points · " +
        (sres.source === "precomputed" ? "precomputed" : sres.source === "remembered" ? "remembered from an earlier visit" : sres.ms + " ms")
      : (isField ? res.columns.position.length : res.columns.x1.length) + " points · " + res.ms + " ms") + " · " +
      (engine.local ? "local Python engine" : "Python running in your browser");
    if (isField) drawFieldFacts(res); else drawFacts(res);
    var notes = (res.notes || []).filter(function (n) { return n.indexOf("azeotrope at") !== 0 && n.indexOf("gas-liquid contact") !== 0; });
    $("notice").hidden = !notes.length; $("notice").textContent = notes.join(" ");
    if (isField) drawFieldTable(res); else drawTable(res);
  }
  function redraw() { if (last) show(last.res, last.req, last.view, last.sres); }

  /* ------------------------------------------------------------ external field */
  function sliderRow(host, label, o) {
    var row = document.createElement("div"); row.className = "param";
    row.innerHTML = "<label>" + label + "</label><input type='range' aria-label='" + label + "'><input type='number' aria-label='" + label + "'>";
    host.appendChild(row);
    var r = row.children[1], n = row.children[2];
    r.min = o.log ? Math.log10(o.min) : o.min; r.max = o.log ? Math.log10(o.max) : o.max;
    r.step = o.log ? 0.005 : o.step; n.step = o.step;
    function show(v) { r.value = o.log ? Math.log10(v) : v; n.value = +(+v).toPrecision(5); }
    show(o.value);
    r.addEventListener("input", function () {
      var v = o.log ? +Math.pow(10, +r.value).toPrecision(4) : +r.value; n.value = v; o.set(v); schedule(150);
    });
    n.addEventListener("change", function () {
      var v = +n.value; if (!isFinite(v) || v < o.min * (o.min > 0 ? 0.999 : 1) || v > o.max * 1.001) { show(o.value); return; }
      show(v); o.set(v); schedule(0);
    });
  }
  function drawField() {
    seg($("field-kind"), [["none", "None"], ["gravity", "Gravity"], ["centrifugal", "Centrifuge"]], state.field, function (k) {
      state.field = k; if (k !== "none") state.view = "2d"; drawField(); drawAxes(); schedule(0);
    });
    var host = $("field-params"); host.innerHTML = "";
    if (state.field === "none") { $("field-note").hidden = false; return; }
    $("field-note").hidden = true;
    if (state.field === "gravity") {
      sliderRow(host, "g (m/s²)", { min: 0.1, max: 100, step: 0.01, log: true, value: state.g, set: function (v) { state.g = v; } });
      sliderRow(host, "Height of the column (m)", { min: 1, max: 5000, step: 1, log: true, value: state.H, set: function (v) { state.H = v; } });
    } else {
      sliderRow(host, "Speed (rpm)", { min: 100, max: 100000, step: 10, log: true, value: state.rpm, set: function (v) { state.rpm = v; } });
      sliderRow(host, "Inner radius (cm)", { min: 0, max: 50, step: 0.1, value: +(state.ri * 100).toFixed(2), set: function (v) { state.ri = v / 100; } });
      sliderRow(host, "Outer radius (cm)", { min: 1, max: 100, step: 0.1, value: +(state.ro * 100).toFixed(2), set: function (v) { state.ro = v / 100; } });
    }
    var h = document.createElement("p"); h.className = "label"; h.style.margin = "0.9rem 0 0.2rem";
    h.textContent = state.field === "gravity" ? "State at the bottom" : "State at the inner radius";
    host.appendChild(h);
    sliderRow(host, "Temperature (K), uniform", { min: pair.ranges.T[0], max: pair.ranges.T[1], step: 1, value: state.T, set: function (v) { state.T = v; } });
    sliderRow(host, "Pressure (bar)", { min: pair.ranges.P[0], max: pair.ranges.P[1] * 10, step: 0.01, log: true, value: state.P, set: function (v) { state.P = v; } });
    sliderRow(host, "Overall mole fraction of " + lower(nameA()), { min: 0.01, max: 0.99, step: 0.01, value: state.z1, set: function (v) { state.z1 = v; } });
  }
  function fieldPos(res) {                                   // m for gravity, cm for a centrifuge
    var g = res.field.kind === "gravity";
    return { g: g, k: g ? 1 : 100, title: g ? "Height (m)" : "Radius (cm)" };
  }
  function plotField(res) {
    var c = res.columns, F = fieldPos(res), Y = state.fy, n = c.position.length;
    var qTitle = { P: "Pressure (bar)", x1: "x of " + lower(nameA()) + " (local phase)", rho: "Density (kg/m³)" }[Y];
    var cfg = { series: [], areas: [], points: [], hlines: [], vlines: [], regions: [], maxHeight: 460, ratio: 0.6,
      aria: qTitle + " along the " + (F.g ? "height of the column" : "radius of the centrifuge") };
    /* one line per stretch of the same phase, coloured by phase */
    var runs = [], cur = null;
    for (var i = 0; i < n; i++) {
      var ph = c.phase[i], v = c[Y][i], pos = c.position[i] * F.k;
      if (v == null) { cur = null; continue; }
      if (!cur || cur.phase !== ph) { cur = { phase: ph, pts: [] }; runs.push(cur); }
      cur.pts.push(F.g ? [v, pos, i] : [pos, v, i]);
    }
    /* join the two phases at each contact: the pressure is continuous there,
       the composition and density jump */
    (res.contacts || []).forEach(function (ct, j) {
      var a = runs[j], b = runs[j + 1], p = ct.position * F.k;
      function val(ph) { return Y === "P" ? ct.P : Y === "x1" ? (ph === "vapor" ? ct.x_vapor : ct.x_liquid) : null; }
      [[a, "end"], [b, "start"]].forEach(function (rr) {
        var r = rr[0]; if (!r) return;
        var v = val(r.phase); if (v == null) return;
        var pt = F.g ? [v, p, null] : [p, v, null];
        if (rr[1] === "end") r.pts.push(pt); else r.pts.unshift(pt);
      });
    });
    var names = { liquid: "Liquid", vapor: "Vapor" };
    runs.forEach(function (r) {
      cfg.series.push({ name: names[r.phase] || r.phase, color: r.phase === "vapor" ? "--vapor" : "--liquid", pts: r.pts, width: 2.5 });
    });
    var posD = [c.position[0] * F.k, c.position[n - 1] * F.k];
    if (Y === "P" && c.P_sat) {                    // the pressure at which the local fluid would split
      var sp = [];
      c.P_sat.forEach(function (v, i) { if (v != null) sp.push(F.g ? [v, c.position[i] * F.k, i] : [c.position[i] * F.k, v, i]); });
      if (sp.length > 1) cfg.series.push({ name: "Saturation pressure of the local fluid", color: "--muted", width: 1.5, dash: "5 4", pts: sp });
    }
    var qs = []; cfg.series.forEach(function (r) { r.pts.forEach(function (p) { qs.push(F.g ? p[0] : p[1]); }); });
    var qD = PEC.niceDomain(qs, Y === "x1" ? { pad: 0.08 } : {});
    if (F.g) { cfg.xDomain = qD; cfg.yDomain = posD; cfg.xTitle = qTitle; cfg.yTitle = F.title; }
    else { cfg.xDomain = posD; cfg.yDomain = qD; cfg.xTitle = F.title; cfg.yTitle = qTitle; }
    var items = [];
    if (runs.some(function (r) { return r.phase !== "vapor"; })) items.push({ name: "Liquid", color: "--liquid" });
    if (runs.some(function (r) { return r.phase === "vapor"; })) items.push({ name: "Vapor", color: "--vapor" });
    if (Y === "P" && cfg.series.some(function (s) { return s.dash; })) items.push({ name: "Bubble P of the local liquid / dew P of the local vapor", color: "--muted" });
    /* phase bands along the field direction */
    if (state.shade) {
      var bands = [], b0 = null;
      for (var k = 0; k < n; k++) {
        var p0 = c.position[k] * F.k;
        if (!b0 || b0.phase !== c.phase[k]) { if (b0) b0.end = p0; b0 = { phase: c.phase[k], start: p0 }; bands.push(b0); }
      }
      if (b0) b0.end = posD[1];
      (res.contacts || []).forEach(function (ct, j) {             // band edges at the exact contacts
        if (bands[j + 1]) { bands[j].end = ct.position * F.k; bands[j + 1].start = ct.position * F.k; }
      });
      bands.forEach(function (bd) {
        var col = bd.phase === "vapor" ? "--vapor" : "--liquid";
        var poly = F.g ? [[qD[0], bd.start], [qD[1], bd.start], [qD[1], bd.end], [qD[0], bd.end]]
                       : [[bd.start, qD[0]], [bd.end, qD[0]], [bd.end, qD[1]], [bd.start, qD[1]]];
        cfg.areas.push({ color: col, opacity: 0.13, poly: poly });
      });
      items.push({ name: "Liquid region", color: "--liquid", area: true });
      if (bands.some(function (bd) { return bd.phase === "vapor"; })) items.push({ name: "Vapor region", color: "--vapor", area: true });
    }
    (res.contacts || []).forEach(function (ct) {
      var p = ct.position * F.k, t = "contact, " + PEC.sig(p, 4) + (F.g ? " m" : " cm");
      if (F.g) cfg.hlines.push({ y: p, text: t }); else cfg.vlines.push({ x: p, text: t });
    });
    cfg.tip = function (s, p) {
      var i = p[2];
      if (i == null) return "";
      if (s.dash) return "<b>" + PEC.sig(c.P_sat[i], 5) + " bar</b><br><span class='dim'>" + (c.phase[i] === "vapor" ? "dew" : "bubble") + " pressure of the fluid here</span>";
      return "<b>" + (F.g ? PEC.sig(c.position[i], 4) + " m" : PEC.sig(c.position[i] * 100, 4) + " cm") + "</b> · " + c.phase[i] +
        "<br><span class='dim'>" + PEC.sig(c.P[i], 5) + " bar · x " + PEC.fmt(c.x1[i], 4) +
        (c.rho[i] != null ? " · " + PEC.sig(c.rho[i], 4) + " kg/m³" : "") + "</span>";
    };
    PEC.legend("legend", items);
    PEC.drawChart($("chart"), cfg);
    $("chart-hint").hidden = false;
    $("chart-hint").textContent = fieldCaption(res, F, Y);
  }
  /* what the profile means, with this case's numbers */
  function fieldCaption(res, F, Y) {
    var c = res.columns, n = c.position.length, along = F.g ? "height" : "radius", a = lower(nameA());
    function grad(phase) {                       // pressure change per 10 m or per cm, within one phase
      var idx = []; for (var i = 0; i < n; i++) if (c.phase[i] === phase && c.P[i] != null) idx.push(i);
      if (idx.length < 2) return null;
      var i0 = idx[0], i1 = idx[idx.length - 1];
      return (c.P[i1] - c.P[i0]) / ((c.position[i1] - c.position[i0]) * F.k) * (F.g ? 10 : 1);
    }
    var unit = F.g ? " bar per 10 m" : " bar per cm", parts = [];
    var gl = grad("liquid"), gv = grad("vapor");
    if (Y === "P") {
      var txt = F.g ? "Each point is the equilibrium state at that height, with T uniform. Pressure falls with height by ρg per metre, so the dense liquid's pressure drops steeply"
                    : "Each point is the equilibrium state at that radius, with T uniform. Pressure rises outward by ρω²r per metre, steeply in the dense liquid";
      if (gl != null) parts.push("liquid " + PEC.sig(gl, 3) + unit);
      if (gv != null) parts.push("vapor " + PEC.sig(gv, 3) + unit);
      return txt + (parts.length ? " (" + parts.join("; ") + ")" : "") + ", while the dilute vapor's pressure hardly changes. " +
        "The dashed line is the pressure at which the local fluid would start to change phase. Where the two meet, the second phase appears: the gas–liquid contact.";
    }
    function span(phase, col, d) {                 // first and last value within one phase
      var v = []; for (var i = 0; i < n; i++) if (c.phase[i] === phase && c[col][i] != null) v.push(c[col][i]);
      return v.length > 1 ? phase + " " + PEC.fmt(v[0], d) + " → " + PEC.fmt(v[v.length - 1], d) : null;
    }
    if (Y === "x1") {
      ["liquid", "vapor"].forEach(function (ph) { var t = span(ph, "x1", 5); if (t) parts.push(t); });
      var within = parts.length && (res.contacts || []).length ? " Within each phase the change is small (" + parts.join("; ") + ", from " + (F.g ? "bottom to top" : "axis side to rim") +
        "), so next to the jump at the contact those lines look " + (F.g ? "vertical" : "flat") + "." : "";
      return xCaption() + within;
    }
    function xCaption() { return "Composition of the phase present at each " + along + ". In a field, ln f + MΦ/RT is the same everywhere for each compound. In the vapor this sends the heavier molecules " +
      (F.g ? "down" : "toward the rim") + "; in the liquid what counts is whether a compound is denser than the mixture around it (its molar mass over the volume it takes up), so a heavy but bulky molecule can move the other way. " +
      "At a contact the composition jumps from the liquid's to the vapor's, as across a tie line of the phase diagram."; }
    var liq = [], vap = [];
    for (var k = 0; k < n; k++) if (c.rho[k] != null) (c.phase[k] === "vapor" ? vap : liq).push(c.rho[k]);
    var jump = liq.length && vap.length ? " At the gas–liquid contact it jumps: the liquid here is about " +
      PEC.sig(Math.min.apply(null, liq) / Math.max.apply(null, vap), 2) + " times denser than the vapor." : "";
    return "Density of the phase present at each " + along + ". It changes slowly within a phase, through compression and the shift in composition." + jump;
  }
  function drawFieldFacts(res) {
    var c = res.columns, n = c.position.length, F = fieldPos(res), out = [], unit = F.g ? " m" : " cm";
    function add(k, v) { out.push("<li><span class='k'>" + k + "</span><span class='v'>" + v + "</span></li>"); }
    var ref = res.reference, where = F.g ? "bottom" : "inner radius";
    add("At the " + where, ref.phase === "two-phase" ? "liquid + vapor, " + PEC.fmt(ref.beta * 100, 0) + "% vapor" : ref.phase);
    (res.contacts || []).forEach(function (ct) { add("Gas–liquid contact", PEC.sig(ct.position * F.k, 4) + unit + ", " + PEC.sig(ct.P, 4) + " bar"); });
    add("Pressure at the " + (F.g ? "top" : "outer radius"), PEC.sig(c.P[n - 1], 4) + " bar");
    var xs = c.x1.filter(function (v) { return v != null; });
    if (xs.length) add("x of " + lower(nameA()) + " from " + where + " to " + (F.g ? "top" : "rim"), PEC.fmt(xs[0], 4) + " → " + PEC.fmt(xs[xs.length - 1], 4));
    $("facts").innerHTML = out.join("");
  }
  function drawFieldTable(res) {
    var c = res.columns, F = fieldPos(res);
    PEC.dataTable("data-table", [F.title, "Phase", "P (bar)", "x₁ (local phase)", "ρ (kg/m³)"],
      c.position.map(function (_, i) {
        return [PEC.sig(c.position[i] * F.k, 5), c.phase[i], PEC.sig(c.P[i], 6), PEC.fmt(c.x1[i], 5), PEC.sig(c.rho[i], 5)]; }));
    $("data-caption").textContent = "Subscript 1 is " + lower(nameA()) + ". The composition is that of the phase present at each position; T is uniform.";
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
  /* 3D: each layer is a diagram at one value of the fixed variable; the surfaces join
     the same curve across layers, matched by liquid composition */
  function plot3d(res, view, sres) {
    var hold = res.hold, X = view.x, Y = view.y, layers = sres.layers.slice().sort(function (a, b) { return a.value - b.value; });
    var base = layers.reduce(function (b, l) { return l.columns.x1.length > b.columns.x1.length ? l : b; }, layers[0]);
    var canon = base.columns.x1.filter(function (v) { return v != null; });
    var first = seriesFor(base.columns, view, hold, "", false);
    var grids = first.map(function () { return []; }), zs = [];
    layers.forEach(function (L) {
      var ser = seriesFor(L.columns, view, hold, "", false), x1 = L.columns.x1;
      first.forEach(function (_, j) {
        grids[j].push(canon.map(function (xc) {
          var i = -1;
          for (var m = 0; m < x1.length; m++) if (x1[m] != null && Math.abs(x1[m] - xc) < 1e-6) { i = m; break; }
          var p = i >= 0 && ser[j] ? ser[j].pts[i] : null;
          if (!p) return null;
          zs.push(p[1]);
          return [p[0], L.value, p[1]];
        }));
      });
    });
    var fixed = hold === "T" ? res.fixed.T : res.fixed.P;
    var slice = seriesFor(res.columns, view, hold, "", false);
    var lines = slice.map(function (s) {
      return { color: "--ink", width: 2.5, pts: s.pts.map(function (p) { if (!p) return null; zs.push(p[1]); return [p[0], fixed, p[1]]; }) };
    });
    var items = first.map(function (s) { return { name: s.name, color: s.color }; });
    items.push({ name: "Diagram at " + (hold === "T" ? PEC.fmt(fixed, 1) + " K" : PEC.sig(fixed, 4) + " bar") + " (slider)", color: "--ink" });
    if (Y === "T" || Y === "P" || Y === "y1") {
      var az = layers.filter(function (l) { return l.azeotrope; }).map(function (l) {
        var a = l.azeotrope; return [a.x1, l.value, Y === "y1" ? a.x1 : (a.P != null ? a.P : a.T)]; });
      if (az.length > 1) { lines.push({ color: "--third", width: 2.5, pts: az }); items.push({ name: "Azeotrope", color: "--third" }); }
    }
    var opts = (Y === "gamma" || Y === "K") ? { include: 1, pad: 0.02 } : { pad: 0.02 };
    var zD = (Y === "y1") ? [0, 1] : PEC.niceDomain(zs, opts), planes = [];
    if (state.shade && first.lens && slice.length > 1) {
      /* the phase regions on the plane of the slider's diagram */
      var Bp = slice[0].pts, Dp = slice[1].pts, liqHigh = Y === "P";
      var zL = liqHigh ? zD[1] : zD[0], zV = liqHigh ? zD[0] : zD[1], qL = [], qV = [], q2 = [];
      var P3 = function (p, z) { return [p[0], fixed, z == null ? p[1] : z]; };
      for (var i = 0; i + 1 < Bp.length; i++) {
        var b0 = Bp[i], b1 = Bp[i + 1], d0 = Dp[i], d1 = Dp[i + 1];
        if (b0 && b1) qL.push([P3(b0), P3(b1), P3(b1, zL), P3(b0, zL)]);
        if (d0 && d1) qV.push([P3(d0), P3(d1), P3(d1, zV), P3(d0, zV)]);
        if (b0 && b1 && d0 && d1) q2.push([P3(b0), P3(d0), P3(d1), P3(b1)]);
      }
      planes = [{ color: "--liquid", alpha: 0.3, quads: qL }, { color: "--vapor", alpha: 0.3, quads: qV },
                { color: "--muted", alpha: 0.35, quads: q2 }];
      items.push({ name: "Liquid", color: "--liquid", area: true }, { name: "Vapor", color: "--vapor", area: true },
        { name: "Liquid + vapor", color: "--muted", area: true });
    }
    var vals = layers.map(function (l) { return l.value; });
    PEC.legend("legend", items);
    PEC.drawSurface($("chart"), {
      surfaces: first.map(function (s, j) { return { name: s.name, color: s.color, grid: grids[j] }; }),
      nested: !!first.lens,
      lines: lines,
      x: { title: X === "xy" ? "x, y of " + lower(nameA()) : axisTitle(X), domain: [0, 1] },
      y: { title: hold === "T" ? "Temperature (K)" : "Pressure (bar)", domain: [Math.min.apply(null, vals), Math.max.apply(null, vals)], log: hold === "P" },
      z: { title: axisTitle(Y), domain: zD },
      planes: planes,
      aria: "3D surface of " + axisTitle(Y) + " against " + axisTitle(X) + " and " + (hold === "T" ? "temperature" : "pressure")
    });
    $("chart-hint").textContent = "Third axis: the fixed " + (hold === "T" ? "temperature, " + PEC.fmt(vals[0], 0) + " to " + PEC.fmt(vals[vals.length - 1], 0) + " K"
      : "pressure, " + PEC.sig(vals[0], 3) + " to " + PEC.sig(vals[vals.length - 1], 3) + " bar") +
      ". Drag to rotate, double-click to reset." + (sres.failed ? " " + sres.failed + " of the diagrams could not be solved and are left out." : "");
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
    var shading = series.lens && state.shade;
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
    /* the state checked below the plot, with its tie line, on a matching diagram */
    var mk = check.mark;
    if (mk && X === "xy" && ((Y === "T" && hold === "P" && Math.abs(res.fixed.P - mk.P) < 1e-6 * mk.P) ||
                             (Y === "P" && hold === "T" && Math.abs(res.fixed.T - mk.T) < 1e-6 * mk.T))) {
      var v = Y === "T" ? mk.T : mk.P;
      if (mk.phase === "two-phase") {
        cfg.series.push({ name: "Tie line", color: "--ink", width: 1.5, dash: "5 4", noHover: true, pts: [[mk.x1, v], [mk.y1, v]] });
        cfg.points.push({ x: mk.x1, y: v, color: "--liquid" }, { x: mk.y1, y: v, color: "--vapor" });
      }
      cfg.points.push({ x: mk.z1, y: v, color: "--ink", text: "your state", dx: 8, dy: -8 });
      items.push({ name: "Your state" + (mk.phase === "two-phase" ? " and its tie line" : ""), color: "--ink", dot: true });
    }
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
    if (shading) {
      /* the same domains drawChart would choose, so the regions can reach the edges */
      var allX = [], allY = [];
      cfg.series.forEach(function (s) { s.pts.forEach(function (p) { if (p) { allX.push(p[0]); allY.push(p[1]); } }); });
      cfg.xDomain = cfg.xDomain || PEC.niceDomain(allX, { pad: 0.02 });
      cfg.yDomain = cfg.yDomain || PEC.niceDomain(allY, cfg.yOpts);
      var rg = phaseRegions(series, X, Y, cfg.xDomain, cfg.yDomain);
      cfg.areas = rg.areas.concat(cfg.areas);
      cfg.regions = (cfg.regions || []).concat(rg.labels);
      items.push({ name: "Liquid", color: "--liquid", area: true }, { name: "Vapor", color: "--vapor", area: true },
        { name: "Liquid + vapor", color: "--muted", area: true });
    }
    PEC.legend("legend", items);
    PEC.drawChart($("chart"), cfg);
  }
  /* liquid, vapor and two-phase regions bounded by the bubble and dew curves.
     The liquid lies at low T or high P; the vapor at high T or low P. */
  function phaseRegions(series, X, Y, xD, yD) {
    var B = series[0].pts.filter(Boolean), D = series[1].pts.filter(Boolean), out = { areas: [], labels: [] };
    if (B.length < 2 || D.length < 2) return out;
    var alongX = !(X === "P" && Y === "T");                   // regions reach the top/bottom edge
    var liqHigh = alongX ? Y === "P" : true;                    // liquid toward high P
    function toEdge(pts, high) {
      var a = pts[0], b = pts[pts.length - 1];
      if (alongX) { var e = high ? yD[1] : yD[0]; return pts.concat([[b[0], e], [a[0], e]]); }
      var ex = high ? xD[1] : xD[0];
      return pts.concat([[ex, b[1]], [ex, a[1]]]);
    }
    out.areas.push({ color: "--liquid", opacity: 0.13, poly: toEdge(B, liqHigh) });
    out.areas.push({ color: "--vapor", opacity: 0.13, poly: toEdge(D, !liqHigh) });
    out.areas.push({ color: "--muted", opacity: 0.16, poly: lensArea(series) });
    if (!alongX) return out;
    var bx = B.map(function (p) { return p[0]; }), by = B.map(function (p) { return p[1]; }),
        dx = D.map(function (p) { return p[0]; }), dy = D.map(function (p) { return p[1]; });
    var best = null;
    for (var k = 1; k < 20; k++) {                              // where the two-phase band is widest
      var xm = xD[0] + (xD[1] - xD[0]) * k / 20, b = PEC.interp(bx, by, xm), d = PEC.interp(dx, dy, xm);
      if (b == null || d == null) continue;
      if (!best || Math.abs(d - b) > best.gap) best = { x: xm, b: b, d: d, gap: Math.abs(d - b) };
    }
    if (!best) return out;
    var lE = liqHigh ? yD[1] : yD[0], vE = liqHigh ? yD[0] : yD[1], span = Math.abs(yD[1] - yD[0]);
    if (Math.abs(lE - best.b) > 0.08 * span) out.labels.push({ x: best.x, y: (best.b + lE) / 2, text: "Liquid", minW: 300 });
    if (Math.abs(vE - best.d) > 0.08 * span) out.labels.push({ x: best.x, y: (best.d + vE) / 2, text: "Vapor", minW: 300 });
    if (best.gap > 0.07 * span) out.labels.push({ x: best.x, y: (best.b + best.d) / 2, text: "Liquid + vapor", minW: 300 });
    return out;
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
    if (last.res.hold === "field")
      use = [["position", "position (m)"], ["phase", "phase"], ["P", "P (bar)"], ["x1", "x1"], ["rho", "rho (kg/m3)"]];
    var lines = ["# " + nameA() + " (1) / " + nameB() + " (2); model " + last.req.model + "; T in K, P and f in bar",
      use.map(function (t) { return t[0]; }).join(",")];
    (c.x1 || c.position).forEach(function (_, i) { lines.push(use.map(function (t) { var v = c[t[0]][i]; return v == null ? "" : v; }).join(",")); });
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
