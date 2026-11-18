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
      var me = document.querySelector("script[src*='assets/app.js']"), v = me && (me.src.split("?v=")[1] || "");
      var worker = new Worker("assets/worker.js" + (v ? "?v=" + v : "")), pending = {}, seq = 0;
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
  var MODEL_SHORT = { PR: "Peng–Robinson", SRK: "SRK", vdW: "van der Waals", ideal: "Ideal", margules: "Margules", nrtl: "NRTL" };
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
    eppr78: "E-PPR78 group contribution" };

  /* the edition (set at build time in site/edition.json) limits what is offered */
  var ED = window.PEC_EDITION || {};
  function offered(kind, v) { return !ED[kind] || ED[kind].indexOf(v) >= 0; }
  var ACT = ["ideal", "margules", "nrtl"].filter(function (m) { return offered("models", m); });
  EOS = EOS.filter(function (m) { return offered("models", m); });
  /* gases and hydrocarbons go with the equations of state, polar liquids with the
     activity models; each family offers its own compounds */
  var SETS = ED.compound_sets || null, DEFPAIRS = ED.default_pairs || {};
  function familyOf(m) { return ["PR", "SRK", "vdW"].indexOf(m) >= 0 ? "eos" : "activity"; }
  function inSet(id, fam) { return !SETS || SETS[fam].indexOf(id) >= 0; }
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
    ["kij", "ppr78", "A12", "A21", "dg12", "dg21", "alpha"].forEach(function (k) { if (h.get(k) != null && isFinite(+h.get(k))) p[k] = +h.get(k); });
    p.__hasT = h.get("T") != null;
    return p;
  }
  function writeHash() {
    if (typeof mode !== "undefined" && mode !== "binary") return;     // the pure-fluid and mixture modes keep their own
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
    if (state.model === "PR" && pair.kij.source === "eppr78") return state.params.ppr78 ? ["ppr78"] : ["ppr78", "kij"];
    if (EOS.indexOf(state.model) >= 0) return ["kij"];
    if (state.model === "margules") return ["A12", "A21"];
    if (state.model === "nrtl") return ["dg12", "dg21", "alpha"];
    return [];
  }
  function defaults() {
    var n = pair.nrtl || { dg12: 0, dg21: 0, alpha: 0.3 };
    var A = pair.A;
    return { kij: state.model === "PR" ? pair.kij.value : 0, ppr78: pair.kij.source === "eppr78" ? 1 : 0,
      A12: A.A12, A21: A.A21, dg12: n.dg12, dg21: n.dg21, alpha: n.alpha };
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
  var PLURAL = { alkane: "Alkanes", alkene: "Alkenes", aromatic: "Aromatics", gas: "Gases", alcohol: "Alcohols", water: "Water", ketone: "Ketones",
    ester: "Esters", ether: "Ethers", chlorinated: "Chlorinated solvents" };
  /* options of a compound list, grouped by family; only the compounds of one model family */
  function compoundOptions(sel, fam, value) {
    var groups = {}, order = [];
    sel.innerHTML = "";
    catalog.compounds.forEach(function (c) {
      if (!inSet(c.id, fam)) return;
      if (!groups[c.family]) { groups[c.family] = []; order.push(c.family); }
      groups[c.family].push(c);
    });
    order.forEach(function (f) {
      var g = document.createElement("optgroup"); g.label = PLURAL[f] || f;
      groups[f].forEach(function (c) { var o = document.createElement("option"); o.value = c.id; o.textContent = c.name + "  (" + c.formula + ")"; g.appendChild(o); });
      sel.appendChild(g);
    });
    sel.value = value;
  }
  function renderCompounds() {
    compoundOptions($("comp-a"), familyOf(state.model), state.a);
    compoundOptions($("comp-b"), familyOf(state.model), state.b);
  }
  /* a pair that belongs to the model's family, keeping the current one if it does */
  function fitPair() {
    var fam = familyOf(state.model);
    if (inSet(state.a, fam) && inSet(state.b, fam)) return false;
    var d = DEFPAIRS[fam] || [state.a, state.b];
    state.a = d[0]; state.b = d[1];
    return true;
  }
  function fillCompounds() {
    renderCompounds();
    [["comp-a", "a"], ["comp-b", "b"]].forEach(function (pr) {
      var sel = $(pr[0]);
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
    function pick(m) {
      state.model = m;
      if (state.compare === m || (state.compare && familyOf(state.compare) !== familyOf(m))) state.compare = "";
      if (fitPair()) {                       // the other family's compounds, on a diagram that suits them
        renderCompounds();
        state.field = "none"; state.x = "xy"; state.view = "2d";
        if (familyOf(m) === "eos") { state.y = "P"; state.hold = "T"; state.T = 263; }
        else { state.y = "T"; state.hold = "P"; state.P = 1.01325; }
        loadPair({ __hasT: true }); return;
      }
      renderCompounds();
      state.params = defaults(); drawModels(); drawParams(); drawCompare(); schedule(0);
    }
    seg($("models-eos"), EOS.map(function (m) { return [m, MODEL_SHORT[m]]; }), state.model, pick);
    seg($("models-act"), ACT.map(function (m) { return [m, MODEL_SHORT[m]]; }), state.model, pick,
      function (m) { return !hasData(m); });
    [].slice.call($("models-act").children).forEach(function (b) {
      if (b.disabled) b.title = "The database has no " + b.textContent + " parameters for this pair";
    });
  }
  function drawCompare() {
    var sel = $("compare");
    sel.innerHTML = "<option value=''>Nothing</option>" + catalog.models.filter(function (m) {
      return m.id !== state.model && offered("models", m.id) && hasData(m.id) && familyOf(m.id) === familyOf(state.model); })
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
  /* NRTL needs parameters from the database for the pair */
  function hasData(m) {
    if (m === "nrtl") return !!pair.nrtl;
    return true;
  }
  function atRef() { return ""; }
  function paramSource() {
    if (state.model === "PR") return pair.kij.source;
    if (state.model === "SRK" || state.model === "vdW") return "none";
    if (state.model === "nrtl") return pair.nrtl ? pair.nrtl.source : "none";
    if (state.model === "ideal") return null;
    return pair.A.source;
  }
  function drawParams() {
    var host = $("params"), keys = activeParams();
    host.innerHTML = "";
    keys.forEach(function (k) {
      if (k === "ppr78") {                                  // a switch, not a slider
        var row0 = document.createElement("div"); row0.className = "param";
        row0.innerHTML = "<label class='check'><input type='checkbox' id='p-ppr78'" + (state.params.ppr78 ? " checked" : "") +
          "> k₁₂ from E-PPR78, varying with temperature</label>";
        host.appendChild(row0);
        row0.querySelector("input").addEventListener("change", function () {
          state.params.ppr78 = this.checked ? 1 : 0; drawParams(); schedule(0);
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
    if (state.model === "PR" && pair.kij.source === "eppr78")
      extra = state.params.ppr78 ? " k₁₂ = " + PEC.fmt(pair.kij.value, 4) + " at " + PEC.fmt(pair.T_ref, 1) + " K." : "";
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
    var p = params();
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
      if (!hasData(state.model)) state.model = familyOf(state.model) === "eos" ? "PR" : "ideal";   // no database values for this pair
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
  function params() {
    var p = {};
    activeParams().forEach(function (k) { if (k !== "ppr78") p[k] = state.params[k]; });
    return p;
  }
  function request() {
    var p = params();
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
  /* one order of each pair is precomputed; the other order is the same surface
     seen from the other side (x1 -> 1 - x1, components 1 and 2 swapped) */
  var SWAP = [["gamma1", "gamma2"], ["K1", "K2"], ["f1", "f2"], ["phiL1", "phiL2"], ["phiV1", "phiV2"]];
  function mirrorParams(p) {
    var q = {};
    Object.keys(p || {}).forEach(function (k) { q[k] = p[k]; });
    if ("dg12" in q) { q.dg12 = p.dg21; q.dg21 = p.dg12; }
    if ("A12" in q) { q.A12 = p.A21; q.A21 = p.A12; }
    return q;
  }
  function mirrorResult(res) {
    var out = {}; Object.keys(res).forEach(function (k) { out[k] = res[k]; });
    out.layers = res.layers.map(function (L) {
      var c = L.columns, m = {};
      Object.keys(c).forEach(function (k) { m[k] = c[k].slice().reverse(); });
      ["x1", "y1"].forEach(function (k) { if (m[k]) m[k] = m[k].map(function (v) { return v == null ? v : +(1 - v).toFixed(10); }); });
      SWAP.forEach(function (pr) { if (m[pr[0]] && m[pr[1]]) { var t = m[pr[0]]; m[pr[0]] = m[pr[1]]; m[pr[1]] = t; } });
      if (m.alpha) m.alpha = m.alpha.map(function (v) { return v == null || v === 0 ? v : 1 / v; });
      var az = L.azeotrope ? Object.assign({}, L.azeotrope, { x1: 1 - L.azeotrope.x1 }) : L.azeotrope;
      return { value: L.value, columns: m, azeotrope: az };
    });
    return out;
  }
  function precomputed(sreq) {
    return surfaceIndex.then(function (ix) {
      if (!ix) return null;
      var c = sreq.components, key = [c[0], c[1], sreq.model, sreq.hold].join("~"), rkey = [c[1], c[0], sreq.model, sreq.hold].join("~");
      var flip = ix.keys.indexOf(key) < 0;
      if (flip && ix.keys.indexOf(rkey) < 0) return null;
      return fetch("surfaces/" + encodeURIComponent(flip ? rkey : key) + ".json?v=" + ix.version).then(function (r) { return r.ok ? r.json() : null; })
        .then(function (f) {
          if (!f || f.version !== ix.version || !f.result || !f.result.ok) return null;
          var q = f.request, comps = flip ? q.components.slice().reverse() : q.components;
          var ok = q.model === sreq.model && q.hold === sreq.hold && q.npts === sreq.npts &&
            comps.join() === sreq.components.join() && sameNumbers(flip ? mirrorParams(q.params) : q.params, sreq.params, 1e-9) &&
            sameNumbers(q.values, sreq.values, 1e-6) && sreq.keep.every(function (k) { return q.keep.indexOf(k) >= 0; });
          if (!ok) return null;
          var res = flip ? mirrorResult(f.result) : f.result;
          res.source = "precomputed";
          return res;
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
  function ant(c) { return c.antoine ? [c.antoine.A, c.antoine.B, c.antoine.C].join(", ") : "none (equations of state only)"; }
  function groups(g) { return g ? Object.keys(g).map(function (k) { return g[k] + " " + k; }).join(", ") : "–"; }
  function drawSources() {
    var A = comp(state.a), B = comp(state.b), used = {};
    function src(k) { used[k] = true; return k; }
    var rows = [
      ["Formula", A.formula, B.formula, ""], ["CAS number", A.cas, B.cas, ""],
      ["Molar mass (g/mol)", PEC.fmt(A.M, 3), PEC.fmt(B.M, 3), src(A.sources.M)],
      ["Critical temperature (K)", PEC.fmt(A.Tc, 2), PEC.fmt(B.Tc, 2), src(A.sources.Tc)],
      ["Critical pressure (bar)", PEC.fmt(A.Pc_bar, 2), PEC.fmt(B.Pc_bar, 2), src(A.sources.Pc)],
      ["Acentric factor", PEC.fmt(A.omega, 3), PEC.fmt(B.omega, 3), src(A.sources.omega)],
      ["Antoine A, B, C", ant(A), ant(B), src("poling")],
      ["Antoine range (K)", A.antoine ? A.antoine.Tmin + " to " + A.antoine.Tmax : "–", B.antoine ? B.antoine.Tmin + " to " + B.antoine.Tmax : "–", "poling"],
      ["Normal boiling point (K)", A.Tb == null ? "–" : PEC.fmt(A.Tb, 2), B.Tb == null ? "–" : PEC.fmt(B.Tb, 2), src("derived")]];
    if (familyOf(state.model) === "eos" && (A.ppr78 || B.ppr78))
      rows.push(["E-PPR78 groups", groups(A.ppr78), groups(B.ppr78), src("eppr78")]);
    var h = "<table><thead><tr><th>Property</th><th>" + A.name + "</th><th>" + B.name + "</th><th>Source</th></tr></thead><tbody>" +
      rows.map(function (r) { return "<tr><td>" + r[0] + "</td><td>" + r[1] + "</td><td>" + r[2] + "</td><td>" + r[3] + "</td></tr>"; }).join("");
    var pr = [];
    if (familyOf(state.model) === "activity") {
      if (pair.nrtl) pr.push(["NRTL Δg₁₂, Δg₂₁ (J/mol), α", pair.nrtl.dg12 + ", " + pair.nrtl.dg21 + ", " + pair.nrtl.alpha, src(pair.nrtl.source)]);
      else pr.push(["NRTL parameters", "not in the database", ""]);
    } else if (pair.kij.source === "eppr78")
      pr.push(["Peng–Robinson k₁₂", PEC.fmt(pair.kij.value, 4) + " at " + PEC.fmt(pair.T_ref, 1) + " K (changes with temperature)", src("eppr78")]);
    else if (pair.kij.source !== "none") pr.push(["Peng–Robinson k₁₂", pair.kij.value, src(pair.kij.source)]);
    else pr.push(["Peng–Robinson k₁₂", "not in the database (0 used)", ""]);
    h += pr.map(function (r) { return "<tr><td>" + r[0] + "</td><td colspan='2'>" + r[1] + "</td><td>" + r[2] + "</td></tr>"; }).join("") + "</tbody></table>";
    $("source-table").innerHTML = h;
    used.libraries = true;
    $("source-notes").innerHTML = "<p class='cite'>Antoine equation: log10(P/Pa) = A − B/(T/K + C).</p>" +
      Object.keys(used).filter(function (k) { return catalog.sources[k]; })
        .map(function (k) { return "<p class='cite'><b>" + k + "</b>: " + catalog.sources[k] + "</p>"; }).join("");
  }

  /* ------------------------------------------------------------ three or four components */
  /* A flash of a multicomponent feed at T and P (web_api.mixture). For three
     components the page also asks for the bubble temperature over a grid of
     liquid compositions at the same P, drawn as a triangle diagram. */
  var mix = { model: null, comps: [], z: [], examples: [], map: null, mapKey: null, seq: 0 };
  function mixModels() { return EOS.concat(ACT.filter(function (m) { return m !== "margules"; })); }
  function bindMixture() {
    engine.call({ op: "mixture_examples" }).then(function (r) {
      mix.examples = (r.examples || []).filter(function (ex) {
        return offered("models", ex.model) && ex.components.every(function (c) { return comp(c) && inSet(c, familyOf(ex.model)); });
      });
      seg($("mix-examples"), mix.examples.map(function (ex, i) { return [String(i), ex.label]; }), null,
        function (i) { loadMixExample(mix.examples[+i]); });
      var start = mix.examples.filter(function (ex) { return ex.components.length === 3 && ex.model === "nrtl"; })[0] || mix.examples[0];
      if (start) loadMixExample(start, true);
    });
    $("mix-add").addEventListener("click", function () {
      readMix();
      if (mix.comps.length >= 4) return;
      var fam = familyOf(mix.model), next = catalog.compounds.filter(function (c) { return inSet(c.id, fam) && mix.comps.indexOf(c.id) < 0; })[0];
      if (!next) return;
      mix.comps.push(next.id); mix.z.push(0.2); drawMix();
    });
    $("mix-run").addEventListener("click", function () { readMix(); runMixture(); });
    ["mix-T", "mix-P"].forEach(function (id) { $(id).addEventListener("keydown", function (e) { if (e.key === "Enter") { readMix(); runMixture(); } }); });
  }
  function loadMixExample(ex, quiet) {
    mix.model = ex.model; mix.comps = ex.components.slice(); mix.z = ex.z.slice();
    $("mix-T").value = ex.T; $("mix-P").value = ex.P;
    seg($("mix-examples"), mix.examples.map(function (e, i) { return [String(i), e.label]; }), String(mix.examples.indexOf(ex)),
      function (i) { loadMixExample(mix.examples[+i]); });
    drawMix();
    runMixture(quiet);
  }
  function readMix() {
    var rows = [].slice.call($("mix-rows").querySelectorAll(".mix-row"));
    mix.comps = rows.map(function (r) { return r.querySelector("select").value; });
    mix.z = rows.map(function (r) { return +r.querySelector("input").value; });
  }
  function drawMix() {
    seg($("mix-model"), mixModels().map(function (m) { return [m, MODEL_SHORT[m]]; }), mix.model, function (m) {
      readMix();
      var fam = familyOf(m);
      mix.model = m;
      if (!mix.comps.every(function (c) { return inSet(c, fam); })) {    // the other family: start from its example
        var ex = mix.examples.filter(function (e) { return familyOf(e.model) === fam && e.components.length === mix.comps.length; })[0] ||
          mix.examples.filter(function (e) { return familyOf(e.model) === fam; })[0];
        if (ex) { loadMixExample(ex); mix.model = m; drawMix(); return; }
      }
      drawMix();
    });
    var host = $("mix-rows"), fam = familyOf(mix.model);
    host.innerHTML = "";
    mix.comps.forEach(function (c, i) {
      var row = document.createElement("div"); row.className = "mix-row";
      row.innerHTML = "<label><span class='sr-only'>Component " + (i + 1) + "</span><select aria-label='Component " + (i + 1) + "'></select></label>" +
        "<label class='mix-z'>Amount<input type='number' step='any' min='0' value='" + mix.z[i] + "'></label>" +
        (mix.comps.length > 3 ? "<button type='button' class='iconbtn' title='Remove' aria-label='Remove this component'>×</button>" : "<span></span>");
      host.appendChild(row);
      compoundOptions(row.querySelector("select"), fam, c);
      var rm = row.querySelector("button");
      if (rm) rm.addEventListener("click", function () { readMix(); mix.comps.splice(i, 1); mix.z.splice(i, 1); drawMix(); });
    });
    $("mix-add").hidden = mix.comps.length >= 4;
  }
  function runMixture(quiet) {
    var box = $("mix-result"), T = +$("mix-T").value, P = +$("mix-P").value;
    var uniq = mix.comps.filter(function (c, i) { return mix.comps.indexOf(c) === i; });
    if (uniq.length !== mix.comps.length) { box.innerHTML = "<p class='notice'>Choose different compounds.</p>"; return; }
    if (!(T > 0) || !(P > 0) || !mix.z.every(function (v) { return v > 0; })) {
      box.innerHTML = "<p class='notice'>Give T > 0, P > 0 and a positive amount of every component.</p>"; return;
    }
    var key = JSON.stringify([mix.comps, mix.model, familyOf(mix.model) === "eos" ? "T" + T : "P" + P]), wantMap = mix.comps.length === 3 && key !== mix.mapKey;
    var id = ++mix.seq;
    box.innerHTML = "<p class='source'>Calculating" + (wantMap ? " the flash and the triangle diagram" : "") + "…</p>";
    engine.call({ op: "mixture", request: { components: mix.comps, model: mix.model, T: T, P: P, z: mix.z, map: wantMap, divisions: 10 } }).then(function (r) {
      if (id !== mix.seq) return;
      if (!r.ok) { box.innerHTML = "<p class='notice'>" + (r.error || "This mixture could not be calculated.") + "</p>"; $("mix-figure").hidden = true; return; }
      if (r.map) { mix.map = r.map; mix.mapKey = key; }
      if (mix.comps.length !== 3) { mix.map = null; mix.mapKey = null; }
      showMixture(r, T, P);
    });
  }
  function showMixture(r, T, P) {
    var names = r.names.map(function (c) { return comp(c).name; });
    var s = function (v) { return v == null ? "–" : PEC.sig(v, 4); }, f = function (v, d) { return v == null ? "–" : PEC.fmt(v, d); };
    var b = r.bounds || {}, two = r.phase === "two-phase";
    var h = "<p class='check-verdict ph-" + r.phase + "'>" + phaseWord(r.phase) + "</p>";
    h += two ? "<p>" + f(r.beta * 100, 1) + "% of the feed is vapor at " + f(T, 2) + " K and " + s(P) + " bar.</p>"
             : "<p>No two-phase equilibrium at this state: the mixture is all " + r.phase + ".</p>";
    var cols = ["Component", "Feed z", "Liquid x", "Vapor y"];
    if (two) cols.push("K = y/x");
    if (r.gamma) cols.push("γ");
    if (r.phiL) cols.push("φ liquid", "φ vapor");
    h += "<div class='table-wrap'><table class='mix-table'><thead><tr>" + cols.map(function (c) { return "<th>" + c + "</th>"; }).join("") + "</tr></thead><tbody>";
    names.forEach(function (n, i) {
      var bar = function (v, cls) { return v == null ? "–" : "<span class='mbar " + cls + "' style='width:" + Math.round(60 * v) + "px'></span>" + PEC.fmt(v, 4); };
      h += "<tr><td>" + n + "</td><td>" + bar(r.z[i], "mz") + "</td><td>" + bar(r.x[i], "mx") + "</td><td>" + bar(r.y[i], "my") + "</td>" +
        (two ? "<td>" + s(r.K[i]) + "</td>" : "") + (r.gamma ? "<td>" + s(r.gamma[i]) + "</td>" : "") +
        (r.phiL ? "<td>" + s(r.phiL[i]) + "</td><td>" + s(r.phiV[i]) + "</td>" : "") + "</tr>";
    });
    h += "</tbody></table></div>";
    h += "<p><b>Where equilibrium occurs for this feed:</b> at " + f(T, 2) + " K, between the dew pressure " + s(b.P_dew) + " bar and the bubble pressure " + s(b.P_bubble) +
      " bar; at " + s(P) + " bar, between the bubble temperature " + f(b.T_bubble, 2) + " K and the dew temperature " + f(b.T_dew, 2) + " K." +
      (b.P_bubble == null || b.T_bubble == null || b.P_dew == null || b.T_dew == null ? " (–: not found, usually because the state is near or above the mixture critical point.)" : "") + "</p>";
    if (r.kij) h += "<p class='source'>k_ij at " + f(T, 2) + " K: " + Object.keys(r.kij).map(function (k) {
      var pr = k.split("|"); return comp(pr[0]).name + "–" + comp(pr[1]).name + " " + PEC.fmt(r.kij[k], 4); }).join(", ") +
      (/PPR78/.test(r.input_text) ? " (E-PPR78)" : "") + "</p>";
    h += "<p class='source'>" + r.route + " · " + r.ms + " ms</p>";
    var esc = function (t) { return t.replace(/&/g, "&amp;").replace(/</g, "&lt;"); };
    h += "<details><summary>Input file and checker report</summary><pre>" + esc(r.input_text) + "</pre><pre>" + esc(r.report) + "</pre></details>";
    $("mix-result").innerHTML = h;
    if (mix.map && r.names.length === 3) drawTriangle(r, names); else $("mix-figure").hidden = true;
  }

  /* triangle diagram: composition of a ternary, the bubble point of each liquid
     shaded from light (low) to dark (high) (bubble temperature at P for an
     activity model, bubble pressure at T, on a log scale, for an equation of
     state), the vapor in equilibrium with each liquid, and the feed's tie line */
  function drawTriangle(r, names) {
    var fig = $("mix-figure"), host = $("mix-chart"), map = mix.map, kind = map.kind || "T";
    var pts = map.points.map(function (p) { return { x: p.x, y: p.y, T: kind === "T" ? p.T : p.P }; });
    fig.hidden = false;
    var n = map.n || Math.round((Math.sqrt(8 * pts.length + 1) - 3) / 2);
    var W = 560, H = 500, top = [W / 2, 34], left = [48, 448], right = [W - 48, 448];
    function xy(c) { return [c[0] * top[0] + c[1] * left[0] + c[2] * right[0], c[0] * top[1] + c[1] * left[1] + c[2] * right[1]]; }
    var at = {};
    pts.forEach(function (p) { at[Math.round(p.x[0] * n) + "," + Math.round(p.x[1] * n)] = p; });
    var Ts = pts.map(function (p) { return p.T; }).filter(function (v) { return v != null; });
    var lo = Math.min.apply(null, Ts), hi = Math.max.apply(null, Ts);
    var dark = (function () {
      var m = getComputedStyle(document.body).backgroundColor.match(/\d+/g);
      return m ? (0.299 * m[0] + 0.587 * m[1] + 0.114 * m[2]) < 128 : false;
    })();
    var ramp = dark ? [[44, 36, 60], [196, 170, 238]] : [[240, 235, 248], [92, 60, 146]];   // one hue, light to dark (dark mode: dim to bright)
    function col(T) {
      var u = hi > lo ? (kind === "T" ? (T - lo) / (hi - lo) : Math.log(T / lo) / Math.log(hi / lo)) : 0.5;
      return "rgb(" + ramp[0].map(function (c, k) { return Math.round(c + (ramp[1][k] - c) * u); }).join(",") + ")";
    }
    var out = [];
    function poly(ps, fill, title) {
      out.push("<polygon points='" + ps.map(function (q) { return xy(q.x).map(function (v) { return v.toFixed(1); }).join(","); }).join(" ") +
        "' fill='" + fill + "' stroke='" + fill + "' stroke-width='0.6'><title>" + title + "</title></polygon>");
    }
    function tip(p) {
      return names.map(function (nm, k) { return nm + " " + PEC.fmt(p.x[k], 2); }).join(", ") +
        (kind === "T" ? "\nbubble T = " + PEC.fmt(p.T, 1) + " K" : "\nbubble P = " + PEC.sig(p.T, 4) + " bar") +
        (p.y ? "\nvapor: " + p.y.map(function (v) { return PEC.fmt(v, 3); }).join(", ") : "");
    }
    for (var i = 0; i < n; i++) for (var j = 0; j < n - i; j++) {
      var a = at[i + "," + j], b2 = at[(i + 1) + "," + j], c = at[i + "," + (j + 1)], d = at[(i + 1) + "," + (j + 1)];
      if (a && b2 && c && a.T != null && b2.T != null && c.T != null) poly([a, b2, c], col((a.T + b2.T + c.T) / 3), tip(a));
      if (j < n - i - 1 && b2 && c && d && b2.T != null && c.T != null && d.T != null) poly([b2, c, d], col((b2.T + c.T + d.T) / 3), tip(d));
    }
    /* gridlines every 0.2 */
    for (var g = 1; g < 5; g++) {
      var u = g / 5;
      [[[u, 1 - u, 0], [u, 0, 1 - u]], [[1 - u, u, 0], [0, u, 1 - u]], [[1 - u, 0, u], [0, 1 - u, u]]].forEach(function (seg2) {
        var p1 = xy(seg2[0]), p2 = xy(seg2[1]);
        out.push("<line x1='" + p1[0] + "' y1='" + p1[1] + "' x2='" + p2[0] + "' y2='" + p2[1] + "' class='tri-grid'/>");
      });
    }
    /* liquid -> vapor at every other grid point */
    pts.forEach(function (p) {
      var ii = Math.round(p.x[0] * n), jj = Math.round(p.x[1] * n);
      if (!p.y || ii % 2 || jj % 2) return;
      var p1 = xy(p.x), p2 = xy(p.y);
      out.push("<line x1='" + p1[0] + "' y1='" + p1[1] + "' x2='" + p2[0] + "' y2='" + p2[1] + "' class='tri-vec'/>" +
        "<circle cx='" + p2[0] + "' cy='" + p2[1] + "' r='2' class='tri-vec-end'/>");
    });
    out.push("<polygon points='" + [top, left, right].map(function (q) { return q.join(","); }).join(" ") + "' class='tri-frame'/>");
    /* the feed and its tie line */
    if (r.phase === "two-phase") {
      var px = xy(r.x), py = xy(r.y);
      out.push("<line x1='" + px[0] + "' y1='" + px[1] + "' x2='" + py[0] + "' y2='" + py[1] + "' class='tri-tie'/>");
      out.push("<circle cx='" + px[0] + "' cy='" + px[1] + "' r='6' class='tri-pt' fill='var(--liquid)'><title>liquid x</title></circle>");
      out.push("<circle cx='" + py[0] + "' cy='" + py[1] + "' r='6' class='tri-pt' fill='var(--vapor)'><title>vapor y</title></circle>");
    }
    var pz = xy(r.z);
    out.push("<rect x='" + (pz[0] - 6) + "' y='" + (pz[1] - 6) + "' width='12' height='12' class='tri-pt' style='fill:var(--ink)'><title>feed z</title></rect>");
    out.push("<text x='" + top[0] + "' y='" + (top[1] - 12) + "' text-anchor='middle' class='tri-label'>" + names[0] + "</text>");
    out.push("<text x='" + left[0] + "' y='" + (left[1] + 26) + "' text-anchor='start' class='tri-label'>" + names[1] + "</text>");
    out.push("<text x='" + right[0] + "' y='" + (right[1] + 26) + "' text-anchor='end' class='tri-label'>" + names[2] + "</text>");
    var cond = kind === "T" ? PEC.sig(map.P, 4) + " bar" : PEC.fmt(map.T, 1) + " K";
    host.innerHTML = "<svg viewBox='0 0 " + W + " " + H + "' role='img' aria-label='Triangle diagram of " + names.join(", ") +
      " at " + cond + "'>" + out.join("") + "</svg>";
    var grad = "linear-gradient(90deg, " + col(lo) + ", " + col(hi) + ")";
    $("mix-legend").innerHTML =
      "<li><span class='swatch-sq' style='background:var(--ink)'></span>Feed z</li>" +
      (r.phase === "two-phase" ? "<li><span class='dot' style='background:var(--liquid)'></span>Liquid x</li><li><span class='dot' style='background:var(--vapor)'></span>Vapor y</li>" : "") +
      "<li><span class='tri-vec-key'></span>Liquid to the vapor it boils into</li>" +
      (kind === "T" ? "<li>Bubble temperature at " + cond + ": " + PEC.fmt(lo, 1) + " K <span class='tri-ramp' style='background:" + grad + "'></span> " + PEC.fmt(hi, 1) + " K</li>"
                    : "<li>Bubble pressure at " + cond + ": " + PEC.sig(lo, 3) + " bar <span class='tri-ramp' style='background:" + grad + "'></span> " + PEC.sig(hi, 3) + " bar (log scale)</li>") +
      (Ts.length < pts.length ? "<li>Blank: no bubble point found (near or above the mixture critical point)</li>" : "");
  }

  /* ------------------------------------------------------------ pure fluid */
  /* One compound with a cubic equation of state (web_api.pure): vapor
     pressure to the critical point, isotherms, coexisting densities, the
     enthalpy of vaporization, and, with a field, the column, the band of
     reference pressures that holds both phases, and the interface. */
  var pure = { compound: "water", model: "PR", T: 350, P: 1.01325, field: "none", g: 9.81, H: 100, rpm: 20000, ri: 0.01, ro: 0.1,
    plot: "PT", res: null, seq: 0, all: [] };
  var PURE_PLOTS = [["PT", "Vapor pressure (P–T)"], ["PV", "Isotherms (P–V)"], ["Trho", "Densities (T–ρ)"], ["dH", "Enthalpy of vaporization"]];
  var PURE_FIELD_PLOTS = [["col", "Pressure along the field"], ["colrho", "Density along the field"], ["band", "Phase boundary shift (P–T)"],
    ["iface", "Interface against field strength"]];
  function bindPure(all) {
    pure.all = all;
    var sel = $("pure-comp"), groups = {}, order = [];
    all.forEach(function (c) { if (!groups[c.family]) { groups[c.family] = []; order.push(c.family); } groups[c.family].push(c); });
    order.forEach(function (f) {
      var g = document.createElement("optgroup"); g.label = PLURAL[f] || f;
      groups[f].forEach(function (c) { var o = document.createElement("option"); o.value = c.id; o.textContent = c.name + "  (" + c.formula + ")"; g.appendChild(o); });
      sel.appendChild(g);
    });
    sel.value = pure.compound;
    sel.addEventListener("change", function () { pure.compound = sel.value; pureDefaults(); runPure(); });
    $("pure-run").addEventListener("click", function () { readPure(); runPure(); });
    ["pure-T", "pure-P"].forEach(function (id) { $(id).addEventListener("keydown", function (e) { if (e.key === "Enter") { readPure(); runPure(); } }); });
    $("pure-download").addEventListener("click", pureCsv);
    drawPureControls();
    pureDefaults();
  }
  function pureDefaults() {
    var c = pure.all.filter(function (x) { return x.id === pure.compound; })[0];
    if (!c) return;
    pure.T = +(c.Tb ? Math.min(c.Tb + 20, 0.9 * c.Tc) : 0.8 * c.Tc).toFixed(1);
    pure.P = 1.01325;
    $("pure-T").value = pure.T; $("pure-P").value = pure.P;
  }
  function readPure() {
    pure.T = +$("pure-T").value; pure.P = +$("pure-P").value;
    var q = function (id) { var e = $(id); return e ? +e.value : null; };
    if (pure.field === "gravity") { pure.g = q("pf-g"); pure.H = q("pf-H"); }
    if (pure.field === "centrifugal") { pure.rpm = q("pf-rpm"); pure.ri = q("pf-ri"); pure.ro = q("pf-ro"); }
  }
  function drawPureControls() {
    seg($("pure-model"), EOS.map(function (m) { return [m, MODEL_SHORT[m]]; }), pure.model, function (m) { readPure(); pure.model = m; drawPureControls(); runPure(); });
    seg($("pure-field"), [["none", "None"], ["gravity", "Gravity"], ["centrifugal", "Centrifuge"]], pure.field, function (f) {
      readPure(); pure.field = f;
      if (f === "none" && PURE_FIELD_PLOTS.some(function (p) { return p[0] === pure.plot; })) pure.plot = "PT";
      if (f !== "none" && !PURE_FIELD_PLOTS.some(function (p) { return p[0] === pure.plot; })) pure.plot = "band";
      drawPureControls(); runPure();
    });
    var fp = $("pure-field-params");
    var inp = function (id, label, v) { return "<label>" + label + "<input type='number' step='any' id='" + id + "' value='" + v + "'></label>"; };
    fp.innerHTML = pure.field === "gravity" ? inp("pf-g", "g (m/s²)", pure.g) + inp("pf-H", "Column height (m)", pure.H)
      : pure.field === "centrifugal" ? inp("pf-rpm", "Speed (rpm)", pure.rpm) + inp("pf-ri", "Inner radius (m)", pure.ri) + inp("pf-ro", "Outer radius (m)", pure.ro) : "";
    $("pure-P-label").textContent = pure.field === "gravity" ? "Pressure at the bottom (bar)" : pure.field === "centrifugal" ? "Pressure at the inner radius (bar)" : "Pressure (bar)";
    var plots = PURE_PLOTS.concat(pure.field === "none" ? [] : PURE_FIELD_PLOTS);
    seg($("pure-plots"), plots, pure.plot, function (k) { pure.plot = k; drawPureControls(); if (pure.res) drawPure(); });
  }
  function runPure() {
    readPure();
    if (!(pure.T > 0) || !(pure.P > 0)) { $("pure-status").textContent = "Give a temperature and a pressure above zero."; return; }
    var req = { compound: pure.compound, model: pure.model, T: pure.T, P: pure.P };
    if (pure.field === "gravity") req.field = { kind: "gravity", g: pure.g, height: pure.H };
    if (pure.field === "centrifugal") req.field = { kind: "centrifugal", rpm: pure.rpm, r_in: pure.ri, r_out: pure.ro };
    var id = ++pure.seq;
    $("pure-status").textContent = "Calculating…";
    engine.call({ op: "pure", request: req }).then(function (r) {
      if (id !== pure.seq) return;
      if (!r.ok) { $("pure-status").textContent = r.error || "This could not be calculated."; return; }
      pure.res = r;
      $("pure-status").textContent = MODEL_SHORT[r.model] + " · " + r.ms + " ms";
      drawPure();
    });
  }
  function zip2(a, b, fa, fb) { return a.map(function (v, i) { return v == null || b[i] == null ? null : [fa ? fa(v) : v, fb ? fb(b[i]) : b[i]]; }); }
  function drawPure() {
    var r = pure.res, c = r.curve, st = r.state, host = $("pure-chart"), cfg, leg = [], cap = "", nm = r.name;
    var f = function (v, d) { return v == null ? "–" : PEC.fmt(v, d); }, sg = function (v) { return v == null ? "–" : PEC.sig(v, 4); };
    var crit = { x: r.Tc, y: r.Pc, color: "--ink", text: "critical point", dx: -8, anchor: "end", dy: -8 };
    if (pure.plot === "PT") {
      var ser = [{ name: "Vapor pressure, " + MODEL_SHORT[r.model], color: "--liquid", pts: zip2(c.T, c.P) }];
      leg.push({ name: "Vapor pressure, " + MODEL_SHORT[r.model], color: "--liquid" });
      if (c.P_antoine.some(function (v) { return v != null; })) {
        ser.push({ name: "Antoine equation (measured data)", color: "--third", dash: "5 4", pts: zip2(c.T, c.P_antoine) });
        leg.push({ name: "Antoine equation (measured data, Poling et al.)", color: "--third" });
      }
      leg.push({ name: "Critical point and your state", color: "--ink", dot: true });
      cfg = { series: ser, points: [crit, { x: st.T, y: st.P, color: "--vapor", text: "T, P", dx: 8 }], xTitle: "Temperature (K)", yTitle: "Pressure (bar)",
        regions: [{ x: c.T[Math.floor(c.T.length * 0.35)], y: r.Pc * 0.75, text: "liquid", minW: 300 }, { x: c.T[Math.floor(c.T.length * 0.75)], y: r.Pc * 0.12, text: "vapor", minW: 300 }],
        fx: function (v) { return "T = " + PEC.fmt(v, 2) + " K"; }, fy: function (v) { return "P = " + PEC.sig(v, 4) + " bar"; } };
      cap = "<b>" + nm + ": vapor pressure</b> from equal fugacities of liquid and vapor, up to the critical point (" + f(r.Tc, 2) + " K, " + sg(r.Pc) + " bar). Above the curve the stable phase is liquid, below it vapor.";
    } else if (pure.plot === "PV") {
      var series = [];
      r.isotherms.forEach(function (it, k) {
        var col = it.Tr < 1 ? "--liquid" : it.Tr === 1 ? "--ink" : "--vapor";
        series.push({ name: "T/Tc = " + it.Tr, color: col, width: it.Tr === 1 ? 2 : 1.5, opacity: it.Tr < 1 ? 0.45 + 0.1 * k : 1, pts: zip2(it.V, it.P) });
        if (it.Psat != null) series.push({ name: "equal-area line, T/Tc = " + it.Tr, color: "--muted", dash: "4 3", width: 1.5, pts: [[it.vL, it.Psat], [it.vV, it.Psat]] });
      });
      series.push({ name: "Two-phase dome", color: "--third", width: 2.5, pts: zip2(r.dome.vL, r.dome.P).concat(zip2(r.dome.vV, r.dome.P).reverse()) });
      leg = [{ name: "Isotherms at T/Tc = 0.80, 0.85, 0.90, 0.95", color: "--liquid" }, { name: "Critical isotherm", color: "--ink" }, { name: "T/Tc = 1.1", color: "--vapor" },
        { name: "Equal-area (Maxwell) line", color: "--muted" }, { name: "Two-phase dome", color: "--third" }];
      cfg = { series: series, xDomain: [0.3, 6], yDomain: [-0.6, 1.8], xTitle: "Molar volume / critical volume", yTitle: "Pressure / critical pressure", xDec: 1, yDec: 1,
        hlines: [{ y: 0 }], fx: function (v) { return "V/Vc = " + PEC.fmt(v, 3); }, fy: function (v) { return "P/Pc = " + PEC.fmt(v, 3); } };
      var a9 = r.isotherms.filter(function (it) { return it.Tr === 0.9; })[0];
      cap = "<b>" + nm + ": isotherms of " + MODEL_SHORT[r.model] + "</b> in reduced units. Below Tc each isotherm has a loop; the stable states follow the dashed line at the vapor pressure, placed so the areas above and below it are equal (at T/Tc = 0.9 they differ by " +
        (a9 && a9.area != null ? PEC.sig(Math.abs(a9.area) * 100, 2) + "% of P_sat(v_V − v_L)" : "–") + "). Parts of the loop can fall below zero pressure: a liquid under tension.";
    } else if (pure.plot === "Trho") {
      cfg = { series: [{ name: "Saturated liquid", color: "--liquid", pts: zip2(c.rhoL, c.T) }, { name: "Saturated vapor", color: "--vapor", pts: zip2(c.rhoV, c.T) }],
        points: [{ x: r.rhoc, y: r.Tc, color: "--ink", text: "critical point", dx: 8, dy: -6 }], xTitle: "Density (kg/m³)", yTitle: "Temperature (K)",
        fx: function (v) { return "ρ = " + PEC.sig(v, 4) + " kg/m³"; }, fy: function (v) { return "T = " + PEC.fmt(v, 2) + " K"; } };
      leg = [{ name: "Saturated liquid", color: "--liquid" }, { name: "Saturated vapor", color: "--vapor" }, { name: "Critical point", color: "--ink", dot: true }];
      cap = "<b>" + nm + ": densities of the coexisting liquid and vapor.</b> They approach each other and meet at the critical point, where the two phases become one.";
    } else if (pure.plot === "dH") {
      cfg = { series: [{ name: "Enthalpy of vaporization", color: "--liquid", pts: zip2(c.T, c.dHvap) }], xTitle: "Temperature (K)", yTitle: "ΔH vap (kJ/mol)",
        fx: function (v) { return "T = " + PEC.fmt(v, 2) + " K"; }, fy: function (v) { return "ΔH = " + PEC.sig(v, 4) + " kJ/mol"; } };
      leg = [{ name: "Clapeyron: T (v_V − v_L) dP_sat/dT", color: "--liquid" }];
      cap = "<b>" + nm + ": enthalpy of vaporization</b> from the Clapeyron equation and the calculated vapor-pressure curve. It falls to zero at the critical point.";
    } else {
      var F = r.field, cl = F.column, posT = F.kind === "gravity" ? "Height above the bottom (m)" : "Radius (m)";
      if (pure.plot === "col" || pure.plot === "colrho") {
        var key = pure.plot === "col" ? "P" : "rho", liq = [], vap = [];
        cl.position.forEach(function (x, i) {
          var p = [x, cl[key][i]], isL = cl.phase[i] === "liquid";
          (isL ? liq : vap).push(p); (isL ? vap : liq).push(null);
        });
        cfg = { series: [{ name: "Liquid", color: "--liquid", pts: liq }, { name: "Vapor or gas", color: "--vapor", pts: vap }],
          vlines: cl.interface != null ? [{ x: cl.interface, text: "interface" }] : [], hlines: key === "P" && cl.Psat != null ? [{ y: cl.Psat, text: "P_sat at " + f(st.T, 1) + " K" }] : [],
          xTitle: posT, yTitle: key === "P" ? "Pressure (bar)" : "Density (kg/m³)", yOpts: { pad: 0.08 },
          fx: function (v) { return "position = " + PEC.sig(v, 4) + " m"; }, fy: function (v) { return key === "P" ? "P = " + PEC.sig(v, 5) + " bar" : "ρ = " + PEC.sig(v, 4) + " kg/m³"; } };
        leg = [{ name: "Liquid", color: "--liquid" }, { name: "Vapor (or one fluid above Tc)", color: "--vapor" }];
        cap = "<b>" + nm + " at " + f(st.T, 1) + " K in " + (F.kind === "gravity" ? "a " + pure.H + " m column, g = " + pure.g + " m/s²" : "a centrifuge at " + pure.rpm + " rpm") +
          ".</b> The fugacity changes along the field as ln f = ln f_ref − M ΔΦ/RT; where it passes f_sat the fluid changes phase. " +
          (cl.interface != null ? "The interface is at " + PEC.sig(cl.interface, 4) + " m." : (cl.notes[0] || ""));
      } else if (pure.plot === "band") {
        var full = F.bands[F.bands.length - 1], series2 = [], shift = function (b, key) {
          return b.T.map(function (t, i) { return b[key][i] == null || b.Psat[i] == null ? null : [t, b[key][i] - b.Psat[i]]; }); };
        F.bands.forEach(function (b) {
          var lab = F.kind === "gravity" ? "g × " + b.fraction : "rpm × " + PEC.fmt(Math.sqrt(b.fraction), 2);
          var w = b.fraction === 1 ? 2 : 1.2, dash = b.fraction === 1 ? "none" : "4 3";
          series2.push({ name: "first liquid, " + lab, color: "--liquid", width: w, dash: dash, pts: shift(b, "P_first") });
          series2.push({ name: "last vapor, " + lab, color: "--vapor", width: w, dash: dash, pts: shift(b, "P_last") });
        });
        var poly = shift(full, "P_first").concat(shift(full, "P_last").reverse()).filter(Boolean);
        var psT = st.Psat != null ? st.P - st.Psat : null;
        cfg = { series: series2, areas: [{ color: "--third", poly: poly, opacity: 0.18 }], hlines: [{ y: 0, text: "no field: P = P_sat(T)" }],
          xTitle: "Temperature (K)", yTitle: "P at the " + (F.kind === "gravity" ? "bottom" : "inner radius") + " − P_sat (bar)", yOpts: { pad: 0.08 },
          points: psT != null ? [{ x: st.T, y: psT, color: "--ink", text: "your state", dx: 8 }] : [],
          fx: function (v) { return "T = " + PEC.fmt(v, 2) + " K"; }, fy: function (v) { return "P − P_sat = " + PEC.sig(v, 4) + " bar"; } };
        leg = [{ name: "First liquid appears", color: "--liquid" }, { name: "Last vapor disappears", color: "--vapor" },
          { name: "Both phases in the " + (F.kind === "gravity" ? "column" : "rotor") + " (full strength)", color: "--third", area: true }];
        cap = "<b>How the field shifts the phase boundary of " + lower(nm) + ".</b> Without a field liquid and vapor coexist only on the line P = P_sat(T) (zero here). In the field the pressure differs along the " +
          (F.kind === "gravity" ? "column" : "rotor") + ", so, read at the reference point, coexistence spreads into a band: from the pressure at which the first liquid forms (where the potential is lowest) to the one at which the last vapor disappears. The plot shows that band as a shift from P_sat. Solid: the chosen field; dashed: a quarter and half of its strength" +
          (F.kind === "gravity" ? " (in g)." : " (in ω²).");
      } else {
        var ic = F.interface_curve;
        if (!ic) { host.innerHTML = "<p class='notice'>Above the critical temperature there is no interface: the fluid is one phase whose density changes along the field.</p>"; $("pure-legend").innerHTML = ""; $("pure-caption").innerHTML = ""; pureFacts(); return; }
        cfg = { series: [{ name: "Interface position", color: "--third", pts: zip2(ic.strength, ic.position) }],
          xTitle: F.kind === "gravity" ? "g (m/s²)" : "Speed (rpm)", yTitle: F.kind === "gravity" ? "Interface height (m)" : "Interface radius (m)",
          fx: function (v) { return (F.kind === "gravity" ? "g = " : "speed = ") + PEC.sig(v, 4); }, fy: function (v) { return "position = " + PEC.sig(v, 4) + " m"; } };
        leg = [{ name: "Interface position", color: "--third" }];
        cap = "<b>Where the liquid–vapor interface sits</b> at " + f(st.T, 1) + " K and " + sg(st.P) + " bar at the reference point, as the field gets stronger. Gaps: the interface lies outside the " + (F.kind === "gravity" ? "column" : "rotor") + " there.";
      }
    }
    cfg.aria = $("pure-plots").querySelector("[aria-pressed=true]") ? $("pure-plots").querySelector("[aria-pressed=true]").textContent + " for " + nm : nm;
    PEC.legend("pure-legend", leg);
    PEC.drawChart(host, cfg);
    $("pure-caption").innerHTML = cap;
    pureFacts(); pureTable(cfg);
  }
  function pureFacts() {
    var r = pure.res, st = r.state, f = function (v, d) { return v == null ? "–" : PEC.fmt(v, d); }, sg = function (v) { return v == null ? "–" : PEC.sig(v, 4); };
    var items = [["Critical point", f(r.Tc, 2) + " K, " + sg(r.Pc) + " bar, " + sg(r.rhoc) + " kg/m³"],
      ["Normal boiling point", f(r.Tb_eos, 2) + " K (" + MODEL_SHORT[r.model] + ")" + (r.Tb_data ? ", " + f(r.Tb_data, 2) + " K from the Antoine data" : "")],
      ["State at " + f(st.T, 2) + " K, " + sg(st.P) + " bar", st.phase + ", ρ = " + sg(st.rho) + " kg/m³, Z = " + PEC.fmt(st.Z, 4)]];
    if (st.Psat != null) items.push(["Vapor pressure at " + f(st.T, 2) + " K", sg(st.Psat) + " bar; ρ_L = " + sg(st.rhoL_sat) + ", ρ_V = " + sg(st.rhoV_sat) + " kg/m³"]);
    if (st.Tsat != null) items.push(["Boiling temperature at " + sg(st.P) + " bar", f(st.Tsat, 2) + " K"]);
    if (r.field) {
      var cl = r.field.column;
      items.push(["In the field", cl.interface != null ? "interface at " + PEC.sig(cl.interface, 4) + " m; " + cl.phase[0] + " at the start, " + cl.phase[cl.phase.length - 1] + " at the end" : (cl.notes[0] || "")]);
      items.push(["Pressure along the field", sg(cl.P[0]) + " to " + sg(cl.P[cl.P.length - 1]) + " bar"]);
    }
    $("pure-facts").innerHTML = items.map(function (it) { return "<li><span class='k'>" + it[0] + "</span><span class='v'>" + it[1] + "</span></li>"; }).join("");
  }
  var pureRows = null;
  function pureTable(cfg) {
    var r = pure.res, c = r.curve, head, rows;
    if (["PT", "Trho", "dH"].indexOf(pure.plot) >= 0) {
      head = ["T (K)", "P_sat (bar)", "ρ liquid (kg/m³)", "ρ vapor (kg/m³)", "ΔH vap (kJ/mol)", "P Antoine (bar)"];
      rows = c.T.map(function (t, i) { return [t, c.P[i], c.rhoL[i], c.rhoV[i], c.dHvap[i], c.P_antoine[i]]; });
    } else if (pure.plot === "PV") {
      head = ["T/Tc", "V/Vc", "P/Pc"]; rows = [];
      r.isotherms.forEach(function (it) { it.V.forEach(function (v, i) { if (i % 4 === 0) rows.push([it.Tr, v, it.P[i]]); }); });
    } else if (pure.plot === "band") {
      head = [r.field.kind === "gravity" ? "g (m/s²)" : "ω (rad/s)", "T (K)", "P_sat (bar)", "P first liquid (bar)", "P last vapor (bar)"]; rows = [];
      r.field.bands.forEach(function (b) { b.T.forEach(function (t, i) { rows.push([b.strength, t, b.Psat[i], b.P_first[i], b.P_last[i]]); }); });
    } else if (pure.plot === "iface" && r.field.interface_curve) {
      head = [r.field.kind === "gravity" ? "g (m/s²)" : "speed (rpm)", "interface position (m)"];
      rows = r.field.interface_curve.strength.map(function (s, i) { return [s, r.field.interface_curve.position[i]]; });
    } else {
      var cl = r.field.column; head = ["position (m)", "phase", "P (bar)", "ρ (kg/m³)"];
      rows = cl.position.map(function (x, i) { return [x, cl.phase[i], cl.P[i], cl.rho[i]]; });
    }
    pureRows = [head, rows];
    var fm = function (v) { return v == null ? "–" : typeof v === "number" ? PEC.sig(v, 5) : v; };
    PEC.dataTable("pure-table", head, rows.map(function (rw) { return rw.map(fm); }));
    $("pure-data-caption").textContent = r.name + ", " + MODEL_SHORT[r.model] + " (" + rows.length + " rows)";
  }
  function pureCsv() {
    if (!pureRows) return;
    var lines = ["# " + pure.res.name + "; " + pure.res.model].concat([pureRows[0].join(",")]).concat(pureRows[1].map(function (rw) { return rw.map(function (v) { return v == null ? "" : v; }).join(","); }));
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([lines.join("\n") + "\n"], { type: "text/csv" }));
    a.download = pure.compound + "_" + pure.model + "_" + pure.plot + ".csv";
    document.body.appendChild(a); a.click(); a.remove();
  }

  /* ------------------------------------------------------------ any input file */
  /* The command-line program in the page: the checker's report, the printed
     results, and the curves the solver computed (web_api.run_input). */
  var files = { list: [], loaded: false };
  function bindFiles() {
    $("file-run").addEventListener("click", runAnyFile);
    $("file-example").addEventListener("change", function () {
      var ex = files.list[+this.value];
      if (ex) { $("file-text").value = ex.text; runAnyFile(); }
    });
  }
  function loadFiles() {
    if (files.loaded) return;
    files.loaded = true;
    engine.call({ op: "examples" }).then(function (r) {
      files.list = r.examples || [];
      $("file-example").innerHTML = files.list.map(function (ex, i) {
        return "<option value='" + i + "'>" + ex.file + "  [" + ex.status + "]</option>"; }).join("");
      var first = files.list.filter(function (e) { return /pure_water/.test(e.file); })[0] || files.list[0];
      if (first) { $("file-example").value = String(files.list.indexOf(first)); $("file-text").value = first.text; runAnyFile(); }
    });
  }
  function runAnyFile() {
    var box = $("file-result"), text = $("file-text").value, host = $("file-charts");
    if (!text.trim()) { box.innerHTML = "<p class='notice'>Write an input file first.</p>"; return; }
    box.innerHTML = "<p class='source'>Checking and solving…</p>"; host.innerHTML = "";
    engine.call({ op: "run_input", text: text }).then(function (r) {
      if (!r.ok) { box.innerHTML = "<p class='notice'>" + (r.error || "Failed") + "</p>"; return; }
      var esc = function (t) { return t.replace(/&/g, "&amp;").replace(/</g, "&lt;"); };
      box.innerHTML = "<p class='check-verdict st-" + r.status.toLowerCase() + "'>" + r.status + "</p>" +
        "<h3>Checker report</h3><pre>" + esc(r.report) + "</pre>" + (r.results ? "<h3>Results</h3><pre>" + esc(r.results) + "</pre>" : "") +
        "<p class='source'>" + r.ms + " ms</p>";
      var ROLE = { liquid: "--liquid", vapor: "--vapor", third: "--third", solid2: "--ink-2", muted: "--muted" };
      (r.charts || []).forEach(function (ch) {
        var fig = document.createElement("figure"); fig.className = "chart-figure";
        fig.innerHTML = "<h3>" + esc(ch.title) + "</h3><ul class='legend'></ul><div class='chart'></div>";
        host.appendChild(fig);
        PEC.legend(fig.querySelector("ul"), ch.series.filter(function (se, i, a) {
          return a.findIndex(function (o) { return o.name === se.name; }) === i; }).map(function (se) { return { name: se.name, color: ROLE[se.role] || "--ink-2" }; })
          .concat(ch.points.length ? [{ name: ch.points.map(function (p) { return p.text; }).join(", "), color: "--ink", dot: true }] : []));
        PEC.drawChart(fig.querySelector(".chart"), {
          aria: ch.title, xTitle: ch.xTitle, yTitle: ch.yTitle, yOpts: { pad: 0.06 },
          series: ch.series.map(function (se) { return { name: se.name, color: ROLE[se.role] || "--ink-2", width: se.role === "muted" ? 1.2 : 2,
            dash: se.role === "muted" ? "4 3" : "none", pts: se.pts.map(function (q) { return q[0] == null || q[1] == null ? null : q; }) }; }),
          points: ch.points.map(function (p) { return { x: p.x, y: p.y, color: "--ink", text: p.text, dx: 8, dy: -6 }; }),
          fx: function (v) { return PEC.sig(v, 4); }, fy: function (v) { return PEC.sig(v, 4); }
        });
      });
    });
  }

  /* ------------------------------------------------------------ modes */
  var mode = "pure", started = false;
  function showShells(m) {
    $("pure-shell").hidden = m !== "pure"; $("main").hidden = m !== "binary"; $("mix-shell").hidden = m !== "mix";
    $("file-shell").hidden = m !== "file";
    [].forEach.call($("mode").querySelectorAll("button"), function (b) { b.setAttribute("aria-pressed", String(b.getAttribute("data-mode") === m)); });
  }
  function setMode(m, quiet) {
    mode = m;
    showShells(m);
    if (!started) return;                  // Python is still starting: the switch takes effect when it is ready
    if (m === "file") { history.replaceState(null, "", "#mode=file"); loadFiles(); }
    if (m === "pure") { history.replaceState(null, "", "#mode=pure"); if (pure.res) drawPure(); else if (!quiet) runPure(); }
    if (m === "mix") history.replaceState(null, "", "#mode=mix");
    if (m === "binary") { writeHash(); if (last) redraw(); }
  }
  function initialMode() {
    var h = location.hash.slice(1), p = new URLSearchParams(h);
    if (h === "mix-panel" || p.get("mode") === "mix") return "mix";
    if (p.get("mode") === "file") return "file";
    if (p.get("mode") === "pure" || !h) return "pure";
    return "binary";
  }

  /* ------------------------------------------------------------ start */
  /* the mode buttons work at once; until Python is ready a line under them says what is happening */
  mode = initialMode();
  showShells(mode);
  [].forEach.call($("mode").querySelectorAll("button"), function (b) { b.addEventListener("click", function () { setMode(b.getAttribute("data-mode")); }); });
  function boot(text, isError) {
    var el = $("boot");
    el.hidden = !text; el.textContent = text || ""; el.classList.toggle("error", !!isError);
  }
  boot("Starting Python in your browser… The first visit downloads about 12 MB and can take up to half a minute; later visits start in a few seconds.");
  engine = new Engine(function (text) {
      overlay(text, "The first visit downloads Python (about 12 MB) and can take up to half a minute. Later visits start in a few seconds.");
      boot(text + "… The first visit downloads about 12 MB and can take up to half a minute.");
    },
    function (err) {
      overlay("Python could not start in this browser", "Check your connection and reload the page. (" + err + ")", true);
      boot("Python could not start in this browser. Check your connection and reload the page. (" + err + ")", true);
    });
  engine.call({ op: "catalog" }).then(function (cat) {
    catalog = cat;
    var allCompounds = catalog.compounds.slice();
    catalog.compounds = catalog.compounds.filter(function (c) { return offered("compounds", c.id); });
    var urlParams = mode === "binary" ? readHash() : {};
    if (!comp(state.a)) state.a = DEF.a;
    if (!comp(state.b) || state.b === state.a) state.b = state.a === DEF.b ? DEF.a : DEF.b;
    if (!MODEL_SHORT[state.model] || !offered("models", state.model)) state.model = DEF.model;
    fitPair();
    if (state.compare && !offered("models", state.compare)) state.compare = "";
    if (!Y_OPTS.some(function (o) { return o[0] === state.y; })) { state.y = "T"; state.x = "xy"; }
    if (!X_OPTS.some(function (o) { return o[0] === state.x; })) state.x = "xy";
    fillCompounds(); bindAxes(); bindMixture(); bindPure(allCompounds); bindFiles();
    started = true;
    boot("");
    setMode(mode, true);
    if (mode === "pure") runPure();
    urlParams.__first = true;
    return loadPair(urlParams);
  });
})();
