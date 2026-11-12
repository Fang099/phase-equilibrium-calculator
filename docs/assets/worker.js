/* Runs the project's Python modules in the browser with Pyodide.
   The page sends JSON requests; web_api.handle() answers them. */
"use strict";
var PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/";
var handle = null;

function status(text) { postMessage({ type: "status", text: text }); }

var ready = (async function () {
  importScripts(PYODIDE + "pyodide.js");
  status("Starting Python in your browser");
  var py = await loadPyodide({ indexURL: PYODIDE });
  status("Loading NumPy");
  await py.loadPackage("numpy");
  status("Loading the calculator's Python modules");
  var manifest = await (await fetch("../py/manifest.json", { cache: "no-cache" })).json();
  py.FS.mkdirTree("/home/pyodide/data");
  await Promise.all(manifest.files.map(async function (f) {
    var r = await fetch("../py/" + f + "?v=" + manifest.version);
    if (!r.ok) throw new Error("could not load " + f);
    py.FS.writeFile("/home/pyodide/" + f, await r.text());
  }));
  py.runPython("import web_api");
  handle = py.globals.get("web_api").handle;
  postMessage({ type: "ready", version: manifest.version });
})().catch(function (err) {
  postMessage({ type: "fatal", error: String(err && err.message || err) });
  throw err;
});

onmessage = async function (e) {
  var id = e.data.id;
  try {
    await ready;
    var out = handle(JSON.stringify(e.data.msg));
    postMessage({ type: "result", id: id, result: JSON.parse(out) });
  } catch (err) {
    postMessage({ type: "result", id: id, result: { ok: false, error: String(err && err.message || err) } });
  }
};
