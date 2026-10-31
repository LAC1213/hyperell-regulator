/* Static viewer: computations and full precision stay in the exported JSON. */
"use strict";
const DIGITS = 10;
const DATA_VERSION = "symmetric-basis-1";
const COLORS = ["#bd5034", "#267f91", "#9267af", "#45956b", "#c28b24", "#3865a4", "#b65787"];
const $ = (id) => document.getElementById(id);
let catalog = null, current = null, selectedCurve = null, loadVersion = 0;
let edgeIndex = 0, sheetIndex = 0;
let selectionFrame = null;
let currentMesh = null;
const cache = new Map();
const meshCache = new Map();

function decimal(value) {
  const n = Number(value);
  if (!Number.isFinite(n)) throw new Error("Non-finite numerical value");
  return n === 0 ? "0.000000000" : n.toPrecision(DIGITS);
}

function status(message, error = false) {
  $("status").textContent = message;
  $("status").classList.toggle("error", error);
}

function cell(row, value, tag = "td") {
  const node = document.createElement(tag);
  node.textContent = value;
  row.append(node);
  return node;
}

function option(select, value, label) {
  const node = document.createElement("option");
  node.value = value;
  node.textContent = label;
  select.append(node);
}

function complexCell(row, value, reverse) {
  const node = cell(row, "");
  node.className = "complex";
  const sign = reverse ? -1 : 1;
  const re = document.createElement("span"), im = document.createElement("span");
  re.textContent = decimal(sign * Number(value[0]));
  im.textContent = (sign * Number(value[1]) < 0 ? "− " : "+ ") + decimal(Math.abs(Number(value[1]))) + " i";
  node.append(re, im);
}

function validate(data) {
  const number = (v) => (typeof v === "string" || typeof v === "number") && String(v).length < 400 && Number.isFinite(Number(v)) && String(v).trim() !== "";
  const complex = (v) => Array.isArray(v) && v.length === 2 && v.every(number);
  if (!data || data.schema_version !== 1 || !Number.isInteger(data.genus) || data.genus < 2 || data.genus > 5 ||
      !Number.isInteger(data.degree) || data.degree < 2 || data.degree > 24 || !data.input ||
      !["F", "A", "B"].every(k => typeof data.input[k] === "string" && data.input[k].length <= 4000) ||
      !complex(data.jump) || !Array.isArray(data.edges) || !data.edges.length || data.edges.length > 100 ||
      !Array.isArray(data.pairings) || data.pairings.length !== data.genus * (data.genus + 1) / 2) {
    throw new Error("This is not a supported cut-system result (schema version 1)");
  }
  data.edges.forEach((edge, index) => {
    if (edge.index !== index + 1 || !Array.isArray(edge.base) || edge.base.length < 2 || edge.base.length > 2000 ||
        !edge.base.every(complex) || !Array.isArray(edge.lifts) || edge.lifts.length !== data.degree) throw new Error("Invalid cut geometry");
    edge.lifts.forEach((lift, sheet) => {
      if (lift.sheet !== sheet + 1 || !Number.isInteger(lift.minus_sheet) || lift.minus_sheet < 1 || lift.minus_sheet > data.degree ||
          !["ordinary", "log_plus", "log_minus"].every(k => Array.isArray(lift[k]) && lift[k].length === data.genus && lift[k].every(complex)) ||
          !Array.isArray(lift.xyz) || lift.xyz.length > 100000 ||
          !lift.xyz.every(p => p === null || (Array.isArray(p) && p.length === 3 && p.every(number)))) throw new Error("Invalid lifted path or integrals");
    });
  });
  data.pairings.forEach(p => {
    if (!Array.isArray(p.forms) || p.forms.length !== 2 || !p.forms.every(n => Number.isInteger(n) && n >= 1 && n <= data.genus) ||
        ![p.period, p.regulator, p.pairing].every(number)) throw new Error("Invalid regulator pairing");
  });
  return data;
}

function validateMesh(mesh, data) {
  if (!mesh || mesh.schema_version !== 1 || mesh.degree !== data.degree ||
      !Array.isArray(mesh.sheets) || mesh.sheets.length !== data.degree ||
      !Array.isArray(mesh.triangles) || !mesh.triangles.length || mesh.triangles.length > 30000) {
    throw new Error("Unsupported surface mesh");
  }
  const count = mesh.sheets[0]?.vertices?.length;
  if (!count || count > 100000 || !mesh.triangles.every(t => Array.isArray(t) && t.length === 3 &&
      t.every(i => Number.isInteger(i) && i >= 0 && i < count))) throw new Error("Invalid mesh triangles");
  mesh.sheets.forEach((sheet, i) => {
    if (sheet.sheet !== i + 1 || !Array.isArray(sheet.vertices) || sheet.vertices.length !== count ||
        !sheet.vertices.every(p => Array.isArray(p) && p.length === 3 && p.every(Number.isFinite))) {
      throw new Error("Invalid mesh vertices or sheet labels");
    }
  });
  return mesh;
}

async function fetchJSON(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`Could not load ${path} (HTTP ${response.status})`);
  return response.json();
}

function renderCatalog() {
  $("curve-select").replaceChildren();
  $("catalog-rows").replaceChildren();
  catalog.curves.forEach(curve => {
    option($("curve-select"), curve.label, curve.label);
    const row = document.createElement("tr");
    row.dataset.label = curve.label;
    const first = cell(row, "");
    const button = document.createElement("button");
    button.textContent = curve.label;
    button.type = "button";
    button.addEventListener("click", () => chooseCurve(curve.label));
    first.append(button);
    cell(row, `(${curve.points[0].join(", ")})`);
    cell(row, `(${curve.points[1].join(", ")})`);
    cell(row, curve.orders[0]);
    cell(row, curve.orders[1]);
    cell(row, decimal(curve.determinant));
    cell(row, curve.q + (curve.integrality_verified ? "" : " †")); // Exact rational torsion/comparison data is never rounded.
    $("catalog-rows").append(row);
  });
  $("curve-select").disabled = false;
}

function chooseCurve(label, functionName = "P") {
  if (functionName === "Q") functionName = "R"; // Preserve older saved links.
  selectedCurve = catalog.curves.find(c => c.label === label);
  if (!selectedCurve) { status("This curve is not in the catalog.", true); return; }
  $("curve-select").value = label;
  $("function-select").replaceChildren();
  selectedCurve.functions.forEach(f => option($("function-select"), f.name, `φ_${f.name} of order ${f.order}`));
  $("function-select").value = selectedCurve.functions.some(f => f.name === functionName) ? functionName : "P";
  $("function-select").disabled = false;
  for (const row of $("catalog-rows").children) row.classList.toggle("selected", row.dataset.label === label);
  return loadFunction();
}

async function loadFunction() {
  const curve = selectedCurve;
  const f = curve.functions.find(f => f.name === $("function-select").value);
  const version = ++loadVersion;
  status(`Loading ${curve.label}, φ_${f.name}…`);
  try {
    let data = cache.get(f.file);
    if (!data) { data = validate(await fetchJSON(`${f.file}?v=${DATA_VERSION}`)); cache.set(f.file, data); }
    let mesh = null, meshError = null;
    if (f.mesh) {
      try {
        mesh = meshCache.get(f.mesh);
        if (!mesh) {
          mesh = validateMesh(await fetchJSON(`${f.mesh}?v=${DATA_VERSION}`), data);
          meshCache.set(f.mesh, mesh);
          while (meshCache.size > 3) meshCache.delete(meshCache.keys().next().value);
        }
      } catch (error) { meshError = error; }
    }
    if (version !== loadVersion) return;
    await show(data, curve, f, mesh);
    if (version !== loadVersion) return;
    if (meshError) status(`Integrals loaded. Surface mesh unavailable: ${meshError.message}`, true);
    const url = new URL(location.href);
    url.hash = `${curve.label}/${f.name}`;
    history.replaceState(null, "", url);
  } catch (error) {
    if (version === loadVersion) status(error.message, true);
  }
}

async function show(data, curve = null, f = null, mesh = null) {
  if (selectionFrame !== null) cancelAnimationFrame(selectionFrame);
  selectionFrame = null;
  current = data;
  currentMesh = mesh;
  edgeIndex = sheetIndex = 0;
  $("results").hidden = false;
  $("export").disabled = false;
  $("result-source").textContent = curve ? `${curve.label} / φ_${f.name}` : "IMPORTED LOCAL RESULT";
  $("result-model").textContent = `y² = ${data.input.F}`;
  $("metadata").textContent = `genus ${data.genus} · ${data.degree} sheets · ${data.edges.length} edges · 10 digits shown`;
  $("exact-inputs").textContent = `F(x) = ${data.input.F}\nA(x) = ${data.input.A}\nB(x) = ${data.input.B}\nφ = A(x) + B(x)y`;
  $("point-info").textContent = curve ? `(φ_${f.name}) = ${f.order} (${f.name} − Q), with ${f.name} = (${f.point.join(", ")}) and Q the common pole at infinity. Normalized determinant: ${decimal(curve.determinant)}.` : "An imported result computed locally with hyperell-regulator.";
  $("edge").replaceChildren();
  data.edges.forEach((e, i) => option($("edge"), i, `σ_${e.index}${e.infinite ? " → ∞" : ""}`));
  $("sheet").replaceChildren();
  for (let s = 0; s < data.degree; s++) option($("sheet"), s, `S_${s + 1}`);
  $("visible-sheet").replaceChildren();
  option($("visible-sheet"), 0, "All sheets");
  for (let s = 1; s <= data.degree; s++) option($("visible-sheet"), s, `Only S_${s}`);
  $("show-mesh").disabled = !mesh;
  $("mesh-caption").textContent = mesh ?
    "Select a colored lift. A single-sheet view includes both banks of its boundary. Projection crossings can occur; surface triangles approach cut endpoints." :
    "Select a colored lift. A single-sheet view includes both banks of its boundary. This result has no surface mesh.";
  $("pairings").replaceChildren();
  data.pairings.forEach(p => {
    const row = document.createElement("tr");
    cell(row, `(${p.forms.join(", ")})`);
    [p.pairing, p.regulator, p.period].forEach(v => cell(row, decimal(v)));
    $("pairings").append(row);
  });
  $("geometry").textContent = data.conventions?.geometry || "The ray is truncated for display. Its integral includes the complete ray to infinity.";
  $("divisor").textContent = data.divisor ? `Verified divisor: N = ${data.divisor.N}; zero: ${data.divisor.P}; pole Q: ${data.divisor.Q}.` : "";
  $("warnings").textContent = (data.warnings || []).join(" ");
  renderIntegrals();
  if (typeof Plotly === "undefined") {
    $("base-plot").textContent = $("lift-plot").textContent = "The plotting library could not load. Use the path selectors below to inspect every integral.";
    status("Data loaded. Plotting library unavailable; path selectors and tables remain usable.", true);
    return;
  }
  await renderPlots();
  status(data.warnings?.length ? "Values loaded. Some display samples are unavailable; see display limits below." : "Ready. Select a cut or lifted path, or use the path selectors below.", Boolean(data.warnings?.length));
}

function renderIntegrals() {
  if (!current) return;
  const edge = current.edges[edgeIndex], lift = edge.lifts[sheetIndex];
  const bank = $("bank").value, reverse = $("orientation").value === "reverse";
  $("path-label").textContent = `e${bank === "plus" ? "+" : "−"}_${bank === "plus" ? lift.sheet : lift.minus_sheet},${edge.index}`;
  $("gluing").textContent = `The same geometric lift bounds S_${lift.sheet} on the positive bank and S_${lift.minus_sheet} on the negative bank. ${edge.infinite ? "This edge ends at infinity; all integral values include the entire tail." : ""}`;
  $("integrals").replaceChildren();
  const logarithmic = bank === "plus" ? lift.log_plus : lift.log_minus;
  for (let a = 0; a < current.genus; a++) {
    const row = document.createElement("tr");
    cell(row, `ω_${a + 1}`);
    complexCell(row, lift.ordinary[a], reverse);
    complexCell(row, logarithmic[a], reverse);
    $("integrals").append(row);
  }
  $("jump").textContent = `L⁺ = L⁻ + Δlog(φ) I, with Δlog(φ) = ${decimal(current.jump[1])} i. Orientation: ${reverse ? "toward 0" : "toward ∞"}.`;
}

const plotConfig = {responsive: true, displaylogo: false, scrollZoom: true, modeBarButtonsToRemove: ["toImage"]};
const axes = {zerolinecolor: "#d5dcd2", gridcolor: "#e8ebe3", tickformat: ".4g", titlefont: {size: 11}};

function pathVisible(lift) {
  const sheet = Number($("visible-sheet").value);
  return sheet === 0 || lift.sheet === sheet || lift.minus_sheet === sheet;
}

async function renderPlots() {
  const data = current;
  const base = data.edges.map((edge, i) => ({
    x: edge.base.map(p => p[0]), y: edge.base.map(p => p[1]), type: "scatter", mode: "lines+markers",
    line: {color: i === edgeIndex ? "#bc542e" : "#44736a", width: i === edgeIndex ? 4 : 2}, marker: {size: 5},
    customdata: edge.base.map(() => i), name: `σ_${edge.index}${edge.infinite ? " → ∞" : ""}`,
    hovertemplate: `σ_${edge.index}<br>Re w = %{x:.10g}<br>Im w = %{y:.10g}<extra></extra>`
  }));
  const critical = data.edges.map(e => e.base[0]);
  base.push({x: critical.map(p => p[0]), y: critical.map(p => p[1]), mode: "markers", type: "scatter",
    marker: {size: 7, color: "#183e3b"},
    customdata: data.edges.map((_, i) => i), hovertemplate: critical.map((_, i) => `z_${i + 1}<br>Re w = %{x:.10g}<br>Im w = %{y:.10g}<extra></extra>`), showlegend: false});
  const lifts = [];
  if (currentMesh) currentMesh.sheets.forEach(sheet => {
    lifts.push({type: "mesh3d", x: sheet.vertices.map(p => p[0]), y: sheet.vertices.map(p => p[1]), z: sheet.vertices.map(p => p[2]),
      i: currentMesh.triangles.map(t => t[0]), j: currentMesh.triangles.map(t => t[1]), k: currentMesh.triangles.map(t => t[2]),
      color: COLORS[(sheet.sheet - 1) % COLORS.length], opacity: 0.2, flatshading: true,
      name: `Surface S_${sheet.sheet}`, hoverinfo: "skip", showlegend: false,
      visible: $("show-mesh").checked, meta: {kind: "surface", sheet: sheet.sheet}});
  });
  data.edges.forEach((edge, i) => edge.lifts.forEach((lift, s) => {
    const chosen = i === edgeIndex && s === sheetIndex;
    lifts.push({type: "scatter3d", mode: "lines", x: lift.xyz.map(p => p && p[0]), y: lift.xyz.map(p => p && p[1]), z: lift.xyz.map(p => p && p[2]),
      line: {color: COLORS[s % COLORS.length], width: chosen ? 8 : 3}, opacity: chosen ? 1 : 0.65,
      visible: pathVisible(lift), meta: {kind: "path", edge: i, sheet: lift.sheet, minus_sheet: lift.minus_sheet},
      customdata: lift.xyz.map(() => [i, s]), name: `e⁺_${s + 1},${edge.index}`, showlegend: false,
      hovertemplate: `e⁺_${s + 1},${edge.index} = e⁻_${lift.minus_sheet},${edge.index}<br>Re w = %{x:.10g}<br>Im w = %{y:.10g}<br>height = %{z:.10g}<extra></extra>`});
  }));
  const shared = {paper_bgcolor: "#fffefa", plot_bgcolor: "#fffefa", font: {family: "sans-serif", color: "#183e3b", size: 11}, margin: {l: 50, r: 15, t: 15, b: 45}, showlegend: false, uirevision: `${data.input.F}/${data.input.A}/${data.input.B}`};
  await Promise.all([
    Plotly.react($("base-plot"), base, {...shared, xaxis: {...axes, title: "Re w"}, yaxis: {...axes, title: "Im w", scaleanchor: "x", scaleratio: 1}, hovermode: "closest"}, plotConfig),
    Plotly.react($("lift-plot"), lifts, {...shared, margin: {l: 0, r: 0, t: 0, b: 0}, scene: {xaxis: {...axes, title: "Re w"}, yaxis: {...axes, title: "Im w"}, zaxis: {...axes, title: "height"}, aspectmode: "cube", camera: {eye: {x: 1.5, y: -1.6, z: 1.2}}}}, plotConfig)
  ]);
  for (const id of ["base-plot", "lift-plot"]) $(id).removeAllListeners("plotly_click");
  $("base-plot").on("plotly_click", event => {
    if (typeof event.points[0].customdata === "number") selectPath(event.points[0].customdata, sheetIndex);
  });
}

function pickLift(clientX, clientY) {
  const plot = $("lift-plot"), scene = plot._fullLayout?.scene?._scene;
  if (!scene?.glplot?.cameraParams) return null;
  // Project the actual displayed polylines. GL's surface picking otherwise
  // intercepts clicks on paths visible through a translucent mesh. These
  // matrices belong to the pinned, bundled Plotly version.
  const rect = plot.getBoundingClientRect(), camera = scene.glplot.cameraParams;
  if (![camera.model, camera.view, camera.projection].every(m => m?.length === 16)) return null;
  const multiply = (matrix, vector) => [0, 1, 2, 3].map(i =>
    vector.reduce((sum, value, j) => sum + matrix[j * 4 + i] * value, 0));
  const project = point => {
    if (!point) return null;
    let v = [...point.map((x, i) => x * scene.dataScale[i]), 1];
    for (const matrix of [camera.model, camera.view, camera.projection]) v = multiply(matrix, v);
    if (v[3] <= 0 || Math.abs(v[2] / v[3]) > 1) return null;
    return [rect.left + (v[0] / v[3] + 1) * rect.width / 2,
      rect.top + (1 - v[1] / v[3]) * rect.height / 2, v[2] / v[3]];
  };
  let best = null, distance = 8 * 8, depth = Infinity;
  current.edges.forEach((edge, e) => edge.lifts.forEach((lift, s) => {
    if (!pathVisible(lift)) return;
    let previous = null;
    for (const point of lift.xyz) {
      const next = project(point);
      if (previous && next) {
        const dx = next[0] - previous[0], dy = next[1] - previous[1];
        const length = dx * dx + dy * dy;
        const t = length ? Math.max(0, Math.min(1,
          ((clientX - previous[0]) * dx + (clientY - previous[1]) * dy) / length)) : 0;
        const d = (clientX - previous[0] - t * dx) ** 2 + (clientY - previous[1] - t * dy) ** 2;
        const z = previous[2] + t * (next[2] - previous[2]);
        if (d < distance - .1 || (d <= distance && z < depth)) {
          distance = d; depth = z; best = [e, s];
        }
      }
      previous = next;
    }
  }));
  return best;
}

let pointerStart = null;
$("lift-plot").addEventListener("pointerdown", event => {
  if (event.button !== 0 || event.target.closest(".modebar")) return;
  pointerStart = {id: event.pointerId, x: event.clientX, y: event.clientY, dragged: false};
}, true);
$("lift-plot").addEventListener("pointermove", event => {
  if (pointerStart && Math.hypot(event.clientX - pointerStart.x, event.clientY - pointerStart.y) > 5) pointerStart.dragged = true;
}, true);
$("lift-plot").addEventListener("pointercancel", () => { pointerStart = null; }, true);
$("lift-plot").addEventListener("pointerup", event => {
  const start = pointerStart;
  pointerStart = null;
  if (!start || start.id !== event.pointerId || start.dragged || !current) return;
  const data = current, x = event.clientX, y = event.clientY;
  requestAnimationFrame(() => {
    if (current !== data) return;
    const path = pickLift(x, y);
    if (path) selectPath(path[0], path[1]);
  });
}, true);

function selectPath(edge, sheet) {
  if (!current || !current.edges[edge]?.lifts[sheet]) return;
  edgeIndex = edge; sheetIndex = sheet;
  let visibleSheet = Number($("visible-sheet").value);
  const lift = current.edges[edge].lifts[sheet];
  if (visibleSheet && lift.sheet !== visibleSheet && lift.minus_sheet !== visibleSheet) {
    visibleSheet = lift.sheet;
    $("visible-sheet").value = String(visibleSheet);
  }
  if (visibleSheet) {
    $("bank").value = lift.sheet === visibleSheet ? "plus" : "minus";
  }
  $("edge").value = edge; $("sheet").value = sheet;
  renderIntegrals();
  if (typeof Plotly !== "undefined") {
    if (selectionFrame !== null) cancelAnimationFrame(selectionFrame);
    // Rebuilding a GL scene inside plotly_click interrupts the browser's
    // mouse-release processing. Defer all redraws until the event unwinds.
    selectionFrame = requestAnimationFrame(() => {
      selectionFrame = null;
      Plotly.restyle($("base-plot"), {"line.width": current.edges.map((_, i) => i === edgeIndex ? 4 : 2), "line.color": current.edges.map((_, i) => i === edgeIndex ? "#bc542e" : "#44736a")}, current.edges.map((_, i) => i));
      const plot = $("lift-plot");
      const indices = [], widths = [], opacities = [];
      plot.data.forEach((trace, i) => {
        if (trace.meta?.kind !== "path") return;
        const chosen = trace.meta.edge === edgeIndex && trace.meta.sheet === sheetIndex + 1;
        indices.push(i); widths.push(chosen ? 8 : 3); opacities.push(chosen ? 1 : .65);
      });
      Plotly.restyle(plot, {"line.width": widths, opacity: opacities}, indices);
      updateSheetView();
    });
  }
}

function updateSheetView() {
  if (!current || typeof Plotly === "undefined" || !$("lift-plot").data) return;
  const plot = $("lift-plot"), sheet = Number($("visible-sheet").value);
  const visible = plot.data.map(trace => trace.meta.kind === "surface" ?
    $("show-mesh").checked && (!sheet || trace.meta.sheet === sheet) :
    !sheet || trace.meta.sheet === sheet || trace.meta.minus_sheet === sheet);
  if (visible.some((value, i) => value !== (plot.data[i].visible !== false))) Plotly.restyle(plot, {visible});
}

$("visible-sheet").addEventListener("change", () => {
  if (!current) return;
  const sheet = Number($("visible-sheet").value);
  const lift = current.edges[edgeIndex].lifts[sheetIndex];
  selectPath(edgeIndex, sheet && lift.sheet !== sheet && lift.minus_sheet !== sheet ? sheet - 1 : sheetIndex);
});
$("show-mesh").addEventListener("change", () => requestAnimationFrame(updateSheetView));
$("curve-select").addEventListener("change", () => chooseCurve($("curve-select").value));
$("function-select").addEventListener("change", loadFunction);
$("edge").addEventListener("change", () => {
  const edge = Number($("edge").value), sheet = Number($("visible-sheet").value);
  const lift = current.edges[edge].lifts[sheetIndex];
  selectPath(edge, sheet && lift.sheet !== sheet && lift.minus_sheet !== sheet ? sheet - 1 : sheetIndex);
});
$("sheet").addEventListener("change", () => selectPath(edgeIndex, Number($("sheet").value)));
$("bank").addEventListener("change", renderIntegrals);
$("orientation").addEventListener("change", renderIntegrals);
$("export").addEventListener("click", () => {
  if (!current) return;
  const url = URL.createObjectURL(new Blob([JSON.stringify(current, null, 2)], {type: "application/json"}));
  const a = document.createElement("a");
  a.href = url; a.download = current.collection ? `${current.collection.entry.label}-${current.collection.function}.json` : "cut-system.json";
  a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
});
$("import").addEventListener("change", async event => {
  const file = event.target.files[0];
  if (!file) return;
  const version = ++loadVersion;
  try {
    if (file.size > 12000000) throw new Error("Import a result smaller than 12 MB");
    const data = validate(JSON.parse(await file.text()));
    if (version !== loadVersion) return;
    $("curve-select").value = "";
    $("function-select").disabled = true;
    for (const row of $("catalog-rows").children) row.classList.remove("selected");
    await show(data);
  } catch (error) { status(error.message, true); }
  event.target.value = "";
});

async function initialize() {
  try {
    catalog = await fetchJSON(`catalog.json?v=${DATA_VERSION}`);
    if (catalog.schema_version !== 1 || !Array.isArray(catalog.curves)) throw new Error("Unsupported curve catalog");
    renderCatalog();
    const hash = decodeURIComponent(location.hash.slice(1)).split("/");
    const label = catalog.curves.some(c => c.label === hash[0]) ? hash[0] : "763.a.763.1";
    await chooseCurve(label, hash[1] || "R");
  } catch (error) { status(error.message + ". Serve this folder over HTTP when viewing locally.", true); }
}
initialize();
