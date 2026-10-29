/* Static viewer: computations and full precision stay in the exported JSON. */
"use strict";
const DIGITS = 10;
const COLORS = ["#bd5034", "#267f91", "#9267af", "#45956b", "#c28b24", "#3865a4", "#b65787"];
const $ = (id) => document.getElementById(id);
let catalog = null, current = null, selectedCurve = null, loadVersion = 0;
let edgeIndex = 0, sheetIndex = 0;
const cache = new Map();

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
    cell(row, curve.q); // Exact rational torsion/comparison data is never rounded.
    $("catalog-rows").append(row);
  });
  $("curve-select").disabled = false;
}

function chooseCurve(label, functionName = "P") {
  selectedCurve = catalog.curves.find(c => c.label === label);
  if (!selectedCurve) { status("This curve is not in the paper catalog.", true); return; }
  $("curve-select").value = label;
  $("function-select").replaceChildren();
  selectedCurve.functions.forEach(f => option($("function-select"), f.name, `φ_${f.name} · N = ${f.order} · (${f.point.join(", ")})`));
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
    if (!data) { data = validate(await fetchJSON(f.file)); cache.set(f.file, data); }
    if (version !== loadVersion) return;
    await show(data, curve, f);
    if (version !== loadVersion) return;
    const url = new URL(location.href);
    url.hash = `${curve.label}/${f.name}`;
    history.replaceState(null, "", url);
  } catch (error) {
    if (version === loadVersion) status(error.message, true);
  }
}

async function show(data, curve = null, f = null) {
  current = data;
  edgeIndex = sheetIndex = 0;
  $("results").hidden = false;
  $("export").disabled = false;
  $("result-source").textContent = curve ? `${curve.label} / φ_${f.name}` : "IMPORTED LOCAL RESULT";
  $("result-model").textContent = `y² = ${data.input.F}`;
  $("metadata").textContent = `genus ${data.genus} · ${data.degree} sheets · ${data.edges.length} edges · 10 digits shown`;
  $("exact-inputs").textContent = `F(x) = ${data.input.F}\nA(x) = ${data.input.A}\nB(x) = ${data.input.B}\nφ = A(x) + B(x)y`;
  $("point-info").textContent = curve ? `div(φ_${f.name}) = ${f.order} (${f.name} − ${data.paper.pole}), with ${f.name} = (${f.point.join(", ")}). Fresh normalized determinant: ${decimal(curve.determinant)}; agrees with the printed paper value.` : "An imported result computed locally with hyperell-regulator.";
  $("edge").replaceChildren();
  data.edges.forEach((e, i) => option($("edge"), i, `σ_${e.index}${e.infinite ? " → ∞" : ""}`));
  $("sheet").replaceChildren();
  for (let s = 0; s < data.degree; s++) option($("sheet"), s, `S_${s + 1}`);
  $("pairings").replaceChildren();
  data.pairings.forEach(p => {
    const row = document.createElement("tr");
    cell(row, `(${p.forms.join(", ")})`);
    [p.pairing, p.regulator, p.period].forEach(v => cell(row, decimal(v)));
    $("pairings").append(row);
  });
  $("geometry").textContent = data.conventions?.geometry || "The ray is truncated for display. Its integral includes the complete ray to infinity.";
  $("divisor").textContent = data.divisor ? `Verified divisor: N = ${data.divisor.N}; zero: ${data.divisor.P}; pole: ${data.divisor.Q}.` : "";
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
  data.edges.forEach((edge, i) => edge.lifts.forEach((lift, s) => {
    const chosen = i === edgeIndex && s === sheetIndex;
    lifts.push({type: "scatter3d", mode: "lines", x: lift.xyz.map(p => p && p[0]), y: lift.xyz.map(p => p && p[1]), z: lift.xyz.map(p => p && p[2]),
      line: {color: COLORS[s % COLORS.length], width: chosen ? 8 : 3}, opacity: chosen ? 1 : 0.65,
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
  $("lift-plot").on("plotly_click", event => {
    const path = event.points[0].customdata;
    if (Array.isArray(path)) selectPath(path[0], path[1]);
  });
}

function selectPath(edge, sheet) {
  edgeIndex = edge; sheetIndex = sheet;
  $("edge").value = edge; $("sheet").value = sheet;
  renderIntegrals();
  if (typeof Plotly !== "undefined") {
    const trace = edgeIndex * current.degree + sheetIndex;
    const count = current.degree * current.edges.length;
    Plotly.restyle($("base-plot"), {"line.width": current.edges.map((_, i) => i === edgeIndex ? 4 : 2), "line.color": current.edges.map((_, i) => i === edgeIndex ? "#bc542e" : "#44736a")}, current.edges.map((_, i) => i));
    Plotly.restyle($("lift-plot"), {"line.width": Array.from({length: count}, (_, i) => i === trace ? 8 : 3), opacity: Array.from({length: count}, (_, i) => i === trace ? 1 : 0.65)});
  }
}

$("curve-select").addEventListener("change", () => chooseCurve($("curve-select").value));
$("function-select").addEventListener("change", loadFunction);
$("edge").addEventListener("change", () => selectPath(Number($("edge").value), sheetIndex));
$("sheet").addEventListener("change", () => selectPath(edgeIndex, Number($("sheet").value)));
$("bank").addEventListener("change", renderIntegrals);
$("orientation").addEventListener("change", renderIntegrals);
$("export").addEventListener("click", () => {
  if (!current) return;
  const url = URL.createObjectURL(new Blob([JSON.stringify(current, null, 2)], {type: "application/json"}));
  const a = document.createElement("a");
  a.href = url; a.download = current.paper ? `${current.paper.entry.label}-${current.paper.function}.json` : "cut-system.json";
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
    catalog = await fetchJSON("catalog.json");
    if (catalog.schema_version !== 1 || !Array.isArray(catalog.curves)) throw new Error("Unsupported curve catalog");
    renderCatalog();
    const hash = decodeURIComponent(location.hash.slice(1)).split("/");
    const label = catalog.curves.some(c => c.label === hash[0]) ? hash[0] : "763.a.763.1";
    await chooseCurve(label, hash[1] || "Q");
  } catch (error) { status(error.message + ". Serve this folder over HTTP when viewing locally.", true); }
}
initialize();
