/* Shared tree and two-sheet views per curve; numerical data are precomputed. */
"use strict";
const $ = id => document.getElementById(id);
const COLORS = ["#bd5034", "#267f91"];
const DATA_VERSION = "mst-current-1";
const HOME_CAMERA = {eye:{x:1.5,y:-1.6,z:1.2},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1},projection:{type:"perspective"}};
let catalog, current, mesh = null, version = 0, label = "";
const decimal = v => Number(v) === 0 ? "0" : Number(v).toPrecision(10);
const pointTuple = p => Array.isArray(p)?"("+p.map(String).join(", ")+")":p==="infinity"?"∞":String(p);
function status(message, error=false) { $("status").textContent=message; $("status").classList.toggle("error",error); }
function option(select,value,text) { const n=document.createElement("option"); n.value=value; n.textContent=text; select.append(n); }
function cell(row,text) { const n=document.createElement("td"); n.textContent=text; row.append(n); return n; }
function complex(z,sign=1) { return decimal(sign*Number(z[0]))+(sign*Number(z[1])<0?" − ":" + ")+decimal(Math.abs(Number(z[1])))+" i"; }
function validate(data) {
  const number=v => (typeof v === "string" || typeof v === "number") && String(v).length<400 && String(v).trim()!=="" && Number.isFinite(Number(v));
  const pair=z => Array.isArray(z) && z.length===2 && z.every(number);
  if (!data || data.schema_version!==2 || data.degree!==2 || !Number.isInteger(data.genus) || data.genus<2 || data.genus>5 ||
      !data.input || typeof data.input.F!=="string" || data.input.F.length>4000 || !Array.isArray(data.weights) || !data.weights.length || data.weights.length>8 ||
      !Array.isArray(data.edges) || !data.edges.length || data.edges.length>200) throw Error("Expected a direct two-sheet result (schema version 2)");
  if (!Array.isArray(data.roots) || !data.roots.every(pair)) throw Error("Invalid branch values");
  if (data.boundary_word !== undefined) {
    const word=data.boundary_word;
    if (!Array.isArray(word) || word.length!==2*data.edges.length || !word.every(b=>Array.isArray(b)&&b.length===2&&Number.isInteger(b[0])&&b[0]>=0&&b[0]<data.edges.length&&[1,-1].includes(b[1])) || new Set(word.map(b=>b.join(","))).size!==word.length) throw Error("Invalid tree boundary order");
  }
  data.weights.forEach(w => {
    if (!number(w.center) || typeof w.label!=="string" || !Array.isArray(w.pairings) || w.pairings.length!==data.genus*(data.genus+1)/2) throw Error("Invalid logarithmic weight");
    w.pairings.forEach(p => { if (![p.period,p.regulator,p.pairing].every(number) || !Array.isArray(p.forms) || p.forms.length!==2 || !p.forms.every(v => Number.isInteger(v)&&v>=1&&v<=data.genus)) throw Error("Invalid pairing"); });
  });
  data.edges.forEach((e,i) => {
    if (e.index!==i+1 || !Array.isArray(e.base) || e.base.length<2 || e.base.length>2000 || !e.base.every(pair) || !Array.isArray(e.jumps) || e.jumps.length!==data.weights.length || !e.jumps.every(pair) || !Array.isArray(e.lifts) || e.lifts.length!==2) throw Error("Invalid cut edge");
    e.lifts.forEach((l,s) => {
      if (l.sheet!==s+1 || ![1,2].includes(l.minus_sheet) || !Array.isArray(l.ordinary) || l.ordinary.length!==data.genus || !l.ordinary.every(pair) || !Array.isArray(l.logs) || l.logs.length!==data.weights.length || !Array.isArray(l.xyz) || l.xyz.length>10000 || !l.xyz.every(p => Array.isArray(p)&&p.length===3&&p.every(number))) throw Error("Invalid lift");
      l.logs.forEach(v => { if (!["minus","plus"].every(k=>Array.isArray(v[k])&&v[k].length===data.genus&&v[k].every(pair))) throw Error("Invalid logarithmic integral"); });
    });
  });
  return data;
}
function validateMesh(m) {
  if (!m || m.degree!==2 || !Array.isArray(m.triangles) || m.triangles.length>30000 || !Array.isArray(m.sheets) || m.sheets.length!==2) throw Error("Invalid surface mesh");
  m.sheets.forEach(s=>{ if (!Array.isArray(s.vertices) || !s.vertices.every(p=>Array.isArray(p)&&p.length===3&&p.every(Number.isFinite)) || !m.triangles.every(t=>Array.isArray(t)&&t.length===3&&t.every(i=>Number.isInteger(i)&&i>=0&&i<s.vertices.length))) throw Error("Invalid surface vertices"); });
  return m;
}
function renderIntegrals() {
  if (!current) return;
  const k=Number($("weight-select").value), e=Number($("edge").value), s=Number($("sheet").value);
  const edge=current.edges[e], lift=edge.lifts[s], sign=$("orientation").value==="reverse"?-1:1;
  $("integrals").replaceChildren(); $("pairings").replaceChildren();
  lift.ordinary.forEach((v,a)=>{const row=document.createElement("tr"); cell(row,"ω"+(a+1)); cell(row,complex(v,sign)); cell(row,complex(lift.logs[k][$("bank").value][a],sign)); $("integrals").append(row);});
  current.weights[k].pairings.forEach(p=>{const row=document.createElement("tr"); cell(row,p.forms.join(", ")); [p.pairing,p.regulator,p.period].forEach(v=>cell(row,decimal(v))); $("pairings").append(row);});
  $("path-label").textContent="Edge "+(e+1)+", sheet "+(s+1);
  $("gluing").textContent="Positive bank sheet "+(s+1)+" meets negative bank sheet "+lift.minus_sheet+".";
  $("jump").textContent="Logarithm jump on this edge: "+complex(edge.jumps[k])+". Both banks use the selected geometric orientation.";
}
function saveCamera() {
  const plot=$("curve-plot"), scene=plot._fullLayout?.scene?._scene;
  if (!scene) return structuredClone(HOME_CAMERA);
  scene.saveLayout(plot.layout);
  return structuredClone(scene.getCamera());
}
function restoreCamera(camera) {
  const plot=$("curve-plot"), scene=plot._fullLayout?.scene?._scene;
  if (!scene) return;
  scene.setViewport({...plot._fullLayout.scene,camera});
  scene.saveLayout(plot.layout);
  scene.glplot.redraw();
}
async function resetCamera() {
  const camera=structuredClone(HOME_CAMERA);
  await Plotly.relayout($("curve-plot"),{"scene.camera":camera});
  restoreCamera(camera);
}
// Read the live GL camera before a mode switch can restore an older layout.
let modeCamera=null;
$("curve-plot").addEventListener("click",event=>{
  modeCamera=null;
  const button=event.target.closest(".modebar-btn");
  if (!button) return;
  const camera=saveCamera();
  if (button.getAttribute("data-attr")==="scene.dragmode") modeCamera=camera;
},true);
$("curve-plot").addEventListener("click",()=>{if(modeCamera){restoreCamera(modeCamera);modeCamera=null;}});
$("camera-home").addEventListener("click",resetCamera);
async function highlightSelection() {
  if (!current || !$("base-plot").data) return;
  const selected=Number($("edge").value), plot=$("curve-plot"), camera=saveCamera();
  const paths=(plot.data||[]).map((trace,i)=>trace.meta?.kind==="path"?i:null).filter(i=>i!==null);
  await Promise.all([Plotly.restyle($("base-plot"), {
    "line.color":current.edges.map((_,i)=>i===selected?"#bc542e":"#44736a"),
    "line.width":current.edges.map((_,i)=>i===selected?4:2)
  }, current.edges.map((_,i)=>i)),
  paths.length?Plotly.restyle(plot,{
    "line.width":paths.map(i=>plot.data[i].meta.edge===selected?8:3),
    opacity:paths.map(i=>plot.data[i].meta.edge===selected?1:.4)
  },paths):Promise.resolve()]);
  restoreCamera(camera);
}
async function renderTree() {
  if (!current) return;
  const selected=Number($("edge").value);
  const traces=current.edges.map((e,i)=>({
    type:"scatter",mode:"lines+text",x:e.base.map(p=>Number(p[0])),y:e.base.map(p=>Number(p[1])),
    text:e.base.map((_,j)=>j===Math.floor(e.base.length/2)?String(e.index):""),textposition:"top center",
    line:{color:i===selected?"#bc542e":"#44736a",width:i===selected?4:2,dash:e.infinite?"dash":"solid"},
    customdata:e.base.map(()=>i),name:"Edge "+e.index+(e.infinite?" → ∞":""),showlegend:false,
    hovertemplate:"Edge "+e.index+(e.infinite?" → ∞":"")+"<br>Re(x) = %{x:.10g}<br>Im(x) = %{y:.10g}<extra></extra>"
  }));
  traces.push({type:"scatter",mode:"markers",x:current.roots.map(p=>Number(p[0])),y:current.roots.map(p=>Number(p[1])),
    marker:{size:8,color:"#253d3c"},name:"Branch values",hovertemplate:"Branch value<br>Re(x) = %{x:.10g}<br>Im(x) = %{y:.10g}<extra></extra>"});
  traces.push({type:"scatter",mode:"markers",x:current.weights.map(w=>Number(w.center)),y:current.weights.map(()=>0),
    marker:{size:13,symbol:"diamond-open",color:"#bc542e",line:{width:2}},name:"Logarithmic centers",
    hovertemplate:"Logarithmic center<br>x=%{x:.10g}<extra></extra>"});
  const rays=current.edges.filter(e=>e.infinite).map(e=>e.base[e.base.length-1]);
  traces.push({type:"scatter",mode:"text",x:rays.map(p=>Number(p[0])),y:rays.map(p=>Number(p[1])),
    text:rays.map(()=>"∞"),textposition:"top right",showlegend:false,hoverinfo:"skip"});
  const plot=$("base-plot");
  await Plotly.react(plot,traces,{margin:{l:55,r:25,t:15,b:75},paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)",
    uirevision:label||current.input.F,dragmode:"pan",hovermode:"closest",font:{color:"#183e3b",size:11},
    xaxis:{title:"Re(x)",zeroline:false},yaxis:{title:"Im(x)",scaleanchor:"x",scaleratio:1,zeroline:false},
    legend:{orientation:"h",y:-.22}},{responsive:true,displaylogo:false});
  plot.removeAllListeners("plotly_click");
  plot.on("plotly_click",event=>{const edge=event.points[0].customdata;if(Number.isInteger(edge)){$("edge").value=edge;renderIntegrals();highlightSelection();}});
}
async function renderPlot() {
  if (!current) return;
  const traces=[], visible=Number($("visible-sheet").value), showMesh=$("show-mesh").checked;
  const selected=Number($("edge").value), camera=saveCamera();
  const dragmode=$("curve-plot")._fullLayout?.scene?.dragmode||"turntable";
  if (mesh && showMesh) mesh.sheets.forEach((s,i)=>{
    if (visible && visible!==i+1) return;
    const v=s.vertices;
    traces.push({type:"mesh3d",x:v.map(p=>p[0]),y:v.map(p=>p[1]),z:v.map(p=>p[2]),i:mesh.triangles.map(t=>t[0]),j:mesh.triangles.map(t=>t[1]),k:mesh.triangles.map(t=>t[2]),color:COLORS[i],opacity:.22,name:"Sheet "+(i+1),hoverinfo:"skip",showscale:false});
  });
  current.edges.forEach((e,i)=>e.lifts.forEach((l,s)=>{
    if (visible && s+1!==visible && l.minus_sheet!==visible) return;
    traces.push({type:"scatter3d",mode:"lines",x:l.xyz.map(p=>p[0]),y:l.xyz.map(p=>p[1]),z:l.xyz.map(p=>p[2]),line:{color:COLORS[s],width:i===selected?8:3},opacity:i===selected?1:.4,meta:{kind:"path",edge:i,sheet:s},name:"Edge "+(i+1)+", sheet "+(s+1),showlegend:false,customdata:l.xyz.map(()=>[i,s]),hovertemplate:"Edge "+(i+1)+" · sheet "+(s+1)+"<extra></extra>"});
  }));
  traces.push({type:"scatter3d",mode:"markers",x:current.roots.map(p=>Number(p[0])),y:current.roots.map(p=>Number(p[1])),z:current.roots.map(()=>0),marker:{size:2,color:"#253d3c"},showlegend:false,name:"Branch values",hovertemplate:"Branch value<br>x=%{x} + %{y} i<extra></extra>"});
  const plot=$("curve-plot");
  await Plotly.react(plot,traces,{margin:{l:0,r:0,t:0,b:0},paper_bgcolor:"rgba(0,0,0,0)",uirevision:"keep-camera",scene:{xaxis:{title:"Re(x)"},yaxis:{title:"Im(x)"},zaxis:{title:"Projected y"},aspectmode:"cube",camera:structuredClone(camera),dragmode,uirevision:"keep-camera"},legend:{orientation:"h"}},{responsive:true,displaylogo:false,scrollZoom:true,doubleClick:false,modeBarButtonsToRemove:["resetCameraDefault3d","resetCameraLastSave3d"]});
  restoreCamera(camera);
  plot.removeAllListeners("plotly_click");
  plot.on("plotly_click",event=>{ const c=event.points[0].customdata; if(c){$("edge").value=c[0];$("sheet").value=c[1];renderIntegrals();highlightSelection();} });
}
function display(data,name) {
  current=validate(data); label=name;
  $("weight-select").replaceChildren(); current.weights.forEach((w,k)=>option($("weight-select"),k,w.label)); $("weight-select").disabled=false;
  $("edge").replaceChildren(); current.edges.forEach((e,k)=>option($("edge"),k,"Edge "+e.index+(e.infinite?" → ∞":"")));
  $("result-source").textContent=name||"Imported curve"; $("result-model").textContent="y² = "+current.input.F;
  $("metadata").textContent="Genus "+current.genus+" · 2 sheets · "+current.precision+"-digit path integrals";
  $("exact-inputs").textContent="F(x) = "+current.input.F+"\nLogarithmic centers: "+current.input.centers.join(", ");
  $("point-info").textContent=current.collection?"P = "+pointTuple(current.collection.points[0])+"; R = "+pointTuple(current.collection.points[1])+"; common pole: "+pointTuple(current.collection.pole):"";
  $("geometry").textContent=current.conventions.geometry;
  $("warnings").textContent=(current.warnings||[]).join(" ");
  $("results").hidden=false; $("export").disabled=false;
  renderIntegrals(); return Promise.all([renderTree(),renderPlot()]);
}
async function chooseCurve(index) {
  const token=++version, c=catalog.curves[index]; status("Loading "+c.label+"…");
  try {
    if (!/^data\/[\w.-]+\.json$/.test(c.file)) throw Error("Invalid curve file");
    const response=await fetch(c.file+"?v="+DATA_VERSION); if(!response.ok) throw Error("Curve data could not be loaded");
    const data=await response.json(); let nextMesh=null;
    if (c.mesh && /^meshes\/[\w.-]+\.json$/.test(c.mesh)) { const r=await fetch(c.mesh+"?v="+DATA_VERSION); if(r.ok) nextMesh=validateMesh(await r.json()); }
    if(token!==version)return; mesh=nextMesh; await display(data,c.label); status("Loaded "+c.label+" · shared tree and sheet views for all weights.");
    $("mesh-caption").textContent=mesh?"Select a colored path. The two translucent sheets use the same validated cut. Drag to rotate; projection crossings can occur.":"Select a colored path. The displayed ray is truncated; integrals include infinity.";
  } catch(e) {if(token===version)status(e.message,true);}
}
function renderCatalog() {
  $("curve-select").replaceChildren(); $("catalog-rows").replaceChildren();
  catalog.curves.forEach((c,i)=>{
    option($("curve-select"),i,c.label);
    const row=document.createElement("tr"), link=document.createElement("button"); link.type="button"; link.className="secondary"; link.textContent=c.label; link.addEventListener("click",()=>{$("curve-select").value=i;chooseCurve(i);}); cell(row,"").append(link);
    c.points.forEach(p=>cell(row,pointTuple(p))); c.orders.forEach(v=>cell(row,v));
    ["det","c_plus","c_minus","L_second_derivative"].forEach(k=>cell(row,c.numerical[k]==null?"—":decimal(c.numerical[k])));
    cell(row,c.comparison_verified?c.q:(c.quotient?decimal(c.quotient):"—"));cell(row,c.relative_agreement?Number(c.relative_agreement).toExponential(2):"—");$("catalog-rows").append(row);
  }); $("curve-select").disabled=false;
}
$("curve-select").addEventListener("change",()=>chooseCurve(Number($("curve-select").value)));
["weight-select","edge","sheet","bank","orientation"].forEach(id=>$(id).addEventListener("change",()=>{renderIntegrals();if(id==="edge")highlightSelection();}));
["visible-sheet","show-mesh"].forEach(id=>$(id).addEventListener("change",renderPlot));
$("export").addEventListener("click",()=>{const url=URL.createObjectURL(new Blob([JSON.stringify(current,null,2)],{type:"application/json"}));const a=document.createElement("a");a.href=url;a.download=(label||"curve")+".json";a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
$("import").addEventListener("change",async event=>{const file=event.target.files[0];if(!file)return;try{if(file.size>20e6)throw Error("Result exceeds 20 MB");const data=validate(JSON.parse(await file.text()));++version;mesh=null;await display(data,file.name);status("Imported "+file.name);}catch(e){status(e.message,true);}event.target.value="";});
fetch("catalog.json?v="+DATA_VERSION).then(r=>{if(!r.ok)throw Error("Catalog could not be loaded");return r.json();}).then(data=>{if(data.schema_version!==2||!Array.isArray(data.curves)||!data.curves.length)throw Error("Unsupported catalog");catalog=data;renderCatalog();return chooseCurve(0);}).catch(e=>status(e.message,true));
