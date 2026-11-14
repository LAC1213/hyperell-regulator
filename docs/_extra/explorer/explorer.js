/* Shared tree and two-sheet views per curve; numerical data are precomputed. */
"use strict";
const $ = id => document.getElementById(id);
const COLORS = ["#bd5034", "#267f91"];
const DATA_VERSION = "mst-current-1";
const HOME_CAMERA = {eye:{x:1.5,y:-1.6,z:1.2},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1},projection:{type:"perspective"}};
let catalog, current, mesh = null, version = 0, label = "";
const decimal = v => Number(v) === 0 ? "0" : Number(v).toPrecision(10);
const pointTuple = p => Array.isArray(p)?"("+p.map(String).join(", ")+")":p==="infinity"?"∞":String(p);
const polynomial = source => source.replace(/\^(\d+)/g,(_,power)=>[...power].map(d=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[d]).join("")).replace(/\*/g,"").replace(/-/g,"−");
const rational = value => String(value).split("/").map(Number).reduce((a,b)=>a/b);
function math(node,tex) {
  node.dataset.tex=tex;
  if (typeof katex==="undefined") { node.textContent=tex; return; }
  katex.render(tex,node,{throwOnError:false,trust:false,maxExpand:1000,output:"htmlAndMathml"});
}
function numberTex(value) {
  const [mantissa,exponent]=decimal(value).split("e");
  return exponent===undefined?mantissa:mantissa+"\\times10^{"+Number(exponent)+"}";
}
function mathCell(row,tex) { const node=cell(row,""); math(node,tex); return node; }
function pointTex(point) {
  return Array.isArray(point)?"\\left("+point.map(CurveGeometry.rationalTex).join(",\\,")+"\\right)":"\\infty";
}
function complexCell(row,z,sign) {
  const node=cell(row,""); node.className="complex-value";
  const real=document.createElement("span"), imaginary=document.createElement("span");
  math(real,numberTex(sign*Number(z[0])));
  math(imaginary,(sign*Number(z[1])<0?"-":"+")+numberTex(Math.abs(Number(z[1])))+"\\,i");
  node.append(real,imaginary);
}
function renderPoints() {
  $("point-info").replaceChildren();
  if (!current.collection) {
    $("point-info").textContent="Logarithmic centers: "+current.input.centers.join(", ");
    return;
  }
  current.collection.points.forEach((point,i)=>{
    const node=document.createElement("span"); node.className="point-coordinate point-"+i;
    math(node,(i===0?"P":"Q")+"="+pointTex(point)); $("point-info").append(node);
  });
}
function centerNames(center) {
  const names=(current.collection?.points||[]).flatMap((p,i)=>rational(p[0])===rational(center)?[i===0?"P":"Q"]:[]);
  if (Array.isArray(current.collection?.pole) && rational(current.collection.pole[0])===rational(center)) names.push("O");
  return names.map(name=>"x("+name+")").join(", ");
}
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
    if (!(number(w.center) || typeof w.center==="string" && /^[+-]?\d+\/[1-9]\d*$/.test(w.center) && w.center.length<400 && Number.isFinite(rational(w.center))) || typeof w.label!=="string" || !Array.isArray(w.pairings) || w.pairings.length!==data.genus*(data.genus+1)/2) throw Error("Invalid logarithmic weight");
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
function renderMatrix(curve) {
  const matrix=curve?.normalized_matrix;
  const available=Array.isArray(matrix) && matrix.length===3 && matrix.every(row=>Array.isArray(row)&&row.length===3);
  $("matrix-values").hidden=!available;
  $("matrix-unavailable").hidden=available;
  $("determinant").replaceChildren();
  if (available) math($("determinant"),"D="+numberTex(curve.numerical.det));
  $("regulator-matrix").replaceChildren();
  if (!available) return;
  matrix.forEach((values,i)=>{
    const row=document.createElement("tr"), heading=document.createElement("th");
    heading.scope="row"; math(heading,["1","2\\operatorname{Re}(x)","|x|^2"][i]); row.append(heading);
    values.forEach(value=>mathCell(row,numberTex(value))); $("regulator-matrix").append(row);
  });
}
function renderIntegrals() {
  if (!current) return;
  const k=Number($("weight-select").value), e=Number($("edge").value), s=Number($("sheet").value);
  const edge=current.edges[e], lift=edge.lifts[s], sign=$("orientation").value==="reverse"?-1:1;
  $("integrals").replaceChildren(); $("pairings").replaceChildren();
  lift.ordinary.forEach((v,a)=>{const row=document.createElement("tr"); mathCell(row,"\\omega_{"+(a+1)+"}"); complexCell(row,v,sign); complexCell(row,lift.logs[k][$("bank").value][a],sign); $("integrals").append(row);});
  current.weights[k].pairings.forEach(p=>{const row=document.createElement("tr"); mathCell(row,"("+p.forms.join(",")+")"); [p.pairing,p.regulator,p.period].forEach(v=>mathCell(row,numberTex(v))); $("pairings").append(row);});
  $("path-label").textContent="Edge "+(e+1)+", sheet "+(s+1);
  $("gluing").textContent="This path borders sheet "+(s+1)+" on its positive side and sheet "+lift.minus_sheet+" on its negative side.";
  $("jump").replaceChildren();
  const jump=document.createElement("span");
  math(jump,"\\Delta\\log(x-a)="+(Number(edge.jumps[k][1])===0?numberTex(edge.jumps[k][0]):numberTex(edge.jumps[k][1])+"\\,i"));
  $("jump").append(jump,". Both sides use the selected orientation.");
  $("weight-info").replaceChildren();
  const weight=document.createElement("span");
  math(weight,"a="+CurveGeometry.rationalTex(current.weights[k].center)+(centerNames(current.weights[k].center)?"="+centerNames(current.weights[k].center):""));
  $("weight-info").append(weight,". Pairings over the whole curve.");
}
function inspectPath(edge,sheet) {
  $("edge").value=edge;
  if (sheet!==undefined) $("sheet").value=sheet;
  renderIntegrals(); highlightSelection();
}
async function renderRealPlot() {
  const plot=$("real-plot");
  try {
    const geometry=CurveGeometry.realLocus(current);
    const traces=geometry.branches.map(branch=>({
      type:"scatter",mode:"lines",x:branch.x,y:branch.y,connectgaps:false,
      line:{color:"#315d58",width:2.3},showlegend:false,meta:{kind:"real-branch"},
      hovertemplate:"x = %{x:.6g}<br>y = %{y:.6g}<extra>C(ℝ)</extra>"
    }));
    geometry.points.forEach((point,i)=>traces.push({
      type:"scatter",mode:"markers+text",x:[point.x],y:[point.y],text:[point.name],
      textposition:point.y<0?"bottom right":"top right",textfont:{family:"Source Serif 4, Georgia, serif",size:16},
      marker:{size:9,color:COLORS[i%COLORS.length],line:{color:"#fffefa",width:1.5}},
      showlegend:false,meta:{kind:"marked-point",point:point.name},
      hovertemplate:point.name+" = (%{x:.8g}, %{y:.8g})<extra></extra>"
    }));
    const axis={showgrid:false,zeroline:true,zerolinecolor:"#ccd4cc",zerolinewidth:1,
      tickfont:{family:"Inter, sans-serif",size:10,color:"#61716b"},ticks:"outside",tickcolor:"#ccd4cc",automargin:true};
    await Plotly.react(plot,traces,{
      margin:{l:38,r:20,t:20,b:35},paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)",
      font:{family:"Source Serif 4, Georgia, serif",color:"#183e3b",size:13},uirevision:current.input.F,
      xaxis:{...axis,title:{text:"x",font:{family:"KaTeX_Math",size:15}},range:geometry.xRange},
      yaxis:{...axis,title:{text:"y",font:{family:"KaTeX_Math",size:15}},range:geometry.yRange},
      dragmode:"pan",hovermode:"closest"
    },{responsive:true,displaylogo:false,scrollZoom:false,modeBarButtonsToRemove:["select2d","lasso2d","toImage"]});
  } catch (error) {
    if (plot.data) Plotly.purge(plot);
    plot.textContent="A real-curve plot is unavailable for this imported model.";
  }
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
    marker:{size:8,color:"#253d3c"},name:"Critical values of x",showlegend:false,hovertemplate:"Critical value of x<br>Re(x) = %{x:.10g}<br>Im(x) = %{y:.10g}<extra></extra>"});
  traces.push({type:"scatter",mode:"markers+text",x:current.weights.map(w=>rational(w.center)),y:current.weights.map(()=>0),
    text:current.weights.map(w=>centerNames(w.center)||"a = "+w.center),textposition:"bottom center",textfont:{size:11,color:"#bc542e"},
    marker:{size:13,symbol:"diamond-open",color:"#bc542e",line:{width:2}},name:"Marked x-coordinates",showlegend:false,
    hovertemplate:current.weights.map(w=>(centerNames(w.center)||"Logarithmic center")+"<br>x=%{x:.10g}<extra></extra>")});
  const rays=current.edges.filter(e=>e.infinite).map(e=>e.base[e.base.length-1]);
  traces.push({type:"scatter",mode:"text",x:rays.map(p=>Number(p[0])),y:rays.map(p=>Number(p[1])),
    text:rays.map(()=>"∞"),textposition:"top right",showlegend:false,hoverinfo:"skip"});
  const plot=$("base-plot");
  await Plotly.react(plot,traces,{margin:{l:45,r:20,t:15,b:45},paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)",
    uirevision:label||current.input.F,dragmode:"pan",hovermode:"closest",font:{family:"Source Serif 4, Georgia, serif",color:"#183e3b",size:13},
    xaxis:{title:"Re(x)",zeroline:false},yaxis:{title:"Im(x)",scaleanchor:"x",scaleratio:1,zeroline:false},
    legend:{orientation:"h",y:-.22}},{responsive:true,displaylogo:false});
  plot.removeAllListeners("plotly_click");
  plot.on("plotly_click",event=>{const edge=event.points[0].customdata;if(Number.isInteger(edge))inspectPath(edge);});
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
  traces.push({type:"scatter3d",mode:"markers",x:current.roots.map(p=>Number(p[0])),y:current.roots.map(p=>Number(p[1])),z:current.roots.map(()=>0),marker:{size:2,color:"#253d3c"},showlegend:false,name:"Ramification points",hovertemplate:"Ramification point on C<br>x=%{x} + %{y} i<extra></extra>"});
  const plot=$("curve-plot");
  await Plotly.react(plot,traces,{margin:{l:0,r:0,t:0,b:0},paper_bgcolor:"rgba(0,0,0,0)",font:{family:"Source Serif 4, Georgia, serif",color:"#183e3b",size:12},uirevision:"keep-camera",scene:{xaxis:{title:"Re(x)"},yaxis:{title:"Im(x)"},zaxis:{title:"Projected y"},aspectmode:"cube",camera:structuredClone(camera),dragmode,uirevision:"keep-camera"},legend:{orientation:"h"}},{responsive:true,displaylogo:false,scrollZoom:true,doubleClick:false,modeBarButtonsToRemove:["resetCameraDefault3d","resetCameraLastSave3d"]});
  restoreCamera(camera);
  plot.removeAllListeners("plotly_click");
  plot.on("plotly_click",event=>{ const c=event.points[0].customdata; if(c)inspectPath(c[0],c[1]); });
}
async function display(data,name,curve=null) {
  current=validate(data); label=name;
  $("weight-select").replaceChildren(); current.weights.forEach((w,k)=>option($("weight-select"),k,(centerNames(w.center)?centerNames(w.center)+" · a = "+w.center:w.label))); $("weight-select").disabled=false;
  $("edge").replaceChildren(); current.edges.forEach((e,k)=>option($("edge"),k,"Edge "+e.index+(e.infinite?" → ∞":"")));
  $("result-source").textContent=name||"Imported curve";
  try { math($("result-model"),CurveGeometry.equationTex(current.input.F)); }
  catch (error) { $("result-model").textContent="C : y² = "+polynomial(current.input.F); }
  $("metadata").textContent="Genus "+current.genus+" · 2 sheets · "+current.precision+"-digit path integrals";
  $("exact-inputs").textContent="F(x) = "+current.input.F+"\nLogarithmic centers: "+current.input.centers.join(", ")+(current.collection?"\nP = "+pointTuple(current.collection.points[0])+"\nQ = "+pointTuple(current.collection.points[1])+"\nCommon pole O = "+pointTuple(current.collection.pole):"");
  renderPoints();
  $("show-mesh").disabled=!mesh;
  $("mesh-caption").textContent="Drag to rotate. Click a colored path to inspect its integrals. "+(mesh?"This is a projection; crossings need not be intersections on C.":"Surface mesh unavailable; the lifted paths are shown.");
  $("geometry").textContent=current.conventions.geometry;
  $("warnings").textContent=(current.warnings||[]).join(" ");
  $("results").hidden=false; $("export").disabled=false;
  renderMatrix(curve);
  renderIntegrals();
  if (document.fonts) await document.fonts.ready;
  return Promise.all([renderRealPlot(),renderTree(),renderPlot()]);
}
async function chooseCurve(index) {
  const token=++version, c=catalog.curves[index]; status("Loading "+c.label+"…");
  try {
    if (!/^data\/[\w.-]+\.json$/.test(c.file)) throw Error("Invalid curve file");
    const response=await fetch(c.file+"?v="+DATA_VERSION); if(!response.ok) throw Error("Curve data could not be loaded");
    const data=await response.json(); let nextMesh=null;
    if (c.mesh && /^meshes\/[\w.-]+\.json$/.test(c.mesh)) { const r=await fetch(c.mesh+"?v="+DATA_VERSION); if(r.ok) nextMesh=validateMesh(await r.json()); }
    if(token!==version)return; mesh=nextMesh; await display(data,c.label,c);
    if(token!==version)return;
    for (const row of $("catalog-rows").children) row.classList.toggle("selected",row.dataset.label===c.label);
    history.replaceState(null,"","#"+encodeURIComponent(c.label));
    status((data.warnings||[]).join(" "),Boolean(data.warnings?.length));
  } catch(e) {if(token===version)status(e.message,true);}
}
function renderCatalog() {
  $("curve-select").replaceChildren(); $("catalog-rows").replaceChildren();
  catalog.curves.forEach((c,i)=>{
    option($("curve-select"),i,c.label);
    const row=document.createElement("tr"), link=document.createElement("button"); row.dataset.label=c.label; link.type="button"; link.className="secondary"; link.textContent=c.label; link.addEventListener("click",()=>{$("curve-select").value=i;chooseCurve(i);}); cell(row,"").append(link);
    c.points.forEach(p=>mathCell(row,pointTex(p))); c.orders.forEach(v=>mathCell(row,String(v)));
    ["det","c_plus","c_minus","L_second_derivative"].forEach(k=>c.numerical[k]==null?cell(row,"—"):mathCell(row,numberTex(c.numerical[k])));
    if(c.comparison_verified)mathCell(row,CurveGeometry.rationalTex(c.q));else if(c.quotient)mathCell(row,numberTex(c.quotient));else cell(row,"—");
    if(c.relative_agreement)mathCell(row,Number(c.relative_agreement).toExponential(2).replace(/e([+-]?\d+)/,"\\times10^{$1}"));else cell(row,"—");$("catalog-rows").append(row);
  }); $("curve-select").disabled=false;
}
$("curve-select").addEventListener("change",()=>chooseCurve(Number($("curve-select").value)));
["weight-select","edge","sheet","bank","orientation"].forEach(id=>$(id).addEventListener("change",()=>{renderIntegrals();if(id==="edge")highlightSelection();}));
["visible-sheet","show-mesh"].forEach(id=>$(id).addEventListener("change",renderPlot));
$("export").addEventListener("click",()=>{const url=URL.createObjectURL(new Blob([JSON.stringify(current,null,2)],{type:"application/json"}));const a=document.createElement("a");a.href=url;a.download=(label||"curve")+".json";a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
$("import").addEventListener("change",async event=>{const file=event.target.files[0];if(!file)return;try{if(file.size>20e6)throw Error("Result exceeds 20 MB");const data=validate(JSON.parse(await file.text()));++version;mesh=null;$("curve-select").value="";for(const row of $("catalog-rows").children)row.classList.remove("selected");await display(data,file.name);status("Imported "+file.name);}catch(e){status(e.message,true);}event.target.value="";});
document.querySelectorAll("[data-tex]").forEach(node=>math(node,node.dataset.tex));
fetch("catalog.json?v="+DATA_VERSION).then(r=>{if(!r.ok)throw Error("Catalog could not be loaded");return r.json();}).then(data=>{if(data.schema_version!==2||!Array.isArray(data.curves)||!data.curves.length)throw Error("Unsupported catalog");catalog=data;renderCatalog();const saved=catalog.curves.findIndex(c=>c.label===decodeURIComponent(location.hash.slice(1)).split("/")[0]);$("curve-select").value=Math.max(0,saved);return chooseCurve(Math.max(0,saved));}).catch(e=>status(e.message,true));
