/* Shared tree and two-sheet views per curve; numerical data are precomputed. */
"use strict";
const $ = id => document.getElementById(id);
const COLORS = ["#bd5034", "#267f91"];
const POINT_COLORS = ["#bd5034", "#267f91", "#7651a0", "#488551", "#b07b25", "#b04d84"];
const PI = "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647";
const DATA_VERSION = "transformed-picture-6";
const HOME_CAMERA = {eye:{x:1.5,y:-1.6,z:1.2},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1},projection:{type:"perspective"}};
let catalog, current, currentCurve = null, mesh = null, version = 0, label = "";
// Round the saved decimal strings directly: Number would discard digits after ~16.
function decimal(value,digits=10) {
  const match=/^([+-]?)(\d*)(?:\.(\d*))?(?:e([+-]?\d+))?$/i.exec(String(value).trim());
  if (!match || !(match[2] || match[3])) throw Error("Invalid decimal value");
  const raw=match[2]+(match[3]||""), first=raw.search(/[1-9]/);
  if (first<0) return "0";
  let exponent=Number(match[4]||0)+match[2].length-first-1;
  const significant=raw.slice(first);
  let rounded=significant.slice(0,digits).padEnd(digits,"0");
  if (Number(significant[digits]||0)>=5) {
    rounded=(BigInt(rounded)+1n).toString();
    if (rounded.length>digits) { exponent++; rounded=rounded.slice(0,digits); }
  }
  const sign=match[1]==="-"?"-":"";
  if (exponent < -6 || exponent >= digits) return sign+rounded[0]+(digits>1?"."+rounded.slice(1):"")+"e"+(exponent>=0?"+":"")+exponent;
  if (exponent<0) return sign+"0."+"0".repeat(-exponent-1)+rounded;
  return sign+rounded.slice(0,exponent+1)+(exponent+1<digits?"."+rounded.slice(exponent+1):"");
}
const polynomial = source => source.replace(/\^(\d+)/g,(_,power)=>[...power].map(d=>"⁰¹²³⁴⁵⁶⁷⁸⁹"[d]).join("")).replace(/\*/g,"").replace(/-/g,"−");
const rational = value => String(value).split("/").map(Number).reduce((a,b)=>a/b);
function math(node,tex) {
  node.dataset.tex=tex;
  if (typeof katex==="undefined") { node.textContent=tex; return; }
  katex.render(tex,node,{throwOnError:false,trust:false,maxExpand:1000,output:"htmlAndMathml"});
}
function numberTex(value,digits=10) {
  const [mantissa,exponent]=decimal(value,digits).split("e");
  return exponent===undefined?mantissa:mantissa+"\\times10^{"+Number(exponent)+"}";
}
function mathCell(row,tex) { const node=cell(row,""); math(node,tex); return node; }
function pointTex(point) {
  return Array.isArray(point)?"\\left("+point.map(CurveGeometry.rationalTex).join(",\\,")+"\\right)":"\\infty";
}
function lmfdbLink(label,text=label+" ↗") {
  const parts=/^(\d+)\.([a-z]+)\.(\d+)\.(\d+)$/.exec(label);
  const link=document.createElement(parts?"a":"span"); link.textContent=text;
  if (parts) {
    link.href="https://www.lmfdb.org/Genus2Curve/Q/"+parts.slice(1).join("/");
    link.target="_blank"; link.rel="noopener noreferrer";
    link.setAttribute("aria-label",label+" in the LMFDB (opens in a new tab)");
  }
  return link;
}
function complexCell(row,z,sign) {
  const node=cell(row,""); node.className="complex-value";
  const real=document.createElement("span"), imaginary=document.createElement("span");
  math(real,numberTex(sign*Number(z[0])));
  math(imaginary,(sign*Number(z[1])<0?"-":"+")+numberTex(Math.abs(Number(z[1])))+"\\,i");
  node.append(real,imaginary);
}
function renderPoints(curve) {
  $("point-info").replaceChildren(); $("packet-info").textContent=""; $("support-info").replaceChildren();
  const collection=current.collection;
  if (!collection?.divisors) {
    $("point-info").textContent="Logarithmic centers: "+current.weights.map(w=>w.center).join(", ");
    return;
  }
  const active=$("divisor-select").value;
  collection.divisors.forEach(D=>{
    const node=document.createElement("div"); node.className="divisor-coordinate";
    node.classList.toggle("active",D.id===active); math(node,D.tex); $("point-info").append(node);
  });
  $("packet-info").textContent="Packet torsion order: "+(curve?.packet_torsion_order||collection.packet_torsion_order||"not recorded")+
    " · Integral trace rank over ℚ: "+(curve?.integral_rank??collection.integral_rank??"not recorded")+
    " · Rational endomorphism algebra: "+(curve?.end_alg||collection.end_alg||"not recorded");
  $("support-info").replaceChildren();
  if (collection.lmfdb_model) {
    const model=document.createElement("div"); model.className="support-coordinate";
    math(model,collection.lmfdb_model.equation_tex); $("support-info").append(model);
  }
  (collection.support||[]).filter(P=>rational(P.coefficients[active]||0)!==0).forEach(P=>{
    const node=document.createElement("span"); node.className="support-coordinate";
    math(node,P.definition_tex+"\\quad ("+CurveGeometry.rationalTex(P.coefficients[active])+")");
    $("support-info").append(node);
  });
}

function centerPosition(center) {
  const k=current.weights.findIndex(w=>w.center===center), weight=current.weights[k];
  return weight?.center_position||current.collection?.logarithmic_centers?.[k]?.position||[rational(center),0];
}
function centerTex(center) {
  const k=current.weights.findIndex(w=>w.center===center);
  return current.collection?.logarithmic_centers?.[k]?.tex||
    (current.weights[k]?.center_exact?.polynomial?"a_{"+(k+1)+"}":CurveGeometry.rationalTex(center));
}
function centerReferences(center) {
  const index=current.weights.findIndex(w=>w.center===center);
  return (current.collection?.support||[]).filter(P=>P.center_index===index).map(P=>"u(P_{"+P.id.slice(1)+"})");
}

function centerNames(center) { return centerReferences(center).map(name=>name.replace(/_\{(\d+)\}/g,"$1")).join(", "); }
function centerColor(center) {
  const index=current.weights.findIndex(w=>w.center===center);
  return POINT_COLORS[Math.max(0,index)%POINT_COLORS.length];
}

function selectedWeight() { return Number($("weight-select").value); }
function treeEdgeColor(edge) {
  const k=selectedWeight();
  return edge.jumps[k].some(value=>Number(value)!==0)?centerColor(current.weights[k].center):"#88958f";
}
function status(message, error=false) { $("status").textContent=message; $("status").classList.toggle("error",error); }
function option(select,value,text) { const n=document.createElement("option"); n.value=value; n.textContent=text; select.append(n); }
function cell(row,text) { const n=document.createElement("td"); n.textContent=text; row.append(n); return n; }
function complex(z,sign=1) { return decimal(sign*Number(z[0]))+(sign*Number(z[1])<0?" − ":" + ")+decimal(Math.abs(Number(z[1])))+" i"; }
function validate(data) {
  const number=v => (typeof v === "string" || typeof v === "number") && String(v).length<400 && String(v).trim()!=="" && Number.isFinite(Number(v));
  const pair=z => Array.isArray(z) && z.length===2 && z.every(number);
  if (!data || data.schema_version!==2 || data.degree!==2 || !Number.isInteger(data.genus) || data.genus<2 || data.genus>5 ||
      !data.input || typeof data.input.F!=="string" || data.input.F.length>4000 || !Array.isArray(data.weights) || !data.weights.length || data.weights.length>128 ||
      !Array.isArray(data.edges) || !data.edges.length || data.edges.length>400) throw Error("Expected a direct two-sheet result (schema version 2)");
  if (!Array.isArray(data.roots) || !data.roots.every(pair)) throw Error("Invalid branch values");
  if (data.boundary_word !== undefined) {
    const word=data.boundary_word;
    if (!Array.isArray(word) || word.length!==2*data.edges.length || !word.every(b=>Array.isArray(b)&&b.length===2&&Number.isInteger(b[0])&&b[0]>=0&&b[0]<data.edges.length&&[1,-1].includes(b[1])) || new Set(word.map(b=>b.join(","))).size!==word.length) throw Error("Invalid tree boundary order");
  }
  data.weights.forEach(w => {
    if (!(pair(w.center_position) || number(w.center) || typeof w.center==="string" && /^[+-]?\d+\/[1-9]\d*$/.test(w.center) && w.center.length<400 && Number.isFinite(rational(w.center))) || typeof w.label!=="string" || !Array.isArray(w.pairings) || w.pairings.length!==data.genus*(data.genus+1)/2) throw Error("Invalid logarithmic weight");
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
  $("matrix-caption").hidden=!available;
  $("import-pairings").hidden=available;
  $("matrix-name").textContent=available?"The regulator matrix":"Logarithmic pairings";
  $("determinant").hidden=!available;
  $("determinant").replaceChildren();
  if (available) {
    const heading=document.createElement("span"), value=document.createElement("span");
    heading.className="determinant-label"; heading.textContent="Determinant:";
    math(value,"D="+numberTex(curve.numerical.det)); $("determinant").append(heading,value);
  }
  $("regulator-matrix").replaceChildren();
  if (!available) return;
  const names=curve.matrix_divisor_ids||["D1","D2"];
  math($("matrix-left"),"q_*B(D_{"+names[0].slice(1)+"})");
  math($("matrix-right"),"q_*B(D_{"+names[1].slice(1)+"})");
  matrix.forEach((values,i)=>{
    const row=document.createElement("tr"), heading=document.createElement("th");
    heading.scope="row"; math(heading,["\\omega_1\\wedge\\overline{\\omega_1}","\\omega_1\\wedge\\overline{\\omega_2}+\\omega_2\\wedge\\overline{\\omega_1}","\\omega_2\\wedge\\overline{\\omega_2}"][i]); row.append(heading);
    values.forEach(value=>mathCell(row,numberTex(value))); $("regulator-matrix").append(row);
  });
}
function linearRelationTex(curve) {
  const terms=curve?.linear_comparison?.terms||[];
  if(terms.length<2)return null;
  return terms.map(term=>CurveGeometry.rationalTex(term.coefficient)+"\\,D_{"+term.basis_indices.map(i=>i+1).join(",")+"}").join("+").replace(/\+\-/g,"-");
}
function renderBeilinson(curve) {
  const values=curve?.numerical, relation=linearRelationTex(curve);
  const available=values && ["det","c_plus","c_minus","L_second_derivative"].every(key=>values[key]!=null&&Number.isFinite(Number(values[key]))) && (relation||curve.comparison_verified||curve.quotient!=null);
  $("beilinson-panel").hidden=!available;
  $("beilinson-values").replaceChildren();
  if (!available) return;
  conjectureFormula="\\pi\\cdot c^+(C)\\cdot c^-(C)\\cdot L''(M_{\\mathrm{std}}(C),1)\\approx "+(relation||"q\\cdot D");
  const maximum=Math.min(100,curve.regulator_digits||100,curve.L_digits||100);
  $("beilinson-digits").max=maximum;
  if (Number($("beilinson-digits").value)>maximum) $("beilinson-digits").value=maximum;
  const digits=Number($("beilinson-digits").value);
  $("beilinson-digits-value").value=String(digits);
  $("beilinson-values").classList.toggle("precise",digits>30);
  const items=[
    ["\\pi",PI],
    ["c^+(C)",values.c_plus],
    ["c^-(C)",values.c_minus],
    ["L''",values.L_second_derivative],
    ["q",curve.comparison_verified?curve.q:curve.quotient],
    ["D",values.det]
  ];
  items.filter(([label])=>!relation||label!=="q").forEach(([label,value])=>{
    const item=document.createElement("div"), term=document.createElement("dt"), definition=document.createElement("dd");
    math(term,label);
    if (label==="q" && curve.comparison_verified) math(definition,"="+CurveGeometry.rationalTex(value));
    else {
      definition.className="conjecture-number";
      const approximation=document.createElement("span"), number=document.createElement("span");
      math(approximation,"\\approx");
      const [mantissa,exponent]=decimal(value,digits).split("e");
      number.className="significant-digits"; number.textContent=mantissa;
      definition.append(approximation,number);
      if (exponent!==undefined) {
        const power=document.createElement("span"); math(power,"\\times10^{"+Number(exponent)+"}"); definition.append(power);
      }
    }
    item.append(term,definition); $("beilinson-values").append(item);
  });
}
function displayEquationTex(multiline=false) {
  return CurveGeometry.equationTex(current.input.F,multiline).replace(/x/g,"u").replace(/y/g,"v");
}
function fitCurveEquation() {
  if (!current) return;
  const node=$("result-model"), available=node.parentElement.clientWidth;
  if (!available) return;
  node.style.fontSize="";
  try {
    math(node,displayEquationTex());
    const width=node.querySelector(".katex").getBoundingClientRect().width;
    const size=parseFloat(getComputedStyle(node).fontSize);
    if (width<=available) return;
    const fitted=size*(available-4)/width;
    if (fitted>=17) node.style.fontSize=fitted+"px";
    else math(node,displayEquationTex(true));
  } catch (error) { node.textContent="C : v² = "+polynomial(current.input.F).replace(/x/g,"u"); }
}
let conjectureFormula=$("beilinson-equation").dataset.tex;
function fitConjectureEquation() {
  if ($("beilinson-panel").hidden) return;
  const node=$("beilinson-equation"), available=node.clientWidth;
  if (!available) return;
  node.style.fontSize="";
  math(node,conjectureFormula);
  const width=node.querySelector(".katex").getBoundingClientRect().width;
  const size=parseFloat(getComputedStyle(node).fontSize);
  if (width<=available) return;
  const fitted=size*(available-4)/width;
  if (fitted>=14) node.style.fontSize=fitted+"px";
  else math(node,"\\begin{gathered}"+conjectureFormula.replace("\\cdot L''","\\\\ {}\\cdot L''")+"\\end{gathered}");
}
let mathResizeFrame=null;
window.addEventListener("resize",()=>{
  if (mathResizeFrame!==null) cancelAnimationFrame(mathResizeFrame);
  mathResizeFrame=requestAnimationFrame(()=>{mathResizeFrame=null;fitCurveEquation();fitConjectureEquation();});
});
function renderIntegrals() {
  if (!current) return;
  const k=Number($("weight-select").value), e=Number($("edge").value), s=Number($("sheet").value);
  const edge=current.edges[e], lift=edge.lifts[s], sign=$("orientation").value==="reverse"?-1:1;
  $("integrals").replaceChildren(); $("pairings").replaceChildren();
  lift.ordinary.forEach((v,a)=>{
    const row=document.createElement("tr"), factor=a===0?"":a===1?"u\\,":"u^{"+a+"}\\,";
    mathCell(row,"\\omega_{"+(a+1)+"}=\\dfrac{"+factor+"\\mathrm{d}u}{v}");
    complexCell(row,v,sign); complexCell(row,lift.logs[k][$("bank").value][a],sign); $("integrals").append(row);
  });
  current.weights[k].pairings.forEach(p=>{const row=document.createElement("tr"); mathCell(row,"("+p.forms.join(",")+")"); [p.pairing,p.regulator,p.period].forEach(v=>mathCell(row,numberTex(v))); $("pairings").append(row);});
  $("path-label").textContent="Edge "+(e+1)+", sheet "+(s+1);
  $("gluing").textContent="This path borders sheet "+(s+1)+" on its positive side and sheet "+lift.minus_sheet+" on its negative side.";
  $("jump").replaceChildren();
  const jump=document.createElement("span");
  math(jump,"\\Delta\\log(u-a)="+(Number(edge.jumps[k][1])===0?numberTex(edge.jumps[k][0]):numberTex(edge.jumps[k][1])+"\\,i"));
  $("jump").append(jump,". Both sides use the selected orientation.");
  $("weight-info").replaceChildren();
  const weight=document.createElement("span");
  math(weight,"a="+[...centerReferences(current.weights[k].center),centerTex(current.weights[k].center)].join("="));
  $("weight-info").append(weight,". Pairings over the whole curve.");
}
function inspectPath(edge,sheet) {
  $("edge").value=edge;
  if (sheet!==undefined) $("sheet").value=sheet;
  renderIntegrals(); highlightSelection();
}
async function renderRealPlot() {
  const plot=$("real-plot"), active=$("divisor-select").value;
  try {
    const geometry=CurveGeometry.realLocus(current,active);
    const traces=geometry.branches.map(branch=>({
      type:"scatter",mode:"lines",x:branch.x,y:branch.y,connectgaps:false,
      line:{color:"#315d58",width:2.3},showlegend:false,meta:{kind:"real-branch"},
      hovertemplate:"u = %{x:.6g}<br>v = %{y:.6g}<extra>C(ℝ)</extra>"
    }));
    geometry.points.forEach(point=>{
      const coefficient=rational(point.coefficient||1), sign=coefficient<0?"−":"+";
      traces.push({type:"scatter",mode:"markers+text",x:[point.x],y:[point.y],
        text:[sign+Math.abs(coefficient)+" "+point.name],textposition:point.y<0?"bottom right":"top right",
        textfont:{family:"Source Serif 4, Georgia, serif",size:14},
        marker:{size:10+Math.min(6,Math.log2(1+Math.abs(coefficient))),symbol:coefficient<0?"circle-open":"circle",
          color:coefficient<0?COLORS[1]:COLORS[0],line:{width:2}},
        showlegend:false,meta:{kind:"divisor-support",point:point.name,coefficient:point.coefficient,divisor:active},
        hovertemplate:point.name+" · coefficient "+point.coefficient+"<br>(%{x:.8g}, %{y:.8g})<extra></extra>"});
    });
    const annotations=geometry.infinity.map((P,i)=>({xref:"paper",yref:"paper",x:1,y:1-.1*i,xanchor:"right",
      text:P.coefficients[active]+" · ∞ (v/u³ = "+decimal(P.y[0],6)+"; "+P.id+")",showarrow:false,font:{size:14,color:COLORS[1]}}));
    const axis={showgrid:false,zeroline:true,zerolinecolor:"#ccd4cc",zerolinewidth:1,
      tickfont:{family:"Inter, sans-serif",size:11,color:"#61716b"},ticks:"outside",tickcolor:"#ccd4cc",automargin:true};
    await Plotly.react(plot,traces,{annotations,
      margin:{l:38,r:20,t:20,b:35},paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)",
      font:{family:"Source Serif 4, Georgia, serif",color:"#183e3b",size:14},uirevision:current.input.F+"/"+active,
      xaxis:{...axis,title:{text:"u"},range:geometry.xRange},yaxis:{...axis,title:{text:"v"},range:geometry.yRange},
      dragmode:"pan",hovermode:"closest"
    },{responsive:true,displaylogo:false,scrollZoom:false,modeBarButtonsToRemove:["select2d","lasso2d","toImage"]});
    $("real-support-caption").textContent="The picture uses transformed coordinates (u, v); divisor formulas use LMFDB coordinates. Filled markers have positive coefficients; open markers have negative coefficients."+
      (geometry.infinity.length?" Points at infinity are listed in the plot annotation.":"")+
      (geometry.nonreal.length?" "+geometry.nonreal.length+" nonreal support points are listed below; their images under the cut coordinate u appear in the tree.":"");
  } catch (error) {
    if (plot.data) Plotly.purge(plot);
    plot.textContent="A real-curve plot is unavailable: "+error.message;
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
    "line.color":current.edges.map(treeEdgeColor),
    "line.width":current.edges.map((_,i)=>i===selected?4:2)
  }, current.edges.map((_,i)=>i)),
  paths.length?Plotly.restyle(plot,{
    "line.width":paths.map(i=>plot.data[i].meta.edge===selected?8:3),
    opacity:paths.map(i=>plot.data[i].meta.edge===selected?1:.4)
  },paths):Promise.resolve()]);
  restoreCamera(camera);
}
async function updateTreeSelection() {
  const plot=$("base-plot");
  if (!current || !plot.data) return;
  const k=selectedWeight();
  // Update the branch cut and active diamond without changing the axes or view.
  await Promise.all([
    Plotly.restyle(plot,{"line.color":current.edges.map(treeEdgeColor)},current.edges.map((_,i)=>i)),
    Plotly.restyle(plot,{
      "marker.size":current.weights.map((_,i)=>i===k?16:13),
      "marker.line.width":current.weights.map((_,i)=>i===k?3:2)
    },current.weights.map((_,i)=>current.edges.length+1+i))
  ]);
}
async function renderTree() {
  if (!current) return;
  const data=current;
  const selected=Number($("edge").value), k=selectedWeight(), center=centerPosition(current.weights[k].center);
  const traces=current.edges.map((e,i)=>({
    type:"scatter",mode:"lines+text",x:e.base.map(p=>Number(p[0])),y:e.base.map(p=>Number(p[1])),
    text:e.base.map((_,j)=>j===Math.floor(e.base.length/2)?String(e.index):""),textposition:"top center",
    line:{color:treeEdgeColor(e),width:i===selected?4:2,dash:e.infinite?"dash":"solid"},
    customdata:e.base.map(()=>i),name:"Edge "+e.index+(e.infinite?" → ∞":""),showlegend:false,
    hovertemplate:"Edge "+e.index+(e.infinite?" → ∞":"")+"<br>Re(u) = %{x:.10g}<br>Im(u) = %{y:.10g}<extra></extra>"
  }));
  traces.push({type:"scatter",mode:"markers",x:current.roots.map(p=>Number(p[0])),y:current.roots.map(p=>Number(p[1])),
    marker:{size:8,color:"#253d3c"},name:"Critical values of u",showlegend:false,hovertemplate:"Critical value of u<br>Re(u) = %{x:.10g}<br>Im(u) = %{y:.10g}<extra></extra>"});
  current.weights.forEach((weight,i)=>{
    const color=centerColor(weight.center), name=centerNames(weight.center)||"a = "+weight.center;
    traces.push({type:"scatter",mode:"markers+text",x:[centerPosition(weight.center)[0]],y:[centerPosition(weight.center)[1]],
      text:[current.weights.length>4?"a"+(i+1):name],textposition:i%2?"bottom center":"top center",textfont:{size:11,color},
      marker:{size:i===k?16:13,symbol:"diamond-open",color,line:{width:i===k?3:2}},name,showlegend:false,
      hovertemplate:"a"+(i+1)+": "+name+"<br>Re(u)=%{x:.10g}<br>Im(u)=%{y:.10g}<extra></extra>"});
  });
  const rays=current.edges.filter(e=>e.infinite).map(e=>e.base[e.base.length-1]);
  traces.push({type:"scatter",mode:"text",x:rays.map(p=>Number(p[0])),y:rays.map(p=>Number(p[1])),
    text:rays.map(()=>"∞"),textposition:"top right",showlegend:false,hoverinfo:"skip"});
  const plot=$("base-plot");
  // Symmetric ranges keep a at the centre while retaining the whole finite tree.
  const points=[...current.edges.flatMap(e=>e.base),...current.roots,...current.weights.map(w=>centerPosition(w.center))];
  const aspect=Math.max(1,plot.clientWidth-65)/Math.max(1,plot.clientHeight-60);
  const halfX=1.15*points.reduce((extent,p)=>Math.max(extent,Math.abs(Number(p[0])-center[0]),Math.abs(Number(p[1])-center[1])*aspect),.1);
  const halfY=halfX/aspect;
  await Plotly.react(plot,traces,{margin:{l:45,r:20,t:15,b:45},paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)",
    uirevision:(label||current.input.F)+"/"+current.weights[k].center,dragmode:"pan",hovermode:"closest",font:{family:"Source Serif 4, Georgia, serif",color:"#183e3b",size:13},
    xaxis:{title:"Re(u)",range:[center[0]-halfX,center[0]+halfX],autorange:false,zeroline:false},yaxis:{title:"Im(u)",range:[center[1]-halfY,center[1]+halfY],autorange:false,scaleanchor:"x",scaleratio:1,zeroline:false},
    legend:{orientation:"h",y:-.22}},{responsive:true,displaylogo:false});
  if (data!==current) { renderIntegrals(); return renderTree(); }
  // Apply the selected centre even if Plotly retained a prior pan or autorange.
  await Plotly.relayout(plot,{
    "xaxis.range":[center[0]-halfX,center[0]+halfX],"xaxis.autorange":false,
    "yaxis.range":[center[1]-halfY,center[1]+halfY],"yaxis.autorange":false
  });
  if (k!==selectedWeight()) { renderIntegrals(); await updateTreeSelection(); }
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
  await Plotly.react(plot,traces,{margin:{l:0,r:0,t:0,b:0},paper_bgcolor:"rgba(0,0,0,0)",font:{family:"Source Serif 4, Georgia, serif",color:"#183e3b",size:12},uirevision:"keep-camera",scene:{xaxis:{title:"Re(u)"},yaxis:{title:"Im(u)"},zaxis:{title:"Projected v"},aspectmode:"cube",camera:structuredClone(camera),dragmode,uirevision:"keep-camera"},legend:{orientation:"h"}},{responsive:true,displaylogo:false,scrollZoom:true,doubleClick:false,modeBarButtonsToRemove:["resetCameraDefault3d","resetCameraLastSave3d"]});
  restoreCamera(camera);
  plot.removeAllListeners("plotly_click");
  plot.on("plotly_click",event=>{ const c=event.points[0].customdata; if(c)inspectPath(c[0],c[1]); });
}
async function display(data,name,curve=null) {
  current=validate(data);
  curve=curve||(data.collection?.experiment?{...data.collection.experiment,...data.collection}:null);
  currentCurve=curve; label=name;
  if (curve?.divisors) current.collection={...(current.collection||{}),
    ...Object.fromEntries(["divisors","support","logarithmic_centers","matrix_divisor_ids","packet_torsion_order","end_alg","coordinate_system","lmfdb_model","integral_rank","integral_rank_field","integral_rank_space","trace_rank_certificate"].map(k=>[k,curve[k]]))};
  $("divisor-select").replaceChildren();
  (current.collection?.divisors||[]).forEach(D=>option($("divisor-select"),D.id,D.id));
  $("divisor-select").disabled=!current.collection?.divisors?.length;
  $("weight-select").replaceChildren(); current.weights.forEach((w,k)=>option($("weight-select"),k,"a = "+[centerNames(w.center),w.center].filter(Boolean).join(" = ")));
  $("weight-select").value="0"; $("weight-select").disabled=false;
  $("edge").replaceChildren(); current.edges.forEach((e,k)=>option($("edge"),k,"Edge "+e.index+(e.infinite?" → ∞":"")));
  $("result-source").replaceChildren();
  if (curve) {
    const curveName=document.createElement("strong"), link=lmfdbLink(curve.label,"LMFDB ↗");
    curveName.textContent=curve.label; link.className="lmfdb-link";
    $("result-source").append(curveName,link);
  }
  else $("result-source").textContent=name||"Imported curve";
  try { math($("result-model"),displayEquationTex()); }
  catch (error) { $("result-model").textContent="C : v² = "+polynomial(current.input.F).replace(/x/g,"u"); }
  math($("cut-model"),"v^2="+CurveGeometry.equationTex(current.input.F).split("y^2=")[1].replace(/x/g,"u"));
  renderPoints(curve);
  $("show-mesh").disabled=!mesh;
  $("mesh-caption").textContent=mesh?"Drag to rotate. Use the sheet selector to view either sheet.":"Drag to rotate the lifted paths. Surface mesh unavailable.";
  $("results").hidden=false; $("export").disabled=false;
  renderMatrix(curve);
  renderBeilinson(curve);
  if (document.fonts) await document.fonts.ready;
  fitCurveEquation();
  fitConjectureEquation();
  renderIntegrals();
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
    option($("curve-select"),i,c.label+" · "+(c.end_alg||"unknown"));
    const row=document.createElement("tr"), button=document.createElement("button");
    row.dataset.label=c.label; button.type="button"; button.className="secondary"; button.textContent=c.label;
    button.addEventListener("click",()=>{$("curve-select").value=i;chooseCurve(i);});
    const link=lmfdbLink(c.label,"LMFDB ↗"); link.className="lmfdb-link"; cell(row,"").append(button,link);
    cell(row,c.end_alg||"—"); cell(row,c.packet_torsion_order||"—"); cell(row,c.integral_rank??"—");
    const divisors=cell(row,""); divisors.className="divisor-cell";
    (c.divisors||[]).forEach(D=>{const node=document.createElement("div");math(node,D.tex);divisors.append(node);});
    ["det","c_plus","c_minus","L_second_derivative"].forEach(k=>c.numerical[k]==null?cell(row,"—"):mathCell(row,numberTex(c.numerical[k])));
    if(linearRelationTex(c))mathCell(row,linearRelationTex(c));else if(c.comparison_verified)mathCell(row,CurveGeometry.rationalTex(c.q));else if(c.quotient)mathCell(row,numberTex(c.quotient));else cell(row,"—");
    if(c.relative_agreement)mathCell(row,numberTex(c.relative_agreement,2));else cell(row,"—");
    $("catalog-rows").append(row);
  }); $("curve-select").disabled=false;
}
$("curve-select").addEventListener("change",()=>chooseCurve(Number($("curve-select").value)));
$("divisor-select").addEventListener("change",()=>{renderPoints(currentCurve);renderRealPlot();});
["weight-select","edge","sheet","bank","orientation"].forEach(id=>$(id).addEventListener("change",()=>{renderIntegrals();if(id==="edge")highlightSelection();if(id==="weight-select")updateTreeSelection();}));
$("beilinson-digits").addEventListener("input",()=>renderBeilinson(currentCurve));
window.addEventListener("pageshow",()=>{if(current){renderIntegrals();updateTreeSelection();}});
["visible-sheet","show-mesh"].forEach(id=>$(id).addEventListener("change",renderPlot));
$("export").addEventListener("click",()=>{const url=URL.createObjectURL(new Blob([JSON.stringify(current,null,2)],{type:"application/json"}));const a=document.createElement("a");a.href=url;a.download=(label||"curve")+".json";a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
$("import").addEventListener("change",async event=>{const file=event.target.files[0];if(!file)return;try{if(file.size>20e6)throw Error("Result exceeds 20 MB");const data=validate(JSON.parse(await file.text()));++version;mesh=null;$("curve-select").value="";for(const row of $("catalog-rows").children)row.classList.remove("selected");await display(data,file.name);status("Imported "+file.name);}catch(e){status(e.message,true);}event.target.value="";});
document.querySelectorAll("[data-tex]").forEach(node=>math(node,node.dataset.tex));
fetch("catalog.json?v="+DATA_VERSION).then(r=>{if(!r.ok)throw Error("Catalog could not be loaded");return r.json();}).then(data=>{if(![2,3].includes(data.schema_version)||!Array.isArray(data.curves)||!data.curves.length)throw Error("Unsupported catalog");catalog=data;renderCatalog();const saved=catalog.curves.findIndex(c=>c.label===decodeURIComponent(location.hash.slice(1)).split("/")[0]);$("curve-select").value=Math.max(0,saved);return chooseCurve(Math.max(0,saved));}).catch(e=>status(e.message,true));
