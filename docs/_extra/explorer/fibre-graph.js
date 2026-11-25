/* Exact graph arithmetic and the marked special-fibre panel. */
"use strict";
window.FibreGraphPanel = (() => {
  const $ = id => document.getElementById(id);
  const gcd = (a,b) => { a=a<0n?-a:a; b=b<0n?-b:b; while(b) [a,b]=[b,a%b]; return a; };
  class Q {
    constructor(n=0n,d=1n) {
      n=BigInt(n); d=BigInt(d);
      if (!d) throw Error("A charge cannot have denominator zero.");
      if (d<0n) { n=-n; d=-d; }
      const g=gcd(n,d); this.n=n/g; this.d=d/g;
    }
    static from(value) {
      if (value instanceof Q) return value;
      const s=String(value).trim();
      if (s.length>120) throw Error("Use at most 120 characters per charge.");
      if (/^[+-]?\d+\/[+-]?\d+$/.test(s)) { const [n,d]=s.split("/"); return new Q(n,d); }
      if (/^[+-]?\d+$/.test(s)) return new Q(s);
      const decimal=/^([+-]?)(\d*)\.(\d+)$/.exec(s);
      if (decimal) return new Q((decimal[1]==="-"?"-":"")+(decimal[2]||"0")+decimal[3],10n**BigInt(decimal[3].length));
      throw Error("Enter an integer, a fraction such as -3/2, or a decimal.");
    }
    add(b) { b=Q.from(b); return new Q(this.n*b.d+b.n*this.d,this.d*b.d); }
    sub(b) { b=Q.from(b); return new Q(this.n*b.d-b.n*this.d,this.d*b.d); }
    mul(b) { b=Q.from(b); return new Q(this.n*b.n,this.d*b.d); }
    div(b) { b=Q.from(b); return new Q(this.n*b.d,this.d*b.n); }
    get zero() { return this.n===0n; }
    toString() { return String(this.n)+(this.d===1n?"":"/"+this.d); }
    valueOf() { return Number(this.n)/Number(this.d); }
    tex() { return this.d===1n?String(this.n):(this.n<0n?"-":"")+"\\frac{"+String(this.n<0n?-this.n:this.n)+"}{"+this.d+"}"; }
  }
  const zero = () => new Q();
  function solve(graph, charges) {
    const n=graph.genera.length, f=charges.map(Q.from);
    if (f.length!==n) throw Error("Enter one charge for every component.");
    if (!f.reduce((a,b)=>a.add(b),zero()).zero) throw Error("The charges must sum to zero: otherwise Lu=f has no solution.");
    const A=graph.laplacian.slice(1).map((row,i)=>[...row.slice(1).map(Q.from),f[i+1]]);
    for (let j=0;j<n-1;j++) {
      const pivot=A.findIndex((row,i)=>i>=j && !row[j].zero);
      if (pivot<0) throw Error("The fibre graph is disconnected.");
      [A[j],A[pivot]]=[A[pivot],A[j]];
      for (let i=j+1;i<n-1;i++) {
        if (A[i][j].zero) continue;
        const factor=A[i][j].div(A[j][j]);
        for (let k=j;k<n;k++) A[i][k]=A[i][k].sub(factor.mul(A[j][k]));
      }
    }
    const u=Array.from({length:n},zero);
    for (let i=n-2;i>=0;i--) {
      let rhs=A[i][n-1];
      for (let j=i+1;j<n-1;j++) rhs=rhs.sub(A[i][j].mul(u[j+1]));
      u[i+1]=rhs.div(A[i][i]);
    }
    return {f,u,...obstruction(graph,u)};
  }
  function obstruction(graph, values) {
    const u=values.map(Q.from), cycles=graph.cycles.map(row=>row.map(Q.from));
    const form=cycles.map(a=>cycles.map(b=>graph.edges.reduce((sum,[i,j],k)=>sum.add(u[i].add(u[j]).mul(a[k]).mul(b[k])),zero())));
    const positive=graph.genera.flatMap((g,i)=>g?[{vertex:i,value:u[i]}]:[]);
    const constantForm=cycles.map(a=>cycles.map(b=>a.reduce((sum,x,k)=>sum.add(x.mul(b[k]).mul(2)),zero())));
    const signature=[...form.flat(),...positive.map(x=>x.value)];
    const tate=[...constantForm.flat(),...positive.map(()=>Q.from(1))];
    const pivot=tate.findIndex(x=>!x.zero);
    const multiple=pivot<0?zero():signature[pivot].div(tate[pivot]);
    const residual=signature.map((x,i)=>x.sub(multiple.mul(tate[i])));
    return {form,positive,constantForm,multiple,residual,inTateSpan:residual.every(x=>x.zero)};
  }
  function tex(node, value) {
    node.dataset.tex=value;
    if (window.katex) katex.render(value,node,{throwOnError:false,trust:false,output:"htmlAndMathml"});
    else node.textContent=value;
  }
  function cell(row,text,tag="td") { const e=document.createElement(tag); e.textContent=text; row.append(e); return e; }
  function matrix(values) { return values.length?"\\begin{pmatrix}"+values.map(row=>row.map(x=>Q.from(x).tex()).join("&")).join("\\\\")+"\\end{pmatrix}":"()"; }
  function vector(values) { return "("+values.map(x=>Q.from(x).tex()).join(",")+")"; }
  function isPrime(p) {
    if (!Number.isSafeInteger(p) || p<2 || p>2147483647) return false;
    for(let k=2;k*k<=p;k+=(k===2?1:2)) if(p%k===0) return false;
    return true;
  }
  let curve=null, collection=null, record=null, graph=null, selected="", currentPrime=0, request=0, result=null;
  const cache=new Map();
  function message(text,error=false) { $("fibre-status").textContent=text; $("fibre-status").classList.toggle("error",error); }
  function support() { return (collection?.support||[]).filter(p=>!Q.from(p.coefficients[selected]||0).zero); }
  function unavailableGraph(entry) {
    const reason=entry?.reason||entry?.places?.map(x=>x.reason).filter(Boolean).join("; ")||"The exact local model is not available.";
    return {reason, label:/time budget|timeout|exceeded/i.test(reason)?"Graph computation unavailable":"Graph not implemented"};
  }
  function setDivisor(id) {
    selected=id; $("fibre-divisor").value=id;
    if (graph) { renderSupport(); useDivisor(); }
  }
  async function setCurve(next, data) {
    const token=++request;
    curve=next; collection=data?.collection||null; record=null; graph=null; result=null;
    $("fibre-workspace").hidden=true; $("fibre-place-label").hidden=true;
    $("fibre-panel").hidden=false;
    $("fibre-divisor").replaceChildren();
    for (const D of collection?.divisors||[]) {
      const option=document.createElement("option"); option.value=D.id; option.textContent=D.id; $("fibre-divisor").append(option);
    }
    selected=$("divisor-select").value||collection?.divisors?.[0]?.id||"";
    $("fibre-divisor").value=selected;
    if (!curve?.label || !collection) { message("Special-fibre data are available for catalog curves."); return; }
    message("Loading the marked special fibres…");
    try {
      if (!cache.has(curve.label)) {
        const response=await fetch("reduction/"+encodeURIComponent(curve.label)+".json?v=semistable-1");
        if (!response.ok) throw Error("No special-fibre data have been exported for this curve.");
        cache.set(curve.label,await response.json());
      }
      if(token!==request) return;
      record=cache.get(curve.label);
      if(record.schema_version!==1 || record.label!==curve.label) throw Error("Invalid special-fibre data.");
      $("fibre-prime").replaceChildren();
      for(const p of [...record.bad_primes].sort((a,b)=>a-b)) {
        const option=document.createElement("option"); option.value=p;
        const fibre=record.fibres.find(f=>f.prime===p);
        option.textContent=String(p)+(fibre?.status==="available"?"":" — "+unavailableGraph(fibre).label.toLowerCase());
        $("fibre-prime").append(option);
      }
      if (!record.bad_primes.length) { message("No bad reduction primes are recorded for this curve."); return; }
      $("fibre-prime").value=record.fibres.find(f=>f.status==="available")?.prime||record.bad_primes[0];
      showPrime();
    } catch(error) { if(token===request) message(error.message,true); }
  }
  function showPrime() {
    if (!record) return;
    const p=Number($("fibre-prime").value.trim());
    if (!isPrime(p) || !record.bad_primes.includes(p)) {
      graph=null; result=null; $("fibre-workspace").hidden=true;
      message("Choose a bad reduction prime for this curve.",true); return;
    }
    currentPrime=p; result=null; graph=null; $("fibre-workspace").hidden=true;
    $("fibre-place").replaceChildren(); $("fibre-place-label").hidden=true;
    const entry=record.fibres.find(x=>x.prime===p);
    if(!entry || entry.status!=="available") {
      const issue=unavailableGraph(entry);
      message(issue.label+" for "+curve.label+" at p = "+p+". "+issue.reason.replace(/[.\s]+$/,"")+". Select another bad prime to inspect an available fibre.",true); return;
    }
    entry.places.forEach((place,i)=>{
      const option=document.createElement("option"); option.value=i;
      option.textContent="Place "+(i+1)+(place.status==="available"?"":" — "+unavailableGraph(place).label.toLowerCase());
      $("fibre-place").append(option);
    });
    $("fibre-place").value=entry.places.findIndex(x=>x.status==="available");
    $("fibre-place-label").hidden=entry.places.length<2;
    showPlace();
  }
  function showPlace() {
    const entry=record.fibres.find(x=>x.prime===currentPrime), place=entry.places[Number($("fibre-place").value)];
    if (place.status!=="available") {
      const issue=unavailableGraph(place);
      graph=null; result=null; $("fibre-workspace").hidden=true;
      message(issue.label+" for this place above p = "+currentPrime+". "+issue.reason.replace(/[.\s]+$/,"")+".",true); return;
    }
    showGraph(place);
  }
  function showGraph(next) {
    graph=next; result=null;
    $("fibre-workspace").hidden=false; $("fibre-solution").hidden=true;
    const positive=graph.genera.filter(g=>g>0).length;
    message("p = "+currentPrime+" · "+graph.genera.length+" components · cycle rank "+graph.cycles.length+" · "+positive+" positive-genus components.");
    $("fibre-model-note").textContent=
      "Geometric fibre after finite extension. All displayed edges have unit thickness; rational chains resolve the nodes."+
      (graph.regularization_scale>1?" The verified metric graph is rescaled by "+graph.regularization_scale+" before resolution.":"");
    $("fibre-extension").textContent=
      "Local model: "+graph.prime_ideal+"; ramification index "+graph.ramification_index+", residue degree "+graph.residue_degree+". Normalization: "+graph.valuation_normalization+".";
    const rows=$("fibre-charges"); rows.replaceChildren();
    graph.genera.forEach((g,i)=>{
      const row=document.createElement("tr"); row.dataset.vertex=i;
      tex(cell(row,"","th"),"Z_{"+i+"}"); cell(row,g);
      const input=document.createElement("input"); input.type="text"; input.inputMode="text"; input.value="0"; input.dataset.vertex=i;
      input.setAttribute("aria-label","Charge on component Z"+i); cell(row,"").append(input);
      cell(row,"—").className="fibre-potential"; rows.append(row);
    });
    renderSupport(); useDivisor();
  }
  function renderSupport() {
    const rows=$("fibre-support"); rows.replaceChildren(); let unresolved=false;
    for (const P of support()) {
      const mark=graph.specializations.find(x=>x.point_id===P.id), row=document.createElement("tr");
      tex(cell(row,""),P.definition_tex); tex(cell(row,""),Q.from(P.coefficients[selected]).tex());
      const target=cell(row,"");
      if(mark?.status==="known") tex(target,"Z_{"+mark.components[0]+"}");
      else if(mark?.components?.length) { target.textContent="Unresolved between "+mark.components.map(i=>"Z"+i).join(" and "); unresolved=true; }
      else { target.textContent="Specialization not yet determined"; unresolved=true; }
      rows.append(row);
    }
    $("fibre-divisor-charge").disabled=unresolved;
    $("fibre-charge-note").textContent=unresolved?"Some point specializations are unresolved. You can still enter a degree-zero charge below and solve on this verified graph.":
      "Use the divisor charge sets f = −Σ aₚ Zₚ, following the sign in the paper’s integrality criterion. You can edit each charge.";
  }
  function useDivisor() {
    if (!graph) return;
    if($("fibre-divisor-charge").disabled) {
      result=null; $("fibre-solution").hidden=true; $("fibre-error").textContent="";
      $("fibre-charges").querySelectorAll("input").forEach(input=>input.value="0");
      $("fibre-charges").querySelectorAll(".fibre-potential").forEach(cell=>cell.textContent="—");
      drawGraph(); return;
    }
    const f=Array.from({length:graph.genera.length},zero);
    if(!$("fibre-divisor-charge").disabled) for(const P of support()) {
      const mark=graph.specializations.find(x=>x.point_id===P.id);
      f[mark.components[0]]=f[mark.components[0]].sub(P.coefficients[selected]);
    }
    [...$("fibre-charges").querySelectorAll("input")].forEach((input,i)=>input.value=f[i].toString());
    solveInputs();
  }
  function solveInputs() {
    if(!graph) return;
    try {
      result=solve(graph,[...$("fibre-charges").querySelectorAll("input")].map(input=>input.value));
      $("fibre-error").textContent=""; $("fibre-solution").hidden=false;
      [...$("fibre-charges").querySelectorAll(".fibre-potential")].forEach((cell,i)=>tex(cell,result.u[i].tex()));
      tex($("fibre-potential"),"u="+vector(result.u)+",\\qquad u_{Z_0}=0");
      tex($("fibre-cycle-form"),"Q_u="+matrix(result.form));
      tex($("fibre-positive"),"(u_Z)_{g(Z)>0}="+vector(result.positive.map(x=>x.value)));
      tex($("fibre-tate-form"),"Q_{\\mathbf1}="+matrix(result.constantForm)+",\\qquad (1)_{g(Z)>0}="+vector(result.positive.map(()=>1)));
      $("fibre-verdict").textContent=result.inTateSpan?"The obstruction lies in the diagonal Tate line.":"The obstruction does not lie in the diagonal Tate line.";
      tex($("fibre-verdict-formula"),result.inTateSpan?
        "\\operatorname{ob}_{\\mathfrak p}(u)="+result.multiple.tex()+"\\,\\operatorname{ob}_{\\mathfrak p}(\\mathbf1)":
        "\\operatorname{ob}_{\\mathfrak p}(u)\\notin\\mathbb Q\\,\\operatorname{ob}_{\\mathfrak p}(\\mathbf1)");
      $("fibre-verdict").classList.toggle("nonzero",!result.inTateSpan);
      renderCycles(); drawGraph();
    } catch(error) {
      result=null; $("fibre-error").textContent=error.message; $("fibre-solution").hidden=true;
      $("fibre-charges").querySelectorAll(".fibre-potential").forEach(cell=>cell.textContent="—"); drawGraph();
    }
  }
  function renderCycles() {
    const body=$("fibre-cycle-basis"); body.replaceChildren();
    graph.cycles.forEach((cycle,i)=>{
      const row=document.createElement("tr"); tex(cell(row,"","th"),"\\gamma_{"+(i+1)+"}");
      tex(cell(row,""),cycle.map((a,k)=>Q.from(a).zero?"":Q.from(a).tex()+"e_{"+(k+1)+"}").filter(Boolean).join("+").replace(/\+\-/g,"-")||"0"); body.append(row);
    });
    $("fibre-no-cycles").hidden=graph.cycles.length>0;
  }
  function positions() {
    const n=graph.genera.length;
    if(n===1) return [[0,0]];
    const pts=Array.from({length:n},(_,i)=>[Math.cos(2*Math.PI*i/n),Math.sin(2*Math.PI*i/n)]);
    for(let step=0;step<150;step++) {
      const forces=pts.map(()=>[0,0]);
      for(let i=0;i<n;i++) for(let j=0;j<i;j++) {
        const dx=pts[i][0]-pts[j][0],dy=pts[i][1]-pts[j][1],d=Math.max(.005,dx*dx+dy*dy),f=.2/d;
        forces[i][0]+=dx*f;forces[i][1]+=dy*f;forces[j][0]-=dx*f;forces[j][1]-=dy*f;
      }
      for(const [a,b] of graph.edges) {
        const dx=pts[a][0]-pts[b][0],dy=pts[a][1]-pts[b][1],d=Math.hypot(dx,dy),f=.8*(d-.55)/Math.max(.01,d);
        forces[a][0]-=dx*f;forces[a][1]-=dy*f;forces[b][0]+=dx*f;forces[b][1]+=dy*f;
      }
      const rate=.1*(1-step/170);
      for(let i=0;i<n;i++) for(let j=0;j<2;j++) pts[i][j]+=Math.max(-.15,Math.min(.15,forces[i][j]*rate));
    }
    return pts;
  }
  function drawGraph() {
    if(!graph || !window.Plotly) return;
    const pts=positions(), traces=[], pairs=new Map();
    graph.edges.forEach(([a,b],k)=>{const key=[Math.min(a,b),Math.max(a,b)].join(",");if(!pairs.has(key))pairs.set(key,[]);pairs.get(key).push(k);});
    for(const indices of pairs.values()) indices.forEach((k,index)=>{
      const [a,b]=graph.edges[k],p=pts[a],q=pts[b],dx=q[0]-p[0],dy=q[1]-p[1],d=Math.max(.01,Math.hypot(dx,dy));
      const bend=(index-(indices.length-1)/2)*.35;
      const x=[],y=[];
      for(let i=0;i<=18;i++){const t=i/18,s=4*t*(1-t)*bend;x.push(p[0]+t*dx-s*dy/d);y.push(p[1]+t*dy+s*dx/d);}
      traces.push({x,y,type:"scatter",mode:"lines",line:{color:"#82968a",width:2},text:"e"+(k+1)+": Z"+a+" → Z"+b,hoverinfo:"text",showlegend:false});
    });
    const marks=graph.genera.map(()=>[]);
    for(const P of support()) {const mark=graph.specializations.find(x=>x.point_id===P.id);if(mark?.status==="known")marks[mark.components[0]].push(P.id);}
    const potentials=result?result.u.map(v=>{const x=Number(v);return Number.isFinite(x)?x:0;}):null;
    const limit=potentials?Math.max(1,...potentials.map(Math.abs)):1;
    traces.push({x:pts.map(p=>p[0]),y:pts.map(p=>p[1]),type:"scatter",mode:"markers+text",text:graph.genera.map((g,i)=>"Z"+i),textposition:"top center",textfont:{family:"Inter, sans-serif",size:12},
      marker:{size:graph.genera.map((g,i)=>g?26:marks[i].length?21:15),color:potentials||graph.genera.map(g=>g?"#267f91":"#fff8ec"),cmin:-limit,cmax:limit,colorscale:[[0,"#267f91"],[.5,"#fff8ec"],[1,"#bd5034"]],line:{color:"#183e3b",width:1.5}},
      customdata:graph.genera.map((g,i)=>i),hovertext:graph.genera.map((g,i)=>"Z"+i+" · genus "+g+"<br>"+(marks[i].join(", ")||"No selected support points")+(result?"<br>f = "+result.f[i]+"; u = "+result.u[i]:"")),hoverinfo:"text",showlegend:false});
    Plotly.react($("fibre-plot"),traces,{margin:{l:20,r:20,t:35,b:20},paper_bgcolor:"rgba(0,0,0,0)",plot_bgcolor:"rgba(0,0,0,0)",xaxis:{visible:false},yaxis:{visible:false,scaleanchor:"x",scaleratio:1},hovermode:"closest",uirevision:curve?.label+":"+currentPrime+":"+$("fibre-place").value},
      {responsive:true,displayModeBar:false}).then(()=>{
        $("fibre-plot").removeAllListeners("plotly_click");
        $("fibre-plot").on("plotly_click",event=>{const i=event.points[0]?.customdata;if(Number.isInteger(i))$("fibre-charges").querySelector('input[data-vertex="'+i+'"]').focus();});
      });
  }
  function init() {
    $("fibre-prime-form").addEventListener("submit",event=>{event.preventDefault();showPrime();});
    $("fibre-prime").addEventListener("change",showPrime);
    $("fibre-place").addEventListener("change",showPlace);
    $("fibre-divisor").addEventListener("change",()=>{const selector=$("divisor-select");selector.value=$("fibre-divisor").value;selector.dispatchEvent(new Event("change"));});
    $("fibre-solve").addEventListener("submit",event=>{event.preventDefault();solveInputs();});
    $("fibre-divisor-charge").addEventListener("click",useDivisor);
    $("fibre-zero-charge").addEventListener("click",()=>{$("fibre-charges").querySelectorAll("input").forEach(input=>input.value="0");solveInputs();});
  }
  return {init,setCurve,setDivisor,math:{Q,solve,obstruction,isPrime}};
})();
