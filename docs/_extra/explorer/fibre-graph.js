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
    if (!f.reduce((a,b)=>a.add(b),zero()).zero) throw Error("The charges must sum to zero.");
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
  let markedPoints=[], diagramWidth=0;
  const cache=new Map();
  function message(text,error=false) { $("fibre-status").textContent=text; $("fibre-status").classList.toggle("error",error); }
  function support() { return (collection?.support||[]).filter(p=>!Q.from(p.coefficients[selected]||0).zero); }
  function unavailableGraph(entry) {
    const reason=entry?.reason||entry?.places?.map(x=>x.reason).filter(Boolean).join("; ")||"The exact local model is not available.";
    return {reason, label:/time budget|timeout|exceeded/i.test(reason)?"Graph computation unavailable":"Graph not implemented"};
  }
  function setDivisor(id) {
    selected=id;
    if (graph) { renderSupport(); useDivisor(); }
  }
  async function setCurve(next, data) {
    const token=++request;
    curve=next; collection=data?.collection||null; record=null; graph=null; result=null;
    $("fibre-workspace").hidden=true; $("fibre-place-label").hidden=true;
    $("fibre-panel").hidden=false;
    $("fibre-divisor").replaceChildren();
    for (const id of CurveGeometry.divisorIds(collection)) {
      const option=document.createElement("option");
      option.value=id; option.textContent=CurveGeometry.divisorLabel(collection,id);
      $("fibre-divisor").append(option);
    }
    selected=$("fibre-divisor").value;
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
    message("");
    const rows=$("fibre-charges"); rows.replaceChildren();
    graph.genera.forEach((_,i)=>{
      const row=document.createElement("tr");
      const component=cell(row,"","th"); component.scope="row"; tex(component,"Z_{"+i+"}");
      const input=document.createElement("input"); input.type="text"; input.inputMode="text"; input.value="0"; input.dataset.vertex=i;
      input.setAttribute("aria-label","Charge on component Z"+i); cell(row,"").append(input);
      cell(row,"—").className="fibre-potential"; rows.append(row);
    });
    renderSupport(); useDivisor();
    if (document.fonts?.status==="loading") {
      document.fonts.ready.then(()=>{if (graph===next) drawGraph();});
    }
  }
  function renderSupport() {
    const specializations=new Map(graph.specializations.map(mark=>[mark.point_id,mark]));
    markedPoints=support().map(point=>{
      const mark=specializations.get(point.id), coefficient=Q.from(point.coefficients[selected]);
      const multiple=coefficient.toString();
      const prefix=multiple==="1"?"":multiple==="-1"?"-":coefficient.tex()+"\\,";
      return {id:point.id, mark, coefficient, formula:prefix+CurveGeometry.pointTex(point)};
    });
    const unresolved=markedPoints.filter(({mark})=>mark?.status!=="known");
    $("fibre-divisor-charge").disabled=unresolved.length>0;
    const note=$("fibre-unresolved");
    note.replaceChildren();
    note.hidden=!unresolved.length;
    if (unresolved.length) {
      const items=document.createElement("ul");
      for (const {mark,formula} of unresolved) {
        const item=document.createElement("li"), point=document.createElement("span");
        tex(point,formula);
        item.append(point);
        if(mark?.components?.length) item.append(" · possible components "+mark.components.map(i=>"Z"+i).join(", "));
        items.append(item);
      }
      note.append("Specialization unresolved for these points. You can enter a custom charge.",items);
    }
  }
  function divisorCharges(model) {
    if (model.status!=="available") return null;
    const specializations=new Map(model.specializations.map(mark=>[mark.point_id,mark]));
    const charges=model.genera.map(zero);
    for (const {id,coefficient} of markedPoints) {
      const mark=specializations.get(id);
      if (mark?.status!=="known") return null;
      const vertex=mark.components[0];
      charges[vertex]=charges[vertex].add(coefficient);
    }
    return charges;
  }
  function useDivisor() {
    if (!graph) return;
    const charges=divisorCharges(graph);
    if (!charges) {
      result=null; $("fibre-solution").hidden=true; $("fibre-error").textContent="";
      $("fibre-charges").querySelectorAll("input").forEach(input=>input.value="0");
      $("fibre-charges").querySelectorAll(".fibre-potential").forEach(cell=>cell.textContent="—");
      drawGraph(); return;
    }
    [...$("fibre-charges").querySelectorAll("input")].forEach((input,i)=>input.value=charges[i].toString());
    solveInputs();
  }
  function solveCharges(model, charges) {
    // The editor displays f in Lu + f = 0. The solver takes the right-hand side.
    return solve(model,charges.map(charge=>zero().sub(charge)));
  }
  function renderVerdict(charges) {
    const verdict=$("fibre-verdict");
    verdict.classList.toggle("nonzero",!result.inTateSpan);
    if (!result.inTateSpan) {
      verdict.textContent="The obstruction does not lie in the diagonal Tate line.";
      return;
    }
    const membership="The obstruction lies in the diagonal Tate line";
    const divisor=divisorCharges(graph);
    if (!divisor?.every((charge,i)=>charge.sub(charges[i]).zero)) {
      verdict.textContent=membership+" for this charge.";
      return;
    }
    const fibre=record.fibres.find(entry=>entry.prime===currentPrime);
    const passesAtPrime=fibre.places.every(place=>{
      if (place===graph) return true;
      const localCharges=divisorCharges(place);
      return localCharges && solveCharges(place,localCharges).inTateSpan;
    });
    const prime=document.createElement("span");
    prime.className="fibre-verdict-prime";
    const value=document.createElement("span");
    tex(value,"p="+currentPrime);
    prime.append(value,".");
    verdict.textContent=membership+(passesAtPrime
      ? ", and hence the class is integral at "
      : " at the selected place above ");
    verdict.append(prime);
  }
  function solveInputs() {
    if(!graph) return;
    try {
      const charges=[...$("fibre-charges").querySelectorAll("input")].map(input=>Q.from(input.value));
      result=solveCharges(graph,charges);
      $("fibre-error").textContent=""; $("fibre-solution").hidden=false;
      [...$("fibre-charges").querySelectorAll(".fibre-potential")].forEach((cell,i)=>tex(cell,result.u[i].tex()));
      tex($("fibre-cycle-form"),"Q_u="+matrix(result.form)+",");
      tex($("fibre-positive"),"(u_Z)_{g(Z)\\geq 1}="+vector(result.positive.map(x=>x.value))+".");
      renderVerdict(charges);
      drawGraph();
    } catch(error) {
      result=null; $("fibre-error").textContent=error.message; $("fibre-solution").hidden=true;
      $("fibre-charges").querySelectorAll(".fibre-potential").forEach(cell=>cell.textContent="—"); drawGraph();
    }
  }
  function drawGraph() {
    const host=$("fibre-plot");
    diagramWidth=host.clientWidth;
    const marks=graph.genera.map(()=>[]);
    for (const {mark,formula} of markedPoints) {
      if (mark?.status==="known") marks[mark.components[0]].push(formula);
    }
    const inputs=[...$("fibre-charges").querySelectorAll("input")];
    FibreDiagram.render(host,graph,{
      marks, charges:inputs.map(input=>input.value), potentials:result?.u, tex,
      onSelect:index=>inputs[index].focus()
    });
  }
  function init() {
    if (typeof ResizeObserver !== "undefined") {
      const host=$("fibre-plot");
      new ResizeObserver(()=>{
        if (graph && host.clientWidth && host.clientWidth!==diagramWidth) drawGraph();
      }).observe(host);
    }
    $("fibre-prime-form").addEventListener("submit",event=>{event.preventDefault();showPrime();});
    $("fibre-prime").addEventListener("change",showPrime);
    $("fibre-place").addEventListener("change",showPlace);
    $("fibre-divisor").addEventListener("change",()=>setDivisor($("fibre-divisor").value));
    $("fibre-solve").addEventListener("submit",event=>{event.preventDefault();solveInputs();});
    $("fibre-divisor-charge").addEventListener("click",useDivisor);
    $("fibre-zero-charge").addEventListener("click",()=>{$("fibre-charges").querySelectorAll("input").forEach(input=>input.value="0");solveInputs();});
  }
  return {init,setCurve,math:{Q,solve,obstruction,isPrime}};
})();
