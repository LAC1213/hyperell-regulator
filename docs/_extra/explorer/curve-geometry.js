/* Display geometry for y^2 = F(x). Numerical regulator data are never changed. */
"use strict";
const CurveGeometry = (() => {
  function rational(value) {
    const parts = String(value).trim().split("/");
    if (parts.length > 2 || !parts.every(v => /^[+-]?\d+$/.test(v))) throw Error("Invalid rational coefficient");
    const result = Number(parts[0]) / (parts.length === 2 ? Number(parts[1]) : 1);
    if (!Number.isFinite(result)) throw Error("Non-finite coefficient");
    return result;
  }

  function rationalTex(value) {
    const parts = String(value).trim().split("/");
    rational(value);
    if (parts.length === 1 || parts[1] === "1") return parts[0];
    const sign = parts[0].startsWith("-") ? "-" : "";
    return sign + "\\frac{" + parts[0].replace(/^[+-]/, "") + "}{" + parts[1] + "}";
  }

  function parsePolynomial(source) {
    const compact = source.replace(/\s+/g, "");
    const pieces = compact.match(/[+-]?[^+-]+/g);
    if (!pieces || pieces.join("") !== compact) throw Error("Expected an expanded polynomial");
    const coefficients = Array(13).fill(0);
    const terms = pieces.map(piece => {
      const negative = piece[0] === "-";
      const unsigned = piece.replace(/^[+-]/, "");
      const match = unsigned.match(/^(?:(\d+(?:\/\d+)?)\*?)?x(?:\^(\d+))?$/);
      let coefficient, power;
      if (match) {
        coefficient = match[1] || "1";
        power = match[2] === undefined ? 1 : Number(match[2]);
      } else {
        if (!/^\d+(?:\/\d+)?$/.test(unsigned)) throw Error("Expected rational coefficients");
        coefficient = unsigned;
        power = 0;
      }
      if (power > 12) throw Error("Display degree exceeds 12");
      const value = rational(coefficient) * (negative ? -1 : 1);
      coefficients[power] += value;
      const variable = power === 0 ? "" : power === 1 ? "x" : "x^{" + power + "}";
      const magnitude = rational(coefficient) === 1 && power ? "" : rationalTex(coefficient);
      return {coefficient: value, power, negative, tex: magnitude + variable};
    });
    while (coefficients.length > 1 && coefficients.at(-1) === 0) coefficients.pop();
    return {coefficients, terms};
  }

  function equationTex(source) {
    const terms = parsePolynomial(source).terms;
    const signed = terms.map((term, i) => (term.negative ? "-" : i ? "+" : "") + term.tex);
    if (terms.length <= 5 && signed.join("").length < 72) return "C:\\quad y^2=" + signed.join("");
    const rows = [];
    for (let i = 0; i < signed.length; i += 3) {
      rows.push((i ? "&\\quad " : "C:\\quad y^2&=") + signed.slice(i, i + 3).join(""));
    }
    return "\\begin{aligned}" + rows.join("\\\\") + "\\end{aligned}";
  }

  function evaluate(coefficients, x) {
    return coefficients.reduceRight((value, coefficient) => value * x + coefficient, 0);
  }

  function realLocus(data) {
    const {coefficients} = parsePolynomial(data.input.F);
    const roots = data.roots
      .filter(([re, im]) => Math.abs(Number(im)) <= 1e-10 * (1 + Math.abs(Number(re))))
      .map(([re]) => Number(re)).sort((a, b) => a - b)
      .filter((x, i, list) => !i || x - list[i - 1] > 1e-10 * (1 + Math.abs(x)));
    const points = (data.collection?.points || []).flatMap((point, i) => {
      if (!Array.isArray(point)) return [];
      const [x, y] = point.map(rational);
      return [{name: i === 0 ? "P" : "Q", x, y}];
    });
    const anchors = [...roots, ...points.map(p => p.x)];
    const low = anchors.length ? Math.min(...anchors) : -1;
    const high = anchors.length ? Math.max(...anchors) : 1;
    const span = Math.max(high - low, 1);
    const xRange = [low - .2 * span, high + .2 * span];
    const bounds = [xRange[0], ...roots, xRange[1]];
    const branches = [];
    let boundedHeight = 0;
    for (let interval = 0; interval < bounds.length - 1; interval++) {
      const a = bounds[interval], b = bounds[interval + 1];
      if (evaluate(coefficients, (a + b) / 2) <= 0) continue;
      const count = Math.max(80, Math.ceil(600 * (b - a) / (xRange[1] - xRange[0])));
      const x = [], positive = [];
      for (let i = 0; i <= count; i++) {
        // Cosine spacing resolves vertical tangents and very short components.
        const value = a + (b - a) * (1 - Math.cos(Math.PI * i / count)) / 2;
        const endpoint = i === 0 && interval > 0 || i === count && interval < bounds.length - 2;
        let square = endpoint ? 0 : evaluate(coefficients, value);
        const tolerance = 64 * Number.EPSILON * coefficients.reduceRight((v, c) => v * Math.abs(value) + Math.abs(c), 0);
        if (square < 0 && square >= -tolerance) square = 0;
        x.push(value);
        positive.push(square < 0 ? null : Math.sqrt(square));
      }
      if (interval > 0 && interval < bounds.length - 2) boundedHeight = Math.max(boundedHeight, ...positive.filter(y => y !== null));
      branches.push({x, y: positive}, {x, y: positive.map(y => y === null ? null : -y)});
    }
    const markedHeight = Math.max(0, ...points.map(p => Math.abs(p.y)));
    const height = Math.max(1.8, markedHeight * 1.6, boundedHeight * 1.18, span * .4);
    return {coefficients, roots, branches, points, xRange, yRange: [-height, height]};
  }

  return {rational, rationalTex, parsePolynomial, equationTex, evaluate, realLocus};
})();
if (typeof module !== "undefined") module.exports = CurveGeometry;
