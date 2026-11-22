/* Polynomial formatting for y^2 = F(x). Numerical regulator data are never changed. */
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

  function equationTex(source, multiline = false) {
    const terms = parsePolynomial(source).terms;
    const signed = terms.map((term, i) => (term.negative ? "-" : i ? "+" : "") + term.tex);
    if (!multiline || terms.length <= 3) return "C:\\quad y^2=" + signed.join("");
    const rows = [];
    for (let i = 0; i < signed.length; i += 3) {
      rows.push((i ? "&\\quad " : "C:\\quad y^2&=") + signed.slice(i, i + 3).join(""));
    }
    return "\\begin{aligned}" + rows.join("\\\\") + "\\end{aligned}";
  }

  return {rational, rationalTex, parsePolynomial, equationTex};
})();
if (typeof module !== "undefined") module.exports = CurveGeometry;
