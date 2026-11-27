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
    return pieces.map(piece => {
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
      const variable = power === 0 ? "" : power === 1 ? "x" : "x^{" + power + "}";
      const magnitude = rational(coefficient) === 1 && power ? "" : rationalTex(coefficient);
      return {negative, tex: magnitude + variable};
    });
  }

  function equationTex(source, multiline = false) {
    const terms = parsePolynomial(source);
    const signed = terms.map((term, i) => (term.negative ? "-" : i ? "+" : "") + term.tex);
    if (!multiline || terms.length <= 3) return "C:\\quad y^2=" + signed.join("");
    const rows = [];
    for (let i = 0; i < signed.length; i += 3) {
      rows.push((i ? "&\\quad " : "C:\\quad y^2&=") + signed.slice(i, i + 3).join(""));
    }
    return "\\begin{aligned}" + rows.join("\\\\") + "\\end{aligned}";
  }

  function divisorIds(collection) {
    const ids = (collection?.divisors || []).map(divisor => divisor.id);
    return [...new Set([...(collection?.matrix_divisor_ids || []), ...ids])].filter(id => ids.includes(id));
  }

  function divisorLabel(collection, id) {
    const index = divisorIds(collection).indexOf(id);
    return index === -1 ? id : "D" + (index + 1);
  }

  function rationalParts(value = "0") {
    return String(value).includes("/") ? String(value).split("/").map(BigInt) : [BigInt(value), 1n];
  }

  function gcd(a, b) {
    while (b !== 0n) [a, b] = [b, a % b];
    return a < 0n ? -a : a;
  }

  function quadraticCoordinates(polynomial, index, coordinates) {
    const fractions = polynomial.map(rationalParts);
    const scale = fractions.reduce((product, [, denominator]) => product * denominator, 1n);
    const [c, b, a] = fractions.map(([numerator, denominator]) => numerator * (scale / denominator));
    const centered = a === 1n && b === 0n;
    const radicand = centered ? -c : b * b - 4n * a * c;
    const rootDenominator = centered ? 1n : 2n * (a < 0n ? -a : a);
    const offset = a < 0n ? b : -b;
    const magnitude = radicand < 0n ? -radicand : radicand;
    const radical = (radicand < 0n ? "i" : "") + (magnitude === 1n ? "" : "\\sqrt{" + magnitude + "}");
    return coordinates.map(([constant, linear]) => {
      const [m, n] = rationalParts(constant), [p, q] = rationalParts(linear);
      let real = m * q * rootDenominator + p * n * offset;
      let surd = (index === 0 ? -p : p) * n;
      let denominator = n * q * rootDenominator;
      const factor = gcd(gcd(real, surd), denominator);
      real /= factor;
      surd /= factor;
      denominator /= factor;
      if (!surd) return rationalTex(real + "/" + denominator);
      const sign = surd < 0n ? "-" : real !== 0n ? "+" : "";
      const coefficient = surd < 0n ? -surd : surd;
      const numerator = (real === 0n ? "" : String(real)) + sign + (coefficient === 1n ? "" : coefficient) + radical;
      return denominator === 1n ? numerator : "\\frac{" + numerator + "}{" + denominator + "}";
    });
  }

  function polynomialTex(coefficients) {
    const terms = [];
    for (let power = coefficients.length - 1; power >= 0; power--) {
      const coefficient = coefficients[power];
      if (coefficient === "0") continue;
      const unsigned = coefficient.replace(/^-/, "");
      const sign = coefficient.startsWith("-") ? "-" : terms.length ? "+" : "";
      const variable = power === 0 ? "" : power === 1 ? "t" : "t^{" + power + "}";
      terms.push(sign + (variable && unsigned === "1" ? "" : rationalTex(unsigned)) + variable);
    }
    return terms.join("") || "0";
  }

  function infinityLabels(tex) {
    return tex.replace(/\\infty\s*_\s*\{\s*Y\s*=\s*/g, "\\infty_{");
  }

  function pointTex(point) {
    const record = point.lmfdb?.point;
    if (!record?.field_polynomial) return infinityLabels(point.tex);
    const {field_polynomial: polynomial, coordinates, root_index: index} = record;
    const infinite = coordinates[0] === "infinity";
    const values = infinite ? coordinates.slice(1) : coordinates;
    const quadratic = polynomial.length === 3;
    const tex = quadratic ? quadraticCoordinates(polynomial, index, values) : values.map(polynomialTex);
    const pair = infinite ? "\\infty_{" + tex[0] + "}" : "\\left(" + tex.join(",") + "\\right)";
    return quadratic ? pair : pair + "\\big|_{t=\\operatorname{Root}_{" + index + "}(" + polynomialTex(polynomial) + ")}";
  }

  function divisorTex(collection, divisor) {
    const points = new Map((collection.support || []).map(point => [point.tex, point]));
    const terms = collection.compact_divisors?.[divisor.id] ?? divisor.tex.slice(divisor.tex.indexOf("=") + 1)
      .replace(/P_\{\d+\}/g, name => points.has(name) ? pointTex(points.get(name)) : name);
    return "D_{" + divisorLabel(collection, divisor.id).slice(1) + "}=" + infinityLabels(terms);
  }

  return {rationalTex, equationTex, divisorIds, divisorLabel, divisorTex, pointTex};
})();
if (typeof module !== "undefined") module.exports = CurveGeometry;
