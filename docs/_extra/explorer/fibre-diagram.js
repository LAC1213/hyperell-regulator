/* Draw the verified dual graph, with its point labels outside the edges. */
"use strict";
const FibreDiagram = (() => {
  function positions(graph, width) {
    const n = graph.genera.length;
    if (n === 1) return [[width / 2, 0]];
    const neighbours = Array.from({length: n}, () => new Set());
    for (const [a, b] of graph.edges) {
      neighbours[a].add(b);
      neighbours[b].add(a);
    }
    const links = neighbours.map(entries => [...entries]);
    const points = Array(n);
    const circle = order => {
      const radius = Math.min(112, width * .26);
      order.forEach((vertex, i) => {
        const angle = 2 * Math.PI * i / order.length - Math.PI / 2;
        points[vertex] = [width / 2 + radius * Math.cos(angle), radius * Math.sin(angle)];
      });
      return points;
    };

    if (links.every(row => row.length <= 2)) {
      const end = links.findIndex(row => row.length === 1);
      const order = [];
      let previous = -1, current = Math.max(0, end);
      while (current !== undefined && order.length < n) {
        order.push(current);
        const next = links[current].find(vertex => vertex !== previous);
        previous = current;
        current = next;
      }
      if (end === -1) return circle(order);
      const columns = Math.min(n, Math.floor((width - 68) / 38) + 1);
      const step = Math.min(112, (width - 68) / (columns - 1));
      const left = (width - (columns - 1) * step) / 2;
      order.forEach((vertex, i) => {
        const row = Math.floor(i / columns), column = i % columns;
        points[vertex] = [left + step * (row % 2 ? columns - 1 - column : column), row * 46];
      });
      return points;
    }

    // Place the leaves of a branched tree in order, then centre each parent.
    if (links.reduce((total, row) => total + row.length, 0) === 2 * (n - 1)) {
      const root = links.reduce((best, row, i) => row.length > links[best].length ? i : best, 0);
      let leaves = 0;
      function place(vertex, parent, depth) {
        const children = links[vertex].filter(child => child !== parent);
        const start = leaves;
        if (children.length) children.forEach(child => place(child, vertex, depth + 1));
        else leaves++;
        points[vertex] = [(start + leaves - 1) / 2, depth * 46];
      }
      place(root, -1, 0);
      return points.map(([x, y]) => [34 + x * (width - 68) / Math.max(1, leaves - 1), y]);
    }
    return circle(graph.genera.map((_, i) => i));
  }

  function svgElement(name, attributes) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", name);
    for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, value);
    return node;
  }

  function drawEdges(svg, graph, points) {
    const pairs = new Map();
    graph.edges.forEach(([a, b], index) => {
      const left = Math.min(a, b), right = Math.max(a, b), key = left + ":" + right;
      if (!pairs.has(key)) pairs.set(key, {left, right, edges: []});
      pairs.get(key).edges.push(index);
    });
    const occupied = [];
    for (const {left, right, edges} of pairs.values()) {
      const [x1, y1] = points[left], [x2, y2] = points[right];
      const dx = x2 - x1, dy = y2 - y1, length = Math.hypot(dx, dy) || 1;
      edges.forEach((index, i) => {
        const offset = (i - (edges.length - 1) / 2) * Math.min(34, length * .8);
        const cx = (x1 + x2) / 2 - dy * offset / length;
        const cy = (y1 + y2) / 2 + dx * offset / length;
        const path = svgElement("path", {
          class: "fibre-edge",
          d: ["M", x1, y1, "Q", cx, cy, x2, y2].join(" ")
        });
        const [a, b] = graph.edges[index];
        const title = svgElement("title", {});
        title.textContent = "Edge e" + (index + 1) + ": Z" + a + " → Z" + b;
        path.append(title);
        svg.append(path);
        // Keep labels clear of the edge, including its curved portion.
        const steps = Math.max(2, Math.ceil((length + Math.abs(offset)) / 7));
        for (let k = 0; k <= steps; k++) {
          const t = k / steps, s = 1 - t;
          occupied.push([s * s * x1 + 2 * s * t * cx + t * t * x2,
            s * s * y1 + 2 * s * t * cy + t * t * y2]);
        }
      });
    }
    return occupied;
  }

  const directions = [[0, -1], [0, 1], [-1, 0], [1, 0], [-1, -1], [1, -1], [-1, 1], [1, 1]];
  function positionLabel(node, point, width, height, centre, occupied, placed) {
    const [x, y] = point, w = node.offsetWidth, h = node.offsetHeight;
    const outward = [x - width / 2, y - centre];
    const priority = ([dx, dy]) => (dx * outward[0] + dy * outward[1]) / Math.hypot(dx, dy);
    const order = [...directions].sort((a, b) => priority(b) - priority(a));
    const clamp = left => Math.max(4, Math.min(width - w - 4, left));
    const fits = (left, top) => top >= 4 && top + h <= height - 4 &&
      !occupied.some(([px, py]) => px >= left - 5 && px <= left + w + 5 && py >= top - 5 && py <= top + h + 5) &&
      !placed.some(box => left < box.left + box.w + 5 && left + w + 5 > box.left &&
        top < box.top + box.h + 5 && top + h + 5 > box.top);
    for (const gap of [14, 26, 42, 64, 92, 130]) {
      for (const [dx, dy] of order) {
        const left = clamp(x + (dx === 0 ? -w / 2 : dx > 0 ? gap : -gap - w));
        const top = y + (dy === 0 ? -h / 2 : dy > 0 ? gap : -gap - h);
        if (fits(left, top)) return {left, top, w, h};
      }
    }
    return {left: clamp(x - w / 2), top: height + 6, w, h};
  }

  function render(host, graph, {marks, charges, potentials, tex, onSelect}) {
    const width = Math.max(300, host.clientWidth);
    const points = positions(graph, width);
    const top = Math.min(...points.map(([, y]) => y));
    const bottom = Math.max(...points.map(([, y]) => y));
    points.forEach(point => point[1] += 74 - top);
    let height = Math.max(168, bottom - top + 148);
    const diagram = document.createElement("div");
    diagram.className = "fibre-diagram";
    diagram.style.width = width + "px";
    const svg = svgElement("svg", {width: "100%", height: "100%", "aria-hidden": true});
    diagram.append(svg);
    host.replaceChildren(diagram);
    const occupied = [...points, ...drawEdges(svg, graph, points)], placed = [];

    // Fit the marked components first so their point lists stay near their vertices.
    const order = graph.genera.map((_, i) => i).sort((a, b) => marks[b].length - marks[a].length);
    for (const i of order) {
      const genus = graph.genera[i];
      if (points.length > 12 && !genus && !marks[i].length) continue;
      const note = document.createElement("div");
      note.id = "fibre-node-note-" + i;
      note.className = "fibre-node-note";
      const title = document.createElement("div"), name = document.createElement("span");
      title.className = "fibre-node-name";
      tex(name, "Z_{" + i + "}");
      title.append(name, " · genus " + genus);
      note.append(title);
      for (const point of marks[i]) {
        const row = document.createElement("div");
        row.className = "fibre-point";
        tex(row, point);
        note.append(row);
      }
      diagram.append(note);
      const [x, y] = points[i];
      const box = positionLabel(note, points[i], width, height, (bottom - top) / 2 + 74, occupied, placed);
      const {left, top: y0, w, h} = box;
      note.style.left = left + "px";
      note.style.top = y0 + "px";
      const endX = Math.max(left, Math.min(left + w, x)), endY = Math.max(y0, Math.min(y0 + h, y));
      svg.append(svgElement("line", {class: "fibre-label-guide", x1: x, y1: y, x2: endX, y2: endY}));
      placed.push(box);
      height = Math.max(height, y0 + h + 6);
    }

    const values = potentials?.map(value => Number.isFinite(Number(value)) ? Number(value) : 0);
    const limit = Math.max(1, ...(values || []).map(Math.abs));
    points.forEach(([x, y], i) => {
      const genus = graph.genera[i], value = values?.[i], size = genus ? 19 : marks[i].length ? 15 : 11;
      const button = document.createElement("button");
      button.className = "fibre-node";
      button.type = "button";
      button.style.width = button.style.height = size + "px";
      button.style.left = x + "px";
      button.style.top = y + "px";
      const color = value < 0 ? "#267f91" : "#bd5034";
      button.style.backgroundColor = value ? "color-mix(in srgb, " + color + " " +
        (25 + 75 * Math.abs(value) / limit) + "%, var(--panel))" : "var(--panel)";
      button.title = "Z" + i + " · genus " + genus +
        (potentials ? " · charge " + charges[i] + " · potential " + potentials[i] : "");
      button.setAttribute("aria-label", button.title);
      if (diagram.querySelector("#fibre-node-note-" + i)) button.setAttribute("aria-describedby", "fibre-node-note-" + i);
      else button.tabIndex = -1;
      button.addEventListener("click", () => onSelect(i));
      diagram.append(button);
    });
    diagram.style.height = height + "px";
  }

  return {render, positions};
})();
if (typeof module !== "undefined") module.exports = FibreDiagram;
