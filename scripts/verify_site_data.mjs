import fs from "node:fs";

const data = JSON.parse(fs.readFileSync("assets/data/site-data.json", "utf8"));

const schema = {
  library_id: 13,
};

const projectionMap = {
  spatial: { x: 0, y: 1 },
  umap: { x: 2, y: 3 },
};

function libraryBounds() {
  const libs = data.annotations.library_id.map(() => ({
    minX: Infinity,
    maxX: -Infinity,
    minY: Infinity,
    maxY: -Infinity,
  }));
  for (const cell of data.cells) {
    const bounds = libs[cell[schema.library_id]];
    if (!bounds) continue;
    const [x, y] = cell;
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    bounds.minX = Math.min(bounds.minX, x);
    bounds.maxX = Math.max(bounds.maxX, x);
    bounds.minY = Math.min(bounds.minY, y);
    bounds.maxY = Math.max(bounds.maxY, y);
  }
  return libs;
}

const libs = libraryBounds();

function coordinateFor(cell, projection) {
  if (projection === "slices") {
    const code = cell[schema.library_id];
    const bounds = libs[code];
    const xSpan = bounds.maxX - bounds.minX || 1;
    const ySpan = bounds.maxY - bounds.minY || 1;
    return [
      (cell[0] - bounds.minX) / xSpan + code * 1.18,
      (cell[1] - bounds.minY) / ySpan,
    ];
  }
  const map = projectionMap[projection];
  return [cell[map.x], cell[map.y]];
}

function boundsFor(projection) {
  if (projection === "slices") {
    return { minX: 0, maxX: libs.length * 1.18 - 0.18, minY: 0, maxY: 1 };
  }
  let minX = Infinity;
  let maxX = -Infinity;
  let minY = Infinity;
  let maxY = -Infinity;
  for (const cell of data.cells) {
    const [x, y] = coordinateFor(cell, projection);
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    minX = Math.min(minX, x);
    maxX = Math.max(maxX, x);
    minY = Math.min(minY, y);
    maxY = Math.max(maxY, y);
  }
  return { minX, maxX, minY, maxY };
}

function projectedSummary(projection, width = 1200, height = 760) {
  const bounds = boundsFor(projection);
  const pad = 26;
  const usableW = width - pad * 2;
  const usableH = height - pad * 2;
  const xSpan = bounds.maxX - bounds.minX || 1;
  const ySpan = bounds.maxY - bounds.minY || 1;
  const scale = Math.min(usableW / xSpan, usableH / ySpan);
  let drawable = 0;
  let inBounds = 0;
  for (const cell of data.cells) {
    const point = coordinateFor(cell, projection);
    if (!point || !Number.isFinite(point[0]) || !Number.isFinite(point[1])) continue;
    drawable += 1;
    const x = pad + (usableW - xSpan * scale) / 2 + (point[0] - bounds.minX) * scale;
    const y = height - (pad + (usableH - ySpan * scale) / 2) - (point[1] - bounds.minY) * scale;
    if (x >= 0 && x <= width && y >= 0 && y <= height) inBounds += 1;
  }
  return { projection, drawable, inBounds, bounds };
}

const summaries = ["spatial", "umap", "slices"].map((projection) => projectedSummary(projection));
console.log(JSON.stringify({
  cells: data.metadata.n_cells,
  genes: data.metadata.n_genes,
  libraries: data.annotations.library_id,
  summaries,
}, null, 2));
