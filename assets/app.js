const DATA_URL = "assets/data/site-data.json";

const projectionMap = {
  spatial: { label: "Spatial map", x: 0, y: 1 },
  umap: { label: "UMAP", x: 2, y: 3 },
  slices: { label: "Library split", split: true },
};

const categoryLabels = {
  cell_types: "Cell type",
  finer_cell_types: "Fine cell type",
  leiden: "Leiden",
  hybrid_leiden: "Hybrid Leiden",
  library_id: "Library",
  timepoint: "Timepoint",
  excitatory: "Excitatory program",
  inhibitory: "Inhibitory program",
};

const schema = {
  cell_types: 9,
  finer_cell_types: 10,
  leiden: 11,
  hybrid_leiden: 12,
  library_id: 13,
  timepoint: 14,
  excitatory: 15,
  inhibitory: 16,
  n_counts: 17,
  n_genes_by_counts: 18,
  total_counts: 19,
};

const els = {
  canvas: document.querySelector("#atlasCanvas"),
  tooltip: document.querySelector("#tooltip"),
  legend: document.querySelector("#legend"),
  bars: document.querySelector("#bars"),
  colorBy: document.querySelector("#colorBy"),
  pointSize: document.querySelector("#pointSize"),
  pointSizeValue: document.querySelector("#pointSizeValue"),
  projectionControls: document.querySelector("#projectionControls"),
  resetFilters: document.querySelector("#resetFilters"),
  plotTitle: document.querySelector("#plotTitle"),
  datasetSource: document.querySelector("#datasetSource"),
  geneSearch: document.querySelector("#geneSearch"),
  geneTable: document.querySelector("#geneTable"),
  qcGrid: document.querySelector("#qcGrid"),
  downloadPng: document.querySelector("#downloadPng"),
  markerChips: document.querySelectorAll(".marker-chip"),
  heroCanvas: document.querySelector("#heroCanvas"),
  stageCanvases: document.querySelectorAll(".stage-card canvas"),
  statCells: document.querySelector("#statCells"),
  statGenes: document.querySelector("#statGenes"),
  statLibraries: document.querySelector("#statLibraries"),
  statVisible: document.querySelector("#statVisible"),
  countLabel: document.querySelector("#countLabel"),
};

const state = {
  data: null,
  projection: "spatial",
  colorBy: "finer_cell_types",
  selectedCodes: new Set(),
  bounds: new Map(),
  libraryBounds: null,
  screenX: new Float32Array(),
  screenY: new Float32Array(),
  visible: new Uint8Array(),
  width: 0,
  height: 0,
};

const ctx = els.canvas.getContext("2d", { alpha: true });
const fmt = new Intl.NumberFormat("en-US");

init();

async function init() {
  const res = await fetch(DATA_URL);
  if (!res.ok) {
    throw new Error(`Could not load ${DATA_URL}`);
  }
  state.data = await res.json();
  state.screenX = new Float32Array(state.data.cells.length);
  state.screenY = new Float32Array(state.data.cells.length);
  state.visible = new Uint8Array(state.data.cells.length);

  hydrateSummary();
  bindEvents();
  resizeCanvas();
  drawPreviewCanvases();
  renderAll();
}

function hydrateSummary() {
  const { metadata } = state.data;
  els.statCells.textContent = fmt.format(metadata.n_cells);
  els.statGenes.textContent = fmt.format(metadata.n_genes);
  els.statLibraries.textContent = fmt.format(state.data.annotations.library_id.length);
  els.datasetSource.textContent = metadata.source_file;
  els.pointSizeValue.textContent = Number(els.pointSize.value).toFixed(1);
  renderQcSummary();
}

function bindEvents() {
  window.addEventListener("resize", () => {
    resizeCanvas();
    drawPreviewCanvases();
  });

  els.projectionControls.addEventListener("click", (event) => {
    const button = event.target.closest("button[data-projection]");
    if (!button) return;
    state.projection = button.dataset.projection;
    els.projectionControls.querySelectorAll("button").forEach((item) => {
      item.classList.toggle("active", item === button);
    });
    renderAll();
  });

  els.colorBy.addEventListener("change", () => {
    state.colorBy = els.colorBy.value;
    state.selectedCodes.clear();
    renderAll();
  });

  els.pointSize.addEventListener("input", () => {
    els.pointSizeValue.textContent = Number(els.pointSize.value).toFixed(1);
    drawPlot();
  });

  els.resetFilters.addEventListener("click", () => {
    state.selectedCodes.clear();
    renderAll();
  });

  els.geneSearch.addEventListener("input", renderGeneTable);
  els.markerChips.forEach((chip) => {
    chip.addEventListener("click", () => {
      els.geneSearch.value = chip.dataset.gene || chip.textContent.trim();
      renderGeneTable();
      document.querySelector("#explorer").scrollIntoView({ block: "start" });
      els.geneSearch.focus({ preventScroll: true });
    });
  });
  els.downloadPng.addEventListener("click", downloadCanvas);
  els.canvas.addEventListener("mousemove", showTooltip);
  els.canvas.addEventListener("mouseleave", () => {
    els.tooltip.hidden = true;
  });
}

function resizeCanvas() {
  const rect = els.canvas.parentElement.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  state.width = Math.max(320, Math.floor(rect.width));
  state.height = Math.max(240, Math.floor(rect.height));
  els.canvas.width = Math.floor(state.width * dpr);
  els.canvas.height = Math.floor(state.height * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  if (state.data) drawPlot();
}

function renderAll() {
  els.plotTitle.textContent = projectionMap[state.projection].label;
  els.countLabel.textContent = categoryLabels[state.colorBy];
  renderLegend();
  renderBars();
  renderGeneTable();
  drawPlot();
}

function drawPreviewCanvases() {
  if (!state.data) return;
  if (els.heroCanvas) {
    drawStaticPreview(els.heroCanvas, "spatial", "finer_cell_types", 2.1, true);
  }
  els.stageCanvases.forEach((canvas) => {
    const card = canvas.closest(".stage-card");
    const live = card?.classList.contains("live");
    drawStaticPreview(canvas, live ? "spatial" : "slices", live ? "finer_cell_types" : "cell_types", live ? 1.35 : 1.6, !live);
  });
}

function drawStaticPreview(canvas, projection, colorBy, radius, muted) {
  const rect = canvas.getBoundingClientRect();
  const width = Math.max(220, Math.floor(rect.width || 420));
  const height = Math.max(140, Math.floor(rect.height || 240));
  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.floor(width * dpr);
  canvas.height = Math.floor(height * dpr);
  const previewCtx = canvas.getContext("2d");
  previewCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
  previewCtx.clearRect(0, 0, width, height);
  previewCtx.fillStyle = "#ffffff";
  previewCtx.fillRect(0, 0, width, height);

  if (muted) {
    previewCtx.strokeStyle = "#dce5df";
    previewCtx.lineWidth = 1;
    for (let x = 18; x < width; x += 42) {
      previewCtx.beginPath();
      previewCtx.moveTo(x, 12);
      previewCtx.lineTo(x, height - 12);
      previewCtx.stroke();
    }
    for (let y = 18; y < height; y += 36) {
      previewCtx.beginPath();
      previewCtx.moveTo(12, y);
      previewCtx.lineTo(width - 12, y);
      previewCtx.stroke();
    }
  }

  const bounds = computeBounds(projection);
  const categories = state.data.annotations[colorBy] || [];
  const pad = muted ? 16 : 22;
  const usableW = Math.max(1, width - pad * 2);
  const usableH = Math.max(1, height - pad * 2);
  const xSpan = bounds.maxX - bounds.minX || 1;
  const ySpan = bounds.maxY - bounds.minY || 1;
  const scale = Math.min(usableW / xSpan, usableH / ySpan);
  const plotW = xSpan * scale;
  const plotH = ySpan * scale;
  const offsetX = pad + (usableW - plotW) / 2;
  const offsetY = pad + (usableH - plotH) / 2;
  const step = Math.max(1, Math.floor(state.data.cells.length / 14000));

  for (let i = 0; i < state.data.cells.length; i += step) {
    const cell = state.data.cells[i];
    const point = coordinateFor(cell, projection);
    if (!point) continue;
    const [xVal, yVal] = point;
    if (!Number.isFinite(xVal) || !Number.isFinite(yVal)) continue;
    const x = offsetX + (xVal - bounds.minX) * scale;
    const y = height - offsetY - (yVal - bounds.minY) * scale;
    const code = cell[schema[colorBy]];
    const color = categories[code]?.color || "#64748b";
    previewCtx.globalAlpha = muted ? 0.2 : 0.78;
    previewCtx.fillStyle = muted ? "#788590" : color;
    previewCtx.beginPath();
    previewCtx.arc(x, y, radius, 0, Math.PI * 2);
    previewCtx.fill();
  }
  previewCtx.globalAlpha = 1;
}

function getCategories(field = state.colorBy) {
  return state.data.annotations[field] || [];
}

function codeFor(cell, field = state.colorBy) {
  return cell[schema[field]];
}

function isVisible(cell) {
  if (state.selectedCodes.size === 0) return true;
  return state.selectedCodes.has(codeFor(cell));
}

function computeBounds(projection) {
  if (state.bounds.has(projection)) return state.bounds.get(projection);
  if (projection === "slices") {
    ensureLibraryBounds();
    const libraryCount = state.libraryBounds.length || 1;
    const bounds = { minX: 0, maxX: libraryCount * 1.18 - 0.18, minY: 0, maxY: 1 };
    state.bounds.set(projection, bounds);
    return bounds;
  }

  const map = projectionMap[projection];
  let minX = Infinity;
  let maxX = -Infinity;
  let minY = Infinity;
  let maxY = -Infinity;

  for (const cell of state.data.cells) {
    const x = cell[map.x];
    const y = cell[map.y];
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    if (x < minX) minX = x;
    if (x > maxX) maxX = x;
    if (y < minY) minY = y;
    if (y > maxY) maxY = y;
  }

  const bounds = Number.isFinite(minX)
    ? { minX, maxX, minY, maxY }
    : { minX: 0, maxX: 1, minY: 0, maxY: 1 };
  state.bounds.set(projection, bounds);
  return bounds;
}

function ensureLibraryBounds() {
  if (state.libraryBounds) return;
  const libraries = state.data.annotations.library_id || [];
  state.libraryBounds = libraries.map(() => ({
    minX: Infinity,
    maxX: -Infinity,
    minY: Infinity,
    maxY: -Infinity,
  }));

  for (const cell of state.data.cells) {
    const code = cell[schema.library_id];
    const bounds = state.libraryBounds[code];
    if (!bounds) continue;
    const x = cell[0];
    const y = cell[1];
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    if (x < bounds.minX) bounds.minX = x;
    if (x > bounds.maxX) bounds.maxX = x;
    if (y < bounds.minY) bounds.minY = y;
    if (y > bounds.maxY) bounds.maxY = y;
  }
}

function coordinateFor(cell, projection = state.projection) {
  if (projection === "slices") {
    ensureLibraryBounds();
    const code = cell[schema.library_id];
    const bounds = state.libraryBounds[code];
    if (!bounds) return null;
    const xSpan = bounds.maxX - bounds.minX || 1;
    const ySpan = bounds.maxY - bounds.minY || 1;
    return [
      (cell[0] - bounds.minX) / xSpan + code * 1.18,
      (cell[1] - bounds.minY) / ySpan,
    ];
  }

  const map = projectionMap[projection];
  if (!map) return null;
  return [cell[map.x], cell[map.y]];
}

function projectPoint(cell, bounds) {
  const pad = 26;
  const usableW = Math.max(1, state.width - pad * 2);
  const usableH = Math.max(1, state.height - pad * 2);
  const xSpan = bounds.maxX - bounds.minX || 1;
  const ySpan = bounds.maxY - bounds.minY || 1;
  const scale = Math.min(usableW / xSpan, usableH / ySpan);
  const plotW = xSpan * scale;
  const plotH = ySpan * scale;
  const offsetX = pad + (usableW - plotW) / 2;
  const offsetY = pad + (usableH - plotH) / 2;
  const point = coordinateFor(cell);
  const x = offsetX + (point[0] - bounds.minX) * scale;
  const y = state.height - offsetY - (point[1] - bounds.minY) * scale;
  return [x, y];
}

function drawPlot() {
  const cells = state.data.cells;
  const bounds = computeBounds(state.projection);
  const categories = getCategories();
  const radius = Number(els.pointSize.value);
  const activeCount = state.selectedCodes.size;
  let visibleCount = 0;

  ctx.clearRect(0, 0, state.width, state.height);
  ctx.save();
  ctx.globalCompositeOperation = "source-over";

  for (let i = 0; i < cells.length; i += 1) {
    const cell = cells[i];
    const visible = isVisible(cell);
    const point = coordinateFor(cell);
    const drawable = visible && point && Number.isFinite(point[0]) && Number.isFinite(point[1]);
    state.visible[i] = drawable ? 1 : 0;
    if (!drawable) continue;

    const [x, y] = projectPoint(cell, bounds);
    state.screenX[i] = x;
    state.screenY[i] = y;
    const code = codeFor(cell);
    const category = categories[code] || {};
    ctx.fillStyle = category.color || "#64748b";
    ctx.globalAlpha = activeCount ? 0.92 : 0.78;
    ctx.beginPath();
    ctx.arc(x, y, radius, 0, Math.PI * 2);
    ctx.fill();
    visibleCount += 1;
  }

  ctx.restore();
  els.statVisible.textContent = fmt.format(visibleCount);
}

function renderLegend() {
  els.legend.replaceChildren();
  const fragment = document.createDocumentFragment();
  getCategories().forEach((category, code) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "legend-item";
    const selected = state.selectedCodes.has(code);
    const hasFilter = state.selectedCodes.size > 0;
    button.classList.toggle("active", selected || !hasFilter);
    button.classList.toggle("dimmed", hasFilter && !selected);
    button.innerHTML = `
      <span class="swatch" style="background:${category.color}"></span>
      <span class="legend-name" title="${escapeHtml(category.label)}">${escapeHtml(category.label)}</span>
      <span class="legend-count">${fmt.format(category.count)}</span>
    `;
    button.addEventListener("click", () => {
      if (state.selectedCodes.has(code)) {
        state.selectedCodes.delete(code);
      } else {
        state.selectedCodes.add(code);
      }
      renderAll();
    });
    fragment.append(button);
  });
  els.legend.append(fragment);
}

function renderBars() {
  els.bars.replaceChildren();
  const categories = [...getCategories()].sort((a, b) => b.count - a.count);
  const maxCount = Math.max(...categories.map((item) => item.count), 1);
  const fragment = document.createDocumentFragment();

  categories.forEach((category) => {
    const row = document.createElement("div");
    row.className = "bar-row";
    const width = Math.max(0.5, (category.count / maxCount) * 100);
    row.innerHTML = `
      <span class="bar-name" title="${escapeHtml(category.label)}">${escapeHtml(category.label)}</span>
      <span class="bar-track"><span class="bar-fill" style="width:${width}%;background:${category.color}"></span></span>
      <span class="bar-count">${fmt.format(category.count)}</span>
    `;
    fragment.append(row);
  });

  els.bars.append(fragment);
}

function renderGeneTable() {
  const query = els.geneSearch.value.trim().toLowerCase();
  const genes = state.data.genes
    .filter((gene) => !query || gene.gene.toLowerCase().includes(query))
    .slice(0, 90);
  const fragment = document.createDocumentFragment();

  genes.forEach((gene) => {
    const row = document.createElement("tr");
    row.innerHTML = `
      <td>${escapeHtml(gene.gene)}</td>
      <td>${fmt.format(gene.n_cells)}</td>
      <td>${Number(gene.mean_counts ?? gene.mean).toFixed(2)}</td>
      <td>${Number(gene.pct_dropout_by_counts).toFixed(1)}%</td>
    `;
    fragment.append(row);
  });

  els.geneTable.replaceChildren(fragment);
}

function renderQcSummary() {
  els.qcGrid.replaceChildren();
  const labels = {
    n_counts: "Molecule counts",
    n_genes_by_counts: "Detected genes",
    total_counts: "Total counts",
  };
  const fragment = document.createDocumentFragment();

  Object.entries(state.data.qc || {}).forEach(([key, values]) => {
    const card = document.createElement("article");
    card.className = "qc-card";
    const median = values.Median ?? values.median ?? 0;
    const mean = values.Mean ?? values.mean ?? 0;
    const min = values["Min."] ?? values.min ?? 0;
    const max = values["Max."] ?? values.max ?? 0;
    card.innerHTML = `
      <strong>${escapeHtml(labels[key] || key)}</strong>
      <div class="qc-stats">
        <div><span>Median</span><br>${fmt.format(Math.round(median))}</div>
        <div><span>Mean</span><br>${fmt.format(Math.round(mean))}</div>
        <div><span>Min</span><br>${fmt.format(Math.round(min))}</div>
        <div><span>Max</span><br>${fmt.format(Math.round(max))}</div>
      </div>
    `;
    fragment.append(card);
  });

  els.qcGrid.append(fragment);
}

function showTooltip(event) {
  const rect = els.canvas.getBoundingClientRect();
  const mx = event.clientX - rect.left;
  const my = event.clientY - rect.top;
  let best = -1;
  let bestD = 100;

  for (let i = 0; i < state.data.cells.length; i += 1) {
    if (!state.visible[i]) continue;
    const dx = state.screenX[i] - mx;
    const dy = state.screenY[i] - my;
    const d = dx * dx + dy * dy;
    if (d < bestD) {
      bestD = d;
      best = i;
    }
  }

  if (best < 0) {
    els.tooltip.hidden = true;
    return;
  }

  const cell = state.data.cells[best];
  const categories = getCategories();
  const finer = state.data.annotations.finer_cell_types[cell[schema.finer_cell_types]];
  const current = categories[codeFor(cell)];
  const library = state.data.annotations.library_id[cell[schema.library_id]];
  els.tooltip.innerHTML = `
    <strong>${escapeHtml(finer?.label || "Cell")}</strong>
    ${escapeHtml(categoryLabels[state.colorBy])}: ${escapeHtml(current?.label || "NA")}<br>
    Library: ${escapeHtml(library?.label || "NA")}<br>
    Genes: ${fmt.format(cell[schema.n_genes_by_counts])}<br>
    Counts: ${fmt.format(Math.round(cell[schema.n_counts]))}
  `;
  els.tooltip.style.left = `${Math.min(mx + 14, state.width - 230)}px`;
  els.tooltip.style.top = `${Math.max(10, my - 14)}px`;
  els.tooltip.hidden = false;
}

function downloadCanvas() {
  const link = document.createElement("a");
  link.download = `e12-cerebellum-${state.projection}-${state.colorBy}.png`;
  link.href = els.canvas.toDataURL("image/png");
  link.click();
}

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}
