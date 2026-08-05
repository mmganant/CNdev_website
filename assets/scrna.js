const SCRNA_MANIFEST_URL = "assets/data/scrna/manifest.json";

const scrnaEls = {
  dataset: document.querySelector("#scrnaDataset"),
  colorBy: document.querySelector("#scrnaColorBy"),
  gene: document.querySelector("#scrnaGene"),
  reset: document.querySelector("#scrnaReset"),
  download: document.querySelector("#scrnaDownload"),
  title: document.querySelector("#scrnaTitle"),
  stats: document.querySelector("#scrnaStats"),
  canvas: document.querySelector("#singleCellCanvas"),
  tooltip: document.querySelector("#scrnaTooltip"),
  legend: document.querySelector("#scrnaLegend"),
};

const scrnaState = {
  manifest: null,
  data: null,
  schema: {},
  colorBy: null,
  geneIndex: -1,
  selected: new Set(),
  screenX: new Float32Array(),
  screenY: new Float32Array(),
  width: 0,
  height: 0,
};

const scrnaCtx = scrnaEls.canvas?.getContext("2d");
const scrnaFmt = new Intl.NumberFormat("en-US");

if (scrnaCtx) initScrnaExplorer();

async function initScrnaExplorer() {
  try {
    const response = await fetch(SCRNA_MANIFEST_URL);
    if (!response.ok) throw new Error("scRNA-seq manifest is unavailable");
    scrnaState.manifest = await response.json();
    for (const dataset of scrnaState.manifest.datasets) {
      scrnaEls.dataset.add(new Option(dataset.title, dataset.id));
    }
    bindScrnaEvents();
    await loadScrnaDataset(scrnaState.manifest.datasets[0].id);
  } catch (error) {
    scrnaEls.title.textContent = "scRNA-seq data could not be loaded";
    scrnaEls.stats.textContent = error.message;
  }
}

function bindScrnaEvents() {
  scrnaEls.dataset.addEventListener("change", () => loadScrnaDataset(scrnaEls.dataset.value));
  scrnaEls.colorBy.addEventListener("change", () => {
    scrnaState.colorBy = scrnaEls.colorBy.value;
    scrnaState.selected.clear();
    renderScrna();
  });
  scrnaEls.gene.addEventListener("change", () => {
    scrnaState.geneIndex = Number(scrnaEls.gene.value);
    renderScrna();
  });
  scrnaEls.reset.addEventListener("click", () => {
    scrnaState.selected.clear();
    renderScrna();
  });
  scrnaEls.download.addEventListener("click", () => {
    const link = document.createElement("a");
    link.download = `${scrnaEls.dataset.value}-${scrnaState.colorBy}-umap.png`;
    link.href = scrnaEls.canvas.toDataURL("image/png");
    link.click();
  });
  scrnaEls.canvas.addEventListener("mousemove", showScrnaTooltip);
  scrnaEls.canvas.addEventListener("mouseleave", () => { scrnaEls.tooltip.hidden = true; });
  window.addEventListener("resize", resizeScrnaCanvas);
}

async function loadScrnaDataset(id) {
  const entry = scrnaState.manifest.datasets.find((item) => item.id === id);
  if (!entry) return;
  scrnaEls.title.textContent = `Loading ${entry.title}`;
  scrnaEls.stats.textContent = "Reading compact browser data…";
  const response = await fetch(entry.data_url);
  if (!response.ok) throw new Error(`Could not load ${entry.data_url}`);
  scrnaState.data = await response.json();
  scrnaState.schema = Object.fromEntries(scrnaState.data.schema.map((name, index) => [name, index]));
  scrnaState.screenX = new Float32Array(scrnaState.data.cells.length);
  scrnaState.screenY = new Float32Array(scrnaState.data.cells.length);
  scrnaState.selected.clear();
  scrnaEls.colorBy.replaceChildren();
  scrnaEls.gene.replaceChildren(new Option("Annotation colors", "-1"));
  for (const field of Object.keys(scrnaState.data.annotations)) {
    scrnaEls.colorBy.add(new Option(humanizeScrnaField(field), field));
  }
  scrnaState.colorBy = scrnaEls.colorBy.value;
  for (const [index, gene] of (scrnaState.data.genes || []).entries()) {
    scrnaEls.gene.add(new Option(gene.gene, String(index)));
  }
  scrnaState.geneIndex = -1;
  scrnaEls.title.textContent = scrnaState.data.metadata.title;
  scrnaEls.stats.textContent = `${scrnaFmt.format(scrnaState.data.metadata.n_cells)} cells · ${scrnaState.data.metadata.source_file}`;
  resizeScrnaCanvas();
  renderScrna();
}

function humanizeScrnaField(field) {
  return field.replace(/[._]+/g, " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
}

function resizeScrnaCanvas() {
  if (!scrnaState.data) return;
  const parent = scrnaEls.canvas.parentElement;
  const dpr = window.devicePixelRatio || 1;
  scrnaState.width = Math.max(320, Math.floor(parent.clientWidth));
  scrnaState.height = Math.max(360, Math.floor(parent.clientHeight));
  scrnaEls.canvas.width = Math.floor(scrnaState.width * dpr);
  scrnaEls.canvas.height = Math.floor(scrnaState.height * dpr);
  scrnaCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
  drawScrna();
}

function renderScrna() {
  renderScrnaLegend();
  drawScrna();
}

function renderScrnaLegend() {
  if (scrnaState.geneIndex >= 0) {
    const gene = scrnaState.data.genes[scrnaState.geneIndex];
    scrnaEls.legend.innerHTML = `<div class="gene-legend"><strong>${gene.gene}</strong><div class="gene-gradient"></div><div><span>0</span><span>${gene.max.toFixed(2)}</span></div><p>Normalized expression, capped at the 99th percentile.</p></div>`;
    return;
  }
  const categories = scrnaState.data.annotations[scrnaState.colorBy] || [];
  scrnaEls.legend.replaceChildren();
  categories.forEach((category, code) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "legend-item";
    button.classList.toggle("muted", scrnaState.selected.size > 0 && !scrnaState.selected.has(code));
    button.innerHTML = `<span class="swatch" style="background:${category.color}"></span><span>${category.label}</span><strong>${scrnaFmt.format(category.count)}</strong>`;
    button.addEventListener("click", () => {
      if (scrnaState.selected.has(code)) scrnaState.selected.delete(code);
      else scrnaState.selected.add(code);
      renderScrna();
    });
    scrnaEls.legend.append(button);
  });
}

function scrnaBounds() {
  const xi = scrnaState.schema.umap_x;
  const yi = scrnaState.schema.umap_y;
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (const cell of scrnaState.data.cells) {
    const x = cell[xi], y = cell[yi];
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    minX = Math.min(minX, x); maxX = Math.max(maxX, x);
    minY = Math.min(minY, y); maxY = Math.max(maxY, y);
  }
  return { minX, maxX, minY, maxY };
}

function drawScrna() {
  if (!scrnaState.data || !scrnaState.width) return;
  const { width, height } = scrnaState;
  scrnaCtx.clearRect(0, 0, width, height);
  scrnaCtx.fillStyle = "#fff";
  scrnaCtx.fillRect(0, 0, width, height);
  const bounds = scrnaBounds();
  const pad = 28;
  const xSpan = bounds.maxX - bounds.minX || 1;
  const ySpan = bounds.maxY - bounds.minY || 1;
  const scale = Math.min((width - pad * 2) / xSpan, (height - pad * 2) / ySpan);
  const offsetX = pad + (width - pad * 2 - xSpan * scale) / 2;
  const offsetY = pad + (height - pad * 2 - ySpan * scale) / 2;
  const xi = scrnaState.schema.umap_x;
  const yi = scrnaState.schema.umap_y;
  const ci = scrnaState.schema[scrnaState.colorBy];
  const categories = scrnaState.data.annotations[scrnaState.colorBy];
  const gene = scrnaState.geneIndex >= 0 ? scrnaState.data.genes[scrnaState.geneIndex] : null;
  const radius = scrnaState.data.cells.length > 45000 ? 1.15 : 1.45;
  scrnaCtx.globalAlpha = 0.76;
  for (let i = 0; i < scrnaState.data.cells.length; i += 1) {
    const cell = scrnaState.data.cells[i];
    const code = cell[ci];
    const visible = scrnaState.selected.size === 0 || scrnaState.selected.has(code);
    const x = offsetX + (cell[xi] - bounds.minX) * scale;
    const y = height - offsetY - (cell[yi] - bounds.minY) * scale;
    scrnaState.screenX[i] = x;
    scrnaState.screenY[i] = y;
    if (!visible) continue;
    scrnaCtx.fillStyle = gene ? expressionColor(gene.values[i], gene.max) : (categories[code]?.color || "#64748b");
    scrnaCtx.beginPath();
    scrnaCtx.arc(x, y, radius, 0, Math.PI * 2);
    scrnaCtx.fill();
  }
  scrnaCtx.globalAlpha = 1;
}

function expressionColor(value, max) {
  const ratio = Math.max(0, Math.min(1, value / (max || 1)));
  const hue = 220 - ratio * 210;
  const lightness = 92 - ratio * 48;
  return `hsl(${hue} 82% ${lightness}%)`;
}

function showScrnaTooltip(event) {
  if (!scrnaState.data) return;
  const rect = scrnaEls.canvas.getBoundingClientRect();
  const x = event.clientX - rect.left;
  const y = event.clientY - rect.top;
  let nearest = -1;
  let distance = 64;
  for (let i = 0; i < scrnaState.data.cells.length; i += 1) {
    const dx = scrnaState.screenX[i] - x;
    const dy = scrnaState.screenY[i] - y;
    const candidate = dx * dx + dy * dy;
    if (candidate < distance) { distance = candidate; nearest = i; }
  }
  if (nearest < 0) { scrnaEls.tooltip.hidden = true; return; }
  const cell = scrnaState.data.cells[nearest];
  const code = cell[scrnaState.schema[scrnaState.colorBy]];
  const category = scrnaState.data.annotations[scrnaState.colorBy][code];
  if (scrnaState.geneIndex >= 0) {
    const gene = scrnaState.data.genes[scrnaState.geneIndex];
    scrnaEls.tooltip.textContent = `${gene.gene}: ${gene.values[nearest].toFixed(3)}`;
  } else {
    scrnaEls.tooltip.textContent = `${humanizeScrnaField(scrnaState.colorBy)}: ${category?.label || "Not available"}`;
  }
  scrnaEls.tooltip.style.left = `${x + 12}px`;
  scrnaEls.tooltip.style.top = `${y + 12}px`;
  scrnaEls.tooltip.hidden = false;
}
