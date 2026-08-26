// Data configuration and per-dataset metadata visibility
const SCRNA_MANIFEST_URL = "assets/data/scrna/manifest.json?v=20260825-2";

const SCRNA_SPECIES_GROUPS = [
  {
    id: "mouse",
    label: "Mouse datasets",
    scientificName: "Mus musculus",
    image: "assets/images/mouse.png",
    datasets: ["url", "all-inhib", "combined-alltp08-2026"],
  },
  {
    id: "chicken",
    label: "Chicken datasets",
    scientificName: "Gallus gallus",
    image: "assets/images/chicken.png",
    datasets: ["chicken-excitatory", "chicken-inhibitory"],
  },
];

const SCRNA_VISIBLE_ANNOTATIONS = {
  url: ["cell_type_coarse", "Sample", "branch", "TimePoint"],
  "all-inhib": [
    "orig.ident",
    "DCN",
    "timepoint",
    "final.clusters",
    "collector",
    "ident",
  ],
  "chicken-excitatory": [
    "cellType",
    "sub.cluster",
    "clusters",
    "ident",
    "sample",
    "collector",
    "orig.ident",
    "species",
  ],
  "chicken-inhibitory": [
    "cellType",
    "sub.cluster",
    "ident",
    "sample",
    "collector",
    "orig.ident",
    "species",
  ],
  "combined-alltp08-2026": [
    "Sample",
    "TimePoint",
    "med_2p_branch_new",
    "med_2p_branch_coarse",
    "DCN",
    "subareas",
    "orig.ident",
  ],
};

function visibleScrnaAnnotations(datasetId, annotations) {
  const requested = SCRNA_VISIBLE_ANNOTATIONS[datasetId];
  const available = Object.keys(annotations);
  return requested ? requested.filter((field) => annotations[field]) : available;
}

const SCRNA_COLOR_OVERRIDES = new Map([
  ["extracerebellar-fated", "#111111"],
  ["intp", "#2ca25f"],
  ["inta/lat", "#f28e2b"],
  ["inta", "#e76f9a"],
  ["lat", "#d62728"],
  ["medearly", "#8c564b"],
  ["med early", "#8c564b"],
  ["early medial", "#8c564b"],
  ["medlate", "#377eb8"],
  ["med late", "#377eb8"],
  ["late medial", "#377eb8"],
  ["rl", "#78cbe6"],
  ["vz", "#78cbe6"],
  ["int/lat prog", "#f2c94c"],
  ["int+latprog", "#f2c94c"],
  ["i1", "#f28e2b"],
  ["i2/3", "#2ca25f"],
  ["i2", "#2ca25f"],
  ["i3", "#2ca25f"],
  ["other", "#d9dedb"],
  ["others", "#d9dedb"],
]);

function applyScrnaColorOverrides(annotations) {
  for (const rows of Object.values(annotations || {})) {
    for (const row of rows) {
      const key = String(row.label).trim().toLowerCase().replace(/_/g, " ").replace(/\s+/g, " ");
      row.color = SCRNA_COLOR_OVERRIDES.get(key) || row.color;
    }
  }
}

// DOM references and explorer state
const scrnaEls = {
  dataset: document.querySelector("#scrnaDataset"),
  colorBy: document.querySelector("#scrnaColorBy"),
  gene: document.querySelector("#scrnaGene"),
  geneResults: document.querySelector("#scrnaGeneResults"),
  embeddingGrid: document.querySelector("#scrnaEmbeddingGrid"),
  explorer: document.querySelector("#scrnaExplorer"),
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
  countIndex: null,
  dataCache: new Map(),
  geneLookup: new Map(),
  shardCache: new Map(),
  activeGene: null,
  schema: {},
  embedding: "umap",
  colorBy: null,
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
    for (const dataset of scrnaState.manifest.datasets) scrnaEls.dataset.add(new Option(dataset.title, dataset.id));
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
  scrnaEls.gene.addEventListener("input", () => renderGeneResults(scrnaEls.gene.value));
  scrnaEls.gene.addEventListener("focus", () => renderGeneResults(scrnaEls.gene.value));
  scrnaEls.gene.addEventListener("keydown", (event) => {
    if (event.key === "Enter") {
      const exact = scrnaState.geneLookup.get(scrnaEls.gene.value.trim().toLowerCase());
      const first = scrnaEls.geneResults.querySelector("button[data-gene]")?.dataset.gene;
      loadSparseGene(exact === undefined ? (first || scrnaEls.gene.value) : scrnaEls.gene.value);
    }
  });
  scrnaEls.reset.addEventListener("click", () => {
    scrnaState.selected.clear();
    scrnaState.activeGene = null;
    scrnaEls.gene.value = "";
    renderScrna();
  });
  scrnaEls.download.addEventListener("click", () => {
    const link = document.createElement("a");
    const color = scrnaState.activeGene?.gene || scrnaState.colorBy;
    link.download = `${scrnaEls.dataset.value}-${color}-${scrnaState.embedding}.png`;
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
  scrnaEls.explorer.setAttribute("aria-busy", "true");
  scrnaEls.title.textContent = `Loading ${entry.title}`;
  scrnaEls.stats.textContent = "Reading compact browser data…";
  scrnaState.data = scrnaState.dataCache.get(id);
  if (!scrnaState.data) {
    const response = await fetch(entry.data_url);
    if (!response.ok) throw new Error(`Could not load ${entry.data_url}`);
    scrnaState.data = await response.json();
    applyScrnaColorOverrides(scrnaState.data.annotations);
    scrnaState.dataCache.set(id, scrnaState.data);
  }
  const countResponse = await fetch(scrnaState.data.metadata.count_index_url);
  if (!countResponse.ok) throw new Error("Sparse count index is unavailable");
  scrnaState.countIndex = await countResponse.json();
  scrnaState.geneLookup = new Map(scrnaState.countIndex.genes.map((gene, index) => [gene.toLowerCase(), index]));
  scrnaState.shardCache.clear();
  scrnaState.activeGene = null;
  scrnaState.embedding = scrnaState.data.metadata.embeddings.includes("umap") ? "umap" : scrnaState.data.metadata.embeddings[0];
  scrnaState.schema = Object.fromEntries(scrnaState.data.schema.map((name, index) => [name, index]));
  scrnaState.screenX = new Float32Array(scrnaState.data.cells.length);
  scrnaState.screenY = new Float32Array(scrnaState.data.cells.length);
  scrnaState.selected.clear();
  scrnaEls.colorBy.replaceChildren();
  const annotationFields = visibleScrnaAnnotations(id, scrnaState.data.annotations);
  for (const field of annotationFields) scrnaEls.colorBy.add(new Option(humanizeScrnaField(field), field));
  scrnaState.colorBy = scrnaEls.colorBy.value;
  scrnaEls.gene.value = "";
  scrnaEls.geneResults.replaceChildren();
  scrnaEls.geneResults.hidden = true;
  scrnaEls.title.textContent = scrnaState.data.metadata.title;
  scrnaEls.stats.textContent = `${scrnaFmt.format(scrnaState.data.metadata.n_cells)} cells · ${scrnaFmt.format(scrnaState.countIndex.n_genes)} genes · ${scrnaState.data.metadata.source_file}`;
  await renderEmbeddingCards();
  resizeScrnaCanvas();
  renderScrna();
  scrnaEls.explorer.setAttribute("aria-busy", "false");
}

function renderGeneResults(query) {
  if (!scrnaState.countIndex) return;
  const normalized = query.trim().toLowerCase();
  const matches = scrnaState.countIndex.genes.filter((gene) => !normalized || gene.toLowerCase().includes(normalized));
  const visible = matches.slice(0, 100);
  scrnaEls.geneResults.replaceChildren();
  const summary = document.createElement("div");
  summary.className = "gene-search-summary";
  summary.textContent = normalized
    ? `${scrnaFmt.format(matches.length)} matching genes${matches.length > 100 ? " · refine to narrow" : ""}`
    : `${scrnaFmt.format(matches.length)} genes · type to search`;
  scrnaEls.geneResults.append(summary);
  for (const gene of visible) {
    const button = document.createElement("button");
    button.type = "button";
    button.dataset.gene = gene;
    button.textContent = gene;
    button.addEventListener("click", () => loadSparseGene(gene));
    scrnaEls.geneResults.append(button);
  }
  scrnaEls.geneResults.hidden = false;
}

function humanizeScrnaField(field) {
  return field.replace(/[._]+/g, " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
}

function embeddingIndices(name = scrnaState.embedding, schema = scrnaState.schema) {
  return name === "umap"
    ? [schema.umap_x, schema.umap_y]
    : [schema[`embedding_${name}_x`], schema[`embedding_${name}_y`]];
}

async function renderEmbeddingCards() {
  scrnaEls.embeddingGrid.replaceChildren();
  const datasetsById = new Map(scrnaState.manifest.datasets.map((dataset) => [dataset.id, dataset]));
  for (const species of SCRNA_SPECIES_GROUPS) {
    const row = document.createElement("section");
    row.className = "scrna-species-row";
    row.setAttribute("aria-labelledby", `scrna-${species.id}-heading`);

    const identity = document.createElement("div");
    identity.className = "scrna-species-identity";
    const speciesImage = document.createElement("img");
    speciesImage.src = species.image;
    speciesImage.alt = "";
    speciesImage.loading = "lazy";
    const speciesCopy = document.createElement("div");
    const heading = document.createElement("h3");
    heading.id = `scrna-${species.id}-heading`;
    heading.textContent = species.label;
    const scientificName = document.createElement("p");
    scientificName.textContent = species.scientificName;
    speciesCopy.append(heading, scientificName);
    identity.append(speciesImage, speciesCopy);

    const cards = document.createElement("div");
    cards.className = "scrna-species-cards";
    row.append(identity, cards);
    scrnaEls.embeddingGrid.append(row);

    for (const datasetId of species.datasets) {
      const dataset = datasetsById.get(datasetId);
      if (!dataset) continue;
      let data = scrnaState.dataCache.get(dataset.id);
      if (!data) {
        const response = await fetch(dataset.data_url);
        if (!response.ok) continue;
        data = await response.json();
        applyScrnaColorOverrides(data.annotations);
        scrnaState.dataCache.set(dataset.id, data);
      }
      if (!data.metadata.embeddings.includes("umap")) continue;
      const button = document.createElement("button");
      button.type = "button";
      button.className = "scrna-embedding-card";
      const active = dataset.id === scrnaEls.dataset.value;
      button.classList.toggle("active", active);
      button.setAttribute("aria-pressed", String(active));
      const canvas = document.createElement("canvas");
      canvas.setAttribute("aria-label", `${dataset.title} UMAP preview`);
      const label = document.createElement("span");
      label.textContent = dataset.title;
      button.append(canvas, label);
      button.addEventListener("click", async () => {
        scrnaEls.dataset.value = dataset.id;
        await loadScrnaDataset(dataset.id);
        scrnaEls.explorer.scrollIntoView({ block: "start" });
      });
      cards.append(button);
      drawEmbeddingPreview(canvas, data, dataset.id);
    }
  }
}

function boundsForEmbedding(name = scrnaState.embedding, data = scrnaState.data, schema = scrnaState.schema) {
  const [xi, yi] = embeddingIndices(name, schema);
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (const cell of data.cells) {
    const x = cell[xi], y = cell[yi];
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    minX = Math.min(minX, x); maxX = Math.max(maxX, x);
    minY = Math.min(minY, y); maxY = Math.max(maxY, y);
  }
  return { minX, maxX, minY, maxY };
}

function drawEmbeddingPreview(canvas, data, datasetId) {
  const name = "umap";
  const schema = Object.fromEntries(data.schema.map((field, index) => [field, index]));
  const width = 260, height = 150, dpr = window.devicePixelRatio || 1;
  canvas.width = width * dpr; canvas.height = height * dpr;
  const context = canvas.getContext("2d");
  context.setTransform(dpr, 0, 0, dpr, 0, 0);
  context.fillStyle = "#fff"; context.fillRect(0, 0, width, height);
  const bounds = boundsForEmbedding(name, data, schema);
  const [xi, yi] = embeddingIndices(name, schema);
  const colorBy = visibleScrnaAnnotations(datasetId, data.annotations)[0];
  const ci = schema[colorBy];
  const categories = data.annotations[colorBy];
  const xSpan = bounds.maxX - bounds.minX || 1, ySpan = bounds.maxY - bounds.minY || 1;
  const scale = Math.min(230 / xSpan, 120 / ySpan);
  const step = Math.max(1, Math.floor(data.cells.length / 12000));
  context.globalAlpha = 0.7;
  for (let i = 0; i < data.cells.length; i += step) {
    const cell = data.cells[i];
    const x = 15 + (230 - xSpan * scale) / 2 + (cell[xi] - bounds.minX) * scale;
    const y = height - 15 - (120 - ySpan * scale) / 2 - (cell[yi] - bounds.minY) * scale;
    context.fillStyle = categories[cell[ci]]?.color || "#64748b";
    context.fillRect(x, y, 1.5, 1.5);
  }
  context.globalAlpha = 1;
}

async function loadSparseGene(requestedGene) {
  const geneIndex = scrnaState.geneLookup.get(requestedGene.trim().toLowerCase());
  if (geneIndex === undefined) {
    if (!requestedGene.trim()) { scrnaState.activeGene = null; renderScrna(); }
    return;
  }
  const gene = scrnaState.countIndex.genes[geneIndex];
  scrnaEls.stats.textContent = `Loading ${gene} sparse counts…`;
  const shard = scrnaState.countIndex.shards.find((item) => geneIndex >= item.start && geneIndex < item.start + item.count);
  const base = scrnaState.data.metadata.count_index_url.replace(/index\.json$/, "");
  const cacheKey = `${scrnaEls.dataset.value}/${shard.file}`;
  let buffer = scrnaState.shardCache.get(cacheKey);
  if (!buffer) {
    const response = await fetch(base + shard.file);
    if (!response.ok) throw new Error(`Could not load count shard for ${gene}`);
    const compressed = await response.arrayBuffer();
    const stream = new Blob([compressed]).stream().pipeThrough(new DecompressionStream("gzip"));
    buffer = await new Response(stream).arrayBuffer();
    scrnaState.shardCache.set(cacheKey, buffer);
  }
  const values = decodeGeneRecord(buffer, geneIndex - shard.start, scrnaState.countIndex.n_cells);
  const nonzero = Array.from(values).filter((value) => value > 0).sort((a, b) => a - b);
  const max = nonzero[Math.min(nonzero.length - 1, Math.floor(nonzero.length * 0.99))] || 1;
  scrnaState.activeGene = { gene, values, max };
  scrnaEls.gene.value = gene;
  scrnaEls.geneResults.hidden = true;
  scrnaEls.stats.textContent = `${scrnaFmt.format(scrnaState.data.metadata.n_cells)} cells · ${scrnaFmt.format(scrnaState.countIndex.n_genes)} genes · ${gene} counts`;
  renderScrna();
}

function decodeGeneRecord(buffer, targetRecord, cellCount) {
  const view = new DataView(buffer);
  const recordCount = view.getUint32(0, true);
  if (targetRecord >= recordCount) throw new Error("Gene record is outside its shard");
  let offset = 4;
  for (let record = 0; record < recordCount; record += 1) {
    const nnz = view.getUint32(offset, true); offset += 4;
    if (record === targetRecord) {
      const values = new Float32Array(cellCount);
      const valueOffset = offset + nnz * 4;
      for (let index = 0; index < nnz; index += 1) {
        values[view.getUint32(offset + index * 4, true)] = view.getFloat32(valueOffset + index * 4, true);
      }
      return values;
    }
    offset += nnz * 8;
  }
  throw new Error("Gene record was not found");
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

function renderScrna() { renderScrnaLegend(); drawScrna(); }

function renderScrnaLegend() {
  if (scrnaState.activeGene) {
    const gene = scrnaState.activeGene;
    scrnaEls.legend.innerHTML = `<div class="gene-legend"><strong>${gene.gene}</strong><div class="gene-gradient"></div><div><span>0</span><span>${gene.max.toFixed(2)}</span></div><p>Sparse counts, colored to the 99th percentile.</p></div>`;
    return;
  }
  const categories = scrnaState.data.annotations[scrnaState.colorBy] || [];
  scrnaEls.legend.replaceChildren();
  categories.forEach((category, code) => {
    const button = document.createElement("button");
    button.type = "button"; button.className = "legend-item";
    const selected = scrnaState.selected.has(code);
    const hasFilter = scrnaState.selected.size > 0;
    button.classList.toggle("muted", hasFilter && !selected);
    button.setAttribute("aria-pressed", String(selected || !hasFilter));
    button.innerHTML = `<span class="swatch" style="background:${category.color}"></span><span>${category.label}</span><strong>${scrnaFmt.format(category.count)}</strong>`;
    button.addEventListener("click", () => {
      if (scrnaState.selected.has(code)) scrnaState.selected.delete(code); else scrnaState.selected.add(code);
      renderScrna();
    });
    scrnaEls.legend.append(button);
  });
}

function drawScrna() {
  if (!scrnaState.data || !scrnaState.width) return;
  const { width, height } = scrnaState;
  scrnaCtx.clearRect(0, 0, width, height); scrnaCtx.fillStyle = "#fff"; scrnaCtx.fillRect(0, 0, width, height);
  const bounds = boundsForEmbedding();
  const pad = 28, xSpan = bounds.maxX - bounds.minX || 1, ySpan = bounds.maxY - bounds.minY || 1;
  const scale = Math.min((width - pad * 2) / xSpan, (height - pad * 2) / ySpan);
  const offsetX = pad + (width - pad * 2 - xSpan * scale) / 2;
  const offsetY = pad + (height - pad * 2 - ySpan * scale) / 2;
  const [xi, yi] = embeddingIndices();
  const ci = scrnaState.schema[scrnaState.colorBy];
  const categories = scrnaState.data.annotations[scrnaState.colorBy];
  const radius = scrnaState.data.cells.length > 45000 ? 1.15 : 1.45;
  scrnaCtx.globalAlpha = 0.76;
  for (let i = 0; i < scrnaState.data.cells.length; i += 1) {
    const cell = scrnaState.data.cells[i], code = cell[ci];
    const x = offsetX + (cell[xi] - bounds.minX) * scale;
    const y = height - offsetY - (cell[yi] - bounds.minY) * scale;
    scrnaState.screenX[i] = x; scrnaState.screenY[i] = y;
    if (scrnaState.selected.size > 0 && !scrnaState.selected.has(code)) continue;
    scrnaCtx.fillStyle = scrnaState.activeGene ? expressionColor(scrnaState.activeGene.values[i], scrnaState.activeGene.max) : (categories[code]?.color || "#64748b");
    scrnaCtx.beginPath(); scrnaCtx.arc(x, y, radius, 0, Math.PI * 2); scrnaCtx.fill();
  }
  scrnaCtx.globalAlpha = 1;
}

function expressionColor(value, max) {
  if (value <= 0) return "hsl(220 18% 91%)";
  const ratio = Math.max(0, Math.min(1, value / (max || 1)));
  return `hsl(${220 - ratio * 210} 82% ${92 - ratio * 48}%)`;
}

function showScrnaTooltip(event) {
  if (!scrnaState.data) return;
  const rect = scrnaEls.canvas.getBoundingClientRect(), x = event.clientX - rect.left, y = event.clientY - rect.top;
  let nearest = -1, distance = 64;
  for (let i = 0; i < scrnaState.data.cells.length; i += 1) {
    const dx = scrnaState.screenX[i] - x, dy = scrnaState.screenY[i] - y, candidate = dx * dx + dy * dy;
    if (candidate < distance) { distance = candidate; nearest = i; }
  }
  if (nearest < 0) { scrnaEls.tooltip.hidden = true; return; }
  if (scrnaState.activeGene) scrnaEls.tooltip.textContent = `${scrnaState.activeGene.gene}: ${scrnaState.activeGene.values[nearest].toFixed(2)} counts`;
  else {
    const code = scrnaState.data.cells[nearest][scrnaState.schema[scrnaState.colorBy]];
    const category = scrnaState.data.annotations[scrnaState.colorBy][code];
    scrnaEls.tooltip.textContent = `${humanizeScrnaField(scrnaState.colorBy)}: ${category?.label || "Not available"}`;
  }
  scrnaEls.tooltip.style.left = `${x + 12}px`; scrnaEls.tooltip.style.top = `${y + 12}px`; scrnaEls.tooltip.hidden = false;
}
