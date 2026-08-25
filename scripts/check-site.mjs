import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const projectRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const failures = [];
let checkedFiles = 0;

function projectPath(relativePath) {
  const filePath = relativePath.split(/[?#]/, 1)[0];
  return path.join(projectRoot, filePath.replace(/^\/+/, ""));
}

function requireFile(relativePath, context) {
  checkedFiles += 1;
  if (!fs.existsSync(projectPath(relativePath))) {
    failures.push(`${context}: missing ${relativePath}`);
    return false;
  }
  return true;
}

function readJson(relativePath, context) {
  if (!requireFile(relativePath, context)) return null;
  try {
    return JSON.parse(fs.readFileSync(projectPath(relativePath), "utf8"));
  } catch (error) {
    failures.push(`${context}: invalid JSON in ${relativePath} (${error.message})`);
    return null;
  }
}

function validateCountIndex(relativePath, context) {
  const index = readJson(relativePath, context);
  if (!index) return;
  if (!Array.isArray(index.genes) || index.genes.length === 0) {
    failures.push(`${context}: count index has no genes`);
  }
  if (Number.isFinite(index.n_genes) && index.genes.length !== index.n_genes) {
    failures.push(`${context}: n_genes does not match the gene catalog`);
  }
}

const html = fs.readFileSync(projectPath("index.html"), "utf8");
const localReferences = [...html.matchAll(/(?:src|href)="([^"]+)"/g)]
  .map((match) => match[1])
  .filter((reference) => !/^(?:https?:|mailto:|#)/.test(reference))
  .map((reference) => reference.split(/[?#]/, 1)[0]);

for (const reference of new Set(localReferences)) {
  requireFile(reference, "index.html");
}

const barseqManifest = readJson("assets/data/barseq-manifest.json", "BARseq3 manifest");
if (barseqManifest) {
  const ids = new Set();
  for (const dataset of barseqManifest.datasets ?? []) {
    if (ids.has(dataset.id)) failures.push(`BARseq3 manifest: duplicate id ${dataset.id}`);
    ids.add(dataset.id);
    requireFile(dataset.data_url, `BARseq3 ${dataset.id}`);
    if (dataset.count_index_url) {
      validateCountIndex(dataset.count_index_url, `BARseq3 ${dataset.id}`);
    }
  }
  if (!ids.has(barseqManifest.default)) {
    failures.push(`BARseq3 manifest: default dataset ${barseqManifest.default} is not defined`);
  }
}

const scrnaManifest = readJson("assets/data/scrna/manifest.json", "scRNA-seq manifest");
if (scrnaManifest) {
  const ids = new Set();
  for (const dataset of scrnaManifest.datasets ?? []) {
    if (ids.has(dataset.id)) failures.push(`scRNA-seq manifest: duplicate id ${dataset.id}`);
    ids.add(dataset.id);
    const browserData = readJson(dataset.data_url, `scRNA-seq ${dataset.id}`);
    const countIndexUrl = browserData?.metadata?.count_index_url;
    if (!countIndexUrl) {
      failures.push(`scRNA-seq ${dataset.id}: browser data has no count_index_url`);
    } else {
      validateCountIndex(countIndexUrl, `scRNA-seq ${dataset.id}`);
    }
  }
}

if (failures.length > 0) {
  console.error("Site validation failed:\n");
  for (const failure of failures) console.error(`- ${failure}`);
  process.exitCode = 1;
} else {
  console.log(`Site validation passed (${checkedFiles} referenced files checked).`);
}
