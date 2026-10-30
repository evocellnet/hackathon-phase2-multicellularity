// Build the experiment slide deck: report/slides/multicellularity_experiment.pptx
// Run from anywhere:  node scripts/report/slides/build_deck.js
// Requires: npm install pptxgenjs react react-dom react-icons sharp (inside scripts/report/slides/)
// Diagrams: run make_diagrams.py first (writes report/slides/img/).
// Every number on the slides comes from the source named in that slide's speaker notes.

const path = require("path");
// let the skill's apply_theme.js find pptxgenjs/jszip installed next to this script
process.env.NODE_PATH = path.join(__dirname, "node_modules");
require("module").Module._initPaths();
const fs = require("fs");
const pptxgen = require("pptxgenjs");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const sharp = require("sharp");
const fa = require("react-icons/fa");

const ROOT = path.resolve(__dirname, "../../..");
const IMG = path.join(ROOT, "report/slides/img");
const OUTDIR = path.join(ROOT, "report/slides");
const OUTFILE = path.join(OUTDIR, "multicellularity_experiment.pptx");
const SKILL = "/home/.claude/skills/synced/d9e26939-73cc-4466-a289-fae5e82e2c75_bb964073-749a-4ca3-8415-3c5e37996ced/pptx";
const { applyTheme } = require(path.join(SKILL, "scripts/apply_theme.js"));

const THEME = {
  name: "Multicellularity Teal",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1F2A2A", lt1: "FFFFFF", dk2: "0E3B39", lt2: "E6F1EF",
    accent1: "1B7F79", accent2: "009E73", accent3: "0072B2", accent4: "E4572E", accent5: "E69F00", accent6: "9AA5A3",
    hlink: "0072B2", folHlink: "6C5B7B",
  },
};
const HEX = THEME.colors;

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.333 x 7.5 in
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "From many origins to single events: an AI-run study of multicellularity";
pres.author = "Claude Opus 5.5 (PI agent) for the AI hackathon";
pres.company = "AI hackathon (Beltrao lab, ETH Zurich)";
const C = pres.SchemeColor;

// ------------------------------------------------------------------------------------------------
// helpers
function pngSize(p) {
  const b = fs.readFileSync(p);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}
function img(slide, file, box, alt, name) {
  const p = path.join(IMG, file);
  const s = pngSize(p);
  const ar = s.w / s.h;
  let w = box.w, h = box.w / ar;
  if (h > box.h) { h = box.h; w = h * ar; }
  const x = box.x + (box.align === "left" ? 0 : (box.w - w) / 2);
  const y = box.y + (box.valign === "top" ? 0 : (box.h - h) / 2);
  slide.addImage({ path: p, x, y, w, h, altText: alt || file, objectName: name || file });
  return { x, y, w, h };
}
function txt(slide, text, o) {
  slide.addText(text, Object.assign({ isTextBox: true, margin: 0, fontSize: 15, color: C.text1, valign: "top" }, o));
}
function caption(slide, text, o) {
  txt(slide, text, Object.assign({ fontSize: 11, italic: true, color: C.text2 }, o));
}
function box(slide, x, y, w, h, fill, name, o) {
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, Object.assign({ x, y, w, h, rectRadius: 0.12, fill: { color: fill },
    line: { color: fill, width: 0 }, objectName: name }, o || {}));
}
const ICON = {};
async function renderIcon(name, color) {
  const key = name + "_" + color;
  if (ICON[key]) return ICON[key];
  const comp = fa[name];
  if (!comp) throw new Error("missing icon " + name);
  const svg = ReactDOMServer.renderToStaticMarkup(React.createElement(comp, { color: "#" + color, size: 256 }));
  const buf = await sharp(Buffer.from(svg)).png().toBuffer();
  ICON[key] = "image/png;base64," + buf.toString("base64");
  return ICON[key];
}
function icon(slide, name, x, y, d, circle, fg) {
  // icon in a coloured circle (deck motif)
  const c = circle || C.accent1;
  slide.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color: c }, line: { color: c, width: 0 }, objectName: "icon-circle-" + name });
  const key = name + "_" + (fg || "FFFFFF");
  const pad = d * 0.24;
  slide.addImage({ data: ICON[key], x: x + pad, y: y + pad, w: d - 2 * pad, h: d - 2 * pad, altText: name, objectName: "icon-" + name });
}
function arrow(slide, x1, y1, x2, y2, color, wpt) {
  const flipH = x2 < x1, flipV = y2 < y1;
  slide.addShape(pres.shapes.LINE, { x: Math.min(x1, x2), y: Math.min(y1, y2), w: Math.max(Math.abs(x2 - x1), 0.001),
    h: Math.max(Math.abs(y2 - y1), 0.001), flipH, flipV,
    line: { color: color || C.text2, width: wpt || 2, endArrowType: "triangle" }, objectName: "arrow" });
}
function stat(slide, x, y, w, big, label, color, bigSize) {
  txt(slide, big, { x, y, w, h: 0.75, fontSize: bigSize || 40, bold: true, color: color || C.accent1, valign: "bottom" });
  txt(slide, label, { x, y: y + 0.8, w, h: 0.75, fontSize: 14, color: C.text1 });
}
const VERDICT = {
  supported: { fill: C.accent2, fg: C.background1, label: "SUPPORTED" },
  notsupported: { fill: C.accent4, fg: C.background1, label: "NOT SUPPORTED" },
  null: { fill: C.accent4, fg: C.background1, label: "NULL" },
  overturned: { fill: C.accent4, fg: C.background1, label: "OVERTURNED / INCONCLUSIVE" },
  notsig: { fill: C.accent5, fg: C.text1, label: "NOT SIGNIFICANT" },
  notdetected: { fill: C.accent5, fg: C.text1, label: "NOT DETECTED (LOW POWER)" },
  lead: { fill: C.accent3, fg: C.background1, label: "EXPLORATORY LEAD" },
  nottested: { fill: C.accent6, fg: C.text1, label: "NOT TESTED" },
  partly: { fill: C.accent6, fg: C.text1, label: "PARTLY TESTED" },
  nottestable: { fill: C.accent6, fg: C.text1, label: "NOT TESTABLE HERE" },
};
function badge(slide, x, y, w, kind, h) {
  const v = VERDICT[kind];
  slide.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h: h || 0.42, rectRadius: 0.2, fill: { color: v.fill },
    line: { color: v.fill, width: 0 }, objectName: "verdict-badge" });
  txt(slide, v.label, { x, y, w, h: h || 0.42, fontSize: 13, bold: true, color: v.fg, align: "center", valign: "middle" });
}

// ------------------------------------------------------------------------------------------------
// layouts
const TITLE_Y = 0.32, TITLE_H = 0.85;
pres.defineSlideMaster({
  title: "DARK_TITLE",
  background: { color: C.text2 },
  objects: [
    { image: { path: path.join(IMG, "motif_teal.png"), x: 9.55, y: 1.55, w: 3.0, h: 3.0 } },
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 1.7, w: 8.4, h: 1.9, fontSize: 40, bold: true,
      color: C.background1, valign: "bottom", align: "left" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 3.8, w: 8.4, h: 1.4, fontSize: 20,
      color: C.background2, valign: "top", align: "left" }, text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "SECTION",
  background: { color: C.text2 },
  objects: [
    { image: { path: path.join(IMG, "motif_teal.png"), x: 9.8, y: 2.1, w: 2.6, h: 2.6 } },
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.05, w: 8.6, h: 1.7, fontSize: 40, bold: true,
      color: C.background1, valign: "bottom", align: "left" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 3.9, w: 8.6, h: 1.0, fontSize: 20,
      color: C.background2, valign: "top", align: "left" }, text: "" } },
  ],
  slideNumber: { x: 12.3, y: 6.95, w: 0.6, h: 0.3, fontSize: 10, color: C.background2, align: "right" },
});
pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: C.background1 },
  objects: [
    { image: { path: path.join(IMG, "motif_teal.png"), x: 12.2, y: 0.47, w: 0.55, h: 0.55 } },
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: TITLE_Y, w: 11.4, h: TITLE_H, fontSize: 32, bold: true,
      color: C.text2, valign: "middle", align: "left" }, text: "" } },
  ],
  slideNumber: { x: 12.3, y: 6.98, w: 0.6, h: 0.3, fontSize: 10, color: C.text2, align: "right" },
});
pres.defineSlideMaster({
  title: "DARK_CONTENT",
  background: { color: C.text2 },
  objects: [
    { image: { path: path.join(IMG, "motif_light.png"), x: 12.2, y: 0.47, w: 0.55, h: 0.55 } },
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: TITLE_Y, w: 11.4, h: TITLE_H, fontSize: 32, bold: true,
      color: C.background1, valign: "middle", align: "left" }, text: "" } },
  ],
  slideNumber: { x: 12.3, y: 6.98, w: 0.6, h: 0.3, fontSize: 10, color: C.background2, align: "right" },
});

let SECTION = "";
function section(title) { SECTION = title; pres.addSection({ title }); }
function content(title, notes, dark) {
  const s = pres.addSlide({ masterName: dark ? "DARK_CONTENT" : "CONTENT", sectionTitle: SECTION });
  s.addText(title, { placeholder: "title" });
  if (notes) s.addNotes(notes);
  return s;
}
function divider(letter, title, sub, notes) {
  const s = pres.addSlide({ masterName: "SECTION", sectionTitle: SECTION });
  txt(s, letter, { x: 0.8, y: 1.05, w: 2, h: 0.9, fontSize: 28, bold: true, color: C.accent5, valign: "bottom" });
  s.addText(title, { placeholder: "title" });
  s.addText(sub, { placeholder: "body" });
  if (notes) s.addNotes(notes);
  return s;
}

// ------------------------------------------------------------------------------------------------
const ICON_LIST = [
  ["FaUniversity", "FFFFFF"], ["FaRobot", "FFFFFF"], ["FaUserCheck", "FFFFFF"], ["FaServer", "FFFFFF"],
  ["FaBookOpen", "FFFFFF"], ["FaDna", "FFFFFF"], ["FaClock", "FFFFFF"], ["FaMicrochip", "FFFFFF"],
  ["FaClipboardCheck", "FFFFFF"], ["FaSitemap", "FFFFFF"], ["FaLightbulb", "FFFFFF"], ["FaPenNib", "FFFFFF"],
  ["FaExclamationTriangle", "FFFFFF"], ["FaThumbsUp", "FFFFFF"], ["FaThumbsDown", "FFFFFF"], ["FaSearch", "FFFFFF"],
  ["FaDatabase", "FFFFFF"], ["FaCubes", "FFFFFF"], ["FaLayerGroup", "FFFFFF"], ["FaProjectDiagram", "FFFFFF"],
  ["FaUsers", "FFFFFF"], ["FaCheckDouble", "FFFFFF"], ["FaLock", "FFFFFF"], ["FaBalanceScale", "FFFFFF"],
  ["FaRedo", "FFFFFF"], ["FaImage", "0E3B39"], ["FaGithub", "FFFFFF"], ["FaLink", "FFFFFF"], ["FaShieldAlt", "FFFFFF"],
  ["FaChartBar", "FFFFFF"], ["FaUserTie", "FFFFFF"], ["FaFlask", "FFFFFF"], ["FaCrosshairs", "FFFFFF"],
  ["FaPuzzlePiece", "FFFFFF"], ["FaCodeBranch", "FFFFFF"], ["FaFileAlt", "FFFFFF"], ["FaTachometerAlt", "FFFFFF"],
  ["FaBug", "FFFFFF"], ["FaSkullCrossbones", "FFFFFF"], ["FaHourglassHalf", "FFFFFF"], ["FaUndo", "FFFFFF"],
  ["FaMicroscope", "FFFFFF"], ["FaGlobeEurope", "FFFFFF"], ["FaHandshake", "FFFFFF"], ["FaCompass", "FFFFFF"],
  ["FaUser", "FFFFFF"], ["FaVial", "FFFFFF"], ["FaBullseye", "FFFFFF"], ["FaExchangeAlt", "FFFFFF"],
  ["FaRulerCombined", "FFFFFF"], ["FaMapSigns", "FFFFFF"], ["FaKey", "FFFFFF"], ["FaCogs", "FFFFFF"],
  ["FaArrowRight", "1B7F79"], ["FaCheck", "FFFFFF"], ["FaTimes", "FFFFFF"], ["FaTree", "FFFFFF"], ["FaCut", "FFFFFF"],
];

const SRC = {
  report: "report/REPORT.md",
  supp: "report/REPORT_SUPPLEMENT.md",
  tree: "RESEARCH_TREE.md",
  lit: "LITERATURE.md",
  score: "results/literature_scorecard.tsv",
  hack: "github.com/evocellnet/ai-hackathon README",
  anc: "external/ancestra-mvp-ai/README.md",
};

async function build() {
  for (const [n, c] of ICON_LIST) await renderIcon(n, c);

  // ============================================================================================
  section("A. The hackathon");
  {
    const s = pres.addSlide({ masterName: "DARK_TITLE", sectionTitle: SECTION });
    s.addText("From many origins to single events", { placeholder: "title" });
    s.addText([
      { text: "An AI-run study of the genomic transitions to multicellularity", options: { breakLine: true } },
      { text: "AI hackathon · [dates: ____] · two rounds, about 31 hours of agent work", options: { fontSize: 16 } },
    ], { placeholder: "body" });
    s.addNotes("Title. This deck documents the whole experiment: the hackathon context, the prompts, how the AI team decided on the method, the method itself, the results, every hypothesis considered, and what the experience says about AI-led research. The agent was Claude Opus 5.5 acting as principal investigator (PI) and delegating to Opus 5.5 / Sonnet 5.5 subagents. Round 1 took about 13.6 h, round 2 about 17.4 h (about 31 h in total). The date field is left blank on purpose. Source: report/REPORT.md (version note), RESEARCH_TREE.md.");
  }
  {
    const s = content("A week-long test of agentic bioinformatics",
      "The hackathon (github.com/evocellnet/ai-hackathon) was run by four labs to assess the bioinformatics capabilities of agentic AI. Each project gives Claude Opus 5.5 a 'predict X' problem; it delegates to Opus/Sonnet subagents, builds, benchmarks and writes up a method with minimal human intervention. Humans verify the methodology and benchmarking afterwards. Hackathon projects run on the ETH Euler cluster. Our project is a phase-2 project that ran on a different university cluster, inside a Singularity container in a /scratch folder: 24 CPUs and 150 GB in round 1; 48 CPUs, 500 GB and one A100 GPU in round 2. Blanks for the dates and number of participants are left for the organisers to fill. Source: hackathon README; PROMPT.md; PROMPT_ROUND2.md.");
    const labs = [["Beltrao lab", "ETH Zürich"], ["Colm Ryan lab", "UCD"], ["Emanuel Gonçalves lab", "IST"], ["Christophe Dessimoz lab", "UNIL"]];
    labs.forEach((l, i) => {
      const x = 0.6 + (i % 2) * 2.95, y = 1.45 + Math.floor(i / 2) * 1.55;
      box(s, x, y, 2.75, 1.35, C.background2, "lab-card");
      icon(s, "FaUniversity", x + 0.18, y + 0.35, 0.65);
      txt(s, [{ text: l[0], options: { bold: true, breakLine: true } }, { text: l[1] }], { x: x + 0.95, y: y + 0.25, w: 1.75, h: 0.9, fontSize: 14, valign: "middle" });
    });
    // blanks
    box(s, 0.6, 4.65, 5.7, 1.0, C.background2, "blanks");
    txt(s, [{ text: "[dates: ____]", options: { bold: true, breakLine: true } }, { text: "[__ participants] · 13 projects" }],
      { x: 0.85, y: 4.75, w: 5.3, h: 0.8, fontSize: 16, valign: "middle", color: C.text2 });
    caption(s, "Blanks to be filled by the organisers.", { x: 0.6, y: 5.75, w: 5.7, h: 0.3 });
    // flow
    const steps = [["FaUser", "Human writes a 'predict X' prompt"], ["FaRobot", "Claude Opus 5.5 leads as PI"],
      ["FaUsers", "Delegates to Opus / Sonnet subagents"], ["FaFlask", "Builds, benchmarks, writes up"], ["FaUserCheck", "Humans verify afterwards"]];
    steps.forEach((st, i) => {
      const y = 1.45 + i * 0.95;
      icon(s, st[0], 7.0, y, 0.7, i === 4 ? C.accent4 : C.accent1);
      txt(s, st[1], { x: 7.9, y: y + 0.05, w: 4.6, h: 0.6, fontSize: 16, valign: "middle" });
      if (i < 4) arrow(s, 7.35, y + 0.72, 7.35, y + 0.93, C.accent6, 1.5);
    });
    caption(s, "This project: own cluster, Singularity in /scratch · R1: 24 CPU, 150 GB · R2: 48 CPU, 500 GB, 1 A100", { x: 7.0, y: 6.3, w: 5.8, h: 0.5 });
  }
  {
    const s = content("The hackathon in one picture",
      "The hackathon project flow, as provided by the organisers (image credit: hackathon organisers). Phase 1, Predict X: method projects such as pocket finding, ligand docking, orthology benchmarks, MS data analysis, kinase targets, gene essentiality, protein stability, fine-mapping and colocalisation. Phase 2, Models and datasets: each method is applied at scale (AFDB pocket atlas, pocket similarity, ortholog calls, MS reprocessing, kinase activity, stability atlas, locus-to-gene atlas). Phase 3, Discovery: projects that combine them (protein-universe similarity, protein complex evolution, kinase function, cancer dependencies, variant testing, QTL buffering). Source: report/slides/img/hackathon_project_flow.png (supplied by the organisers).");
    img(s, "hackathon_project_flow.png", { x: 0.6, y: 1.25, w: 12.15, h: 5.5 }, "Hackathon project flow: Phase 1 Predict X, Phase 2 Models and datasets, Phase 3 Discovery", "hackathon-flow-image");
    caption(s, "Image: hackathon organisers.", { x: 0.6, y: 6.85, w: 6, h: 0.3 });
  }
  {
    const s = content("Where this project sits in the flow",
      "Placing this project in the organisers' flow. It is a Phase-2 project (workspace phase2-multicellularity) built on the output of the Phase-1 orthology project: Orthology benchmarks produced Ancestra, which feeds the Phase-2 'Ortholog calls' step; this project applies Ancestra to reconstruct ancestral gene content at the origins of multicellularity. In round 2 it also touched the Phase-3 'Protein complex evolution' idea (orthology plus co-fractionation mass spectrometry), using the published LECA complex map of Cox et al. 2026. The project does not appear as its own box in the official diagram; its exact position is to be confirmed by the organisers. Source: report/slides/img/hackathon_project_flow.png; PROMPT.md; report/REPORT.md §9.");
    const ph = [["Phase 1 · Predict X", 0.6], ["Phase 2 · Models & datasets", 4.75], ["Phase 3 · Discovery", 8.9]];
    ph.forEach((p) => {
      box(s, p[1], 1.35, 3.85, 5.0, C.background2, "phase-column");
      txt(s, p[0], { x: p[1] + 0.2, y: 1.45, w: 3.5, h: 0.45, fontSize: 16, bold: true, color: C.text2 });
    });
    // Phase 1 box
    box(s, 0.8, 2.15, 3.45, 1.15, C.text2, "p1-orthology");
    txt(s, [{ text: "Orthology benchmarks", options: { bold: true, breakLine: true } }, { text: "built Ancestra", options: { fontSize: 14 } }], { x: 1.0, y: 2.15, w: 3.1, h: 1.15, fontSize: 16, color: C.background1, valign: "middle" });
    // Phase 2: ortholog calls + this project
    box(s, 4.95, 2.15, 3.45, 1.15, C.accent1, "p2-orthologcalls");
    txt(s, [{ text: "Ortholog calls", options: { bold: true, breakLine: true } }, { text: "Ancestra at scale", options: { fontSize: 14 } }], { x: 5.15, y: 2.15, w: 3.1, h: 1.15, fontSize: 16, color: C.background1, valign: "middle" });
    arrow(s, 4.27, 2.72, 4.93, 2.72, C.text2, 2);
    box(s, 4.95, 3.85, 3.45, 1.75, C.accent4, "this-project");
    icon(s, "FaTree", 5.1, 4.05, 0.6, C.text2);
    txt(s, [{ text: "THIS PROJECT", options: { bold: true, breakLine: true } }, { text: "phase2-multicellularity: ancestral gene content at origins of multicellularity", options: { fontSize: 13 } }],
      { x: 5.8, y: 3.9, w: 2.5, h: 1.65, fontSize: 16, color: C.background1, valign: "middle" });
    arrow(s, 6.67, 3.33, 6.67, 3.83, C.accent4, 2.5);
    txt(s, "uses Ancestra", { x: 6.8, y: 3.38, w: 1.6, h: 0.4, fontSize: 12, color: C.accent4, bold: true });
    // Phase 3: protein complex evolution
    box(s, 9.1, 3.85, 3.45, 1.75, C.background1, "p3-complex", { line: { color: C.accent4, width: 1.5 } });
    txt(s, [{ text: "Protein complex evolution", options: { bold: true, breakLine: true } }, { text: "orthology + co-frac + MS", options: { fontSize: 14 } }], { x: 9.3, y: 3.85, w: 3.1, h: 1.75, fontSize: 16, color: C.text1, valign: "middle" });
    s.addShape(pres.shapes.LINE, { x: 8.42, y: 4.72, w: 0.66, h: 0.001, line: { color: C.accent4, width: 2, dashType: "dash", endArrowType: "triangle" }, objectName: "touches-arrow" });
    txt(s, "touches (Cox et al. LECA complexes)", { x: 9.1, y: 5.65, w: 3.5, h: 0.5, fontSize: 12, color: C.text2, italic: true });
    // note
    box(s, 0.8, 3.85, 3.45, 1.75, C.background1, "tbc-note", { line: { color: C.accent5, width: 2, dashType: "dash" } });
    txt(s, "[position in official flow: to be confirmed]", { x: 0.95, y: 3.85, w: 3.15, h: 1.15, fontSize: 15, bold: true, color: C.text2, valign: "middle", align: "center" });
    txt(s, "no box of its own in the diagram", { x: 0.95, y: 4.95, w: 3.15, h: 0.5, fontSize: 13, color: C.text2, align: "center" });
    caption(s, "Simplified from the organisers' flow diagram; only the boxes this project connects to are shown.", { x: 0.6, y: 6.55, w: 12, h: 0.35 });
  }
  {
    const s = content("Thirteen projects, one shared recipe",
      "The 13 hackathon projects as listed in the hackathon README. Orthology prediction produced Ancestra, the tool this project builds on: our multicellularity study is a phase-2 project that uses Ancestra and its results (highlighted). Source: hackathon README (Projects list).");
    const P = [["FaDna", "Cancer essentiality"], ["FaCogs", "Cancer metabolic model"], ["FaPuzzlePiece", "Co-complex from co-essentiality"],
      ["FaExchangeAlt", "CNV → trans-transcriptome"], ["FaCrosshairs", "Pocket prediction (BURROW)"], ["FaKey", "Ligand docking (ADock)"],
      ["FaBalanceScale", "Pocket similarity"], ["FaBullseye", "Kinase targets"], ["FaTachometerAlt", "Kinase activity"],
      ["FaRulerCombined", "Protein stability (physDDG)"], ["FaMapSigns", "Coloc analysis"], ["FaVial", "MS analysis"],
      ["FaSitemap", "Orthology: Ancestra"]];
    const cols = 5, w = 2.3, h = 1.05, gx = 0.12, gy = 0.18, x0 = 0.6, y0 = 1.4;
    P.forEach((p, i) => {
      const x = x0 + (i % cols) * (w + gx), y = y0 + Math.floor(i / cols) * (h + gy);
      const hi = p[1].startsWith("Orthology");
      box(s, x, y, w, h, hi ? C.accent1 : C.background2, "project-tile");
      icon(s, p[0], x + 0.12, y + 0.22, 0.6, hi ? C.text2 : C.accent1);
      txt(s, p[1], { x: x + 0.82, y: y + 0.08, w: w - 0.9, h: h - 0.16, fontSize: 14, valign: "middle", bold: hi, color: hi ? C.background1 : C.text1 });
    });
    // this project tile
    const tx = x0 + 3 * (w + gx) + 0.45, ty = y0 + 2 * (h + gy);
    box(s, tx, ty, 2 * w + gx - 0.45, h, C.accent4, "this-project-tile");
    icon(s, "FaTree", tx + 0.12, ty + 0.22, 0.6, C.text2);
    txt(s, [{ text: "This deck: phase-2 project", options: { bold: true, breakLine: true } }, { text: "multicellularity with Ancestra" }],
      { x: tx + 0.85, y: ty + 0.08, w: 3.3, h: h - 0.16, fontSize: 15, valign: "middle", color: C.background1 });
    // arrow from Ancestra tile (index 11: col 1,row 2) to this tile
    const axx = x0 + 2 * (w + gx) + w, ayy = y0 + 2 * (h + gy) + h / 2;
    arrow(s, axx + 0.02, ayy, tx - 0.04, ayy, C.accent4, 2.5);
    caption(s, "Same recipe for all: one prompt, Claude as PI, subagents, human verification after the fact.", { x: 0.6, y: 5.35, w: 12, h: 0.35 });
    txt(s, "builds on", { x: axx - 0.6, y: ay(), w: 1.7, h: 0.3, fontSize: 12, color: C.accent4, align: "center" });
    function ay() { return y0 + 2 * (h + gy) + h + 0.02; }
  }
  {
    const s = content("Every report opens with a human disclaimer",
      "Each project repository and report starts with a human-written disclosure (quoted in shortened form). It says the work is an exercise in assessing agentic AI, that results were generated entirely by AI with minimal human intervention, that it is not a traditional publication, and that no result has yet been verified by a human. Our project-specific line: investigate multicellularity across the tree of life, executed by Claude Opus 5.5 with Opus/Sonnet subagents. Source: hackathon README (AI disclosure header); README.md; report/REPORT.md.");
    box(s, 0.6, 1.45, 7.6, 4.9, C.background2, "quote-card");
    txt(s, "“", { x: 0.85, y: 1.45, w: 1, h: 1, fontSize: 72, bold: true, color: C.accent1 });
    txt(s, [
      { text: "This is an exercise in assessing the bioinformatics capabilities of agentic AI … All results in this project were generated entirely by AI with minimal human intervention.", options: { breakLine: true } },
      { text: " ", options: { breakLine: true, fontSize: 8 } },
      { text: "Currently none of the results have been verified by a human.", options: { bold: true } },
    ], { x: 1.0, y: 2.35, w: 6.9, h: 3.0, fontSize: 18, color: C.text2, italic: true });
    caption(s, "Shortened from the hackathon's human-written header.", { x: 1.0, y: 5.75, w: 6.9, h: 0.35 });
    const items = [["FaRobot", C.accent1, "Generated by AI", "minimal human intervention"], ["FaFileAlt", C.accent3, "Not a publication", "an exercise in agentic science"],
      ["FaUserCheck", C.accent4, "Not yet verified", "human checks pending"]];
    items.forEach((it, i) => {
      const y = 1.55 + i * 1.6;
      icon(s, it[0], 8.7, y, 0.9, it[1]);
      txt(s, [{ text: it[2], options: { bold: true, breakLine: true } }, { text: it[3] }], { x: 9.8, y: y, w: 3.0, h: 0.9, fontSize: 16, valign: "middle" });
    });
  }

  // ============================================================================================
  section("B. The prompt");
  divider("B", "The initial prompt", "What the human asked for, which parts are generic, and what changed in round 2",
    "Section B: the prompts. Source: PROMPT.md (round 1, about 1,470 words, saved verbatim), PROMPT_ROUND2.md (round 2).");
  {
    const s = content("Anatomy of the round-1 prompt",
      "The round-1 prompt (PROMPT.md) splits cleanly into two kinds of instructions. GENERAL parts (teal) are a reusable recipe for any agent-run project: the agent is the PI and delegates; a time budget (48 h allocation) with milestones and START.txt; hardware and sandbox (24 CPUs, 150 GB, no GPU); literature review and design critique in at most three review rounds; reproducibility rules (never change metrics to look better, fixed seeds, DATA.md, verification by an agent that did not produce the result); a research tree with a status legend; a LESSONS file; reporting and writing standards; and the human disclaimer. PROJECT-SPECIFIC parts (orange) define the science: kinds of multicellularity (clonal, aggregative, facultative, grey zones), the use of Ancestra to reconstruct ancestral genomes before and after the transitions, datasets (reference proteomes, EukProt, LukProt and 'Eukaryoma'), a menu of analyses (gain, duplication, loss, enrichment, ancestral structure, domain expansion), and the question of detecting a shift in selection pressure in sequences or structures. Source: PROMPT.md.");
    // centre doc icon
    icon(s, "FaFileAlt", 6.05, 1.5, 1.2, C.text2);
    txt(s, "PROMPT.md\n≈1,470 words", { x: 5.4, y: 2.8, w: 2.5, h: 0.7, fontSize: 14, align: "center", color: C.text2, bold: true });
    txt(s, "GENERAL · reusable", { x: 0.6, y: 1.35, w: 4.8, h: 0.4, fontSize: 18, bold: true, color: C.accent1 });
    txt(s, "PROJECT-SPECIFIC", { x: 8.0, y: 1.35, w: 4.8, h: 0.4, fontSize: 18, bold: true, color: C.accent4, align: "right" });
    const gen = ["Role: PI agent who delegates", "Time budget + milestones", "Hardware and sandbox", "Literature + ≤3 design critiques",
      "Reproducibility rules", "Research tree", "AI lessons file", "Reporting + writing standards", "Human disclaimer"];
    const spec = ["Kinds of multicellularity", "Ancestra: ancestral genomes", "EukProt, LukProt, 'Eukaryoma'", "Menu of analyses",
      "Selection-pressure question"];
    gen.forEach((g, i) => {
      const y = 1.85 + i * 0.55;
      box(s, 0.6, y, 4.6, 0.46, C.accent1, "general-chip");
      txt(s, g, { x: 0.8, y, w: 4.3, h: 0.46, fontSize: 15, color: C.background1, valign: "middle" });
    });
    spec.forEach((g, i) => {
      const y = 1.85 + i * 0.99;
      box(s, 8.2, y, 4.55, 0.85, C.accent4, "specific-chip");
      txt(s, g, { x: 8.4, y, w: 4.2, h: 0.85, fontSize: 16, bold: true, color: C.background1, valign: "middle" });
    });
    // proportion bar
    txt(s, "9 generic blocks", { x: 5.3, y: 4.05, w: 2.7, h: 0.35, fontSize: 14, color: C.accent1, align: "center", bold: true });
    txt(s, "5 science blocks", { x: 5.3, y: 4.45, w: 2.7, h: 0.35, fontSize: 14, color: C.accent4, align: "center", bold: true });
    caption(s, "Most of the prompt is a reusable research protocol; the science fits in one paragraph.", { x: 5.35, y: 5.0, w: 2.6, h: 1.2, align: "center" });
  }
  {
    const s = content("Round 2: seven new scientific asks",
      "The round-2 prompt (PROMPT_ROUND2.md) added resources (48 CPUs, 500 GB, one A100 GPU, about 30 h soft deadline), asked for a department-head review (Fable agent) of round 1 to be published, and listed seven scientific points: (1) 'Eukaryoma' was a misunderstanding: it is EukaryOMA, a precomputed orthology database; (2) 3D structures: use Ancestra's structure mode, predicting structure letters with ProstT5 on the GPU; (3) the step from simple to complex multicellularity; (4) increased resolution: one Ancestra run per multicellularity event; (5) linear-motif binding domains, cell-death genes and cell-type transcription factors; (6) specific Pfam families and their domain architectures; (7) proteomics: protein complexes from the Marcotte co-fractionation study. Source: PROMPT_ROUND2.md.");
    box(s, 0.6, 1.35, 12.15, 0.75, C.text2, "resources-strip");
    icon(s, "FaMicrochip", 0.8, 1.42, 0.6, C.accent5, "FFFFFF");
    txt(s, "New resources: 48 CPU · 500 GB · one A100 GPU · ~30 h  —  and a department-head review first", { x: 1.6, y: 1.35, w: 11.0, h: 0.75, fontSize: 17, bold: true, color: C.background1, valign: "middle" });
    const pts = [["FaDatabase", "EukaryOMA", "the database round 1 'could not find'"], ["FaCubes", "3D structure", "ProstT5 structure letters on GPU"],
      ["FaLayerGroup", "Simple → complex", "not only origins of multicellularity"], ["FaSearch", "One run per event", "dense sampling around each origin"],
      ["FaPuzzlePiece", "Motifs, death, TFs", "motif-binding, cell-death, cell-type TF genes"], ["FaProjectDiagram", "Specific Pfam families", "domain architectures, functions"],
      ["FaUsers", "Protein complexes", "co-fractionation map (Marcotte lab)"]];
    pts.forEach((p, i) => {
      const col = i % 4, row = Math.floor(i / 4);
      const x = 0.6 + col * 3.07, y = 2.45 + row * 2.2;
      box(s, x, y, 2.9, 1.95, C.background2, "ask-card");
      icon(s, p[0], x + 0.18, y + 0.2, 0.7, C.accent1);
      txt(s, String(i + 1), { x: x + 2.3, y: y + 0.15, w: 0.45, h: 0.45, fontSize: 20, bold: true, color: C.accent4, align: "right" });
      txt(s, [{ text: p[1], options: { bold: true, breakLine: true, fontSize: 16 } }, { text: p[2], options: { fontSize: 14 } }],
        { x: x + 0.18, y: y + 1.0, w: 2.6, h: 0.9 });
    });
    box(s, 0.6 + 3 * 3.07, 4.65, 2.9, 1.95, C.accent1, "review-card");
    icon(s, "FaUserTie", 0.6 + 3 * 3.07 + 0.18, 4.85, 0.7, C.text2);
    txt(s, [{ text: "Head review", options: { bold: true, breakLine: true, fontSize: 16 } }, { text: "Fable agent reads all of it, published", options: { fontSize: 14 } }],
      { x: 0.6 + 3 * 3.07 + 0.18, y: 5.65, w: 2.6, h: 0.9, color: C.background1 });
  }
  {
    const s = content("An external review set the round-2 agenda",
      "Bridge between rounds. A department-head reviewer (Fable agent) read the round-1 report and research tree and wrote a one-page review, published at https://claude.ai/artifact/NfrexbqRBFmS3A97VvTqcR. Wrong or weak: averaging seven very different origins with 13-35% power at a 1.5-fold effect; overstated validation (deep-node positive controls failed 0 of 3); EukaryOMA wrongly declared not found. Good: hashed pre-registration, a calibrated null, simulated gain inflation, independent re-derivation, and the diagnosis of positive-control failures; two leads (simple-to-complex novelty, new domain architectures on the animal stem). Not evaluated: EukaryOMA cross-check, per-event dense runs, simple-to-complex step, specific Pfam families. Source: reviews/round2_department_head_review.txt; reviews/review_round1.html.");
    const cols = [["FaExclamationTriangle", C.accent4, "What was wrong", ["7 heterogeneous origins averaged: 13–35% power", "Validation overstated: deep-node controls 0 of 3", "EukaryOMA wrongly 'not found'"]],
      ["FaThumbsUp", C.accent2, "What was good", ["Hashed pre-registration, calibrated null", "Diagnosis of why controls failed", "Lead: novelty at simple → complex"]],
      ["FaCompass", C.accent3, "What to do next", ["Cross-check with EukaryOMA", "Dense per-event Ancestra runs", "Simple → complex, specific families"]]];
    cols.forEach((c, i) => {
      const x = 0.6 + i * 4.1;
      box(s, x, 1.45, 3.9, 4.7, C.background2, "review-col");
      icon(s, c[0], x + 0.25, 1.65, 0.75, c[1]);
      txt(s, c[2], { x: x + 1.15, y: 1.65, w: 2.6, h: 0.75, fontSize: 19, bold: true, color: C.text2, valign: "middle" });
      txt(s, c[3].map((t, j) => ({ text: t, options: { bullet: true, breakLine: j < c[3].length - 1 } })),
        { x: x + 0.25, y: 2.65, w: 3.45, h: 3.3, fontSize: 16, paraSpaceAfter: 10 });
    });
    caption(s, "Fable agent as department head · reviews/round2_department_head_review.txt · public page linked at the end", { x: 0.6, y: 6.35, w: 12, h: 0.35 });
  }

  // ============================================================================================
  section("C. Overview");
  divider("C", "Overview", "Two rounds, one AI team, about 31 hours",
    "Section C: timeline, team and key numbers. Source: RESEARCH_TREE.md; report/REPORT.md.");
  {
    const s = content("Two rounds in about 31 hours",
      "Timeline of both rounds from the research tree milestones. Round 1 started 2026-10-28 12:21 and ended with the final report at 01:55 the next night (13.6 h): literature at t=0.4 h, pre-registration frozen at 14:25, the large eukaryote Ancestra run finished 22:36, all tests and the report by 01:55. Round 2 started 2026-10-29 15:07 with the department-head review and a new plan (15:41), eight dense panels built by 17:33, the confirmatory family closed at 22:03, structure analyses finished 05:20, verification and the final report by 08:31 on 30 October (17.4 h). Between the rounds the human read the report and wrote the round-2 prompt. Source: RESEARCH_TREE.md (M0-M5, R2-M0-R2-M4); START.txt; START_ROUND2.txt.");
    img(s, "gantt.png", { x: 0.6, y: 1.45, w: 12.15, h: 4.6 }, "Gantt chart of both rounds");
    caption(s, "Milestones from RESEARCH_TREE.md. Wall-clock time; the agent team worked in parallel inside each round.", { x: 0.6, y: 6.35, w: 12, h: 0.35 });
  }
  {
    const s = content("One PI agent, many specialised agents",
      "Team structure. The human gave two prompts and will verify the results. A PI agent (Claude Opus 5.5) planned, delegated, reviewed and integrated, and kept the research record (research tree, lessons, decisions). Subagents (Opus 5.5 or Sonnet 5.5) did literature review, data building, analysis, design review (external PI, evolutionary biology, statistics, data/compute, AI failure modes), independent adversarial verification, citation and numbers audits, and writing and figures. A separate Fable agent acted as department head and reviewed round 1 and the final report. Subagents wrote plain-text summaries so that the PI could keep its context lean. Source: report/REPORT.md §11; RESEARCH_TREE.md; reviews/.");
    // human
    box(s, 4.9, 1.35, 3.5, 0.85, C.background2, "human");
    icon(s, "FaUser", 5.05, 1.43, 0.68, C.accent4);
    txt(s, [{ text: "Human", options: { bold: true, breakLine: true } }, { text: "2 prompts, verifies later", options: { fontSize: 13 } }], { x: 5.85, y: 1.38, w: 2.5, h: 0.8, fontSize: 16, valign: "middle" });
    arrow(s, 6.65, 2.22, 6.65, 2.6, C.text2, 2);
    // PI
    box(s, 4.3, 2.65, 4.7, 1.15, C.text2, "pi");
    icon(s, "FaRobot", 4.5, 2.83, 0.8, C.accent1);
    txt(s, [{ text: "PI agent · Claude Opus 5.5", options: { bold: true, breakLine: true } }, { text: "plans, delegates, decides, keeps the record", options: { fontSize: 14 } }],
      { x: 5.45, y: 2.7, w: 3.5, h: 1.05, fontSize: 17, color: C.background1, valign: "middle" });
    // Fable
    box(s, 9.6, 2.65, 3.15, 1.15, C.accent4, "fable");
    icon(s, "FaUserTie", 9.75, 2.85, 0.75, C.text2);
    txt(s, [{ text: "Department head", options: { bold: true, breakLine: true } }, { text: "Fable agent, external", options: { fontSize: 13 } }],
      { x: 10.6, y: 2.7, w: 2.1, h: 1.05, fontSize: 15, color: C.background1, valign: "middle" });
    arrow(s, 9.58, 3.22, 9.05, 3.22, C.accent4, 2);
    const subs = [["FaBookOpen", "Literature"], ["FaDatabase", "Data"], ["FaChartBar", "Analysis"], ["FaBalanceScale", "Design reviewers"],
      ["FaCheckDouble", "Verifiers"], ["FaClipboardCheck", "Audits"], ["FaPenNib", "Writing, figures"]];
    subs.forEach((sb, i) => {
      const x = 0.6 + i * 1.75;
      icon(s, sb[0], x + 0.42, 4.75, 0.75, i === 4 ? C.accent4 : C.accent1);
      txt(s, sb[1], { x: x, y: 5.6, w: 1.6, h: 0.6, fontSize: 14, align: "center", bold: true });
      arrow(s, 6.65, 3.85, x + 0.8, 4.7, C.accent6, 1.25);
    });
    caption(s, "Subagents: Opus 5.5 or Sonnet 5.5, chosen by the PI · verifiers never checked their own work", { x: 0.6, y: 6.4, w: 12, h: 0.35 });
  }
  {
    const s = content("The experiment in numbers",
      "Key numbers. Round 1: 187 eukaryotic and 82 prokaryotic proteomes, 2.97 million eukaryotic proteins (36 from transcriptomes). Round 2: eight dense panels with 227 proteomes (3.06 million proteins), 9 Ancestra reconstructions (8 panels plus one sensitivity panel without de novo proteomes), in addition to the two round-1 runs; 3 frozen pre-registrations as git tags (prereg-v1, prereg-r2, prereg-r2-p1); 16 literature hypotheses from about 250 references; every confirmatory number re-derived by an agent that did not produce it, and 2 claims overturned by that verification; about 31 h of wall time. Source: report/REPORT.md §2, §4, Methods; RESEARCH_TREE.md (M1, D11-D15).");
    const S = [["187 + 82", "eukaryotic + prokaryotic proteomes (round 1)"], ["2.97 M", "eukaryotic proteins in round 1"],
      ["8 · 227", "dense panels · proteomes (round 2)"], ["2 + 9", "Ancestra runs, round 1 + round 2"],
      ["16", "literature hypotheses (~250 references)"], ["3", "frozen pre-registrations (git tags)"],
      ["2", "own claims overturned by verification"], ["≈31 h", "total wall time, two rounds"]];
    S.forEach((st, i) => {
      const x = 0.6 + (i % 4) * 3.07, y = 1.5 + Math.floor(i / 4) * 2.55;
      box(s, x, y, 2.9, 2.3, C.background2, "stat-card");
      stat(s, x + 0.2, y + 0.2, 2.6, st[0], st[1], i === 6 ? C.accent4 : C.accent1, 36);
    });
  }

  // ============================================================================================
  section("D. Decision process");
  divider("D", "How the method was decided", "Literature, critique, pre-registration, calibration, verification",
    "Section D: decision process. Source: RESEARCH_TREE.md (decision log D0-D16); DESIGN.md; DESIGN_ROUND2.md; PREREG.md.");
  {
    const s = content("A pipeline of checks before any claim",
      "Decision flow used in both rounds. A literature review (LITERATURE.md, 16 hypotheses) fed a design (DESIGN.md), critiqued by independent agent reviewers in three rounds in round 1 (5 reviewers, then 3, then a pre-registration audit) and one round of 3 reviewers in round 2. The design was then frozen with hashed inputs and a git tag before any reconstruction output existed. Pilots and calibration (simulated genomes, null calibration on the real tree, GPU pilot) came before long runs. After the runs, an agent that did not produce a result re-derived it; then an external senior review and a numbers audit preceded the report. The department-head review of round 1 closed the loop into round 2. Source: RESEARCH_TREE.md M1-M2, R2-M0; report/REPORT.md §3, §11.");
    const row1 = [["FaBookOpen", "Literature", "16 hypotheses"], ["FaPenNib", "Design", "DESIGN.md"], ["FaBalanceScale", "Critique", "3 rounds R1 · 1 round R2"],
      ["FaLock", "Pre-register", "hash + git tag"], ["FaVial", "Pilot + calibrate", "simulations, nulls, GPU"]];
    const row2 = [["FaCogs", "Run", "Ancestra + tests"], ["FaCheckDouble", "Verify", "independent agent"], ["FaUserTie", "External review", "department head"],
      ["FaClipboardCheck", "Audit + report", "numbers, citations"]];
    const W = 2.25, G = 0.22;
    row1.forEach((r, i) => {
      const x = 0.6 + i * (W + G), y = 1.55;
      box(s, x, y, W, 1.85, C.background2, "flow-step");
      icon(s, r[0], x + 0.75, y + 0.15, 0.72, C.accent1);
      txt(s, [{ text: r[1], options: { bold: true, breakLine: true } }, { text: r[2], options: { fontSize: 13 } }], { x: x + 0.1, y: y + 0.95, w: W - 0.2, h: 0.85, fontSize: 16, align: "center" });
      if (i < 4) arrow(s, x + W + 0.02, y + 0.92, x + W + G - 0.02, y + 0.92, C.text2, 2);
    });
    arrow(s, 0.6 + 4 * (W + G) + W / 2, 3.45, 0.6 + 4 * (W + G) + W / 2, 3.95, C.text2, 2);
    row2.forEach((r, i) => {
      const x = 0.6 + (4 - i) * (W + G), y = 4.0;
      const hi = i === 1 || i === 2;
      box(s, x, y, W, 1.85, hi ? C.accent4 : C.background2, "flow-step");
      icon(s, r[0], x + 0.75, y + 0.15, 0.72, hi ? C.text2 : C.accent1);
      txt(s, [{ text: r[1], options: { bold: true, breakLine: true } }, { text: r[2], options: { fontSize: 13 } }], { x: x + 0.1, y: y + 0.95, w: W - 0.2, h: 0.85, fontSize: 16, align: "center", color: hi ? C.background1 : C.text1 });
      if (i < 3) arrow(s, x - 0.02, y + 0.92, x - G + 0.02, y + 0.92, C.text2, 2);
    });
    // loop back
    box(s, 0.6, 4.0, W, 1.85, C.text2, "round2-loop");
    icon(s, "FaRedo", 0.6 + 0.75, 4.15, 0.72, C.accent5);
    txt(s, [{ text: "Round 2", options: { bold: true, breakLine: true } }, { text: "review → new design", options: { fontSize: 13 } }], { x: 0.7, y: 4.95, w: W - 0.2, h: 0.85, fontSize: 16, align: "center", color: C.background1 });
    arrow(s, 0.6 + W / 2, 3.97, 0.6 + W / 2, 3.45, C.accent5, 2.5);
    caption(s, "Red: steps done by agents that did not produce the result.", { x: 0.6, y: 6.25, w: 12, h: 0.35 });
  }
  function decisionSlide(title, notes, items) {
    const s = content(title, notes);
    items.forEach((d, i) => {
      const x = 0.6 + i * 4.1, w = 3.85;
      txt(s, d[0], { x, y: 1.35, w, h: 0.4, fontSize: 14, bold: true, color: C.accent4 });
      box(s, x, 1.8, w, 1.45, C.background2, "observed");
      txt(s, [{ text: "Observed  ", options: { bold: true, color: C.text2 } }, { text: d[1] }], { x: x + 0.2, y: 1.85, w: w - 0.4, h: 1.35, fontSize: 15, valign: "middle" });
      arrow(s, x + w / 2, 3.28, x + w / 2, 3.68, C.text2, 2.5);
      box(s, x, 3.72, w, 1.45, C.accent1, "decided");
      txt(s, [{ text: "Decided  ", options: { bold: true } }, { text: d[2] }], { x: x + 0.2, y: 3.77, w: w - 0.4, h: 1.35, fontSize: 15, valign: "middle", color: C.background1 });
      arrow(s, x + w / 2, 5.2, x + w / 2, 5.6, C.text2, 2.5);
      txt(s, [{ text: "Effect  ", options: { bold: true, color: C.accent1 } }, { text: d[3] }], { x: x + 0.1, y: 5.65, w: w - 0.2, h: 1.1, fontSize: 15 });
    });
    return s;
  }
  decisionSlide("Key decisions, round 1",
    "Three round-1 decisions. D0 (t=0): eukaryote and prokaryote reconstructions kept separate, because orthology across domains adds cost and no transition spans domains; the 1.4 h prokaryote run became the pipeline dry run, recovering both prokaryotic positive controls. Method calibration before results: on simulated genomes Ancestra inflated gene gains 2.1-2.6-fold (an old family split into several groups makes old genes look new), so gains were filtered to sequence-novel ones (no DIAMOND hit outside the clade), and a burst of novelty would be claimed only above 2.2-fold (Dev-4), fixed before eukaryote output existed. D6: on the real tree the planned relabelling null was anti-conservative (14.5% false positives at alpha 0.05), so a calibrated permutation scheme (4.9%) became primary. Source: RESEARCH_TREE.md D0, D6, M2, M3; report/REPORT.md §2.",
    [["1 · Separate runs", "Orthology across bacteria and eukaryotes adds cost; no transition spans both", "Separate eukaryote and prokaryote Ancestra runs", "The 1.4 h prokaryote run became the dry run"],
     ["2 · Novelty ladder", "Simulations: raw gene gains inflated 2.1–2.6×", "Count only sequence-novel gains; claim a burst only above 2.2×", "Rule fixed before any eukaryote output"],
     ["3 · Calibrated null", "Planned null: 14.5% false positives at α = 0.05 on the real tree", "Calibrated permutation scheme as the primary test", "4.9% false positives: p-values mean what they say"]]);
  decisionSlide("Key decisions, round 2",
    "Three round-2 decisions. After the review, eight dense per-event panels replaced one giant run: small runs avoid the giant similarity component that joined 73% of round-1 proteins, and each took 15-45 min; all used the fast preset, on which the 2.2-fold bias rule was calibrated. The confirmatory simple-to-complex novelty test (P1) was restricted to events with no round-1 novelty value and frozen (tag prereg-r2-p1) before their reconstructions existed, to avoid testing a lead on the data that suggested it; the cost was low power (18-41%). The rule for Ancestra's structure mode (adopt only if deep-node controls improve and EukaryOMA agreement does not fall) was fixed before any structure-mode output; it failed (23.3% to 20.7%), so the default was kept. Source: report/REPORT.md §3, §4, §8.2; RESEARCH_TREE.md D11-D15.",
    [["4 · Dense panels", "One average over 7 origins: 13–35% power; giant component of 73% of proteins", "8 small per-event panels, fast preset, 15–45 min each", "Sharper events, but only 5–9 controls per panel"],
     ["5 · Unseen events only", "Round-1 nested events already looked novel (p = 0.016)", "Confirmatory test on 4 new events, frozen before their runs", "No double dipping; power only 18–41%"],
     ["6 · Rule before data", "Deep nodes were Ancestra's weak point", "Use structure mode only if controls improve and agreement holds", "Rule failed (23.3% → 20.7%): default kept"]]);

  // ============================================================================================
  section("E. Method");
  divider("E", "The method, explained", "Simple pictures for the parts that non-experts find hard",
    "Section E: method explainers. Source: report/REPORT.md §1.3, §4-§9; external/ancestra-mvp-ai/README.md.");
  {
    const s = content("What is ancestral gene content?",
      "Ancestral gene content: given the proteins of living species and their family tree, infer which genes their common ancestors had. Genes are grouped into families (descendants of one ancestral gene). Reading the tree branch by branch gives gains (a family first appears), duplications (one ancestral gene becomes two) and losses. A hierarchical orthologous group (HOG) is a family seen at one level of the tree: each HOG at a node is one inferred ancestral gene. Schematic: the duplication of family 2 and the gain of family 3 on one branch; the gain of family 4 later; the loss of family 3 in species E. Source: report/REPORT.md §1.3.");
    img(s, "explain_ancestral.png", { x: 0.5, y: 1.35, w: 7.9, h: 5.4 }, "Schematic tree with gene families gained, duplicated and lost");
    const defs = [["FaLayerGroup", "Gene family", "genes descended from one ancestral gene"], ["FaSitemap", "HOG", "a family at one tree level = one ancestral gene"],
      ["FaCodeBranch", "Events", "gain, duplication or loss, read branch by branch"]];
    defs.forEach((d, i) => {
      const y = 1.6 + i * 1.65;
      icon(s, d[0], 8.75, y, 0.8, C.accent1);
      txt(s, [{ text: d[1], options: { bold: true, breakLine: true, fontSize: 18 } }, { text: d[2] }], { x: 9.75, y: y - 0.05, w: 3.0, h: 1.3, fontSize: 15 });
    });
    caption(s, "Schematic, not real data.", { x: 8.75, y: 6.45, w: 4, h: 0.3 });
  }
  {
    const s = content("Ancestra in three steps",
      "Ancestra (github.com/evocellnet/ancestra-mvp-ai) was itself built by AI agents in another hackathon project. Step 1: a reference-free family step clusters proteins with DIAMOND (no all-against-all search, no reference database), with a Leiden cap on giant components. Step 2: similarity is computed only inside families. Step 3: a level-synchronous merge walks the guide species tree from the tips to the root; at each node the ancestral genes of the child lineages are merged by mutual best matches and in-paralog absorption; ambiguous clusters get a small local tree (accurate preset: FastTree; fast preset: neighbour joining). Each merged ancestral gene is a HOG. The lower figure is Ancestra's own method figure, reused with credit. Source: external/ancestra-mvp-ai/README.md; docs/figures/fig1_method.png.");
    img(s, "explain_ancestra_steps.png", { x: 0.6, y: 1.3, w: 12.1, h: 3.1 }, "Ancestra simplified in three steps");
    img(s, "../../../external/ancestra-mvp-ai/docs/figures/fig1_method.png", { x: 0.6, y: 4.55, w: 8.4, h: 2.5, align: "left" }, "Ancestra method figure (from the Ancestra repository)");
    caption(s, "Lower figure: Ancestra's own method overview, reused with credit (evocellnet/ancestra-mvp-ai, docs/figures/fig1_method.png).", { x: 9.3, y: 5.1, w: 3.4, h: 1.2 });
  }
  {
    const s = content("Ancestra among orthology methods",
      "Comparison with other approaches. Common to all orthology methods: homology inferred from sequence similarity, usually as a graph; OMA/FastOMA, Hieranoid and Ancestra report HOGs at every level of a species tree (Hieranoid and Ancestra build them bottom-up). OrthoFinder gives orthogroups and gene trees that need a separate step to be mapped to ancestors. Pfam presence/absence with Dollo or Count parsimony skips orthology and reconstructs copy numbers on the tree. Ancestra's difference: a reference-free family step feeding a level-synchronous merge over all families at once. Benchmarks from the Ancestra README: OrthoBench test F (36 held-out reference groups) and QfO 2022 SwissTree recall (TPR). Caveats: on OrthoBench non-inferiority to OrthoFinder and SonicParanoid2 is not shown (wide intervals); on QfO Ancestra's recall is clearly below OrthoFinder (0.721 vs 0.825); OMA comparators are not independent of the benchmark. In this project, round-2 family tests used Pfam counts with phylogenetic regression, and EukaryOMA (OMA/FastOMA) served as the cross-check. Source: external/ancestra-mvp-ai/README.md (headline results); LITERATURE.md §5.1.");
    const hdr = ["Method", "How families are built", "Output", "OrthoBench F", "QfO SwissTree recall"];
    const rows = [
      ["Ancestra (this project)", "reference-free DIAMOND families, then level-by-level merge", "HOGs at every tree node", "76.2 / 74.8 (accurate / fast)", "0.721 / 0.731"],
      ["OrthoFinder 3.2", "all-vs-all similarity graph + gene trees", "orthogroups, gene trees", "69.7", "0.825"],
      ["OMA · FastOMA (EukaryOMA)", "pairwise orthologs; FastOMA places genes on a reference", "HOGs at every level", "60.4 / 33.1 ‡", "0.696 / 0.699"],
      ["SonicParanoid2 · InParanoid · Hieranoid", "graph / pairwise; Hieranoid bottom-up", "orthogroups, pairs, HOGs", "70.4 (SonicParanoid2)", "0.717 · 0.701 · 0.726"],
      ["Pfam counts + Dollo / Count", "domain annotation, no orthology", "copy numbers on the tree", "n/a", "n/a"],
    ];
    const head = hdr.map((h) => ({ text: h, options: { bold: true, color: C.background1, fill: { color: C.text2 }, fontSize: 14 } }));
    const body = rows.map((r, i) => r.map((c, j) => ({ text: c, options: { fontSize: 13, bold: j === 0, color: C.text1,
      fill: { color: i === 0 ? C.background2 : C.background1 } } })));
    s.addTable([head, ...body], { x: 0.6, y: 1.35, w: 12.15, colW: [2.75, 3.6, 2.2, 1.9, 1.7], border: { type: "solid", pt: 0.75, color: "C9D6D4" },
      margin: 0.06, valign: "middle", rowH: 0.62, objectName: "method-table" });
    box(s, 0.6, 5.35, 5.9, 1.2, C.background2, "common");
    icon(s, "FaHandshake", 0.75, 5.55, 0.8, C.accent1);
    txt(s, [{ text: "Shared  ", options: { bold: true, color: C.accent1 } }, { text: "homology from similarity graphs; a HOG = all genes from one ancestral gene" }], { x: 1.75, y: 5.4, w: 4.6, h: 1.1, fontSize: 15, valign: "middle" });
    box(s, 6.85, 5.35, 5.9, 1.2, C.background2, "different");
    icon(s, "FaExchangeAlt", 7.0, 5.55, 0.8, C.accent4);
    txt(s, [{ text: "Different  ", options: { bold: true, color: C.accent4 } }, { text: "no reference database; all families merged level by level at once" }], { x: 8.0, y: 5.4, w: 4.6, h: 1.1, fontSize: 15, valign: "middle" });
    caption(s, "Numbers from the Ancestra README. ‡ OMA comparators overlap the benchmark. Non-inferiority to OrthoFinder not shown; QfO recall below OrthoFinder.", { x: 0.6, y: 6.65, w: 11.5, h: 0.4 });
  }
  {
    const s = content("Why 'new gene' counts are inflated",
      "Gains are the least certain quantity. If a method splits one old family into two groups, one of them appears to arise on an inner branch: a false gain. On simulated genomes Ancestra inflated gains 2.1-2.6-fold. The project therefore used a novelty ladder: a gain counts as sequence-novel only if none of up to three representative proteins has a DIAMOND hit (E <= 1e-5, >= 30% query coverage) to any protein outside the clade. In the round-2 panels 41-80% of Ancestra's gains were sequence-novel. In round 2 a structural search (Foldseek) found that about half of sequence-novel families have a structural relative outside their clade, so sequence novelty still overstates true novelty. Source: report/REPORT.md §2, §4, §8.3.");
    img(s, "explain_split.png", { x: 0.5, y: 1.35, w: 4.9, h: 5.0 }, "A split family creates a false gain");
    img(s, "explain_ladder.png", { x: 5.6, y: 1.5, w: 7.2, h: 4.4 }, "Novelty ladder: all gains, sequence-novel, structurally unmatched");
    caption(s, "Each rung removes 'new' genes that have detectable relatives outside the clade.", { x: 5.7, y: 6.1, w: 7, h: 0.35 });
  }
  {
    const s = content("Compare like with like",
      "Transition branches are compared with control branches: internal branches inside unicellular lineages that are not events. A regression fitted on control branches (negative binomial: branch length, number of descendant species, data quality, distance to the outgroup) predicts how many sequence-novel families a branch should have; the rate ratio is observed divided by predicted. The tip-versus-internal artefact: a gain on a tip (single-species) branch is any gene of that species without a detectable relative, including orphan genes and mis-predicted gene models. Against internal controls every kind of tip scores high whatever its biology: unicellular tips with no event 8.1-fold, origin tips 19-fold, reversal tips 20.7-fold, simple-to-complex tips 40.8-fold. Tips must be compared with tips. Source: report/REPORT.md §6.1, §6.3; Figure 3b.");
    img(s, "explain_controls.png", { x: 0.5, y: 1.35, w: 5.4, h: 4.9 }, "Event branch compared with control branches");
    img(s, "explain_tips.png", { x: 6.2, y: 1.35, w: 6.5, h: 4.6 }, "Rate ratios of tips against internal controls");
    caption(s, "Right: geometric mean / pooled rate ratios from report Figure 3b. A tip looks 'novel' just for being a tip.", { x: 6.3, y: 6.15, w: 6.4, h: 0.5 });
  }
  {
    const s = content("Phylogenetic regression in one picture",
      "Phylogenetic generalised least squares (PGLS): related species share ancestry, so they are not independent data points. PGLS is a regression whose residuals are correlated according to the species tree, with the strength of that correlation estimated as Pagel's lambda. The family tests asked whether multicellular species carry more copies of a family set than unicellular species of the same proteome size (covariates: log proteome size and disordered fraction), on 148 genome-based proteomes. The right panel is a schematic with illustrative data. Source: report/REPORT.md §7.1 (Pagel 1999; Ho & Ane 2014).");
    img(s, "explain_pgls.png", { x: 0.6, y: 1.35, w: 12.1, h: 4.75 }, "PGLS schematic");
    txt(s, [{ text: "PGLS: ", options: { bold: true, color: C.accent1 } }, { text: "regression whose errors follow the species tree (Pagel's λ) · 148 genomes · covariates: proteome size, disorder" }],
      { x: 0.6, y: 6.25, w: 12, h: 0.45, fontSize: 15 });
  }
  {
    const s = content("Pre-registration and multiple testing",
      "Pre-registration: the design, test statistics and inputs were written down, hashed (sha256) and tagged in git before any relevant output existed (round 1: prereg-v1 with 15 hashed inputs; round 2: prereg-r2 for design and family lists, prereg-r2-p1 for the four confirmatory events). Every later change was logged as a dated deviation with a 'Seen:' line saying which results were visible. Holm correction: with three confirmatory tests, the smallest p-value must beat 0.05/3, the next 0.05/2, the last 0.05; testing stops at the first failure. Cell death (p = 0.0047) beat 0.0167; simple-to-complex novelty (p = 0.0504) did not beat 0.025; motif-binding (p = 0.72) was null. Source: report/REPORT.md §3, §6.2, §7, Methods; PREREG_HASHES.tsv; data/r2_prereg/PREREG_R2_HASHES.tsv.");
    img(s, "explain_prereg.png", { x: 0.6, y: 1.45, w: 5.9, h: 2.1 }, "Pre-registration steps");
    const tags = [["prereg-v1", "round 1 · 15 hashed inputs"], ["prereg-r2", "round 2 design + family lists"], ["prereg-r2-p1", "4 confirmatory events"]];
    tags.forEach((t, i) => {
      const y = 3.85 + i * 0.85;
      icon(s, "FaLock", 0.7, y, 0.6, C.text2);
      txt(s, [{ text: t[0], options: { bold: true, breakLine: true } }, { text: t[1], options: { fontSize: 14 } }], { x: 1.5, y: y - 0.05, w: 4.8, h: 0.75, fontSize: 15 });
    });
    img(s, "explain_holm.png", { x: 6.8, y: 1.35, w: 5.95, h: 4.4 }, "Holm thresholds for the three confirmatory tests");
    caption(s, "Holm: smallest p must beat α/3, the next α/2, …; stop at the first failure.", { x: 6.9, y: 6.0, w: 5.8, h: 0.4 });
  }
  {
    const s = content("Structure from sequence, on one GPU",
      "Predicting full 3D structures for three million proteins was not feasible in the time available. Foldseek describes each residue of a structure by one of 20 3Di letters (local geometry of contact with its nearest neighbour), so a structure becomes a searchable string. ProstT5, a protein language model, predicts 3Di letters directly from sequence. A pilot showed Foldseek's GPU build silently ran the model on the CPU (15 proteins/s); the team's own implementation of the ProstT5 encoder and prediction head ran at 48 proteins/s on the A100, agreed with Foldseek's ProstT5 output at 97.1% of residues, and matched letters computed from AlphaFold models at 65.5% of residues (chance 16.7%; always the commonest letter 31.4%; published range 41-65%). Source: report/REPORT.md §8.1; RESEARCH_TREE.md R2-3, R2-4 (van Kempen et al. 2024; Heinzinger et al. 2024).");
    img(s, "explain_prostt5.png", { x: 0.6, y: 1.35, w: 12.1, h: 2.9 }, "Sequence to ProstT5 to 3Di letters to Foldseek");
    img(s, "explain_gpu_bars.png", { x: 0.6, y: 4.35, w: 8.6, h: 2.6, align: "left" }, "GPU throughput and accuracy bars");
    box(s, 9.5, 4.55, 3.25, 2.1, C.background2, "gpu-note");
    txt(s, [{ text: "Lesson  ", options: { bold: true, color: C.accent4 } }, { text: "check the device is really used: the pilot caught a silent CPU fallback in the first hour" }],
      { x: 9.7, y: 4.65, w: 2.9, h: 1.9, fontSize: 15, valign: "middle" });
  }
  {
    const s = content("Cross-checking with a second method",
      "EukaryOMA (eukaryoma.unil.ch, release All.Oct2026) is a precomputed orthology database for 1,328 eukaryotic genomes built with OMA/FastOMA. 125 of our 187 species have a EukaryOMA genome; 87.1% of our proteins were mapped by identical sequence or reciprocal best hits. Both sets of HOGs were restricted to the shared species, stripped of their own branch labels, and placed on the same pruned tree by the same parser; a family's origin is the common ancestor of its members. Then: do both methods place the origin of the same family on the same branch? On six transition stems, 29% of Ancestra's stem gains (540 of 1,859) were placed on the same stem by EukaryOMA. It is a cross-method check on largely the same proteins, not independent data. Source: report/REPORT.md §5; Supplement S3.");
    img(s, "explain_xcheck.png", { x: 0.6, y: 1.35, w: 12.1, h: 4.9 }, "Two methods, same proteins, same tree: where did each family arise?");
    caption(s, "Same proteins, same tree, same parser: only the orthology method differs. Agreement measures consistency, not accuracy.", { x: 0.6, y: 6.4, w: 12, h: 0.35 });
  }
  {
    const s = content("Round 2: eight focused reconstructions",
      "Round-2 design. Eight dense panels (18-40 proteomes each, 227 in all), each built around one or a few origins with as many related genomes as were available, with a fixed rooted species tree and a table of coded events (59 in total: origins, nested simple-to-complex events, reversals). Each reconstruction took 15-45 min, against about 9-10 h for round 1's single eukaryote run. Panels sharpen single events but do not add independent origins, and most panels have only 5-9 control branches, which limits every per-event test. Figure: report Figure 1a-b (panels by state and coded events by kind). Source: report/REPORT.md §4, Table 1, Figure 1.");
    img(s, "fig1_ab.png", { x: 0.6, y: 1.35, w: 12.1, h: 3.55 }, "Report Figure 1a-b: proteomes per panel and coded events");
    const S = [["18–40", "proteomes per panel (227 in all)"], ["15–45 min", "per reconstruction (round 1: ~9–10 h)"], ["5–9", "control branches per panel: coarse tests"]];
    S.forEach((st, i) => stat(s, 0.8 + i * 4.1, 5.0, 3.8, st[0], st[1], i === 2 ? C.accent4 : C.accent1, 34));
  }

  // ============================================================================================
  section("F. Main results");
  divider("F", "Main results", "One message per slide",
    "Section F: main results. Master source: report/REPORT.md.");
  {
    const s = content("Round 1: no shared signature of seven origins",
      "Round 1 compared the transition stems of seven independent eukaryotic origins (animals, fungal hyphae, dictyostelids, streptophytes, red algae, brown-algal lineage, volvocine algae) with matched control branches. None of the five pre-registered primary tests (sequence-novel gains, transcription-factor duplications, adhesion/signalling duplications, flagellar loss, protein evolutionary rate) was significant after Holm correction (smallest adjusted p = 0.10). The gene-content tests had only 13-35% power at a 1.5-fold effect, so they are absence of evidence; the rate test, however, excludes a 1.5-fold or larger average shift (pooled rate ratio 1.10, 95% CI 0.96-1.26). Figure: round-1 Figure 3. Source: report/REPORT_round1.md; report/REPORT.md §2; RESEARCH_TREE.md M4.");
    img(s, "../../../figures/fig3_primary_tests.png", { x: 0.5, y: 1.3, w: 8.6, h: 5.5 }, "Round-1 Figure 3: five primary tests");
    stat(s, 9.4, 1.5, 3.4, "0 of 5", "primary tests significant after Holm", C.accent4, 38);
    stat(s, 9.4, 3.2, 3.4, "1.10×", "rate ratio [0.96–1.26]: no ≥1.5-fold shift", C.accent1, 38);
    stat(s, 9.4, 4.9, 3.4, "13–35%", "power of gene-content tests at 1.5-fold", C.accent5, 38);
  }
  {
    const s = content("Errata: the review caught three overstatements",
      "An agent checked the reviewer's points with new scripts. The depth imbalance between stems and controls was 2.06, not 1.77, control standard deviations. Positive controls were overstated as validation of Ancestra: of the 10 controls that need Ancestra to place an event on a branch, 7 passed, but one is a null control; excluding it, 6 of 9 passed, and all three deep-node placement controls failed. The rate test depended on a random seed (p 0.19-0.61); the exact, orientation-averaged version over 1,024 orientations gives p = 0.31. The power table had used placeholder branch lengths (re-running changed no value by more than 0.025), and the round-1 record had wrongly declared EukaryOMA not found. None changes a round-1 conclusion. Source: report/REPORT.md §2; results/round1_errata/ERRATA.txt.");
    const rows = [["Depth imbalance, stems vs controls", "1.77 SD", "2.06 SD"],
      ["Positive controls that test Ancestra", "'11 of 16'", "6 of 9; deep-node 0 of 3"],
      ["Rate-test p-value", "seed-dependent 0.19–0.61", "exact: p = 0.31"],
      ["'Eukaryoma' database", "'not found'", "EukaryOMA: 1,328 genomes"]];
    txt(s, "Claim", { x: 1.6, y: 1.35, w: 4.2, h: 0.4, fontSize: 14, bold: true, color: C.text2 });
    txt(s, "Round 1 said", { x: 6.0, y: 1.35, w: 2.6, h: 0.4, fontSize: 14, bold: true, color: C.accent4 });
    txt(s, "Corrected", { x: 9.4, y: 1.35, w: 3.3, h: 0.4, fontSize: 14, bold: true, color: C.accent2 });
    rows.forEach((r, i) => {
      const y = 1.85 + i * 1.08;
      icon(s, "FaUndo", 0.65, y + 0.12, 0.7, C.accent1);
      box(s, 1.5, y, 11.25, 0.92, C.background2, "erratum-row");
      txt(s, r[0], { x: 1.7, y, w: 4.2, h: 0.92, fontSize: 16, bold: true, valign: "middle" });
      txt(s, r[1], { x: 6.0, y, w: 2.6, h: 0.92, fontSize: 16, valign: "middle", color: C.accent4 });
      arrow(s, 8.65, y + 0.46, 9.25, y + 0.46, C.text2, 2);
      txt(s, r[2], { x: 9.4, y, w: 3.3, h: 0.92, fontSize: 16, valign: "middle", bold: true, color: C.text2 });
    });
    box(s, 1.5, 6.25, 11.25, 0.55, C.accent1, "erratum-footer");
    txt(s, "None of these changes a round-1 conclusion.", { x: 1.7, y: 6.25, w: 10.9, h: 0.55, fontSize: 16, bold: true, color: C.background1, valign: "middle" });
  }
  {
    const s = content("Two orthology methods disagree on where genes arose",
      "Ancestra and EukaryOMA place the origin of only a minority of the same gene families on the same branch: 29% of Ancestra's stem gains (540 of 1,859) are placed on the same stem by EukaryOMA, only 9% in the other direction; duplications agree in 4-8% of cases on stems and controls alike. On the animal stem Ancestra reports about 6,600 duplications and 804 new families, EukaryOMA 196 duplications and 2,502 new families. Stems do not agree better than controls once branch depth is taken into account (depth-adjusted expectation 22-24.7%). A count of new families or duplications on a stem therefore says as much about the method as about the biology. Figure: report Figure 2a-b. Source: report/REPORT.md §5.");
    img(s, "fig2_ab.png", { x: 0.6, y: 1.3, w: 12.1, h: 3.85 }, "Report Figure 2a-b: new families and duplications per stem by method");
    const S = [["29%", "of Ancestra stem gains on the same branch in EukaryOMA"], ["9%", "of EukaryOMA stem gains confirmed by Ancestra"],
      ["4–8%", "agreement on duplications"], ["6,617 vs 196", "animal-stem duplications, Ancestra vs EukaryOMA"]];
    S.forEach((st, i) => stat(s, 0.7 + i * 3.07, 5.15, 2.9, st[0], st[1], i === 0 ? C.accent4 : C.accent1, 30));
  }
  {
    const s = content("Simple to complex: no significant burst of new genes",
      "Confirmatory test P1: four simple-to-complex events never analysed in round 1 (Astrephomene, Volvox carteri, Volvox africanus/reticuliferus, the Fucophycidae brown-algal stem). Combined one-sided p = 0.0504 (Holm-adjusted 0.10): not significant. A variant with an undocumented branch-length floor gives 0.072. The large per-event ratios (13.6 and 14.9) belong to the two single-species tip events, and every kind of tip scores high against internal controls (panel b). Compared tip with tip, Astrephomene scores 0.56 and V. carteri 0.80, combined p = 0.49. Three of the four events share one panel and its controls; treated as one effect per panel, the smallest attainable p (0.031) cannot pass Holm. Figure: report Figure 3a-b. Source: report/REPORT.md §6.2-6.3.");
    img(s, "fig3_a.png", { x: 0.5, y: 1.3, w: 8.7, h: 2.6, align: "left" }, "Report Figure 3a: confirmatory test");
    img(s, "fig3_b.png", { x: 0.5, y: 4.05, w: 8.7, h: 2.6, align: "left" }, "Report Figure 3b: tips score high against internal controls");
    stat(s, 9.5, 1.45, 3.3, "p = 0.050", "combined, Holm-adjusted 0.10", C.accent5, 34);
    stat(s, 9.5, 3.15, 3.3, "8.1×", "unicellular tips with no event, vs internal controls", C.accent4, 34);
    stat(s, 9.5, 4.85, 3.3, "0.56–0.80", "the tip events, compared tip with tip", C.accent1, 34);
  }
  {
    const s = content("Internal branches only: a weak lead",
      "Restricted to internal branches, seven simple-to-complex events (Pezizomycotina, Agaricomycetes, land plants, V. africanus/reticuliferus, Florideophyceae, Fucophycidae, heterocystous Nostocales) pool to a rate ratio of 2.81 (95% CI 1.34-5.89) as specified (2.70 as implemented, plotted), against 0.82 for twelve internal origins from unicellular ancestors. This is the one gene-content lead that survives round 2, but it is weak: three of seven events were already seen in round 1, it was not pre-registered, events lie deeper than controls, and the lower bound does not clear the 2.2-fold method-bias threshold (itself from a single simulation). Figure: report Figure 3d. Source: report/REPORT.md §6.4.");
    img(s, "fig3_d.png", { x: 0.5, y: 1.3, w: 8.4, h: 5.5, align: "left" }, "Report Figure 3d: internal-only estimate");
    stat(s, 9.3, 1.5, 3.5, "2.8×", "simple → complex, internal [1.34–5.89]", C.accent3, 40);
    stat(s, 9.3, 3.2, 3.5, "0.82×", "origins of simple multicellularity", C.accent2, 40);
    box(s, 9.3, 4.95, 3.45, 1.6, C.background2, "bias-note");
    txt(s, "Lower bound below the 2.2× bias line: exploratory only", { x: 9.5, y: 5.0, w: 3.1, h: 1.5, fontSize: 15, valign: "middle", bold: true, color: C.text2 });
  }
  {
    const s = content("Cell death: the one confirmed signal",
      "Confirmatory test P3: 22 pre-defined cell-death families (caspases, metacaspases, NB-ARC, NACHT, HET, HeLo, gasdermins, Bcl-2 and others; 12 kept by the pre-set rule). On 148 genome-based proteomes, PGLS with proteome size and disorder as covariates: multicellular species carry +0.167 SD more (95% CI 0.042-0.293), one-sided p = 0.0047, a geometric-mean 1.21-fold excess (1.04-1.42); Holm-adjusted p = 0.014. The excess is concentrated in complex multicellularity (+0.249, p = 0.0005; simple +0.111, p = 0.073). Robust to leaving out any origin or family, reproduced independently; but without the ten largest proteomes the effect shrinks about 28% (p = 0.025) and as a share of the proteome it is not significant (p = 0.14): 'more than expected for proteome size'. Figure: report Figure 4a. Source: report/REPORT.md §7.2.");
    img(s, "fig4_a.png", { x: 0.5, y: 1.3, w: 8.6, h: 4.2, align: "left" }, "Report Figure 4a: PGLS estimates for three family sets");
    stat(s, 9.4, 1.4, 3.4, "1.21×", "more cell-death proteins [1.04–1.42], Holm p = 0.014", C.accent2, 44);
    stat(s, 9.4, 3.2, 3.4, "+0.249 SD", "complex vs unicellular (simple: +0.111)", C.accent3, 32);
    box(s, 0.6, 5.75, 12.15, 0.95, C.background2, "celldeath-caveat");
    icon(s, "FaExclamationTriangle", 0.75, 5.88, 0.7, C.accent5, "FFFFFF");
    txt(s, "Caveat: size-dependent. Without the 10 largest proteomes p = 0.025; as a share of the proteome p = 0.14. Read as 'more than expected for proteome size'.",
      { x: 1.65, y: 5.75, w: 11.0, h: 0.95, fontSize: 15, valign: "middle" });
  }
  {
    const s = content("The extra proteins are lineage-specific receptors",
      "What the extra cell-death proteins are (descriptive follow-up on 9,031 proteins carrying any of the 22 domains). Complex species carry on average 138 cell-death proteins, simple multicellular 18, unicellular 22. Land plants: NB-ARC proteins with leucine-rich repeats (plant NLR immune receptors; +84 per species vs Zygnematophyceae). Red algae: NB-ARC with WD40 repeats (+29). Filamentous fungi: HET domains (+37) and NACHT proteins with WD40 or ankyrin repeats (non-self recognition). Animals: death-fold adaptors (Death, CARD, DED), NACHT, caspases, Bcl-2, BIR. Dictyostelids and volvocines: no excess. Not one shared cell-death program: convergence on a function using different receptor architectures. Figure: report Figure 5a. Source: report/REPORT.md §7.3 (Jones & Dangl 2006; Dyrka et al. 2014).");
    img(s, "fig5_a.png", { x: 0.5, y: 1.3, w: 7.9, h: 5.5, align: "left" }, "Report Figure 5a: cell-death domains per species by lineage");
    const L = [[C.accent2, "Land plants", "NB-ARC–LRR immune receptors (+84 / species)"], [C.accent4, "Red algae", "NB-ARC–WD40 (+29)"],
      [C.accent5, "Filamentous fungi", "HET (+37), NACHT: non-self recognition"], [C.accent3, "Animals", "Death, CARD, DED adaptors; caspases"],
      [C.accent6, "Dictyostelids, volvocines", "no excess"]];
    L.forEach((l, i) => {
      const y = 1.4 + i * 1.07;
      s.addShape(pres.shapes.OVAL, { x: 8.7, y: y + 0.12, w: 0.36, h: 0.36, fill: { color: l[0] }, line: { color: l[0], width: 0 }, objectName: "lineage-dot" });
      txt(s, [{ text: l[1], options: { bold: true, breakLine: true } }, { text: l[2], options: { fontSize: 14 } }], { x: 9.2, y, w: 3.6, h: 0.95, fontSize: 16 });
    });
  }
  {
    const s = content("Motif-binding domains and TF families: null",
      "Confirmatory test P2: 115 Pfam families that bind short linear motifs (SH2, SH3, PDZ, WW, 14-3-3 and others, from the ELM resource) show no excess in multicellular species: -0.028 SD (95% CI -0.124 to 0.067), one-sided p = 0.72, a 0.98-fold ratio; six of seven leave-one-origin-out estimates are negative. A Bayesian count model gave a small positive estimate (sensitivity; read as a null with unresolved model dependence). Eleven cell-type transcription-factor families (exploratory): +0.032 SD (-0.088 to 0.152), p = 0.30. Known expansions reappear within single origins (T-box and Runt in animals, TEA in hyphal fungi, SBP in land plants). Figure: report Figure 4b-c (leave-one-out; proteome-size sensitivity of cell death). Source: report/REPORT.md §7.4.");
    img(s, "fig4_bc.png", { x: 0.5, y: 1.3, w: 8.6, h: 5.4, align: "left" }, "Report Figure 4b-c: leave-one-origin-out and proteome-size sensitivity");
    stat(s, 9.4, 1.5, 3.4, "0.98×", "motif-binding domains, p = 0.72", C.accent4, 44);
    stat(s, 9.4, 3.3, 3.4, "+0.03 SD", "cell-type TF families, p = 0.30", C.accent4, 36);
    caption(s, "Left: motif-binding estimates are negative in 6 of 7 leave-one-out runs; cell death stays positive in all 7.", { x: 9.4, y: 5.1, w: 3.35, h: 1.2 });
  }
  {
    const s = content("Structure did not help Ancestra",
      "Predicted structure letters are accurate enough to search with (65.5% residue agreement with AlphaFold-derived letters; panel a). Ancestra's optional structure mode rescues small families into large ones by structural similarity; the rule fixed before any output: replace the default only if the deep-node positive controls improve and agreement with EukaryOMA does not fall. On the Holozoa panel the rescue attached 25,035 small families to 2,400 large ones (a guard confirmed only 17%); deep-node controls stayed 1 of 3 in all three runs; agreement with EukaryOMA on new families at the two animal stems was 23.3% (default), 20.7% (structure) and 21.5% (with guard). The rule failed; the default was kept. One panel, one preset. Figure: report Figure 7a-b. Source: report/REPORT.md §8.1-8.2.");
    img(s, "fig7_a.png", { x: 0.5, y: 1.3, w: 6.0, h: 4.6 }, "Report Figure 7a: predicted vs AlphaFold-derived letters");
    img(s, "fig7_b.png", { x: 6.7, y: 1.3, w: 6.0, h: 4.6 }, "Report Figure 7b: structure mode inside Ancestra");
    const S = [["65.5%", "letters match AlphaFold (chance 16.7%)"], ["23.3% → 20.7%", "EukaryOMA agreement, default → structure"], ["1 of 3", "deep-node controls, in every run"]];
    S.forEach((st, i) => stat(s, 0.7 + i * 4.1, 5.85, 3.9, st[0], st[1], i === 0 ? C.accent2 : C.accent4, 26));
  }
  {
    const s = content("Two structural signals did not survive controls",
      "Structural novelty: event branches first seemed to carry more structurally unmatched sequence-novel families than controls (43.7% vs 33.3%, odds ratio 1.44, p = 0.004). An adversarial verifier showed unmatched families are much shorter (median 211 vs 519 aa) and more disordered: 'no detectable fold', not 'new fold'. Adjusted for length, disorder, family size, depth and search space, the odds ratio fell from 1.50 to 1.18 (0.93-1.49), and an exhaustive search found homologues for 21% of 'unmatched' families. Structure-per-sequence change: inside multicellular crowns 3Di change per amino-acid change was 7-12% higher, but a shuffled-letter negative control gave a comparable, in animals twice as large, effect (0.222 vs 0.110). Inconclusive. Figure: report Figure 7c-d. Source: report/REPORT.md §8.3-8.4.");
    img(s, "fig7_c.png", { x: 0.5, y: 1.3, w: 6.0, h: 4.6 }, "Report Figure 7c: structurally unmatched families, adjusted");
    img(s, "fig7_d.png", { x: 6.7, y: 1.3, w: 6.0, h: 4.6 }, "Report Figure 7d: structure per sequence change with controls");
    stat(s, 0.7, 5.85, 5.8, "OR 1.44 → 1.18 [0.93–1.49]", "after adjusting for length, disorder, depth", C.accent4, 26);
    stat(s, 6.9, 5.85, 5.8, "0.222 vs 0.110", "shuffled control vs real data (animals)", C.accent4, 26);
  }
  {
    const s = content("Protein complexes add nothing detectable",
      "Cox et al. (2026) reconstructed protein complexes of the last eukaryotic common ancestor (LECA) from co-fractionation mass spectrometry in 31 species. Only this published LECA map was used (raw elution profiles not reanalysed); it contains only complexes conserved across eukaryotes. Among testable duplications (13 species mappable through UniProt), 66.6% on the animal stem and 72.1% on the fungal-hyphae stem involved a complex subunit, against 63.6% on control branches; brown-algal and volvocine stems were depleted; no complex was enriched across stems after correction; co-retention was, if anything, slightly lower in multicellular species. Treated as a null. Source: report/REPORT.md §9; Supplement S7.");
    s.addChart(pres.charts.BAR, [{ name: "duplications involving a LECA complex subunit (%)", labels: ["Animal stem", "Fungal-hyphae stem", "Control branches"], values: [66.6, 72.1, 63.6] }], {
      x: 0.6, y: 1.4, w: 7.4, h: 5.0, barDir: "bar", chartColors: [HEX.accent3, HEX.accent2, HEX.accent6],
      showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: '0.0"%"', dataLabelFontSize: 14, dataLabelFontFace: "+mn-lt", dataLabelColor: HEX.dk1,
      valAxisMinVal: 0, valAxisMaxVal: 100, valAxisLabelFontSize: 12, valAxisLabelColor: HEX.dk2, valAxisLabelFontFace: "+mn-lt",
      catAxisLabelFontSize: 14, catAxisLabelColor: HEX.dk1, catAxisLabelFontFace: "+mn-lt", catAxisOrientation: "maxMin",
      valGridLine: { color: "DDE5E4", size: 0.75 }, catGridLine: { style: "none" }, showLegend: false,
      showTitle: true, title: "Duplications involving a LECA complex subunit", titleFontSize: 15, titleColor: HEX.dk2, titleFontFace: "+mn-lt",
      valAxisTitle: "% of testable duplications", showValAxisTitle: true, valAxisTitleFontSize: 12, valAxisTitleColor: HEX.dk2, valAxisTitleFontFace: "+mn-lt",
      objectName: "complex-chart" });
    icon(s, "FaUsers", 8.6, 1.6, 0.9, C.accent1);
    txt(s, "No complex enriched across stems after correction", { x: 9.7, y: 1.55, w: 3.1, h: 1.0, fontSize: 17, bold: true, valign: "middle" });
    txt(s, [{ text: "Map holds only ancient (LECA) complexes", options: { bullet: true, breakLine: true } },
      { text: "13 species mappable", options: { bullet: true, breakLine: true } },
      { text: "Duplication counts are method-dependent", options: { bullet: true } }], { x: 8.6, y: 2.95, w: 4.1, h: 2.2, fontSize: 15, paraSpaceAfter: 8 });
    caption(s, "Data: Cox et al. 2026, LECA interactome from co-fractionation mass spectrometry.", { x: 8.6, y: 5.6, w: 4.1, h: 0.7 });
  }
  {
    const s = content("Take-home messages",
      "Five take-home messages. (1) No large shared genomic signature across origins; the average protein evolutionary rate did not shift 1.5-fold or more. (2) The one confirmed signal: about 1.2-fold more cell-death proteins than expected for proteome size, mostly in complex lineages, made of lineage-specific immune and self/non-self receptors. (3) Counts of new families and duplications on a branch depend on the orthology method (29% agreement). (4) Tip branches look novel whatever their biology; compare tips with tips. (5) Predicted structure is a good search tool but did not help Ancestra and gave no robust evolutionary signal. Source: report/REPORT.md §10.1.");
    const M = [["FaGlobeEurope", C.accent1, "No shared signature", "across origins; no ≥1.5-fold rate shift"],
      ["FaShieldAlt", C.accent2, "Cell death, ~1.2×", "lineage-specific immune / self–non-self receptors"],
      ["FaExchangeAlt", C.accent4, "Method matters", "two orthology methods agree on 29% of stem gains"],
      ["FaBalanceScale", C.accent5, "Compare like with like", "single-species branches always look 'novel'"],
      ["FaCubes", C.accent3, "Structure: search, not signal", "did not help Ancestra; signals failed controls"]];
    M.forEach((m, i) => {
      const y = 1.4 + i * 1.08;
      icon(s, m[0], 0.7, y, 0.85, m[1]);
      txt(s, String(i + 1), { x: 1.75, y, w: 0.5, h: 0.85, fontSize: 28, bold: true, color: m[1], valign: "middle" });
      txt(s, [{ text: m[2] + "  ", options: { bold: true, fontSize: 20 } }, { text: m[3], options: { fontSize: 16 } }], { x: 2.35, y, w: 10.3, h: 0.85, valign: "middle" });
    });
  }

  // ============================================================================================
  section("G. Hypotheses");
  divider("G", "Every hypothesis, one by one", "16 from the literature review, plus the round-2 tests",
    "Section G: all hypotheses. H1-H16 come from LITERATURE.md §7 and were pre-registered or explored in round 1 (verdicts in results/literature_scorecard.tsv). Round-2 tests come from DESIGN_ROUND2.md §6 and report/REPORT.md. Verdict colours: green supported; red not supported / null / overturned; amber not detected at this power or not significant; blue exploratory lead; grey not tested.");
  {
    const s = content("Scorecard of all hypotheses",
      "Overview of all hypotheses and verdicts. Round 1 (H1-H16, LITERATURE.md §7): five pre-registered primary tests (H1/H2, H4, H5, H9, H12); others exploratory or cut for time. Round 2: three confirmatory tests (P1 simple-to-complex novelty, P2 motif-binding domains, P3 cell death) and exploratory analyses (structure, complexes). Sources: results/literature_scorecard.tsv; report/REPORT.md §10.1.");
    const H = [["H1–H2", "co-option vs burst", "notdetected"], ["H3", "pre-adaptation", "nottested"], ["H4", "TF expansion", "notdetected"],
      ["H5", "adhesion, signalling", "notdetected"], ["H6", "receptor kinases", "partly"], ["H7–H8", "aggregative vs clonal", "nottested"],
      ["H9", "relaxed selection", "notsupported"], ["H10", "disorder", "notdetected"], ["H11", "domain shuffling", "lead"],
      ["H12", "loss of unicellular genes", "notdetected"], ["H13", "cell-cycle co-option", "nottested"], ["H14–H15", "HGT, site convergence", "nottested"],
      ["H16", "regulation, not coding", "nottestable"], ["R2-P1", "simple → complex novelty", "notsig"], ["R2-P2", "motif-binding domains", "null"],
      ["R2-P3", "cell death", "supported"], ["R2", "structural novelty", "overturned"], ["R2", "protein complexes", "null"]];
    H.forEach((h, i) => {
      const col = i % 3, row = Math.floor(i / 3);
      const x = 0.6 + col * 4.1, y = 1.35 + row * 0.8;
      const v = VERDICT[h[2]];
      box(s, x, y, 3.9, 0.68, v.fill, "score-chip");
      txt(s, [{ text: h[0] + "  ", options: { bold: true } }, { text: h[1] }], { x: x + 0.15, y, w: 3.65, h: 0.68, fontSize: 15, color: v.fg, valign: "middle" });
    });
    const legend = [["supported", "supported"], ["notdetected", "not detected / n.s."], ["notsupported", "not supported / null"], ["lead", "exploratory lead"], ["nottested", "not tested"]];
    legend.forEach((l, i) => {
      const v = VERDICT[l[0]];
      s.addShape(pres.shapes.OVAL, { x: 0.6 + i * 2.45, y: 6.32, w: 0.3, h: 0.3, fill: { color: v.fill }, line: { color: v.fill, width: 0 }, objectName: "legend-dot" });
      txt(s, l[1], { x: 1.0 + i * 2.45, y: 6.27, w: 2.0, h: 0.4, fontSize: 13, valign: "middle" });
    });
  }
  function hypSlide(h) {
    const s = content(h.title, h.notes);
    img(s, h.img, { x: 0.45, y: 1.3, w: 6.3, h: 3.55 }, "Schematic of hypothesis " + h.title);
    box(s, 0.6, 5.05, 6.0, 1.75, C.background2, "literature-box");
    icon(s, "FaBookOpen", 0.75, 5.2, 0.55, C.text2);
    txt(s, [{ text: "Literature  ", options: { bold: true, color: C.text2 } }, { text: h.lit }], { x: 1.45, y: 5.12, w: 5.0, h: 1.6, fontSize: 12 });
    txt(s, [{ text: "Predicts  ", options: { bold: true, color: C.accent1 } }, { text: h.pred }], { x: 7.05, y: 1.35, w: 5.7, h: 1.1, fontSize: 15 });
    // pros
    box(s, 7.05, 2.55, 5.7, 1.35, C.background2, "pros");
    icon(s, "FaThumbsUp", 7.2, 2.72, 0.55, C.accent2);
    txt(s, h.pro, { x: 7.9, y: 2.6, w: 4.7, h: 1.25, fontSize: 14, valign: "middle" });
    box(s, 7.05, 4.05, 5.7, 1.35, C.background2, "cons");
    icon(s, "FaThumbsDown", 7.2, 4.22, 0.55, C.accent4);
    txt(s, h.con, { x: 7.9, y: 4.1, w: 4.7, h: 1.25, fontSize: 14, valign: "middle" });
    badge(s, 7.05, 5.6, 3.3, h.verdict);
    txt(s, h.key, { x: 7.05, y: 6.1, w: 5.7, h: 0.75, fontSize: 14, bold: true, color: C.text2 });
    return s;
  }
  const HYP = [
    { title: "H1–H2 · Co-option or a burst of new genes?", img: "hyp_h01_h02.png", verdict: "notdetected",
      pred: "Burst: more new gene families on transition stems than on controls. Co-option: no excess; the genes were already there.",
      pro: "Bursts reported at animal, land-plant and brown-algal origins; choanoflagellates already carry cadherins and kinases.",
      con: "Gain counts inflate 2.1–2.6×; denser outgroup sampling moves 'animal-specific' families earlier.",
      key: "pooled 1.67× [0.64–4.34], p = 0.11 (Holm 0.45); 80–100% of expanded families present before",
      lit: "Paps & Holland 2018; Bowles et al. 2020; Denoeud et al. 2024; King et al. 2008; Sebé-Pedrós et al. 2010; Prochnik et al. 2010; Kiss et al. 2019; Richter et al. 2018",
      notes: "H1 (co-option precedes the transition) and H2 (burst of novelty) were tested as one comparison: sequence-novel gains per stem, observed versus expected from matched controls, with a co-option index as companion. Result: pooled observed/expected 1.67 (0.64-4.34), parametric-bootstrap p = 0.113, Holm 0.45; below the 2.2-fold method-bias threshold. Co-option index 0.80-1.00: 80-100% of expanded Pfam families were present before the transition in all seven origins (descriptive support for co-option). Minimum detectable effect 2.4-5.9-fold, so not detected at this power. Sources: LITERATURE.md §7; results/literature_scorecard.tsv; report/REPORT_round1.md." },
    { title: "H3 · Pre-adaptation before the transition?", img: "hyp_h03.png", verdict: "nottested",
      pred: "Ancestors of lineages that later turned multicellular already had more adhesion, TF and signalling genes than sisters that never did.",
      pro: "'Predisposition' reported for complex fungi; gradual assembly of the toolkit in holozoans.",
      con: "Needs many transitioning and non-transitioning sister pairs with good genomes.",
      key: "not run in either round",
      lit: "Ocaña-Pallarès et al. 2022; Grau-Bové et al. 2017; Merényi et al. 2020",
      notes: "H3, gradual pre-adaptive accumulation: multicellularity-related genes accumulate on branches leading to the transition, not only on the transition branch. It was listed as exploratory in the design but not run (pre-transition accumulation on outgroup ancestors was not analysed). Sources: LITERATURE.md §7; results/literature_scorecard.tsv." },
    { title: "H4 · Do transcription-factor families expand?", img: "hyp_h04.png", verdict: "notdetected",
      pred: "More TF duplications on transition stems than on controls, strongest in complex lineages.",
      pro: "TF repertoires track cell-type number; expansions reappear within origins (T-box, Runt in animals; TEA in fungi; SBP in plants).",
      con: "Much of the TF toolkit predates multicellular lineages.",
      key: "R1: odds ratio 1.16 [0.62–2.19], p = 0.36 · R2 TF families: +0.03 SD, p = 0.30",
      lit: "de Mendoza et al. 2013; Larroux et al. 2008; Lang et al. 2010; Catarino et al. 2016; Knoll 2011; Nagy et al. 2014",
      notes: "H4: TF families expand at the transition, most strongly in complex lineages. Round 1 (primary): share of stem duplications in TF DNA-binding-domain HOGs vs matched controls, T = 0.548, p = 0.361, Holm 0.81, pooled odds ratio 1.16 (0.62-2.19); orthology-free version p = 0.18. Round 2 (exploratory): 11 cell-type TF families, +0.032 SD, p = 0.30. Known expansions reappear within single origins. Sources: results/literature_scorecard.tsv; report/REPORT.md §7.4." },
    { title: "H5 · Convergent adhesion and signalling expansion?", img: "hyp_h05.png", verdict: "notdetected",
      pred: "The same adhesion/ECM and signalling domains expand independently in two or more lineages.",
      pro: "Every multicellular body needs cell adhesion and communication.",
      con: "Classical cadherins are new domain combinations inside old families, invisible to family-level dating.",
      key: "odds ratio 0.83 [0.55–1.27], p = 0.77; Pfam-wide expansions not above chance",
      lit: "Knoll 2011; Abedin & King 2010; Pereira Lobo et al. 2025; Rokas 2008a",
      notes: "H5: adhesion/ECM and communication domains expand convergently. Round 1 primary: T = 0.412, p = 0.772, Holm 0.81; pooled odds ratio 0.83 (0.55-1.27); orthology-free version opposite sign. Pfam-wide convergent expansion 370 vs 643 expected (p = 0.79), top hits transposon domains; in round 2, expansions at two or more origins 42 vs 54 expected (p = 0.98), whereas contractions recur (157 vs 90, p = 0.001). Cadherin and integrin positive controls fail in both Ancestra and EukaryOMA: dating a domain gained inside an older family is a general limit. Sources: results/literature_scorecard.tsv; report/REPORT.md §5, §7.5." },
    { title: "H6 · Receptor kinases in complex lineages?", img: "hyp_h06.png", verdict: "partly",
      pred: "Receptor-kinase families duplicate on brown-algal, plant and animal stems, but not on fungal hyphae.",
      pro: "Receptor kinases expanded in Ectocarpus and Arabidopsis; R2: one new kinase family at 5 of 7 complex events.",
      con: "Receptor-kinase classification not run; kinase recurrence is only suggestive (q = 0.06).",
      key: "fungal 'no expansion' control passed; kinase recurrence p = 0.023, q = 0.06",
      lit: "Cock et al. 2010; Shiu & Bleecker 2001; Manning et al. 2008; Kiss et al. 2019",
      notes: "H6: receptor kinases and signalling expand independently in complex lineages but not in fungal hyphae. The receptor-kinase classification (transmembrane or extracellular domains) was not run. The expected null on the fungal-hyphae stem was recovered as positive control C15 (no kinase, adhesion or signalling expansion). Round 2 found one new protein-kinase family at five of seven internal simple-to-complex events (recurrence p = 0.023, q = 0.06; suggestive). Sources: results/literature_scorecard.tsv; report/REPORT.md §7.5." },
    { title: "H7–H8 · Aggregative versus clonal", img: "hyp_h07_h08.png", verdict: "nottested",
      pred: "H7: aggregative origins gain secreted and membrane signals. H8: clonal origins show more gene turnover.",
      pro: "Dictyostelid development relies on extracellular signals; clonal development allows more cell types.",
      con: "Few aggregative lineages with genomes; most are single species, so the stem is also a tip.",
      key: "one aggregative primary origin (dictyostelids); not run",
      lit: "Glöckner et al. 2016; Grosberg & Strathmann 2007; Fisher et al. 2013; Márquez-Zacarías et al. 2021",
      notes: "H7 (aggregative transitions need new extracellular signals and sensors) and H8 (clonal transitions involve more genomic change than aggregative ones) both concern the aggregative/clonal contrast. H8 was demoted to exploratory at design because aggregative transitions are mostly single species; the signal-peptide/transmembrane predictor for H7 was not run. Only one aggregative primary origin (dictyostelids) was available; the dictyostelid polyketide-synthase control (16 vs 1 proteins per species) was recovered but is not the H7 test. Sources: DESIGN.md §4.2; results/literature_scorecard.tsv." },
    { title: "H9 · Relaxed selection after the transition?", img: "hyp_h09.png", verdict: "notsupported",
      pred: "Smaller populations weaken selection: faster protein evolution and larger proteomes after eukaryotic transitions.",
      pro: "Drift theory of genome complexity, recently extended to multicellularity.",
      con: "Rates differ by lineage (4 stems faster, 3 slower); R2 structure-per-sequence test failed its control.",
      key: "pooled rate ratio 1.10 [0.96–1.26]: no ≥1.5-fold shift; exact p = 0.31",
      lit: "Lynch & Conery 2003; Lynch 2007; Bingham & Ratcliff 2024; Pentz et al. 2023",
      notes: "H9, the non-adaptive drift hypothesis: multicellularity lowers effective population size, leading to faster, less constrained protein evolution. Round 1 primary relative-rate test: pooled rate ratio 1.10 (0.96-1.26), excluding a 1.5-fold or larger average shift; the seed-dependent p was replaced in round 2 by the exact orientation-averaged p = 0.31. Gene-level RER found 0 genes at q < 0.10; radical/conservative substitution tests p = 0.60 and 0.16. Round 2 structure-per-sequence divergence was inconclusive (shuffled-letter control as large as the signal). Sources: results/literature_scorecard.tsv; report/REPORT.md §2, §8.4." },
    { title: "H10 · More disorder with complexity?", img: "hyp_h10.png", verdict: "notdetected",
      pred: "Disorder in transcription factors rises with cell-type number and on transition stems.",
      pro: "Disorder and alternative splicing scale with organismal complexity in earlier surveys.",
      con: "Mixed directions: disorder up in 4 of 7 transitions, down in volvocines and animals.",
      key: "percentiles vs controls 0.34–0.73; regression on cell-type number not run",
      lit: "Schad et al. 2011; Yruela et al. 2017; Chen et al. 2014; Dunker et al. 2015",
      notes: "H10: disorder and splicing scale with complexity. Exploratory round-1 analysis on shared single-copy HOGs: disorder rose in 4 of 7 transitions (+0.019 to +0.037) and fell in volvocines and animals; percentiles relative to controls 0.34-0.73. The phylogenetic regression on cell-type number was not run. Sources: results/literature_scorecard.tsv; results/exploratory/euk/summary.md §4." },
    { title: "H11 · New domain combinations?", img: "hyp_h11.png", verdict: "lead",
      pred: "More new multi-domain architectures on complex-lineage stems than on controls.",
      pro: "Domain shuffling builds new proteins; animal stem: 1,131 new architectures vs control median 34.",
      con: "Fails a completeness filter (BUSCO ≥ 80%: p = 0.16); exploratory only.",
      key: "p = 0.005, survives leave-one-out, not robust to data quality",
      lit: "Ekman et al. 2007; Denoeud et al. 2024; Ocaña-Pallarès et al. 2022",
      notes: "H11: domain rearrangement and gene fusion increase at complex transitions. Exploratory round-1 result: new multi-domain architectures on stems T = 0.80, p = 0.005; Metazoa 1,131 vs matched-control median 34; survives leave-one-out (p <= 0.017) and the genome-only subset, but not a BUSCO >= 80% completeness filter (p = 0.16); pooled observed/expected 1.02 (0.36-2.91). Flagged by the department head as worth follow-up. Sources: results/literature_scorecard.tsv; reviews/round2_department_head_review.txt." },
    { title: "H12 · Shedding unicellular functions?", img: "hyp_h12.png", verdict: "notdetected",
      pred: "Flagellar and phagotrophy genes are lost at transitions; reversals lose multicellularity genes.",
      pro: "Fungi and dictyostelids lost all four flagellar modules, as expected ('shedding of protist genes').",
      con: "Only two loss-expected stems (power ≤ 0.33); reversals do not preferentially shed transition genes.",
      key: "p = 0.020, Holm-adjusted 0.10; reversal odds 0.86–1.39",
      lit: "Merényi et al. 2023; Yang et al. 2026; Nguyen et al. 2017",
      notes: "H12: losses of unicellular functions accompany transitions; reversals lose multicellularity genes. Round 1 primary (flagellar-module loss with pre-registered direction per transition): T = 0.738, p = 0.020, Holm 0.10, not significant after correction; both loss-expected stems (fungal hyphae, dictyostelids) lose all four flagellar modules. The reviewer noted the test is largely a confirmation that spent Holm correction. Reversals (5 yeasts, Zygnematophyceae): odds of losing genes gained at the preceding transition 0.86-1.39 (not supported). Sources: results/literature_scorecard.tsv; reviews/round2_department_head_review.txt." },
    { title: "H13 · Is cell-cycle control co-opted?", img: "hyp_h13.png", verdict: "nottested",
      pred: "RB/E2F, ubiquitin–proteasome and ESCRT/septin genes duplicate or speed up on clonal stems only.",
      pro: "Gonium co-opted the retinoblastoma pathway; incomplete cell division links cells by bridges.",
      con: "Duplication and rate contrasts were not run; only one co-option control.",
      key: "RB pathway already present in the unicellular volvocine ancestor (control recovered)",
      lit: "Hanschen et al. 2016; Featherston et al. 2018; Chaigne & Brunet 2022",
      notes: "H13: cell-cycle and cytokinesis control is co-opted on clonal transition branches. Not tested as specified; the clonal vs aggregative contrast was not run. Positive control C11 was recovered: the RB pathway is present at the unicellular volvocine ancestor, consistent with co-option (Hanschen et al. 2016). Sources: results/literature_scorecard.tsv." },
    { title: "H14–H15 · Gene transfer and site convergence", img: "hyp_h14_h15.png", verdict: "nottested",
      pred: "H14: non-animal origins gain more genes by horizontal transfer. H15: the same amino-acid changes recur at independent origins.",
      pro: "Gene transfer is common in fungi, protists and brown algae; strong molecular convergence in complex fungi.",
      con: "Both need per-family gene trees; cut for time at the design stage.",
      key: "cut for time (DESIGN.md §4.2)",
      lit: "Ocaña-Pallarès et al. 2022; Denoeud et al. 2024; Sheikh et al. 2024; Glöckner et al. 2016; Merényi et al. 2020; Iizuka et al. 2026",
      notes: "H14 (horizontal gene transfer contributes more to non-animal transitions) and H15 (molecular convergence at the level of sites, testable with CSUBST/PCOC) were both cut for time in the round-1 design (DESIGN.md §4.2). They remain open questions for future work. Sources: LITERATURE.md §5.2, §7; DESIGN.md." },
    { title: "H16 · Regulation rather than new proteins?", img: "hyp_h16.png", verdict: "nottestable",
      pred: "Multicellularity mainly rewires gene regulation; protein-coding content changes little.",
      pro: "Dynamic regulatory genome in Capsaspora; our many nulls are compatible with it.",
      con: "Regulatory change is invisible in proteomes; compatibility is not evidence.",
      key: "outside the scope of protein-only data",
      lit: "Sebé-Pedrós et al. 2016; Jiménez-Marín et al. 2022; Nagy et al. 2014",
      notes: "H16: regulatory rather than coding change dominates. Only partly testable with proteomes (expected negative results for novelty in volvocines and fungal hyphae, with TF expansions). The report lists it as a limitation: only protein-coding content was studied; regulatory change, which much of the literature emphasises, is invisible here. Sources: LITERATURE.md §7; report/REPORT.md §10.4." },
    { title: "R2-P1 · New genes at simple-to-complex steps?", img: "hyp_r2_p1.png", verdict: "notsig",
      pred: "Four unseen simple → complex events carry more sequence-novel families than matched controls.",
      pro: "Round-1 lead at nested events (1.9–3.3×); land-plant bursts; convergent expansions in complex fungi.",
      con: "Power 18–41%; two events are single-species tips; three share one panel's controls.",
      key: "combined p = 0.050 (Holm 0.10); internal-only 2.8× [1.3–5.9], exploratory",
      lit: "Knoll 2011; Bowles et al. 2020; Merényi et al. 2020; Nagy et al. 2018; Herron & Michod 2008; Lindsey et al. 2021",
      notes: "Round-2 confirmatory test P1. Events frozen before reconstructions: Astrephomene, Volvox carteri, V. africanus/reticuliferus (Volvox is polyphyletic), Fucophycidae stem. Combined p = 0.0504 as specified (0.072 as implemented), Holm 0.10; tip artefact explains the large per-event ratios; with one effect per panel the design could not have reached significance. Internal-only exploratory estimate 2.81 (1.34-5.89) vs origins 0.82. Sources: report/REPORT.md §3, §6; DESIGN_ROUND2.md §6." },
    { title: "R2-P2 · More motif-binding domains?", img: "hyp_r2_p2.png", verdict: "null",
      pred: "Domains binding short linear motifs (SH2, SH3, PDZ, WW, 14-3-3…) are more abundant in multicellular genomes.",
      pro: "Proposed to expand with regulatory complexity; 115 families fixed from ELM before testing.",
      con: "Expansion only within animals; contraction in hyphal fungi.",
      key: "0.98-fold, p = 0.72; 6 of 7 leave-one-out estimates negative",
      lit: "Kumar et al. 2024 (ELM resource)",
      notes: "Round-2 confirmatory test P2: 115 motif-binding Pfam families from ELM; PGLS on 148 genomes; -0.028 SD (95% CI -0.124 to 0.067), p = 0.72, 0.98-fold. A Bayesian phylogenetic count model gave a small positive estimate (sensitivity; considered anti-conservative for many-family sets). Expansion visible within the animal origin and contraction in hyphal fungi: lineage-specific. Sources: report/REPORT.md §7.4; Supplement S5." },
    { title: "R2-P3 · More cell-death machinery?", img: "hyp_r2_p3.png", verdict: "supported",
      pred: "Multicellular genomes carry more programmed-cell-death proteins than expected for proteome size.",
      pro: "Cell death can benefit relatives; a joint origin with virus arms races has been proposed.",
      con: "Size-dependent; mostly immune receptors; core apoptosis machinery already in the common ancestor of eukaryotes.",
      key: "1.21× [1.04–1.42], Holm p = 0.014; mainly complex lineages",
      lit: "Durand et al. 2016; Iranzo et al. 2014; Klim et al. 2018; Jones & Dangl 2006; Dyrka et al. 2014",
      notes: "Round-2 confirmatory test P3, the one supported test: +0.167 SD, one-sided p = 0.0047, 1.21-fold (1.04-1.42), Holm 0.014; complex vs unicellular +0.249 (p = 0.0005). Robust to leaving out any origin or family; independently re-derived; random prevalence-matched sets reached p < 0.05 in 2.5% of cases. Caveats: without the ten largest proteomes p = 0.025; as share of the proteome p = 0.14; mostly lineage-specific immune and self/non-self receptors; cannot separate developmental cell death from immunity. Sources: report/REPORT.md §7.2-7.3, §10.2." },
    { title: "R2 · Structural novelty and constraint", img: "hyp_r2_struct.png", verdict: "overturned",
      pred: "Families new at events lack structural relatives; structure changes faster per sequence change in multicellular clades.",
      pro: "Folds outlast sequences; ProstT5 + Foldseek make proteome-scale structural search feasible.",
      con: "'No match' tracks short, disordered proteins; a shuffled-letter control gave the same 'effect'.",
      key: "OR 1.44 → 1.18 [0.93–1.49] after adjustment; divergence test inconclusive",
      lit: "Weisman et al. 2020; van Kempen et al. 2024; Heinzinger et al. 2024; Lynch & Conery 2003",
      notes: "Round-2 exploratory structure analyses. Structural novelty: unadjusted OR 1.44 (p = 0.004) for structurally unmatched families on event branches; overturned by adversarial verification (length, disorder, depth; adjusted 1.18, 0.93-1.49; exhaustive search found homologues for 21%). A weak simple-to-complex lead remains (family-level adjusted OR 1.73, not significant with branches as units). Structure per sequence change: +7-12% inside crowns, but shuffled control equal or larger. Sources: report/REPORT.md §8.3-8.4." },
    { title: "R2 · Protein complexes kept or duplicated?", img: "hyp_r2_complex.png", verdict: "null",
      pred: "Duplications on multicellular stems hit subunits of ancient complexes; complexes stay intact more often.",
      pro: "LECA interactome from co-fractionation mass spectrometry across 31 species.",
      con: "Map holds only ancient complexes; 13 mappable species; duplication counts depend on the method.",
      key: "66.6% (animals), 72.1% (fungi) vs 63.6% controls; none enriched after correction",
      lit: "Cox et al. 2026 (Marcotte lab; preprint Cox et al. 2024)",
      notes: "Round-2 exploratory complex analysis using only the published LECA complex map (raw co-fractionation profiles not reanalysed). Stem duplications show no consistent preference for complex subunits; co-retention slightly lower in multicellular species. Treated as a null because duplication counts are method-dependent. Sources: report/REPORT.md §9; RESEARCH_TREE.md W7." },
  ];
  HYP.forEach(hypSlide);

  // ============================================================================================
  section("H. Decision making");
  divider("H", "Decisions that transfer", "What this project did, what went wrong, and what future projects should copy",
    "Section H. Source: LESSONS.md (14 rules); agents/retrospective_R2.txt; report/REPORT.md §11.");
  {
    const s = content("From this project to future projects",
      "Transferable practices. Hashed pre-registration (3 tags) kept confirmatory tests separate from later sensitivity analyses; next time, also label tests with power below about 40% as estimation. Calibration on simulations and control data (2.1-2.6-fold gain inflation; null false-positive rate 14.5% to 4.9%) should be extended to piloting every statistic on control-only data, including tip/internal comparability, before freezing. Independent adversarial verification overturned two claims; authors should never verify their own results. External senior review steered round 2. Many small dense runs (15-45 min) beat one giant run. A second method (EukaryOMA) gave the most sobering result. GPU and long-run pilots with a device check, and a LESSONS file of 14 rules, should be carried into the next project. Source: LESSONS.md; agents/retrospective_R2.txt; report/REPORT.md §11.");
    txt(s, "This project", { x: 0.6, y: 1.3, w: 5.4, h: 0.45, fontSize: 18, bold: true, color: C.accent1 });
    txt(s, "Future projects", { x: 7.3, y: 1.3, w: 5.4, h: 0.45, fontSize: 18, bold: true, color: C.accent4 });
    const R = [["FaLock", "Hashed pre-registration, 3 tags", "Label tests with < 40% power as estimation"],
      ["FaVial", "Calibrated on simulations (2.1–2.6×)", "Pilot every statistic on control-only data"],
      ["FaCheckDouble", "Adversarial verification: 2 claims overturned", "Never let authors verify their own work"],
      ["FaUserTie", "Department-head review steered round 2", "External review at every round boundary"],
      ["FaLayerGroup", "Dense panels, 15–45 min each", "Many small runs beat one giant run"],
      ["FaExchangeAlt", "EukaryOMA cross-check (29%)", "Run a second method on key counts"],
      ["FaMicrochip", "GPU pilot caught a silent CPU fallback", "Check the device before long jobs"]];
    R.forEach((r, i) => {
      const y = 1.85 + i * 0.7;
      icon(s, r[0], 0.6, y, 0.56, C.accent1);
      box(s, 1.3, y, 4.95, 0.56, C.background2, "this-row");
      txt(s, r[1], { x: 1.45, y, w: 4.7, h: 0.56, fontSize: 14, valign: "middle" });
      arrow(s, 6.35, y + 0.28, 7.15, y + 0.28, C.accent4, 2.25);
      box(s, 7.3, y, 5.45, 0.56, C.accent1, "future-row");
      txt(s, r[2], { x: 7.45, y, w: 5.2, h: 0.56, fontSize: 14, valign: "middle", color: C.background1, bold: true });
    });
    caption(s, "Condensed from LESSONS.md (14 rules) and the fresh-agent retrospective.", { x: 0.6, y: 6.85, w: 11, h: 0.3 });
  }
  {
    const s = content("What went wrong",
      "What went wrong, briefly. Errata: round 1 overstated three statements, and a final numbers audit of the round-2 report found eight further mismatches (166 of 174 claims matched). Overturned claims: structural novelty at transitions and better cross-method agreement at stems, both already written into summaries. Tool incidents: twice a pattern-based kill command killed its own caller, once leaving a GPU job unsupervised; a wait loop launched an analysis prematurely on four of eight panels; an undocumented branch-length floor changed the confirmatory p from 0.050 to 0.072. Low power: the simple-to-complex test had 18-41% power and, read with its dependence, could never reach significance. Source: report/REPORT.md §11; RESEARCH_TREE.md R2-M4; LESSONS.md rules 5, 8, 12.");
    const T = [["FaUndo", C.accent5, "Errata", "3 round-1 overstatements; 8 mismatches in the final numbers audit"],
      ["FaTimes", C.accent4, "Overturned claims", "structural novelty and stem agreement, after verification"],
      ["FaBug", C.text2, "Tool incidents", "self-killing pattern kill (×2), orphaned GPU job, premature wait-loop run, undocumented floor"],
      ["FaHourglassHalf", C.accent3, "Low-power test", "simple → complex test: 18–41% power, could never pass Holm"]];
    T.forEach((t, i) => {
      const x = 0.6 + (i % 2) * 6.15, y = 1.45 + Math.floor(i / 2) * 2.65;
      box(s, x, y, 5.95, 2.4, C.background2, "wrong-card");
      icon(s, t[0], x + 0.25, y + 0.3, 0.9, t[1]);
      txt(s, [{ text: t[2], options: { bold: true, breakLine: true, fontSize: 20 } }, { text: t[3], options: { fontSize: 16 } }], { x: x + 1.4, y: y + 0.25, w: 4.4, h: 1.95 });
    });
    caption(s, "All were caught and documented by the process itself; none is hidden in the report.", { x: 0.6, y: 6.75, w: 12, h: 0.35 });
  }

  // ============================================================================================
  section("I. Conclusion");
  divider("I", "Conclusion", "AI-led research: outcome and future possibilities",
    "Section I. Source: report/REPORT.md §10-§11.");
  {
    const s = content("AI-led research: honest, fast, self-correcting", "AI-led research outcome. Speed: about 31 hours of wall time for two full rounds, including literature review, data assembly of hundreds of proteomes, 11 reconstructions, GPU structure prediction and two reports. Strengths: honest nulls with stated power, self-correction (verification overturned two of the team's own claims; analysis agents caught the tip artefact themselves), a complete audit trail (pre-registration, dated deviations with 'Seen:' lines, research tree). Weaknesses: confirmatory tests frozen with known low power, design problems found only after seeing data, undocumented analysis choices, lagging records and tool handling. The outcome is a careful, mostly negative study whose clearest product is methodological: where ancestral gene-content reconstruction cannot yet be trusted. Nothing has been verified by a human yet. Source: report/REPORT.md §11; agents/retrospective_R2.txt.", true);
    txt(s, "≈31 h", { x: 0.6, y: 1.5, w: 3.6, h: 1.3, fontSize: 66, bold: true, color: C.accent5 });
    txt(s, "two complete rounds, from literature to reviewed report", { x: 0.6, y: 2.9, w: 3.4, h: 1.0, fontSize: 16, color: C.background1 });
    txt(s, "Not yet verified by a human", { x: 0.6, y: 4.6, w: 3.4, h: 1.0, fontSize: 18, bold: true, color: C.accent4 });
    const good = ["Honest nulls, with power", "Self-correction: 2 own claims overturned", "Full audit trail: pre-registration, deviations", "Parallel team, cheap dense runs"];
    const bad = ["Low-power confirmatory tests", "Design flaws found after seeing data", "Undocumented choices, stale numbers", "Tool handling, lagging records"];
    [[good, "Did well", "FaCheck", C.accent2], [bad, "Did badly", "FaTimes", C.accent4]].forEach((col, k) => {
      const x = 4.55 + k * 4.2;
      icon(s, col[2], x, 1.5, 0.7, col[3]);
      txt(s, col[1], { x: x + 0.85, y: 1.5, w: 3, h: 0.7, fontSize: 22, bold: true, color: C.background1, valign: "middle" });
      col[0].forEach((g, i) => {
        box(s, x, 2.45 + i * 1.08, 3.95, 0.9, "1B5450", "cell");
        txt(s, g, { x: x + 0.2, y: 2.45 + i * 1.08, w: 3.6, h: 0.9, fontSize: 16, color: C.background1, valign: "middle" });
      });
    });
  }
  {
    const s = content("Future possibilities", "Future possibilities. Human verification of methodology and benchmarking is pending (the disclaimer will be updated). Data: EukaryOMA's 1,328 genomes and new genome sisters for transcriptome-only lineages. Method: better orthology at deep nodes and domain-level dating (cadherin-like innovations are new domain combinations inside old families); calibrating the gain bias with a gene-birth simulator (recommended, not done). Structure at scale: 3Di prediction stopped at 271 of 322 proteomes; full coverage and 3D models would allow real fold analyses. Biology: the cell-death/immune-receptor result invites experimental and population-level follow-up (developmental cell death vs immunity); regulatory change needs non-protein data. Source: report/REPORT.md §10.3-10.4; RESEARCH_TREE.md R2-M4.", true);
    const F = [["FaUserCheck", C.accent4, "Human verification", "methods and benchmarks, then update the disclaimer"],
      ["FaGlobeEurope", C.accent1, "More genomes", "EukaryOMA's 1,328 genomes; genome sisters for transcriptome-only lineages"],
      ["FaSitemap", C.accent3, "Deep-node orthology", "domain-level dating; gene-birth simulator for the bias"],
      ["FaMicrochip", C.accent2, "Structure at scale", "finish 3Di for all proteomes; real 3D models"],
      ["FaMicroscope", C.accent5, "Experimental links", "cell death vs immunity; regulatory data"]];
    F.forEach((f, i) => {
      const x = 0.6 + i * 2.45;
      box(s, x, 1.55, 2.3, 4.2, "1B5450", "future-card");
      icon(s, f[0], x + 0.7, 1.85, 0.9, f[1]);
      txt(s, [{ text: f[2], options: { bold: true, breakLine: true, fontSize: 18 } }, { text: f[3], options: { fontSize: 15 } }],
        { x: x + 0.18, y: 3.0, w: 1.95, h: 2.6, color: C.background1, align: "center" });
    });
    txt(s, "First step: a human check of methods and benchmarks; nothing here has been verified yet.", { x: 0.6, y: 6.1, w: 12, h: 0.5, fontSize: 16, bold: true, color: C.accent5 });
  }
  {
    const s = pres.addSlide({ masterName: "DARK_TITLE", sectionTitle: SECTION });
    s.addText("Read more", { placeholder: "title" });
    s.addText([
      { text: "Final report: claude.ai/artifact/9rkyi2Xwe1VyhCgUcjZuoJ", options: { breakLine: true, hyperlink: { url: "https://claude.ai/artifact/9rkyi2Xwe1VyhCgUcjZuoJ" } } },
      { text: "Department-head review: claude.ai/artifact/NfrexbqRBFmS3A97VvTqcR", options: { breakLine: true, hyperlink: { url: "https://claude.ai/artifact/NfrexbqRBFmS3A97VvTqcR" } } },
      { text: "Repository: report/REPORT.md · RESEARCH_TREE.md · LESSONS.md · reviews/", options: { breakLine: true, fontSize: 16 } },
      { text: "Hackathon: github.com/evocellnet/ai-hackathon · [repository URL: ____]", options: { fontSize: 16 } },
    ], { placeholder: "body", fontSize: 18 });
    s.addNotes("Links. Final report page: https://claude.ai/artifact/9rkyi2Xwe1VyhCgUcjZuoJ. Department-head review of round 1: https://claude.ai/artifact/NfrexbqRBFmS3A97VvTqcR. Repository paths: report/REPORT.md (final report), report/REPORT_SUPPLEMENT.md, report/REPORT_round1.md, RESEARCH_TREE.md, LESSONS.md, LITERATURE.md, reviews/. The public GitHub URL of this project's repository is left blank for the PI. Source: RESEARCH_TREE.md R2-M4; report/REPORT.md Data and code availability.");
  }

  await pres.writeFile({ fileName: OUTFILE });
  await applyTheme(OUTFILE, THEME);
  console.log("wrote", OUTFILE);
}

build().catch((e) => { console.error(e); process.exit(1); });
