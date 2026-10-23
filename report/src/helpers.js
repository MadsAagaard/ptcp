// ── Helpers (germline template primitives) ───────────────────────────────────
function esc(s) {
  if (s == null) return '—';
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}
// keepCase: genotypes such as αα/αα must not be upper-cased to AA/AA
function badge(text, cls, sm, keepCase) { return `<span class="badge badge-${cls}${sm ? ' sm' : ''}${keepCase ? ' kc' : ''}">${esc(text)}</span>`; }
function kvi(label, value, sub, subCls) {
  return `<div class="kv-item"><div class="kv-label">${label}</div><div class="kv-value">${value}</div>
    ${sub ? `<div class="kv-sub${subCls ? ' ' + subCls : ''}">${sub}</div>` : ''}</div>`;
}
function kviL(label, value, sub, subCls, valCls) {
  return `<div class="kv-item"><div class="kv-label">${label}</div>
    <div class="kv-value large${valCls ? ' ' + valCls : ''}">${value}</div>
    ${sub ? `<div class="kv-sub${subCls ? ' ' + subCls : ''}">${sub}</div>` : ''}</div>`;
}
function src(key) {
  const s = DATA._sources && DATA._sources[key];
  return s ? `&#9432;&nbsp;${esc(s)}` : '';
}
function noteRow(kind, source, body) {
  const icon = kind === 'caution' ? '&#9888;' : '&#9432;';
  return `<div class="note-row ${kind}"><div class="note-icon">${icon}</div>
    <div class="note-src">${esc(source)}</div><div class="note-body">${body}</div></div>`;
}
function notesCard(rows, prefix) {
  if (!rows.length) return '';
  const nCaution = rows.filter(r => r.indexOf('note-row caution') !== -1).length;
  return `<div class="card notes-card">
    <button class="reg-toggle" id="${prefix}-toggle-0" onclick="toggleSubsec('${prefix}', 0)">
      <span class="arrow">&#9654;</span>
      <span class="compare-header-title">Notes &amp; Caveats</span>
      <span class="compare-header-sub">${rows.length} item${rows.length > 1 ? 's' : ''}</span>
      <span class="spacer"></span>
      ${nCaution ? badge(`${nCaution} caution`, 'warn') : badge('info only', 'nd')}
    </button>
    <div class="subsec-body" id="${prefix}-body-0">${rows.join('')}</div>
  </div>`;
}
function toggleSubsec(prefix, idx) {
  const toggle = document.getElementById(prefix + '-toggle-' + idx);
  const body   = document.getElementById(prefix + '-body-'   + idx);
  if (!toggle || !body) return;
  const open = body.classList.toggle('visible');
  toggle.classList.toggle('open', open);
}
function interp(lede, ledeCls, body) {
  return `<div class="interp-body"><span class="interp-lede ${ledeCls}">${lede}</span>${body}</div>`;
}
function conclusion(title, call) {
  return `<div class="compare-header" style="border-top:1px solid var(--border)">
      <span class="compare-header-title">Conclusion</span>
      <span class="compare-header-sub">${title}</span><span class="spacer"></span>${badge(call.txt, call.cls, false, call.kc)}
    </div>${interp('Interpretation', call.lede, call.body)}`;
}
function num(v, d) { return v == null ? '—' : Number(v).toFixed(d == null ? 1 : d); }
function grp(v)  { return v == null ? '—' : Number(v).toLocaleString('en'); }
function pct(f, d) { return f == null ? '—' : `${(100 * f).toFixed(d == null ? 1 : d)}%`; }
// ptcp-qc reports read accuracy in percent (99.995); readers think in Q.
function qOf(acc) {
  if (acc == null) return null;
  const err = 1 - acc / 100;
  return err <= 0 ? 60 : Math.min(60, -10 * Math.log10(err));
}
function median(a) {
  const s = a.filter(v => v != null).slice().sort((x, y) => x - y);
  if (!s.length) return null;
  const m = Math.floor(s.length / 2);
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}

// ── Thresholds (configuration; placeholders until approved) ──────────────────
const THRESHOLDS = (DATA.thresholds && typeof DATA.thresholds === 'object') ? DATA.thresholds : {};
const TH_APPROVED = !!(THRESHOLDS.approved_by && THRESHOLDS.approved_date);
function thresholdsFor(section, builtins) { return Object.assign({}, builtins, THRESHOLDS[section] || {}); }
function thresholdNote(criteria) {
  return TH_APPROVED
    ? noteRow('info', 'thresholds', `Acceptance criteria (${criteria}) are this laboratory's own, read from <code>${esc(THRESHOLDS.source_file || 'the thresholds file')}</code> and approved by <b>${esc(THRESHOLDS.approved_by)}</b> on <b>${esc(THRESHOLDS.approved_date)}</b>.`)
    : noteRow('caution', 'thresholds', `The thresholds shown (${criteria}) are <b>placeholders</b>, not this laboratory's acceptance criteria. They become criteria once the collector's thresholds file carries <code>approved_by</code> and <code>approved_date</code>.`);
}
const TQC  = thresholdsFor('qc',  { locus_reads_fail: 50, locus_reads_warn: 150, on_target_frac_warn: 0.10, read_q_warn: 30, share_ratio_warn: 0.5 });
const TSEX = thresholdsFor('sex', { y_xy_min: 0.15, y_xx_max: 0.05, x_xx_min: 0.75, x_xy_max: 0.65 });
const TPP  = thresholdsFor('paraphase', {
  smn1_expected: 2,
  // pseudogene-like when ANY limit is reached (see section 03)
  hap_class: { gba: { min_del: 500, min_patho: 10, min_pseudo: 10 }, cyp21: { min_del: 64, min_patho: 10, min_pseudo: 10 }, hbb: null, rpgr: null },
  total_cn_expected: { gba: 4, cyp21: 4, hbb: 2, rpgr: { M: 1, K: 2 } },
});
// Repeat classes in UNITS of the counted motif(s): count = indexes into TRGT
// motif_count. FMR1 CGG:AGG counts both (clinical sizing), FXN A:GAA counts GAA.
const TTR  = thresholdsFor('trgt', {
  FXS_FMR1:   { gene: 'FMR1', unit: 'CGG', count: [0, 1], interruption: 1, xlinked: true,
                classes: [[0, 44, 'normal', 'pass'], [45, 54, 'intermediate', 'info'], [55, 200, 'premutation', 'warn'], [201, 1e9, 'full mutation', 'fail']] },
  FRDA_FXN:   { gene: 'FXN',  unit: 'GAA', count: [1], xlinked: false,
                classes: [[0, 33, 'normal', 'pass'], [34, 65, 'premutation', 'warn'], [66, 1e9, 'full penetrance', 'fail']] },
  FRAXE_AFF2: { gene: 'AFF2', unit: 'GCC', count: [0], xlinked: true,
                classes: [[0, 39, 'normal', 'pass'], [40, 60, 'intermediate', 'info'], [61, 200, 'premutation', 'warn'], [201, 1e9, 'full mutation', 'fail']] },
  EIEE1_ARX:  { gene: 'ARX', tract: 'PA1', unit: 'GCN', what: 'polyalanine tract 1', count: [0], xlinked: true,
                classes: [[0, 16, 'normal', 'pass'], [17, 1e9, 'expanded', 'warn']] },
  PRTS_ARX:   { gene: 'ARX', tract: 'PA2', unit: 'GCN', what: 'polyalanine tract 2', count: [0], xlinked: true,
                classes: [[0, 12, 'normal', 'pass'], [13, 1e9, 'expanded', 'warn']] },
});
Object.entries(TTR).forEach(([k, c]) => { if (c.tract) c.gene = `${c.gene} ${c.tract}`; });

// Carrier targets (ptcp-qc locus groups), in display order.
const TARGETS = THRESHOLDS.carrier_targets || [
  ['smn1', 'SMN1/2'], ['gba', 'GBA1'], ['cyp21', 'CYP21A2'], ['hba', 'HBA1/2'], ['hbb', 'HBB'], ['rpgr', 'RPGR'],
  ['f8inv1', 'F8 int1'], ['f8inv22', 'F8 int22'], ['FXS_FMR1', 'FMR1'], ['FRDA_FXN', 'FXN'], ['AFF2', 'AFF2'], ['ARX', 'ARX'],
];

// ── Sample QC verdict (shared by both reports) ───────────────────────────────
const VERDICT_TXT = { pass: 'pass', warn: 'review', fail: 'fail', nd: 'no data' };
const VERDICT_CLS = { pass: 'pass', warn: 'warn', fail: 'fail', nd: 'gated' };
function locusState(n) { return n == null || n < TQC.locus_reads_fail ? 'fail' : n < TQC.locus_reads_warn ? 'warn' : 'ok'; }
function sexCheck(ratios, sexReported) {
  const r = ratios || {};
  let data = null;
  if (r.ratio_y_auto != null && r.ratio_x_auto != null) {
    if (r.ratio_y_auto >= TSEX.y_xy_min && r.ratio_x_auto <= TSEX.x_xy_max) data = 'M';
    else if (r.ratio_y_auto <= TSEX.y_xx_max && r.ratio_x_auto >= TSEX.x_xx_min) data = 'K';
    else data = '?';
  }
  const state = data == null ? 'nd' : data === '?' ? 'warn' : data === sexReported ? 'ok' : 'fail';
  return { data, state, r };
}
// qc: { total_reads, on_target{total}, on_target_fraction, read_quality_median, sex_ratios }
// loci: { group: n } or { group: { n } }
// shares / compTargets (optional): the sample's target shares and the comparison's
// { target: { median } }; a share below share_ratio_warn x the median adds a review issue.
function sampleQC(qc, loci, sexReported, shares, compTargets, compLabel) {
  const q = qc || {};
  const nOf = v => v == null ? null : (typeof v === 'object' ? v.n : v);
  const loc = TARGETS.map(([k]) => nOf((loci || {})[k]));
  const sx = sexCheck(q.sex_ratios, sexReported);
  const qv = qOf(q.read_quality_median);
  const onFrac = q.on_target_fraction != null ? q.on_target_fraction
               : (q.total_reads && q.on_target ? q.on_target.total / q.total_reads : null);
  if (!loci || !Object.keys(loci).length) return { verdict: 'nd', issues: ['no ptcp-qc sample JSON'], minLoc: null, loc, sx, qv, onFrac };
  const minLoc = Math.min(...loc.map(v => v == null ? 0 : v));
  const issues = [];
  let verdict = 'pass';
  const lower = v => { if (verdict === 'pass') verdict = v; };
  if (sx.state === 'fail') { verdict = 'fail'; issues.push('sex mismatch'); }
  const low = TARGETS.filter((t, i) => loc[i] == null || loc[i] < TQC.locus_reads_fail).map(t => t[1]);
  const weak = TARGETS.filter((t, i) => loc[i] != null && loc[i] >= TQC.locus_reads_fail && loc[i] < TQC.locus_reads_warn).map(t => t[1]);
  if (low.length) { verdict = 'fail'; issues.push(`< ${TQC.locus_reads_fail} reads: ${low.join(', ')}`); }
  if (weak.length) { lower('warn'); issues.push(`< ${TQC.locus_reads_warn} reads: ${weak.join(', ')}`); }
  if (onFrac != null && onFrac < TQC.on_target_frac_warn) { lower('warn'); issues.push('low on-target'); }
  if (qv != null && qv < TQC.read_q_warn) { lower('warn'); issues.push('low read Q'); }
  if (shares && compTargets) {
    const lowShare = TARGETS.filter(([k]) => compTargets[k] && compTargets[k].median && shares[k] != null
                                         && shares[k] / compTargets[k].median < TQC.share_ratio_warn).map(t => t[1]);
    if (lowShare.length) { lower('warn'); issues.push(`low share vs ${compLabel || 'comparison'}: ${lowShare.join(', ')}`); }
  }
  if (sx.state === 'warn') { lower('warn'); issues.push('sex ambiguous'); }
  if (sx.state === 'nd') issues.push('no sex ratios');
  return { verdict, issues, minLoc, loc, sx, qv, onFrac };
}
