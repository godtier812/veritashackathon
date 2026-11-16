(function () {
  var toolurl = /\/explorer\/?$/.test(location.pathname) ? "../tool/" : "";
  var order = "ACDEFGHIKLMNPQRSTVWY";
  var names = {
    A: "alanine", C: "cysteine", D: "aspartate", E: "glutamate", F: "phenylalanine", G: "glycine",
    H: "histidine", I: "isoleucine", K: "lysine", L: "leucine", M: "methionine", N: "asparagine",
    P: "proline", Q: "glutamine", R: "arginine", S: "serine", T: "threonine", V: "valine",
    W: "tryptophan", Y: "tyrosine"
  };
  var na = "not available";
  var demonote = "At lab-tested positions a counselor would use the lab result. This demo shows how the model does where we can check it.";
  var colors = { harmful: 0xff7a66, cant: 0xf2b94b, harmless: 0x3fcfb8, dim: 0x3b4a66, bright: 0xf2f6ff, off: 0x3a4659, zinc: 0xc8d3e6 };
  var icons = {
    harmful: '<path d="M8 2.5 14.2 13.5H1.8Z"/>',
    harmless: '<circle cx="8" cy="8" r="5.6"/>',
    cant: '<path d="M8 1.8 14.2 8 8 14.2 1.8 8Z"/>',
    untested: '<circle cx="8" cy="8" r="5.6" stroke-dasharray="2.2 2.2"/>',
    agree: '<path d="M2.8 8.6 6.4 12.2 13.2 4.4"/>',
    miss: '<path d="M3.5 3.5 12.5 12.5M12.5 3.5 3.5 12.5"/>',
    alarm: '<path d="M3.5 3.5 12.5 12.5M12.5 3.5 3.5 12.5"/>',
    grey: '<path d="M2.5 8h11"/>'
  };
  var agreeinfo = {
    untested: ["Not lab-tested", "Nobody has measured this exact swap in the lab, so there is no answer key to check the model against."],
    cant: ["Model can't tell", "The score landed in the middle zone, where the model is wrong too often to make a call. Here the lab result is the answer."],
    grey: ["Lab result in between", "The lab measured something between broken and working, so the model's call can't be marked right or wrong."],
    agree: ["Model and lab agree", "The model's call and the lab's answer point the same way. That is one match, not a promise about the next one."],
    miss: ["Model missed this one", "The lab found this change breaks the protein, but the model said it looked fine. This is the kind of mistake that matters most."],
    alarm: ["Model raised a false alarm", "The model flagged this change, but in the lab the protein kept working normally."]
  };
  var verdictplain = {};
  var verdictkey = {};

  var reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  var data, cfg, cut, bypos = {}, positions = [], swaps = [], stats = {};
  var state = { region: "start", mode: "model", filter: "all", text: "", sel: null, auto: false };
  var viewer = null, models = [], loaded = null, pending = null, frame = 0, flextimer = null, spinpause = null, fallback = false;
  var zincpairs = [];
  var structcache = {};
  var mouse = { x: 0, y: 0 };
  var playtimer = null;

  function $(s) { return document.querySelector(s); }
  function esc(t) { return String(t).replace(/[&<>"']/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#x27;" }[c]; }); }
  function minus(t) { return t.charAt(0) === "-" ? "−" + t.slice(1) : t; }
  function fmt(v, d) { return v === null || v === undefined || isNaN(v) ? na : minus(Number(v).toFixed(d)); }
  function icon(k, cls) { return '<svg class="ico ' + (cls || "") + '" viewBox="0 0 16 16" aria-hidden="true">' + icons[k] + "</svg>"; }
  function hex(n) { return "#" + ("000000" + n.toString(16)).slice(-6); }
  function mix(a, b, t) {
    t = Math.max(0, Math.min(1, t));
    var r = ((a >> 16) & 255) + (((b >> 16) & 255) - ((a >> 16) & 255)) * t;
    var g = ((a >> 8) & 255) + (((b >> 8) & 255) - ((a >> 8) & 255)) * t;
    var bl = (a & 255) + ((b & 255) - (a & 255)) * t;
    return (Math.round(r) << 16) | (Math.round(g) << 8) | Math.round(bl);
  }
  function three(lo, hi, v) {
    if (v === null || v === undefined) return colors.off;
    var mid = (lo + hi) / 2;
    if (v <= lo) return colors.harmful;
    if (v >= hi) return colors.harmless;
    if (v < mid) return mix(colors.harmful, colors.cant, (v - lo) / (mid - lo));
    return mix(colors.cant, colors.harmless, (v - mid) / (hi - mid));
  }
  function modelcolor(v) { return three(-cut.b2, -cut.a, v); }
  function labcolor(v) { return three(cfg.lof_below, cfg.func_above, v); }
  function discolor(v) { return v === null ? colors.off : mix(colors.dim, colors.bright, v); }
  function poscolor(p, mode) {
    var s = stats[p];
    if (!s) return colors.off;
    if (mode === "lab") return labcolor(s.meanlab);
    if (mode === "disagree") return discolor(s.frac);
    return modelcolor(s.meanesm);
  }
  function inregion(p, r) { var b = cfg[r]; return p >= b[0] && p <= b[1]; }
  function regionof(p) { return inregion(p, "start") ? "start" : "end"; }

  function prepare() {
    cfg = data.config;
    var c = rules.cutoffs(cfg);
    cut = { a: c.a, b: c.b, b2: c.b === null ? c.a : c.b };
    verdictplain[rules.harmful] = "The model thinks this change probably breaks the protein.";
    verdictplain[rules.harmless] = "The model thinks the protein probably copes with this change.";
    verdictplain[rules.cant] = "The model is not sure about this one, so it does not make a call.";
    verdictkey[rules.harmful] = "harmful"; verdictkey[rules.harmless] = "harmless"; verdictkey[rules.cant] = "cant";
    data.positions.forEach(function (p) {
      bypos[p.p] = p;
      positions.push(p.p);
      var se = 0, sl = 0, n = 0, dis = 0;
      p.swaps.forEach(function (s) {
        var v = rules.verdict(s[1], cfg);
        var item = { p: p.p, wt: p.wt, nw: s[0], esm: s[1], lc: s[2], ls: s[3], bl: s[4], v: v, a: rules.agreement(v, s[2]) };
        item.name = p.wt + p.p + s[0];
        s.item = item;
        if (s[2]) {
          swaps.push(item);
          se += s[1]; sl += s[3]; n += 1;
          if (item.a === "miss" || item.a === "alarm") dis += 1;
        }
      });
      stats[p.p] = { n: n, dis: dis, meanesm: n ? se / n : null, meanlab: n ? sl / n : null, frac: n ? dis / n : null };
    });
  }

  function swapat(p, nw) {
    var e = bypos[p];
    if (!e) return null;
    for (var i = 0; i < e.swaps.length; i++) if (e.swaps[i][0] === nw) return e.swaps[i].item;
    return null;
  }
  function defaultswap(p) {
    var e = bypos[p], best = null;
    e.swaps.forEach(function (s) { if (s[2] && (!best || s[1] < best.esm)) best = s.item; });
    return best || e.swaps[0].item;
  }

  function select(item, opts) {
    opts = opts || {};
    state.sel = item;
    state.auto = !!opts.auto;
    if (toolurl) $("#tool-link").href = toolurl + "?q=" + encodeURIComponent(item.name);
    try { history.replaceState(null, "", location.pathname + "?m=" + encodeURIComponent(item.name)); } catch (err) {}
    var r = regionof(item.p);
    if (r !== state.region) setregion(r, true);
    renderreadout(opts);
    marklist(opts.fromlist);
    drawtrack();
    if (fallback) drawribbon();
    else paint(opts.quick ? 120 : 900);
  }
  function selectpos(p, opts) {
    var cur = state.sel && swapat(p, state.sel.nw);
    if (cur && cur.lc) select(cur, opts);
    else select(defaultswap(p), Object.assign({ auto: true }, opts || {}));
  }

  function countup(el, to, d) {
    if (to === null || to === undefined) { el.textContent = na; return; }
    if (reduced) { el.textContent = fmt(to, d); return; }
    var t0 = performance.now(), dur = 600, done = false;
    setTimeout(function () { done = true; el.textContent = fmt(to, d); }, dur + 80);
    function step(t) {
      if (done) return;
      var k = Math.max(0, Math.min(1, (t - t0) / dur));
      var e = 1 - Math.pow(1 - k, 3);
      el.textContent = fmt(to * e, d);
      if (k < 1) requestAnimationFrame(step); else el.textContent = fmt(to, d);
    }
    requestAnimationFrame(step);
  }

  function pct(v, lo, hi) { return Math.max(0, Math.min(100, (v - lo) / (hi - lo) * 100)); }

  var esmlo, esmhi, lablo, labhi, absmin, rankesm, ranklab;
  function share(arr, v) {
    var lo = 0, hi = arr.length;
    while (lo < hi) { var m = (lo + hi) >> 1; if (arr[m] < v) lo = m + 1; else hi = m; }
    var left = lo; hi = arr.length;
    while (lo < hi) { var n = (lo + hi) >> 1; if (arr[n] <= v) lo = n + 1; else hi = n; }
    return (arr.length - lo + 0.5 * (lo - left)) / arr.length * 100;
  }
  function rgauge(title, sh, raw, c1, c2, z) {
    var a = 100 - c1, b = 100 - c2;
    var h = '<div class="gauge"><div class="lab"><span>' + title + '</span><span class="num">' + (sh === null ? "not lab-tested" : "more harmful than " + Math.round(sh) + "%") + "</span></div>";
    h += '<div class="gtrack">';
    h += '<i class="z" style="width:' + a + "%;background:" + z[0] + '"></i>';
    h += '<i class="z" style="width:' + (b - a) + "%;background:" + z[1] + '"></i>';
    h += '<i class="z" style="width:' + (100 - b) + "%;background:" + z[2] + '"></i>';
    if (sh !== null) h += '<i class="mk" style="left:' + (100 - sh) + '%"></i>';
    return h + "</div></div>";
  }
  function ranges() {
    var e = [], l = [];
    rankesm = []; ranklab = [];
    data.positions.forEach(function (p) { p.swaps.forEach(function (s) { e.push(s[1]); if (s[3] !== null) l.push(s[3]); if (s[2]) { rankesm.push(s[1]); ranklab.push(s[3]); } }); });
    rankesm.sort(function (x, y) { return x - y; }); ranklab.sort(function (x, y) { return x - y; });
    esmlo = Math.floor(Math.min.apply(null, e)) - 1; esmhi = Math.ceil(Math.max.apply(null, e)) + 1;
    lablo = Math.floor(Math.min.apply(null, l)) - 1; labhi = Math.ceil(Math.max.apply(null, l)) + 1;
    absmin = Math.min.apply(null, e);
  }

  function structinfo(p) {
    var s = data.structures[regionof(p)];
    if (!s) return { has: false, id: "" };
    return { has: s.missing.indexOf(p) < 0 && p >= s.first && p <= s.last, id: s.id };
  }
  function zincat(p) { return data.structures.start.zinc_positions.indexOf(p) >= 0; }

  function renderreadout(opts) {
    var it = state.sel, el = $("#readout");
    var vk = verdictkey[it.v];
    var ai = agreeinfo[it.a];
    var si = structinfo(it.p);
    var arm1 = data.arms.filter(function (a) { return a.arm === "arm1"; })[0];
    var bltext = it.bl === null || it.bl === undefined ? na : String(it.bl);
    if (arm1 && arm1.cutoff !== null && it.bl !== null) {
      var harm = 0 - it.bl;
      bltext += " · the BLOSUM62 rule flags harm (−BLOSUM62) at or above " + arm1.cutoff.toFixed(1) + ", so it calls this " + (harm >= arm1.cutoff ? "harmful" : "harmless") + ".";
    }
    var labtext = it.lc ? ({ LOF: "LOF: loss of function in the lab assay", FUNC: "FUNC: the protein kept working in the lab assay", INT: "INT: in between" }[it.lc]) : "not lab-tested";
    var regtext = regionof(it.p) === "start"
      ? "Start region (RING, 1–101). Practice: the two cutoffs were set here."
      : "End region (BRCT, 1631–1855). Sealed exam: opened once, after the cutoffs were frozen.";
    var zinc = zincat(it.p) ? "<dt>Zinc</dt><dd>In 1JM7 this " + names[it.wt] + "'s " + (it.wt === "H" ? "nitrogen" : "sulfur") + " sits within 2.6 Å of a zinc ion (model 1).</dd>" : "";
    var h = "";
    h += '<div class="ro-name"><h1>' + esc(it.name) + "</h1><span>" + cap(names[it.wt]) + " to " + names[it.nw] + " at position " + it.p + "</span></div>";
    if (state.auto) h += '<p class="strip-note">Picked the lab-tested swap the model found most surprising at this position. Choose another letter below.</p>';
    h += '<div class="ro-verdict st-' + vk + '">' + icon(vk) + "<div><b>The model's call: " + esc(it.v) + "</b><p>" + verdictplain[it.v] + "</p></div></div>";
    h += '<div class="bignums"><div class="bignum"><div class="k">Model score</div><div class="n" id="n-esm"></div></div>';
    h += '<div class="bignum"><div class="k">Lab function score</div><div class="n' + (it.lc ? "" : " na") + '" id="n-lab"></div></div></div>';
    h += '<div class="badge st-' + it.a + '">' + icon(it.a) + "<div><b>" + ai[0] + "</b><p>" + ai[1] + "</p></div></div>";
    h += '<p class="strip-note">Both bars use the same scale: where this swap ranks among the ' + rankesm.length.toLocaleString("en-US") + " lab-tested swaps (LOF, INT and FUNC; ties count half). Top: the model's prediction. Bottom: the lab measurement. Further left means more harmful. If the model's ranking matched the lab's, the markers would line up.</p>";
    h += rgauge("Model: how harmful it looks", share(rankesm, it.esm), "score " + fmt(it.esm, 2), share(rankesm, -cut.b2), share(rankesm, -cut.a),
      ["rgba(255,122,102,.55)", "rgba(242,185,75,.45)", "rgba(63,207,184,.45)"]);
    h += rgauge("Lab: how harmful it was", it.lc ? share(ranklab, it.ls) : null, it.lc ? "score " + fmt(it.ls, 3) : "", share(ranklab, cfg.lof_below), share(ranklab, cfg.func_above),
      ["rgba(255,122,102,.55)", "rgba(125,138,161,.45)", "rgba(63,207,184,.45)"]);
    h += "<details><summary>Show the numbers</summary><dl>";
    h += "<dt>Cutoffs</dt><dd>Harmful at or below <span class=num>" + fmt(-cut.b2, 2) + "</span>; harmless above <span class=num>" + fmt(-cut.a, 2) + "</span>; can't tell in between. Set on the start region only.</dd>";
    h += "<dt>Region</dt><dd>" + regtext + "</dd>";
    h += "<dt>Lab class</dt><dd>" + labtext + '<br><span class="strip-note">LOF below ' + fmt(cfg.lof_below, 3) + " · INT " + fmt(cfg.lof_below, 3) + " to " + fmt(cfg.func_above, 3) + " · FUNC above " + fmt(cfg.func_above, 3) + " (Findlay et al. 2018)</span></dd>";
    h += "<dt>BLOSUM62</dt><dd><span class=num>" + bltext + "</span></dd>";
    h += "<dt>Entropy</dt><dd><span class=num>" + fmt(bypos[it.p].entropy, 3) + "</span> · how spread out the model's guesses were at this position. Higher means more spread out.</dd>";
    h += "<dt>3D</dt><dd>" + (si.has ? "Shown in PDB " + si.id + "." : "Not in the 3D structure.") + "</dd>";
    h += zinc + "</dl></details>";
    h += '<p class="strip-head">Position ' + it.p + ": every possible swap, with the normal letter hidden</p>";
    h += '<div class="strip" id="strip"></div>';
    h += '<p class="strip-note">Bar = log P(new letter) − log P(normal letter), from esm_scores.csv. These are log-ratios, not probabilities; longer bars down mean more surprising. A dot marks a swap the lab tested.</p>';
    h += '<p class="note">' + demonote + "</p>";
    el.innerHTML = h;
    countup($("#n-esm"), it.esm, 2);
    if (it.lc) countup($("#n-lab"), it.ls, 3); else $("#n-lab").textContent = "not lab-tested";
    renderstrip(opts);
  }
  function cap(s) { return s.charAt(0).toUpperCase() + s.slice(1); }

  function renderstrip() {
    var it = state.sel, e = bypos[it.p], box = $("#strip");
    var k = 62 / Math.abs(absmin);
    var h = "";
    for (var i = 0; i < order.length; i++) {
      var l = order[i];
      if (l === e.wt) {
        h += '<button type="button" class="wt" disabled aria-label="normal letter ' + l + '"><span class="barwrap"><span class="mask">?</span></span><span class="let" data-wt="' + l + '">?</span><span class="dot none"></span></button>';
        continue;
      }
      var s = swapat(it.p, l);
      var hgt = Math.max(2, Math.abs(Math.min(0, s.esm)) * k);
      var up = s.esm > 0 ? Math.min(10, s.esm * k) : 0;
      var col = hex(modelcolor(s.esm));
      h += '<button type="button" data-nw="' + l + '" class="' + (l === it.nw ? "sel" : "") + '" title="' + esc(s.name) + ": " + fmt(s.esm, 2) + (s.lc ? ", lab " + s.lc : ", not lab-tested") + '">';
      h += '<span class="barwrap"><span class="bar" style="margin-top:' + (10 - up) + "px;height:" + (s.esm > 0 ? up : hgt) + "px;background:" + col + ';transform:scaleY(' + (reduced ? 1 : 0) + ')"></span></span>';
      h += '<span class="let">' + l + '</span><span class="dot' + (s.lc ? "" : " none") + '"></span></button>';
    }
    box.innerHTML = h;
    box.querySelectorAll("button[data-nw]").forEach(function (b) {
      b.addEventListener("click", function () { select(swapat(it.p, b.getAttribute("data-nw")), { quick: true }); });
    });
    var reveal = function () {
      box.querySelectorAll(".bar").forEach(function (b) { b.style.transform = "scaleY(1)"; });
      var w = box.querySelector("[data-wt]");
      if (w) w.textContent = w.getAttribute("data-wt");
      var m = box.querySelector(".mask");
      if (m) m.textContent = "";
    };
    if (reduced) reveal(); else { requestAnimationFrame(function () { requestAnimationFrame(function () { box.querySelectorAll(".bar").forEach(function (b) { b.style.transform = "scaleY(1)"; }); }); }); setTimeout(reveal, 650); }
  }

  function renderlist() {
    var q = state.text.trim().toUpperCase().replace(/^P\./, "");
    var shown = swaps.filter(function (s) {
      if (state.filter !== "all" && s.a !== state.filter) return false;
      if (q && s.name.indexOf(q) < 0) return false;
      return true;
    });
    var ul = $("#list"), frag = document.createDocumentFragment();
    shown.forEach(function (s) {
      var li = document.createElement("li");
      li.innerHTML = '<button type="button" data-k="' + s.name + '"><span class="st-' + s.a + '">' + icon(s.a) + '</span><span class="mut">' + s.name + '</span><span class="v">' + fmt(s.esm, 2) + '</span><span class="v">' + fmt(s.ls, 3) + "</span></button>";
      li.firstChild.addEventListener("click", function () { select(s, { fromlist: true }); });
      li.firstChild.setAttribute("aria-label", s.name + ", " + agreeinfo[s.a][0]);
      frag.appendChild(li);
    });
    ul.innerHTML = "";
    ul.appendChild(frag);
    $("#list-count").textContent = "Showing " + shown.length + " of " + swaps.length + " lab-tested swaps · model score · lab score";
    marklist(true);
  }
  function marklist(noscroll) {
    var cur = document.querySelector(".list button.sel");
    if (cur) cur.classList.remove("sel");
    if (!state.sel) return;
    var b = document.querySelector('.list button[data-k="' + state.sel.name + '"]');
    if (b) { b.classList.add("sel"); if (!noscroll) b.scrollIntoView({ block: "nearest" }); }
  }

  function deeplink() {
    var m = new URLSearchParams(location.search).get("m");
    if (!m) return null;
    $("#search").value = m.slice(0, 40);
    var before = state.sel;
    var r = rules.parse(m, cfg.protein_length);
    if (r.error) { $("#search-msg").textContent = r.error; return null; }
    var e = bypos[r.p];
    if (!e || e.wt !== r.wt) { search(m); state.sel = before; return null; }
    return swapat(r.p, r.nw);
  }

  function search(text) {
    var msg = $("#search-msg");
    var r = rules.parse(text, cfg.protein_length);
    if (r.error) { msg.textContent = r.error; return; }
    var e = bypos[r.p];
    if (!e) { msg.textContent = "Position " + r.p + " is not one of the " + cfg.n_tested_positions + " lab-tested positions (1 to 101 and 1631 to 1855). This demo only answers where its answers can be checked against the lab, so it will not guess here."; return; }
    if (inregion(r.p, "end") && !cfg.end_labels_included) { msg.textContent = "End-region answers (positions 1631 to 1855) are hidden until the one test run is done."; return; }
    if (e.wt !== r.wt) {
      msg.textContent = e.wt === r.nw ? "BRCA1 has " + e.wt + " at position " + r.p + ", not " + r.wt + "."
        : "BRCA1 has " + e.wt + " at position " + r.p + ", not " + r.wt + ". Did you mean " + e.wt + r.p + r.nw + "?";
      return;
    }
    msg.textContent = "";
    select(swapat(r.p, r.nw), { fromlist: false });
  }

  var track, tctx, geo = {};
  function setuptrack() {
    track = $("#track");
    tctx = track.getContext("2d");
    var drag = false;
    function at(ev) {
      var r = track.getBoundingClientRect();
      var x = ev.clientX - r.left;
      var best = null, bd = 1e9;
      positions.forEach(function (p) { var d = Math.abs(geo.x(p) - x); if (d < bd) { bd = d; best = p; } });
      return best;
    }
    track.addEventListener("pointerdown", function (ev) { stopplay(); drag = true; track.setPointerCapture(ev.pointerId); selectpos(at(ev), { quick: true }); });
    track.addEventListener("pointermove", function (ev) {
      if (!drag) return;
      var p = at(ev);
      if (state.sel && p !== state.sel.p) selectpos(p, { quick: true });
    });
    track.addEventListener("pointerup", function () { drag = false; });
    track.addEventListener("keydown", function (ev) {
      if (!state.sel) return;
      var i = positions.indexOf(state.sel.p);
      if (ev.key === "ArrowRight" && i < positions.length - 1) { selectpos(positions[i + 1], { quick: true }); ev.preventDefault(); }
      if (ev.key === "ArrowLeft" && i > 0) { selectpos(positions[i - 1], { quick: true }); ev.preventDefault(); }
    });
  }
  function drawtrack() {
    if (!track) return;
    var dpr = window.devicePixelRatio || 1;
    var w = track.clientWidth, h = track.clientHeight;
    if (track.width !== Math.round(w * dpr) || track.height !== Math.round(h * dpr)) { track.width = Math.round(w * dpr); track.height = Math.round(h * dpr); }
    var c = tctx;
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.clearRect(0, 0, w, h);
    var gap = 34, n1 = 0, n2 = 0;
    positions.forEach(function (p) { if (regionof(p) === "start") n1++; else n2++; });
    var cw = (w - gap) / (n1 + n2);
    var idx = {};
    positions.forEach(function (p, i) { idx[p] = i; });
    geo.x = function (p) { var i = idx[p]; return i < n1 ? (i + 0.5) * cw : gap + (i + 0.5) * cw; };
    var top = 26, rh = 16, rg = 4;
    var rows = [
      function (p) { return modelcolor(stats[p].meanesm); },
      function (p) { return labcolor(stats[p].meanlab); },
      function (p) { return discolor(stats[p].frac); }
    ];
    positions.forEach(function (p) {
      var x = geo.x(p) - cw / 2;
      var dim = regionof(p) !== state.region ? 0.55 : 1;
      c.globalAlpha = dim;
      rows.forEach(function (f, r) { c.fillStyle = hex(f(p)); c.fillRect(x, top + r * (rh + rg), Math.max(1, cw - (cw > 3 ? 0.6 : 0)), rh); });
    });
    c.globalAlpha = 1;
    var bx = n1 * cw;
    c.strokeStyle = "rgba(160,185,225,.5)";
    c.lineWidth = 1;
    for (var k = 0; k < 2; k++) {
      c.beginPath(); c.moveTo(bx + 11 + k * 8, top + 3 * (rh + rg)); c.lineTo(bx + 17 + k * 8, top - 4); c.stroke();
    }
    c.fillStyle = "#8090a8";
    c.font = "10px " + getComputedStyle(document.body).getPropertyValue("--mono");
    c.textAlign = "center";
    (w < 600 ? [1, 101, 1631, 1855] : [1, 25, 50, 75, 101, 1631, 1675, 1725, 1775, 1825, 1855]).forEach(function (p) {
      if (idx[p] === undefined) return;
      c.fillText(String(p), Math.min(Math.max(geo.x(p), 12), w - 14), 18);
    });
    c.fillText(String(cfg.end[0] - cfg.start[1] - 1) + " unscored", bx + gap / 2, top + 3 * (rh + rg) + 10 > h ? h - 2 : top + 3 * (rh + rg) + 10);
    if (state.sel) {
      var x = geo.x(state.sel.p);
      c.fillStyle = "#e6ecf6";
      c.fillRect(Math.round(x) - 1, top - 6, 2, 3 * (rh + rg) + 4);
      c.beginPath(); c.moveTo(x - 5, top - 10); c.lineTo(x + 5, top - 10); c.lineTo(x, top - 4); c.fill();
      var lab = state.sel.wt + state.sel.p;
      c.fillStyle = "#0a1120";
      var tw = c.measureText(lab).width + 10;
      var lx = Math.min(Math.max(x, tw / 2), w - tw / 2);
      c.fillRect(lx - tw / 2, 3, tw, 14);
      c.fillStyle = "#e6ecf6";
      c.fillText(lab, lx, 14);
    }
  }

  function play() {
    if (playtimer) { stopplay(); return; }
    var list = positions.filter(function (p) { return regionof(p) === state.region; });
    var i = state.sel ? list.indexOf(state.sel.p) : 0;
    if (i < 0 || i >= list.length - 1) i = 0;
    $("#play").textContent = "Pause";
    playtimer = setInterval(function () {
      i += 1;
      if (i >= list.length) { stopplay(); return; }
      selectpos(list[i], { quick: true });
    }, 140);
  }
  function stopplay() { if (playtimer) clearInterval(playtimer); playtimer = null; $("#play").textContent = "Play"; }

  function setregion(r, keep) {
    state.region = r;
    document.querySelectorAll("[data-region]").forEach(function (b) { b.classList.toggle("on", b.getAttribute("data-region") === r); });
    stopflex();
    if (fallback) { drawribbon(); captions(); return; }
    loadstructure(r).then(function () {
      if (!keep) {
        var first = positions.filter(function (p) { return regionof(p) === r; });
        var pick = r === "start" ? 61 : first[0];
        if (!state.sel || regionof(state.sel.p) !== r) { selectpos(pick); return; }
      }
      paint(0);
    });
    captions();
    drawtrack();
  }

  function captions() {
    var s = data.structures[state.region];
    var cap = state.region === "start"
      ? "PDB " + s.id + " · solution NMR, BRCA1/BARD1 RING domains · model " + (frame + 1) + " of " + Math.max(1, models.length) + " · BRCA1 chain A, residues " + s.first + "–" + s.last
      : "PDB " + s.id + " · X-ray, 1.85 Å, BRCT domains with a BACH1 peptide · BRCA1 chain A, residues " + s.first + "–" + s.last;
    if (fallback) cap = "3D view unavailable here: scored positions drawn as a ribbon of dots.";
    $("#pdb-caption").textContent = cap;
    var note = "Structures are experimental (NMR for the start region, X-ray for the end region). The model never saw them. ";
    if (state.region === "start") note += "Start-region structure includes the partner protein BARD1, shown faint and not scored. Grey spheres: zinc ions in the file; thin lines join them to the letters that hold them.";
    else {
      var m = s.missing;
      note += "Positions " + m[0] + " to " + m[m.length - 1] + " are not in the crystal structure. Grey sticks: BACH1 peptide, from the crystal structure.";
    }
    if (fallback) note = "Each dot is one scored position, in order along the protein. Click a dot to select it.";
    $("#stage-note").textContent = note;
    $("#flex").hidden = fallback || state.region !== "start" || models.length < 2;
    $("#flex-label").textContent = models.length + " NMR conformers: the range of shapes the experiment allows. Showing model " + (frame + 1) + ".";
    legend();
  }

  function legend() {
    var el = $("#legend");
    var g, ends;
    if (state.mode === "lab") { g = [colors.harmful, colors.cant, colors.harmless]; ends = ["broken", "working"]; }
    else if (state.mode === "disagree") { g = [colors.dim, colors.bright]; ends = ["none", "all"]; }
    else { g = [colors.harmful, colors.cant, colors.harmless]; ends = [fmt(-cut.b2, 2), fmt(-cut.a, 2)]; }
    var title = { model: "Model score, mean of lab-tested swaps", lab: "Lab function score, mean", disagree: "Share of tested swaps where they disagree" }[state.mode];
    el.innerHTML = "<span>" + title + '</span><span class="bar" style="background:linear-gradient(90deg,' + g.map(hex).join(",") + ')"></span><span class="ends"><span>' + ends[0] + "</span><span>" + ends[1] + "</span></span>";
  }

  function fetchtext(url) {
    if (structcache[url]) return Promise.resolve(structcache[url]);
    return fetch(url).then(function (r) { if (!r.ok) throw new Error(url); return r.text(); }).then(function (t) { structcache[url] = t; return t; });
  }

  var loading = null;
  function loadstructure(r) {
    var file = r === "start" ? "structures/1jm7.pdb" : "structures/1t15.pdb";
    loading = file;
    return fetchtext(file).then(function (text) {
      if (loading !== file) return;
      viewer.removeAllModels();
      viewer.removeAllShapes();
      viewer.removeAllLabels();
      models = [];
      loaded = null;
      frame = 0;
      if (r === "start") {
        text.split(/^ENDMDL.*$/m).forEach(function (chunk) {
          if (/^ATOM/m.test(chunk)) models.push(viewer.addModel(chunk.replace(/^MODEL.*$/m, ""), "pdb"));
        });
        models.forEach(function (m, i) { if (i) m.hide(); });
        zincpairs = data.structures.start.zinc.map(function (z) { return z; });
      } else {
        models.push(viewer.addModel(text, "pdb"));
        zincpairs = [];
      }
      loaded = r;
      models.forEach(function (m) {
        viewer.setHoverable({ model: m, chain: "A" }, true, function (atom) { hover(atom); }, function () { $("#tip").hidden = true; });
        viewer.setClickable({ model: m, chain: "A" }, true, function (atom) { if (bypos[atom.resi]) { stopplay(); selectpos(atom.resi); } });
      });
      captions();
      draw(0, true);
      viewer.zoomTo({ chain: "A", model: models[0] });
      viewer.zoom(0.95);
      if (state.sel && regionof(state.sel.p) === r && structinfo(state.sel.p).has) viewer.center({ model: models[0], chain: "A", resi: state.sel.p }, reduced ? 0 : 600);
      viewer.render();
    });
  }

  function hover(atom) {
    var tip = $("#tip");
    var e = bypos[atom.resi];
    var letter = e ? e.wt : atom.resn;
    tip.textContent = (e ? letter + atom.resi : atom.resn + " " + atom.resi + " (not scored)");
    tip.hidden = false;
    var r = $(".stage").getBoundingClientRect();
    tip.style.left = Math.min(mouse.x - r.left + 12, r.width - 90) + "px";
    tip.style.top = (mouse.y - r.top + 12) + "px";
  }

  function paint(dur, nocamera) {
    if (!viewer) return;
    var first = !pending;
    pending = { dur: pending ? Math.max(pending.dur, dur) : dur, nocamera: pending ? pending.nocamera && nocamera : nocamera };
    if (first) requestAnimationFrame(function () { var a = pending; pending = null; draw(a.dur, a.nocamera); });
  }
  function draw(dur, nocamera) {
    if (!viewer || !models.length || !loaded) return;
    var m = models[frame];
    var mode = state.mode;
    var reg = loaded;
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    viewer.setStyle({ model: m }, {});
    viewer.setStyle({ model: m, chain: "A" }, { cartoon: { colorfunc: function (a) { return poscolor(a.resi, mode); }, thickness: 0.5 } });
    if (reg === "start") {
      viewer.setStyle({ model: m, chain: "B" }, { cartoon: { color: 0x56637a, opacity: 0.35 } });
      viewer.setStyle({ model: m, chain: "A", resn: "ZN" }, { sphere: { radius: 0.75, color: colors.zinc } });
      viewer.setStyle({ model: m, chain: "B", resn: "ZN" }, { sphere: { radius: 0.6, color: 0x3a4659 } });
      var znat = m.selectedAtoms({ chain: "A", resn: "ZN" });
      zincpairs.forEach(function (z) {
        var lig = m.selectedAtoms({ chain: "A", resi: z.resi, atom: z.atom })[0];
        var zn = znat.filter(function (a) { return a.resi === z.zn; })[0];
        if (lig && zn) viewer.addCylinder({ start: { x: lig.x, y: lig.y, z: lig.z }, end: { x: zn.x, y: zn.y, z: zn.z }, radius: 0.07, color: colors.zinc, dashed: true, fromCap: 1, toCap: 1 });
      });
      var cb = m.selectedAtoms({ chain: "B", atom: "CA" });
      if (cb.length) {
        var c = centroid(cb);
        viewer.addLabel("BARD1, partner protein, not scored", { position: c, fontSize: 11, fontColor: "#a9b6cc", backgroundColor: 0x0a1120, backgroundOpacity: 0.6, borderThickness: 0, inFront: true });
      }
    } else {
      viewer.setStyle({ model: m, chain: "B" }, { stick: { color: 0x8a96ab, radius: 0.18 } });
      var pb = m.selectedAtoms({ chain: "B", atom: "CA" });
      if (pb.length) viewer.addLabel("BACH1 peptide, from the crystal structure", { position: centroid(pb), fontSize: 11, fontColor: "#a9b6cc", backgroundColor: 0x0a1120, backgroundOpacity: 0.6, borderThickness: 0, inFront: true });
    }
    if (state.sel && regionof(state.sel.p) === reg) {
      var p = state.sel.p;
      var ca = m.selectedAtoms({ chain: "A", resi: p, atom: "CA" })[0];
      if (ca) {
        viewer.addStyle({ model: m, chain: "A", resi: p }, { stick: { color: 0xffffff, radius: 0.28 } });
        viewer.addSphere({ center: { x: ca.x, y: ca.y, z: ca.z }, radius: 2.6, color: 0xffffff, opacity: 0.16 });
        viewer.addLabel(state.sel.wt + p + " → " + state.sel.nw + " (letter only)", { position: { x: ca.x, y: ca.y + 2.2, z: ca.z }, fontSize: 12, fontColor: "#e6ecf6", backgroundColor: 0x111c31, backgroundOpacity: 0.9, borderThickness: 1, borderColor: 0x8fa3c4, inFront: true });
        if (!nocamera) viewer.center({ model: m, chain: "A", resi: p }, reduced ? 0 : dur);
      }
    }
    viewer.render();
  }
  function centroid(list) {
    var x = 0, y = 0, z = 0;
    list.forEach(function (a) { x += a.x; y += a.y; z += a.z; });
    return { x: x / list.length, y: y / list.length, z: z / list.length };
  }

  function stopflex() { if (flextimer) clearInterval(flextimer); flextimer = null; var b = $("#flex-btn"); if (b) b.textContent = "Cycle conformers"; }
  function flex() {
    if (flextimer) { stopflex(); return; }
    $("#flex-btn").textContent = "Stop";
    flextimer = setInterval(function () {
      models[frame].hide();
      frame = (frame + 1) % models.length;
      models[frame].show();
      draw(0, true);
      captions();
    }, 500);
  }

  function startspin() { if (viewer && !reduced) viewer.spin("y", 0.25); }
  function setupviewer() {
    var ok = typeof $3Dmol !== "undefined" && !/fallback/.test(location.search);
    if (ok) {
      try {
        viewer = $3Dmol.createViewer($("#viewer"), { backgroundColor: "black", backgroundAlpha: 0, antialias: true });
        viewer.setBackgroundColor(0x000000, 0);
        viewer.setHoverDuration(40);
      } catch (err) { viewer = null; }
    }
    if (!viewer) {
      fallback = true;
      $(".orient").textContent = "Drawn for orientation only. The model saw sequence alone.";
      $("#viewer").hidden = true;
      $("#fallback").hidden = false;
      setupribbon();
      return;
    }
    var box = $("#viewer");
    box.addEventListener("pointerdown", function () { viewer.spin(false); clearTimeout(spinpause); });
    window.addEventListener("pointerup", function () { clearTimeout(spinpause); spinpause = setTimeout(startspin, 2500); });
    box.addEventListener("wheel", function () { viewer.spin(false); clearTimeout(spinpause); spinpause = setTimeout(startspin, 2500); }, { passive: true });
  }

  var ribbon, rctx, rpts = [];
  function setupribbon() {
    ribbon = $("#ribbon");
    rctx = ribbon.getContext("2d");
    function near(ev) {
      var r = ribbon.getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top, best = null, bd = 18 * 18;
      rpts.forEach(function (q) { var d = (q.x - x) * (q.x - x) + (q.y - y) * (q.y - y); if (d < bd) { bd = d; best = q; } });
      return best;
    }
    ribbon.addEventListener("mousemove", function (ev) {
      var q = near(ev), tip = $("#tip");
      if (!q) { tip.hidden = true; return; }
      tip.textContent = bypos[q.p].wt + q.p; tip.hidden = false;
      tip.style.left = (q.x + 12) + "px"; tip.style.top = (q.y + 12) + "px";
    });
    ribbon.addEventListener("click", function (ev) { var q = near(ev); if (q) selectpos(q.p); });
  }
  function drawribbon() {
    if (!ribbon) return;
    var dpr = window.devicePixelRatio || 1, w = ribbon.clientWidth, h = ribbon.clientHeight;
    ribbon.width = w * dpr; ribbon.height = h * dpr;
    var c = rctx; c.setTransform(dpr, 0, 0, dpr, 0, 0); c.clearRect(0, 0, w, h);
    var list = positions.filter(function (p) { return regionof(p) === state.region; });
    var turns = state.region === "start" ? 2 : 4;
    rpts = list.map(function (p, i) {
      var t = i / (list.length - 1);
      return { p: p, x: 50 + t * (w - 100), y: h * 0.5 + Math.sin(t * Math.PI * 2 * turns) * h * 0.22 };
    });
    c.strokeStyle = "rgba(160,185,225,.2)"; c.lineWidth = 2; c.beginPath();
    rpts.forEach(function (q, i) { if (i) c.lineTo(q.x, q.y); else c.moveTo(q.x, q.y); }); c.stroke();
    rpts.forEach(function (q) {
      var s = state.sel && state.sel.p === q.p;
      c.fillStyle = hex(poscolor(q.p, state.mode));
      c.beginPath(); c.arc(q.x, q.y, s ? 7 : 4, 0, Math.PI * 2); c.fill();
      if (s) { c.strokeStyle = "#ffffff"; c.lineWidth = 2; c.beginPath(); c.arc(q.x, q.y, 11, 0, Math.PI * 2); c.stroke(); }
    });
  }

  function arm(a) {
    var au = data.auroc[a.arm];
    var h = '<div class="arm"><div class="name">' + esc(a.name) + "<span>AUROC " + (au === undefined ? na : au.toFixed(4)) + (a.cutoff === null ? "" : " · cutoff " + a.cutoff.toFixed(2)) + "</span></div>";
    h += '<div><div class="sbar"><i style="width:' + (a.caught / a.n_harmful * 100) + '%;background:var(--agree)"></i><i style="width:' + (a.missed / a.n_harmful * 100) + '%;background:var(--hatch)"></i></div>';
    h += '<div class="slab"><span>caught <b data-c="' + a.caught + '">0</b></span><span>missed <b data-c="' + a.missed + '">0</b></span><span>of ' + a.n_harmful + " harmful</span></div></div>";
    h += '<div><div class="sbar"><i style="width:' + (a.false_alarms / a.n_harmless * 100) + '%;background:var(--hatch)"></i><i style="width:' + (a.cleared / a.n_harmless * 100) + '%;background:var(--agree)"></i></div>';
    h += '<div class="slab"><span>false alarms <b data-c="' + a.false_alarms + '">0</b></span><span>cleared <b data-c="' + a.cleared + '">0</b></span><span>of ' + a.n_harmless + " harmless</span></div></div></div>";
    return h;
  }
  function guess(t) { return t.indexOf("Guess:") === 0 ? "<b>Guess:</b>" + esc(t.slice(6)) : "<b>Guess:</b> " + esc(t); }
  function failcard(kind, title) {
    var f = data.failures.filter(function (x) { return x.kind === kind && x.rank === 1; })[0];
    if (!f) return "<p>" + title + ": " + na + "</p>";
    return '<div class="fail"><p>' + title + '</p><p class="mono">' + esc(f.mutant) + " · model <span class=num>" + fmt(f.esm_score, 2) + "</span> · lab <span class=num>" + fmt(f.lab_score, 3) + "</span> (" + esc(f.lab_class) + ')</p><p class="guess">' + guess(f.guess) + "</p></div>";
  }
  function renderscore() {
    var s = data.sentences, b = data.bands;
    var h = '<h2 id="score-title">Scoreboard: the sealed exam (end region)</h2>';
    h += '<p class="lead">' + esc(s.S6 || na) + "</p>";
    h += "<p>" + esc(s.S1 || na) + " " + esc(s.S3 || "") + "</p>";
    h += '<div class="arms">' + data.arms.map(arm).join("") + "</div>";
    h += "<p class=strip-note>Counts from headline.csv. Each arm flags at one cutoff, so here the model's can't-tell band counts as flagged. The three-way split is below.</p>";
    if (b && b.LOF) {
      h += "<h3>The model's three calls against the lab (band_counts.csv)</h3><dl>";
      [["LOF", "Lab: loss of function"], ["FUNC", "Lab: kept working"], ["INT", "Lab: in between"]].forEach(function (r) {
        var x = b[r[0]];
        if (x) h += "<dt>" + r[1] + "</dt><dd class=num>looks harmful " + x.looks_harmful + " · can't tell " + x.cant_tell + " · looks harmless " + x.looks_harmless + " · total " + x.total + "</dd>";
      });
      h += "</dl>";
    }
    h += "<h3>Where the model was most wrong</h3>";
    h += '<div class="two">' + failcard("miss", "Worst miss") + failcard("false_alarm", "Worst false alarm") + "</div>";
    $("#score-body").innerHTML = h;
    $("#score-body").querySelectorAll("[data-c]").forEach(function (el) { countup(el, Number(el.getAttribute("data-c")), 0); });
  }
  function renderabout() {
    var s = data.structures;
    var h = '<h2 id="about-title">BRCA1 is our exam, not the use case</h2>';
    h += "<p>The real use is a gene where nobody has tested a swap yet. We picked BRCA1 because a lab tested almost every single-letter DNA change in two of its regions, so the model's calls on the protein changes those edits make can be graded against the lab's answer.</p>";
    h += "<p>BRCA1 is a tumour-suppressor protein. A swap marked harmful lost function in a cell-based lab test. That is not a patient diagnosis. The model's “Looks harmful” is a prediction of that lab result.</p>";
    h += "<h3>How it works</h3><ol class=steps>";
    h += "<li>A protein language model (" + esc(cfg.model_id || na) + ") reads BRCA1's letters. It never saw the lab results.</li>";
    h += "<li>At each position we hide the normal letter and ask the model how likely each letter is. A swap's score is log P(new letter) − log P(normal letter).</li>";
    h += "<li>Two cutoffs, set on the start region only, turn the score into a call: looks harmful (at or below " + fmt(-cut.b2, 2) + "), can't tell, or looks harmless (above " + fmt(-cut.a, 2) + ").</li>";
    h += "<li>We compare each call with the lab's answer: loss of function, kept working, or in between.</li></ol>";
    h += "<p>The lab edits DNA; we study the protein letter that edit changes. Every change here is an amino-acid change: C to G means cysteine to glycine.</p>";
    h += '<h3>Practice and sealed exam</h3><div class="two"><div><b>Start region, 1–101 (RING)</b><p>Practice. The cutoffs were tuned here.</p></div><div><b>End region, 1631–1855 (BRCT)</b><p>Sealed exam. Scored once, after the cutoffs were frozen.</p></div></div>';
    h += "<h3>Three limits</h3><ol class=steps><li>Blind to splicing: we could not test this.</li><li>Cells in a dish, not people (a 2025 preprint notes cell-type effects).</li><li>Stronger tools exist (AlphaMissense, free to download), and zero-shot scoring here is not a new method.</li></ol>";
    h += "<h3>About the 3D view</h3><p>Structure shown for orientation only. The model saw sequence alone. Structures are experimental (NMR for the start region, X-ray for the end region). The model never saw them.</p>";
    h += "<p>Letter check: " + s.start.id + " matches our letters at " + s.start.n_letter_match + " of " + s.start.n_with_coords + " scored positions it contains; " + s.end.id + " at " + s.end.n_letter_match + " of " + s.end.n_with_coords + ". Zinc-holding letters found in " + s.start.id + " model 1: " + s.start.zinc_positions.join(", ") + ".</p>";
    h += '<h3 id="tool">The simple tool and the code</h3><p>The official demo is the simple tool in <span class=mono>app/</span> (run <span class=mono>python app/app.py</span>). The notebook is <span class=mono>notebook.ipynb</span>. Both are in the project repository, alongside this explorer in <span class=mono>explorer/</span>.</p>';
    h += "<h3>Credits</h3><p>Lab data: Findlay et al. 2018, <i>Nature</i> 562:217–222, via ProteinGym. Model: ESM-2 (Lin et al. 2023). Structures from the RCSB Protein Data Bank: 1JM7 (Brzovic et al. 2001, <i>Nat. Struct. Biol.</i> 8:833) and 1T15 (Clapperton et al. 2004, <i>Nat. Struct. Mol. Biol.</i> 11:512). 3D drawing: 3Dmol.js 2.4.2 (BSD-3-Clause).</p>";
    $("#about-body").innerHTML = h;
  }
  function openov(id) {
    if (id === "score") renderscore();
    if (id === "about") renderabout();
    var o = document.getElementById(id);
    o.hidden = false;
    o.querySelector("[data-close]").focus();
  }

  function intro(done) {
    var box = $("#intro"), seq = $("#intro-seq"), capel = $("#intro-cap");
    if (reduced) { done(); return; }
    box.hidden = false;
    var lo = 49, hi = 73, spans = [], timers = [], over = false;
    for (var p = lo; p <= hi; p++) {
      var sp = document.createElement("span");
      sp.textContent = bypos[p].wt;
      if (p === 61) sp.innerHTML = bypos[p].wt + '<em class="pos">61</em>';
      seq.appendChild(sp); spans.push(sp);
    }
    function later(ms, f) { timers.push(setTimeout(f, ms)); }
    function finish() {
      if (over) return; over = true;
      timers.forEach(clearTimeout);
      box.classList.add("out");
      setTimeout(function () { box.hidden = true; }, 600);
      document.removeEventListener("keydown", finish);
      done();
    }
    spans.forEach(function (s, i) { later(i * 40, function () { s.classList.add("in"); }); });
    var c61g = swapat(61, "G");
    later(1300, function () {
      var s = spans[61 - lo];
      s.classList.add("hit");
      s.innerHTML = 'G<em class="pos">61</em><svg class="wave" viewBox="0 0 24 8" preserveAspectRatio="none"><path d="M0 4 Q2 0 4 4 T8 4 T12 4 T16 4 T20 4 T24 4" fill="none" stroke="#ff7a66" stroke-width="1.6"/></svg>';
      capel.textContent = "One letter changed. Does it matter?";
    });
    later(2500, function () {
      var t = "BRCA1, residue 61. Cysteine to glycine.";
      if (zincat(61)) t += " Cysteine 61 is one of the letters that hold a zinc ion in place.";
      capel.textContent = t;
    });
    later(3900, function () {
      capel.textContent = "The model's call: " + c61g.v + " (" + fmt(c61g.esm, 2) + "). The lab's answer: " + (c61g.lc === "LOF" ? "loss of function in the lab assay" : c61g.lc) + " (" + fmt(c61g.ls, 3) + ").";
    });
    later(5300, finish);
    box.addEventListener("click", finish);
    document.addEventListener("keydown", finish);
  }

  function wire() {
    document.querySelectorAll("[data-region]").forEach(function (b) { b.addEventListener("click", function () { stopplay(); setregion(b.getAttribute("data-region")); }); });
    document.querySelectorAll("[data-mode]").forEach(function (b) {
      b.addEventListener("click", function () {
        state.mode = b.getAttribute("data-mode");
        document.querySelectorAll("[data-mode]").forEach(function (x) { x.classList.toggle("on", x === b); });
        legend();
        if (fallback) drawribbon(); else paint(0, true);
      });
    });
    document.querySelectorAll("[data-filter]").forEach(function (b) {
      b.addEventListener("click", function () {
        state.filter = b.getAttribute("data-filter");
        document.querySelectorAll("[data-filter]").forEach(function (x) { x.classList.toggle("on", x === b); });
        renderlist();
      });
    });
    $("#search").addEventListener("input", function (ev) { state.text = ev.target.value; $("#search-msg").textContent = ""; renderlist(); });
    $("#search-form").addEventListener("submit", function (ev) { ev.preventDefault(); search($("#search").value); });
    $("#open-score").addEventListener("click", function () { openov("score"); });
    $("#open-about").addEventListener("click", function () { openov("about"); });
    var tl = $("#tool-link");
    if (toolurl) { tl.href = toolurl; tl.target = "_blank"; tl.rel = "noopener"; }
    else tl.addEventListener("click", function (ev) { ev.preventDefault(); openov("about"); var t = document.getElementById("tool"); if (t) t.scrollIntoView(); });
    document.querySelectorAll(".overlay").forEach(function (o) {
      o.addEventListener("click", function (ev) { if (ev.target === o || ev.target.hasAttribute("data-close")) o.hidden = true; });
    });
    document.addEventListener("keydown", function (ev) { if (ev.key === "Escape") document.querySelectorAll(".overlay").forEach(function (o) { o.hidden = true; }); });
    $("#play").addEventListener("click", play);
    $("#flex-btn").addEventListener("click", flex);
    document.addEventListener("mousemove", function (ev) { mouse.x = ev.clientX; mouse.y = ev.clientY; });
    var rt;
    window.addEventListener("resize", function () {
      clearTimeout(rt);
      rt = setTimeout(function () { drawtrack(); if (viewer) { viewer.resize(); viewer.render(); } if (fallback) drawribbon(); }, 120);
    });
  }

  function start() {
    fetch("data.json").then(function (r) { return r.json(); }).then(function (d) {
      data = d;
      prepare();
      ranges();
      wire();
      setuptrack();
      setupviewer();
      renderlist();
      var target = swapat(61, "G");
      var asked = deeplink();
      if (asked) target = asked;
      state.sel = target;
      state.region = regionof(target.p);
      document.querySelectorAll("[data-region]").forEach(function (b) { b.classList.toggle("on", b.getAttribute("data-region") === state.region); });
      captions();
      renderreadout({});
      drawtrack();
      var first = fallback ? Promise.resolve() : loadstructure(state.region);
      if (fallback) drawribbon();
      var go = function () { first.then(function () { select(target, { fromlist: !asked }); startspin(); }); };
      if (asked || new URLSearchParams(location.search).get("m")) go(); else intro(go);
    }).catch(function (err) {
      $("#readout").innerHTML = "<p>Could not load data.json. Open this page through a local web server.</p>";
    });
  }
  start();
})();
