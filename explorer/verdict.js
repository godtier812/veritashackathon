(function (scope) {
  var harmful = "Looks harmful";
  var cant = "Can't tell";
  var harmless = "Looks harmless";

  function cutoffs(config) {
    var a = Number(config.cutoff_a);
    var b = config.cutoff_b === null || config.cutoff_b === undefined ? null : Number(config.cutoff_b);
    if (!config.has_band || b === null || !(b > a)) b = null;
    return { a: a, b: b };
  }

  function verdict(esm, config) {
    var cut = cutoffs(config);
    var harm = 0.0 - Number(esm);
    if (cut.b !== null) {
      if (harm >= cut.b) return harmful;
      if (harm < cut.a) return harmless;
      return cant;
    }
    return harm >= cut.a ? harmful : harmless;
  }

  function agreement(v, labclass) {
    if (!labclass) return "untested";
    if (v === cant) return "cant";
    if (labclass === "INT") return "grey";
    if ((v === harmful && labclass === "LOF") || (v === harmless && labclass === "FUNC")) return "agree";
    if (v === harmless && labclass === "LOF") return "miss";
    return "alarm";
  }

  var three = {
    ALA: "A", CYS: "C", ASP: "D", GLU: "E", PHE: "F", GLY: "G", HIS: "H", ILE: "I", LYS: "K", LEU: "L",
    MET: "M", ASN: "N", PRO: "P", GLN: "Q", ARG: "R", SER: "S", THR: "T", VAL: "V", TRP: "W", TYR: "Y"
  };
  var order = "ACDEFGHIKLMNPQRSTVWY";

  function parse(text, length) {
    if (text === null || text === undefined || String(text).trim() === "") {
      return { error: "Type a mutation such as C61G: the normal letter, the position, then the new letter." };
    }
    var raw = String(text);
    var s = raw.trim();
    var body = s.slice(0, 2).toLowerCase() === "p." ? s.slice(2).trim() : s;
    var bad = 'I could not read "' + raw.slice(0, 40).replace(/[`<>*_\[\]]/g, "") +
      '". Use the form C61G: the normal letter, the position, then the new letter. Three-letter names such as Cys61Gly also work.';
    var wt, pos, nw;
    var m1 = body.match(/^([A-Za-z])\s*(\d{1,5})\s*([A-Za-z])$/);
    if (m1) {
      wt = m1[1].toUpperCase(); nw = m1[3].toUpperCase(); pos = m1[2];
      if (order.indexOf(wt) < 0 || order.indexOf(nw) < 0) return { error: bad };
    } else {
      var m2 = body.match(/^([A-Za-z]{3})\s*(\d{1,5})\s*([A-Za-z]{3})$/);
      if (!m2) return { error: bad };
      wt = three[m2[1].toUpperCase()]; nw = three[m2[3].toUpperCase()]; pos = m2[2];
      if (!wt || !nw) return { error: bad };
    }
    var p = parseInt(pos, 10);
    if (wt === nw) return { error: '"' + wt + p + nw + '" is not a change: both letters are ' + wt + "." };
    if (p < 1 || p > length) return { error: "Position " + p + " is outside BRCA1, which has 1,863 positions." };
    return { wt: wt, p: p, nw: nw };
  }

  var api = { verdict: verdict, agreement: agreement, parse: parse, cutoffs: cutoffs, harmful: harmful, cant: cant, harmless: harmless };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else scope.rules = api;
})(this);
