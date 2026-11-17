/* FBW tetracyclone trapping explorer: energy diagrams, structure tables, path player and the SI energy tables.
   Data come from data/data.js (window.SB_DATA) and data/paths.js (window.SB_PATHS), both written by
   tools/build_site.py; 3D views use 3Dmol.js. */
(function () {
  "use strict";

  var D = window.SB_DATA;
  var PATHS = window.SB_PATHS || [];
  if (!D) { document.body.innerHTML = "<p style='padding:20px'>data/data.js did not load.</p>"; return; }
  var SP = {};
  D.species.forEach(function (s) { SP[s.id] = s; });
  var ST = D.states;
  var M = D.meta;
  var SVGNS = "http://www.w3.org/2000/svg";

  // ------------------------------------------------------------------ helpers
  function $(sel, root) { return (root || document).querySelector(sel); }
  function $$(sel, root) { return Array.prototype.slice.call((root || document).querySelectorAll(sel)); }
  function cssVar(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function esc(s) {
    return String(s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }
  function signed(v, d) {
    if (v === null || v === undefined) return "–";
    var r = Number(v).toFixed(d === undefined ? 1 : d);
    if (Number(r) === 0) return r.replace("-", "");
    return (v > 0 ? "+" : "−") + r.replace("-", "");
  }
  function eh(v) { return (v === null || v === undefined) ? "–" : Number(v).toFixed(6).replace("-", "−"); }
  function kindName(k) { return k === "ts" ? "TS" : (k === "complex" ? "complex" : "minimum"); }
  function pill(k) { return "<span class='pill " + k + "'>" + kindName(k) + "</span>"; }
  function statusText(s) {
    return s.status === "complete" ? "" : s.status === "freq"
      ? (s.group === "irc" ? "frequencies still running" : "no frequencies (job stopped before them)")
      : "optimisation not converged";
  }
  function badge(s) { var t = statusText(s); return t ? " <span class='badge'>" + esc(t) + "</span>" : ""; }
  function svg(tag, attrs, text) {
    var e = document.createElementNS(SVGNS, tag);
    for (var a in attrs) { if (attrs[a] !== undefined && attrs[a] !== null) e.setAttribute(a, attrs[a]); }
    if (text !== undefined) e.textContent = text;
    return e;
  }
  function niceStep(raw) {
    var p = Math.pow(10, Math.floor(Math.log10(raw)));
    var f = raw / p;
    return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p;
  }

  // ------------------------------------------------------------------ xyz handling
  function parseXYZ(txt) {
    var L = txt.split("\n"), n = parseInt(L[0], 10), atoms = [];
    for (var i = 2; i < 2 + n; i++) {
      var t = L[i].trim().split(/\s+/);
      atoms.push({ el: t[0], x: +t[1], y: +t[2], z: +t[3] });
    }
    return atoms;
  }
  function toXYZ(atoms, title) {
    return atoms.length + "\n" + (title || "") + "\n" + atoms.map(function (a) {
      return a.el + " " + a.x.toFixed(6) + " " + a.y.toFixed(6) + " " + a.z.toFixed(6);
    }).join("\n") + "\n";
  }
  function centroid(atoms) {
    var c = { x: 0, y: 0, z: 0 };
    atoms.forEach(function (a) { c.x += a.x; c.y += a.y; c.z += a.z; });
    c.x /= atoms.length; c.y /= atoms.length; c.z /= atoms.length;
    return c;
  }
  // Several molecules side by side along x (first one keeps its atom order, so key-bond indices stay valid)
  function combine(list) {
    var out = [], xmax = null;
    list.forEach(function (atoms) {
      var c = centroid(atoms);
      var moved = atoms.map(function (a) { return { el: a.el, x: a.x - c.x, y: a.y - c.y, z: a.z - c.z }; });
      if (xmax !== null) {
        var xmin = Math.min.apply(null, moved.map(function (a) { return a.x; }));
        var shift = xmax + 3.5 - xmin;
        moved.forEach(function (a) { a.x += shift; });
      }
      xmax = Math.max.apply(null, moved.map(function (a) { return a.x; }));
      out = out.concat(moved);
    });
    return out;
  }
  function dist(a, b) { return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z)); }
  function keyText(atoms, keys) {
    return (keys || []).map(function (k) {
      var a = atoms[k[0]], b = atoms[k[1]];
      if (!a || !b) return "";
      return a.el + k[0] + "···" + b.el + k[1] + " " + dist(a, b).toFixed(2) + " Å";
    }).filter(Boolean).join(", ");
  }
  var parsed = {};
  function speciesAtoms(id) { if (!parsed[id]) parsed[id] = parseXYZ(SP[id].xyz); return parsed[id]; }

  // ------------------------------------------------------------------ 3D viewer
  function MolView(el) { this.el = el; this.v = null; }
  MolView.prototype.ensure = function () {
    if (!window.$3Dmol) { this.el.innerHTML = "<div class='hint'>3Dmol.js did not load.</div>"; return null; }
    if (!this.v) {
      var hint = $(".hint", this.el);
      if (hint) hint.remove();
      this.v = window.$3Dmol.createViewer(this.el, { backgroundColor: cssVar("--mol-bg") || "white", antialias: true });
    }
    return this.v;
  };
  MolView.prototype.show = function (xyz, keys, keepView) {
    var v = this.ensure();
    if (!v) return;
    var view = (keepView && this.hasModel) ? v.getView() : null;
    v.removeAllModels(); v.removeAllShapes(); v.removeAllLabels();
    var m = v.addModel(xyz, "xyz");
    v.setStyle({}, { stick: { radius: 0.13 }, sphere: { scale: 0.23 } });
    var atoms = m.selectedAtoms({});
    var keyColor = cssVar("--key") || "#f08c00";
    (keys || []).forEach(function (k) {
      var a = atoms[k[0]], b = atoms[k[1]];
      if (!a || !b) return;
      v.addCylinder({ start: { x: a.x, y: a.y, z: a.z }, end: { x: b.x, y: b.y, z: b.z }, radius: 0.07,
                      color: keyColor, dashed: true, fromCap: 1, toCap: 1 });
      v.addLabel(dist(a, b).toFixed(2) + " Å", {
        position: { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, z: (a.z + b.z) / 2 },
        fontSize: 12, fontColor: "#1d2330", backgroundColor: "#ffd8a8", backgroundOpacity: 0.9,
        showBackground: true, inFront: true, borderThickness: 0
      });
    });
    if (view) v.setView(view); else v.zoomTo();
    this.hasModel = true;
    v.render();
  };
  MolView.prototype.refresh = function () {
    if (this.v && this.el.offsetParent) { this.v.setBackgroundColor(cssVar("--mol-bg") || "white"); this.v.resize(); this.v.render(); }
  };
  var views = {};
  function view(id) { if (!views[id]) views[id] = new MolView(document.getElementById(id)); return views[id]; }

  // ------------------------------------------------------------------ tabs
  var TABS = ["routes", "species", "paths", "energies", "about"];
  function setTab(name) {
    $$("#tabs button").forEach(function (b) {
      var on = b.getAttribute("data-tab") === name;
      b.classList.toggle("on", on); b.setAttribute("aria-selected", on ? "true" : "false");
    });
    $$(".tab").forEach(function (t) { t.classList.toggle("on", t.id === "tab-" + name); });
    Object.keys(views).forEach(function (k) { views[k].refresh(); });
    if (name === "routes") drawDiagram();
    if (name === "paths" && curPath) drawPathPlot();
    try { history.replaceState(null, "", "#" + name); } catch (e) { /* file:// may refuse */ }
  }
  $$("#tabs button").forEach(function (b) { b.addEventListener("click", function () { setTab(b.getAttribute("data-tab")); }); });

  // ------------------------------------------------------------------ routes + energy diagram
  var curRoute = D.routes[0], mode = "G", selState = null;
  function routeList() {
    var box = $("#routelist");
    box.innerHTML = "<h3>Routes</h3>";
    D.routes.forEach(function (r) {
      var b = document.createElement("button");
      b.innerHTML = esc(r.label) + "<small>" + r.states.length + " states</small>";
      b.addEventListener("click", function () { curRoute = r; selState = null; routeList(); drawDiagram(); routeInfo(); });
      if (r === curRoute) b.classList.add("on");
      box.appendChild(b);
    });
    var lg = document.createElement("div");
    lg.className = "legend";
    lg.innerHTML = "<span class='minimum'>minimum</span><span class='complex'>complex</span><span class='ts'>TS</span>";
    box.appendChild(lg);
  }
  $$("#emode button").forEach(function (b) {
    b.addEventListener("click", function () {
      mode = b.getAttribute("data-mode");
      $$("#emode button").forEach(function (x) {
        var on = x === b; x.classList.toggle("on", on); x.setAttribute("aria-pressed", on ? "true" : "false");
      });
      drawDiagram();
    });
  });
  function wrapLabel(text, maxChars) {
    var words = text.split(" "), lines = [], cur = "";
    words.forEach(function (w) {
      if (cur && (cur + " " + w).length > maxChars) { lines.push(cur); cur = w; } else { cur = cur ? cur + " " + w : w; }
    });
    if (cur) lines.push(cur);
    return lines;
  }
  function overlapArea(a, b) {
    var w = Math.min(a.x2, b.x2) - Math.max(a.x1, b.x1), h = Math.min(a.y2, b.y2) - Math.max(a.y1, b.y1);
    return (w > 0 && h > 0) ? w * h : 0;
  }
  function score(cands, placed) {
    var s = 0;
    cands.forEach(function (a) { placed.forEach(function (b) { s += overlapArea(a, b); }); });
    return s;
  }
  function drawDiagram() {
    var el = $("#diagram");
    if (!el || !el.offsetParent) return;
    $("#routetitle").textContent = curRoute.label;
    $("#routedesc").textContent = curRoute.desc;
    var W = Math.max(520, el.clientWidth || 800), H = 450;
    var m = { l: 60, r: 14, t: 52, b: 46 };
    var ids = curRoute.states.filter(function (id) { return ST[id]; });
    function val(id) { var v = ST[id][mode]; return (v === null || v === undefined) ? null : v; }
    var present = ids.filter(function (id) { return val(id) !== null; });
    var vals = present.map(val);
    var lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
    var span = (hi - lo) || 20;
    lo -= span * 0.08; hi += span * 0.12;
    var n = ids.length, slot = (W - m.l - m.r) / n, bw = Math.min(64, slot * 0.66);
    function X(i) { return m.l + slot * (i + 0.5); }
    function Y(v) { return m.t + (hi - v) / (hi - lo) * (H - m.t - m.b); }
    var s = svg("svg", { viewBox: "0 0 " + W + " " + H, role: "img", "aria-label": "Energy diagram: " + curRoute.label });
    // grid + axis
    var step = niceStep((hi - lo) / 6);
    for (var t = Math.ceil(lo / step) * step; t <= hi; t += step) {
      s.appendChild(svg("line", { x1: m.l, x2: W - m.r, y1: Y(t), y2: Y(t), "class": "grid" }));
      s.appendChild(svg("text", { x: m.l - 8, y: Y(t) + 4, "text-anchor": "end", "class": "tick" }, signed(Math.round(t), 0)));
    }
    s.appendChild(svg("line", { x1: m.l, x2: m.l, y1: m.t - 10, y2: H - m.b, "class": "axis" }));
    s.appendChild(svg("text", { x: 14, y: (m.t + H - m.b) / 2, transform: "rotate(-90 14 " + (m.t + H - m.b) / 2 + ")",
                                "text-anchor": "middle", "class": "axislab" },
                      (mode === "G" ? "G (313.15 K, 1 M)" : "Electronic E") + " (kJ/mol)"));
    s.appendChild(svg("text", { x: W - m.r, y: H - 6, "text-anchor": "end", "class": "axislab" }, "reaction co-ordinate →"));
    // links, bridging over states without a value in this mode
    var links = curRoute.links.map(function (l) { return [l[0], l[1], l[2] || ""]; });
    ids.filter(function (id) { return val(id) === null; }).forEach(function (mid) {
      var ins = links.filter(function (l) { return l[1] === mid; });
      var outs = links.filter(function (l) { return l[0] === mid; });
      links = links.filter(function (l) { return l[0] !== mid && l[1] !== mid; });
      ins.forEach(function (a) { outs.forEach(function (b) { links.push([a[0], b[1], (a[2] || b[2]) ? "bridge conf" : "bridge"]); }); });
    });
    links.forEach(function (l) {
      var ia = ids.indexOf(l[0]), ib = ids.indexOf(l[1]);
      if (ia < 0 || ib < 0 || val(l[0]) === null || val(l[1]) === null) return;
      var left = ia < ib, xa = X(ia) + (left ? bw / 2 : -bw / 2), xb = X(ib) + (left ? -bw / 2 : bw / 2);
      s.appendChild(svg("line", { x1: xa, y1: Y(val(l[0])), x2: xb, y2: Y(val(l[1])), "class": "link " + l[2] }));
    });
    // levels: the label goes above the bar (value below) or below it (value above), whichever overlaps less of the
    // labels, values and bars already drawn; bars are reserved first so no text covers another level
    var maxChars = Math.max(7, Math.floor(slot / 6.4)), small = slot < 78, cw = small ? 5.9 : 6.5, lh = small ? 12 : 13;
    var placed = [];
    ids.forEach(function (id, i) {
      if (val(id) !== null) placed.push({ x1: X(i) - bw / 2, x2: X(i) + bw / 2, y1: Y(val(id)) - 4, y2: Y(val(id)) + 4 });
    });
    function tbox(cx, base, text, w8) { var w = text.length * (w8 || cw) + 2; return { x1: cx - w / 2, x2: cx + w / 2, y1: base - 10, y2: base + 3 }; }
    ids.forEach(function (id, i) {
      var v = val(id);
      if (v === null) return;
      var st = ST[id], y = Y(v);
      var lines = wrapLabel(st.short, maxChars);
      var vtxt = signed(v);
      function layout(place) {
        var labs = lines.map(function (ln, k) {
          return { text: ln, base: place === "above" ? y - 8 - (lines.length - 1 - k) * lh : y + 21 + k * lh };
        });
        var vb = place === "above" ? y + 16 : y - 7;
        return { place: place, labs: labs, vbase: vb,
                 boxes: labs.map(function (l) { return tbox(X(i), l.base, l.text); }).concat([tbox(X(i), vb, vtxt, 7)]) };
      }
      var A = layout("above"), B = layout("below");
      var sa = score(A.boxes, placed), sb = score(B.boxes, placed);
      var L = (sa === 0 || sa <= sb) ? A : B;
      placed = placed.concat(L.boxes);
      var ys = L.boxes.map(function (b) { return b.y1; }).concat([y - 6]), ye = L.boxes.map(function (b) { return b.y2; }).concat([y + 6]);
      var g = svg("g", { "class": "lvl" + (id === selState ? " sel" : ""), "data-id": id, tabindex: 0, role: "button",
                         "aria-label": st.label + ", " + vtxt + " kJ/mol" });
      g.appendChild(svg("rect", { x: X(i) - slot / 2 + 2, y: Math.min.apply(null, ys) - 3, width: slot - 4,
                                  height: Math.max.apply(null, ye) - Math.min.apply(null, ys) + 6, "class": "hit" }));
      g.appendChild(svg("line", { x1: X(i) - bw / 2, x2: X(i) + bw / 2, y1: y, y2: y, "class": "lvlbar k-" + st.kind }));
      L.labs.forEach(function (l) {
        g.appendChild(svg("text", { x: X(i), y: l.base, "text-anchor": "middle", "class": "lab" + (small ? " small" : "") }, l.text));
      });
      g.appendChild(svg("text", { x: X(i), y: L.vbase, "text-anchor": "middle", "class": "val k-" + st.kind }, vtxt));
      g.addEventListener("click", function () { selectState(id); });
      g.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); selectState(id); } });
      s.appendChild(g);
    });
    el.innerHTML = "";
    el.appendChild(s);
    var missing = ids.filter(function (id) { return val(id) === null; }).map(function (id) { return ST[id].label; });
    $("#diagramnote").textContent = (missing.length ? "Not shown (no frequencies yet): " + missing.join("; ") + ". " : "") +
      M.reference + ". " + (mode === "G" ? "M06-2X free energies (GoodVibes, 313.15 K, 1 mol/L)." :
      "M06-2X electronic energies (no ZPE).") + " Dotted links: a conformational change that has not been computed.";
  }
  function selectState(id) {
    selState = id;
    drawDiagram();
    var st = ST[id];
    var frags = st.members.map(speciesAtoms);
    var atoms = frags.length > 1 ? combine(frags) : frags[0];
    var keys = SP[st.members[0]].keys;
    view("routeviewer").show(toXYZ(atoms, st.label), keys, false);
    routeInfo();
  }
  function speciesBlock(sp, withHead) {
    var atoms = speciesAtoms(sp.id);
    var kt = keyText(atoms, sp.keys);
    var c = sp.compare;
    var cmp = "";
    if (c) {
      cmp = "<div>Compared with</div><div>" + esc(SP[c.match] ? SP[c.match].label : c.match) + ": " +
        (c.same ? "<b>the same structure</b>" : "<b>a different conformer</b>") + " (" + esc(c.measure) + " " +
        c.value.toFixed(3) + " Å" + (c.heavy !== null && c.heavy !== undefined ? ", heavy atoms " + c.heavy.toFixed(3) + " Å" : "") +
        "; ΔE " + signed(c.dE, 2) + (c.dG !== null && c.dG !== undefined ? ", ΔG " + signed(c.dG, 2) : "") + " kJ/mol)</div>";
    }
    return "<div class='member'>" + (withHead ? "<b>" + esc(sp.label) + "</b>" + pill(sp.kind) + badge(sp) : "") +
      "<div class='kv'>" +
      (sp.imag && sp.imag.length ? "<div>Imaginary frequency</div><div>" + sp.imag.map(function (f) { return f.toFixed(1); }).join(", ") + " cm⁻¹</div>" : "") +
      (kt ? "<div>" + (sp.kind === "ts" ? "Forming/breaking bonds" : "Key distances") + "</div><div class='mono'>" + esc(kt) + "</div>" : "") +
      "<div>Eel (M06-2X)</div><div class='mono'>" + eh(sp.E_Eh) + " Eh</div>" +
      "<div>G (1 M)</div><div class='mono'>" + (sp.G_Eh === null ? "– (" + esc(statusText(sp)) + ")" : eh(sp.G_Eh) + " Eh") + "</div>" +
      (sp.Gcorr_Eh !== null ? "<div>G − Eel</div><div class='mono'>" + eh(sp.Gcorr_Eh) + " Eh</div>" : "") +
      (sp.Ecc_Eh !== null ? "<div>Eel (DLPNO-CCSD(T))</div><div class='mono'>" + eh(sp.Ecc_Eh) + " Eh (T1 " + sp.T1 + ")</div>" : "") +
      "<div>Atoms</div><div>" + sp.natoms + "</div>" +
      "<div>Source</div><div>" + esc(sp.source) + "</div>" + cmp +
      "<div>File</div><div><a href='" + esc(sp.file) + "' download>" + esc(sp.file) + "</a></div>" +
      "</div></div>";
  }
  function routeInfo() {
    var box = $("#routeinfo");
    if (!selState) { box.innerHTML = ""; return; }
    var st = ST[selState];
    var fb = st.fallback_for && st.fallback_for.length
      ? "<p class='note'>The IRC end point for this state is still waiting for its frequencies, so the structure below " +
        "stands in until they are done.</p>" : "";
    box.innerHTML = "<h2>" + esc(st.label) + pill(st.kind) + "</h2>" +
      "<div class='kv'><div>G (1 M)</div><div>" + signed(st.G) + (st.G === null ? " (no frequencies yet)" : " kJ/mol") + "</div>" +
      "<div>G (1 atm)</div><div>" + signed(st.G_atm) + (st.G_atm === null ? "" : " kJ/mol") + "</div>" +
      "<div>Electronic E</div><div>" + signed(st.E) + " kJ/mol</div></div>" + fb +
      st.members.map(function (mid) { return speciesBlock(SP[mid], true); }).join("");
  }

  // ------------------------------------------------------------------ structures table
  var selSpecies = null;
  function speciesRow(s) {
    var stl = s.states.filter(function (k) { return ST[k].G !== null || ST[k].E !== null; }).map(function (k) {
      return ST[k].short + " (" + signed(ST[k].G !== null ? ST[k].G : ST[k].E) + (ST[k].G === null ? " E" : "") + ")";
    }).join("; ");
    return "<tr data-id='" + s.id + "' tabindex='0'" + (s.id === selSpecies ? " class='on'" : "") + ">" +
      "<td>" + esc(s.label) + pill(s.kind) + badge(s) + "<br><small>" + esc(s.source) + "</small></td>" +
      "<td class='num'>" + s.natoms + "</td>" +
      "<td class='num mono'>" + eh(s.E_Eh) + "</td>" +
      "<td class='num mono'>" + eh(s.G_Eh) + "</td>" +
      "<td class='num'>" + (s.imag.length ? s.imag.map(function (f) { return f.toFixed(1); }).join(", ") : (s.nimag === 0 ? "none" : "–")) + "</td>" +
      "<td><small>" + esc(stl || "–") + "</small></td>" +
      "<td><a href='" + esc(s.file) + "' download>xyz</a></td></tr>";
  }
  function speciesTable() {
    var t = $("#sptable");
    var main = D.species.filter(function (s) { return s.group !== "irc"; });
    var irc = D.species.filter(function (s) { return s.group === "irc"; });
    t.innerHTML = "<thead><tr><th>Structure</th><th>Atoms</th><th>Eel (Eh)</th><th>G, 1 M (Eh)</th>" +
      "<th>Imag. freq. (cm⁻¹)</th><th>In diagrams (kJ/mol)</th><th></th></tr></thead><tbody>" +
      "<tr class='grp'><td colspan='7'>Stationary points</td></tr>" + main.map(speciesRow).join("") +
      "<tr class='grp'><td colspan='7'>Optimised IRC end points (Opt + Freq from each end of each IRC)</td></tr>" +
      irc.map(speciesRow).join("") + "</tbody>";
    $$("tbody tr[data-id]", t).forEach(function (tr) {
      function pick(e) {
        if (e.target.tagName === "A") return;
        selectSpecies(tr.getAttribute("data-id"));
      }
      tr.addEventListener("click", pick);
      tr.addEventListener("keydown", function (e) { if (e.key === "Enter") pick(e); });
    });
  }
  function selectSpecies(id) {
    selSpecies = id;
    $$("#sptable tbody tr").forEach(function (tr) { tr.classList.toggle("on", tr.getAttribute("data-id") === id); });
    var sp = SP[id];
    view("spviewer").show(sp.xyz, sp.keys, false);
    $("#spinfo").innerHTML = "<h2>" + esc(sp.label) + pill(sp.kind) + badge(sp) + "</h2>" + speciesBlock(sp, false);
  }

  // ------------------------------------------------------------------ paths
  var curPath = null, frame = 0, timer = null, pathAtoms = [];
  function pathList() {
    var box = $("#pathlist");
    box.innerHTML = "";
    var groups = [];
    PATHS.forEach(function (p) { if (groups.indexOf(p.group) < 0) groups.push(p.group); });
    groups.forEach(function (gname) {
      var h = document.createElement("h3"); h.textContent = gname; box.appendChild(h);
      PATHS.filter(function (p) { return p.group === gname; }).forEach(function (p) {
        var b = document.createElement("button");
        b.innerHTML = esc(p.label) + "<small>" + p.frames.length + " frames</small>";
        if (curPath === p) b.classList.add("on");
        b.addEventListener("click", function () { selectPath(p); });
        box.appendChild(b);
      });
    });
  }
  function selectPath(p) {
    stop();
    curPath = p; frame = 0;
    pathAtoms = p.frames.map(parseXYZ);
    pathList();
    var sl = $("#pslider");
    sl.max = String(p.frames.length - 1); sl.value = "0";
    view("pathviewer").show(p.frames[0], p.keys, false);
    setFrame(0, true);
  }
  function setFrame(i, fresh) {
    if (!curPath) return;
    frame = Math.max(0, Math.min(curPath.frames.length - 1, i));
    if (!fresh) view("pathviewer").show(curPath.frames[frame], curPath.keys, true);
    $("#pslider").value = String(frame);
    var e = curPath.energies[frame];
    $("#pframe").textContent = "frame " + (frame + 1) + "/" + curPath.frames.length + (e === null ? "" : " · " + signed(e) + " kJ/mol");
    drawPathPlot();
    var kt = keyText(pathAtoms[frame], curPath.keys);
    $("#pathinfo").innerHTML = "<h2>" + esc(curPath.label) + "</h2><p class='desc'>" + esc(curPath.note) + "</p>" +
      "<div class='kv'><div>Frame</div><div>" + (frame + 1) + " of " + curPath.frames.length + "</div>" +
      "<div>Energy</div><div>" + (e === null ? "–" : signed(e) + " kJ/mol (E, relative to " + esc(curPath.ref) + ")") + "</div>" +
      (kt ? "<div>Key distances</div><div class='mono'>" + esc(kt) + "</div>" : "") +
      "<div>File</div><div><a href='" + esc(curPath.file) + "' download>" + esc(curPath.file) + "</a></div></div>";
  }
  function drawPathPlot() {
    var el = $("#pathplot");
    if (!el || !curPath || !el.offsetParent) return;
    var W = Math.max(320, el.clientWidth || 600), H = 230, m = { l: 58, r: 12, t: 14, b: 32 };
    var E = curPath.energies, n = E.length;
    var pts = [];
    E.forEach(function (v, i) { if (v !== null) pts.push([i, v]); });
    if (!pts.length) { el.innerHTML = ""; return; }
    var lo = Math.min.apply(null, pts.map(function (p) { return p[1]; }));
    var hi = Math.max.apply(null, pts.map(function (p) { return p[1]; }));
    var span = (hi - lo) || 10; lo -= span * 0.08; hi += span * 0.1;
    function X(i) { return m.l + (n > 1 ? i / (n - 1) : 0.5) * (W - m.l - m.r); }
    function Y(v) { return m.t + (hi - v) / (hi - lo) * (H - m.t - m.b); }
    var s = svg("svg", { viewBox: "0 0 " + W + " " + H, role: "img", "aria-label": "Energy along " + curPath.label });
    var step = niceStep((hi - lo) / 5);
    for (var t = Math.ceil(lo / step) * step; t <= hi; t += step) {
      s.appendChild(svg("line", { x1: m.l, x2: W - m.r, y1: Y(t), y2: Y(t), "class": "grid" }));
      s.appendChild(svg("text", { x: m.l - 8, y: Y(t) + 4, "text-anchor": "end", "class": "tick" }, signed(Math.round(t), 0)));
    }
    s.appendChild(svg("text", { x: 14, y: (m.t + H - m.b) / 2, transform: "rotate(-90 14 " + (m.t + H - m.b) / 2 + ")",
                                "text-anchor": "middle", "class": "axislab" }, "E (kJ/mol)"));
    s.appendChild(svg("text", { x: W - m.r, y: H - 8, "text-anchor": "end", "class": "axislab" }, "frame →"));
    s.appendChild(svg("polyline", { points: pts.map(function (p) { return X(p[0]) + "," + Y(p[1]); }).join(" "), "class": "eline" }));
    pts.forEach(function (p) { s.appendChild(svg("circle", { cx: X(p[0]), cy: Y(p[1]), r: 2.2, "class": "edot" })); });
    if (E[frame] !== null) {
      s.appendChild(svg("line", { x1: X(frame), x2: X(frame), y1: m.t, y2: H - m.b, "class": "cursor" }));
      s.appendChild(svg("circle", { cx: X(frame), cy: Y(E[frame]), r: 5, "class": "emark" }));
    }
    var hit = svg("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, "class": "hitarea" });
    hit.addEventListener("click", function (ev) {
      var r = s.getBoundingClientRect();
      var x = (ev.clientX - r.left) * (W / r.width);
      setFrame(Math.round((x - m.l) / (W - m.l - m.r) * (n - 1)));
    });
    s.appendChild(hit);
    el.innerHTML = "";
    el.appendChild(s);
  }
  function stop() { if (timer) { clearInterval(timer); timer = null; } $("#pplay").textContent = "▶"; }
  $("#pplay").addEventListener("click", function () {
    if (!curPath) return;
    if (timer) { stop(); return; }
    if (frame >= curPath.frames.length - 1) setFrame(0);
    $("#pplay").textContent = "⏸";
    timer = setInterval(function () {
      if (frame >= curPath.frames.length - 1) { stop(); return; }
      setFrame(frame + 1);
    }, 140);
  });
  $("#pfirst").addEventListener("click", function () { stop(); setFrame(0); });
  $("#pprev").addEventListener("click", function () { stop(); setFrame(frame - 1); });
  $("#pnext").addEventListener("click", function () { stop(); setFrame(frame + 1); });
  $("#pslider").addEventListener("input", function () { stop(); setFrame(parseInt(this.value, 10)); });

  // ------------------------------------------------------------------ energies (the SI tables)
  function table(head, rows, cls) {
    return "<div class='tablewrap'><table class='static " + (cls || "") + "'><thead><tr>" +
      head.map(function (h) { return "<th>" + h + "</th>"; }).join("") + "</tr></thead><tbody>" + rows.join("") +
      "</tbody></table></div>";
  }
  function td(v, cls) { return "<td" + (cls ? " class='" + cls + "'" : "") + ">" + v + "</td>"; }
  function energies() {
    var ccTxt = Object.keys(M.cc_check || {}).map(function (k) { return esc(k) + " " + signed(M.cc_check[k]); }).join(", ");
    var h = "<div class='card'><h2>How these numbers are made</h2><ul class='tight'>" +
      "<li><b>Structures:</b> " + esc(M.level) + ".</li>" +
      "<li><b>Free energies:</b> " + esc(M.thermo) + ". Command: <code>" + esc(M.gv_command) + "</code>. The 1 mol/L values " +
      "add RT ln(Vm/V°) = " + M.ss.toFixed(2) + " kJ/mol to each species' 1 atm value, so bimolecular associations are " +
      M.ss.toFixed(2) + " kJ/mol lower at 1 M than at 1 atm.</li>" +
      "<li><b>Single points:</b> " + esc(M.cc_level) + "; composite G = Eel(DLPNO-CCSD(T)) + [G − Eel](M06-2X) " +
      "(GoodVibes <code>--spc</code>). <b>" + M.n_cc + " of " + M.n_cc_needed + "</b> structures have their DLPNO-CCSD(T) energy " +
      "so far, so every relative energy here is an M06-2X value for now." +
      (ccTxt ? " Check of the pipeline with the four that are done: the composite FBW barrier and alkyne come out at " +
        ccTxt + " kJ/mol, as in the PDF (41.1 and −44.8)." : "") + "</li>" +
      "<li><b>Units:</b> kJ/mol (1 Eh = " + M.eh.toFixed(2) + " kJ/mol), relative to the separated reactants of each step.</li></ul></div>";

    // relative energies by reaction
    var rrows = [];
    D.reactions.forEach(function (rx) {
      if (rx.id === "all") return;
      rrows.push("<tr class='grp'><td colspan='4'>" + esc(rx.label) + "</td></tr>");
      rx.rows.forEach(function (r) {
        rrows.push("<tr>" + td(esc(r.label) + pill(r.kind)) + td(signed(r.dE), "num") + td(signed(r.dG), "num") +
          td(signed(r.dG_atm), "num") + "</tr>");
      });
    });
    h += "<div class='card'><h2>Relative energies by reaction (kJ/mol, M06-2X, 313.15 K)</h2>" +
      "<p class='note'>Each row relative to the first row of its reaction (the separated reactants). – means no frequencies yet.</p>" +
      table(["Structure", "ΔEel", "ΔG (1 M)", "ΔG (1 atm)"], rrows) + "</div>";

    // TS table
    var trows = D.ts_table.map(function (t) {
      return "<tr>" + td(esc(t.label)) + td(t.imag.map(function (x) { return x.toFixed(2); }).join(", ") || "–", "num") +
        td(t.forming.map(function (b) { return b.d.toFixed(3); }).join(", ") || "—", "num") +
        td(t.breaking.map(function (b) { return b.d.toFixed(3); }).join(", ") || "—", "num") +
        td("<small>" + esc(t.forming.concat(t.breaking).map(function (b) { return b.bond; }).join(", ")) + "</small>") + "</tr>";
    });
    h += "<div class='card'><h2>Transition structures</h2>" +
      table(["Transition structure", "Imag. freq. (cm⁻¹)", "Forming (Å)", "Breaking (Å)", "Bonds"], trows) +
      "<p class='note'>FBW shift: distance from the migrating carbon to its destination (forming) and to its origin (breaking). " +
      "Additions: the two C···C distances between the reacting fragments. CO extrusion: the two breaking C–C(O) bonds.</p></div>";

    // IRC table
    var irows = [];
    D.ircs.forEach(function (i) {
      irows.push("<tr class='grp'><td colspan='6'>" + esc(i.label) + " — from " + esc(SP[i.ts].label) +
        (i.done ? " (converged in " + i.converged + " directions)" : " (still running)") + "</td></tr>");
      i.ends.forEach(function (e) {
        var c = e.compare, verdict = "–";
        if (c) {
          verdict = (c.same ? "same as " : "<b>differs from</b> ") + esc(SP[c.match].label) + " <small>(" + esc(c.measure) + " " +
            c.value.toFixed(3) + " Å; ΔE " + signed(c.dE, 2) + (c.dG !== null && c.dG !== undefined ? ", ΔG " + signed(c.dG, 2) : "") +
            " kJ/mol)</small>";
        }
        irows.push("<tr>" + td(e.dir) + td(e.present ? esc(e.label) : "–") +
          td(e.status === "complete" ? "done" : e.status === "freq" ? "frequencies running" : esc(e.status)) +
          td(e.nimag === null ? "–" : e.nimag, "num") + td(signed(e.dE_ts), "num") + td(verdict) + "</tr>");
      });
    });
    h += "<div class='card'><h2>IRCs and their optimised end points</h2>" +
      table(["Direction", "End point", "Opt + Freq", "N<sub>imag</sub>", "ΔE vs TS", "Compared with"], irows) +
      "<p class='note'>ΔE vs TS: electronic energy of the optimised end point relative to its TS (kJ/mol). RMSD after " +
      "best-fit overlay of all atoms; where the atom order differs, the RMS difference of the sorted interatomic distances. " +
      "Same structure: below 0.3 Å (0.03 Å for the distance measure) and within 0.5 kJ/mol.</p></div>";

    // absolute energies
    var arows = [];
    [["Stationary points", "main"], ["Optimised IRC end points", "irc"]].forEach(function (g) {
      arows.push("<tr class='grp'><td colspan='8'>" + g[0] + "</td></tr>");
      D.species.filter(function (s) { return (s.group === "irc") === (g[1] === "irc"); }).forEach(function (s) {
        arows.push("<tr>" + td(esc(s.label) + badge(s)) + td(s.nimag === null ? "–" : s.nimag, "num") +
          td(eh(s.E_Eh), "num mono") + td(eh(s.Ecc_Eh), "num mono") + td(eh(s.Gcorr_Eh), "num mono") +
          td(eh(s.G_Eh), "num mono") + td(eh(s.Gc_Eh), "num mono") + td(s.T1 === null ? "–" : s.T1.toFixed(3), "num") + "</tr>");
      });
    });
    h += "<div class='card'><h2>Absolute energies (Eh)</h2>" +
      table(["Structure", "N<sub>imag</sub>", "Eel M06-2X", "Eel DLPNO-CCSD(T)", "G − Eel", "G (M06-2X)",
             "G (DLPNO-CCSD(T))", "T1"], arows) +
      "<p class='note'>G − Eel is the M06-2X thermal correction at 313.15 K and 1 mol/L; G (DLPNO-CCSD(T)) = Eel(DLPNO-CCSD(T)) " +
      "+ (G − Eel).</p></div>";

    // method check
    var mrows = D.methods.map(function (mm) {
      return "<tr>" + td(esc(SP[mm.id] ? SP[mm.id].label : mm.id) + "<br><small>" + esc(mm.what) + "</small>") +
        td("<span class='mono'>" + esc(mm.method_line) + "</span>") + td(esc(mm.freq)) + td(esc(mm.version), "nowrap") +
        td(mm.ok ? "<span class='ok'>✓ matches</span>" : "<span class='bad'>✗ " + esc(mm.failed.join(", ")) + "</span>") + "</tr>";
    });
    h += "<div class='card'><h2>Method check</h2><p class='note'>Every input checked against Part II of the SI: M06-2X/def2-TZVP " +
      "D3(0), RIJCOSX with def2/J, DefGrid2, TightSCF (TightOPT for optimisations), CPCM ε 9.08 / n 1.424, ORCA 6.1.0; " +
      "DLPNO-CCSD(T)/def2-TZVPP with def2-TZVPP/C and TightPNO. Frequencies: analytic (FREQ) except the original FBW TS " +
      "(numerical).</p>" + table(["Calculation", "Method line", "Frequencies", "ORCA", "Against the SI"], mrows) + "</div>";

    h += "<div class='card'><h2>Downloads</h2><ul class='tight'>" +
      [["thermochemistry/GoodVibes_1M.dat", "GoodVibes output, 1 mol/L"], ["thermochemistry/GoodVibes_1atm.dat", "GoodVibes output, 1 atm"],
       ["thermochemistry/GoodVibes_DLPNO_1M.dat", "GoodVibes output with the DLPNO-CCSD(T) energies done so far"],
       ["thermochemistry/energies_absolute.csv", "absolute energies"], ["thermochemistry/energies_relative.csv", "relative energies by reaction"],
       ["thermochemistry/transition_structures.csv", "transition structures"], ["thermochemistry/irc_end_points.csv", "IRC end points"],
       ["thermochemistry/method_check.csv", "method check"], ["thermochemistry/SI_coordinates.txt", "coordinates with energies, in the SI's format"]]
        .map(function (f) { return "<li><a href='" + f[0] + "' download>" + f[0] + "</a>: " + f[1] + "</li>"; }).join("") +
      "</ul><p class='note'>The ORCA outputs these come from are in <code>calculations/</code> in the repository.</p></div>";
    $("#energies").innerHTML = h;
  }

  // ------------------------------------------------------------------ about
  function about() {
    var pend = (M.pending || []).map(function (p) { return "<li>" + esc(p) + "</li>"; }).join("");
    var probs = (M.problems || []).length ? "<h2>Checks that failed</h2><ul>" +
      M.problems.map(function (p) { return "<li>" + esc(p) + "</li>"; }).join("") + "</ul>" : "";
    $("#about").innerHTML =
      "<h2>What this shows</h2>" +
      "<p>The [2.2.2] carbene (the 6→7 system) can be trapped by tetracyclone in two ways. It can add directly across a " +
      "C=C bond of tetracyclone ([1+2], on either face, giving the S,S,Z and S,S,E cyclopropanes). Or it can first undergo " +
      "the Fritsch–Buttenberg–Wiechell (FBW) rearrangement to the strained [3.2.2] alkyne, which tetracyclone traps by a " +
      "Diels–Alder ([4+2]) reaction; the adduct then loses CO. This is section B of the study's ORCA work summary.</p>" +
      "<h2>Method</h2><p>As in Part II of the Supplementary Information: " + esc(M.level) + ". " + esc(M.thermo) + ". " +
      "The Energies tab has the details, the SI-style tables and a check of every calculation's input against the SI.</p>" +
      "<h2>Where the structures come from</h2><ul>" +
      "<li><b>Minima of the original study</b>: the carbene, tetracyclone, alkyne, the two cyclopropanes, the [4+2] adduct, " +
      "the decarbonylated product and CO.</li>" +
      "<li><b>The four trapping TSs</b> (marked NOT FOUND in the original summary) were located by OptTS with an analytic " +
      "Hessian, each with one imaginary frequency along the expected bonds, and confirmed by relaxed scans and NEB-TS. A " +
      "lower [4+2] TS was found later from the top of a long NEB band; the diagrams use it.</li>" +
      "<li><b>IRCs</b> were run in both directions from each TS, and both ends optimised with frequencies. Where an end " +
      "point's frequencies are done it is used in the diagrams (it is the structure the DLPNO-CCSD(T) single points are run " +
      "on); until then the earlier structure stands in.</li>" +
      "<li><b>The FBW TS</b> is the original study's (numerical frequencies); its IRC connects the carbene and the alkyne.</li></ul>" +
      "<h2>What the IRCs show</h2><ul>" +
      "<li>The S,S,Z and S,S,E TSs connect the carbene·tetracyclone complexes (Z and E faces) with the S,S,Z and S,S,E " +
      "cyclopropanes of the original study.</li>" +
      "<li>The FBW TS connects the [2.2.2] carbene with the [3.2.2] alkyne.</li>" +
      "<li>The lowest [4+2] TS connects the alkyne·tetracyclone complex with <b>a different conformer of the [4+2] adduct</b> " +
      "from the original study's (the alkyne cage is turned). The CO-extrusion IRC, still running, will show which adduct " +
      "conformer the CO-extrusion TS starts from.</li></ul>" +
      "<h2>Reading the views</h2><ul>" +
      "<li>Dashed orange lines mark the bonds that form or break in each step, labelled with their length.</li>" +
      "<li>Atom numbers are ORCA's, counting from 0 (C15 is the carbene carbon in the trapping steps).</li>" +
      "<li>States made of separate molecules (for example carbene + tetracyclone) show them side by side.</li>" +
      "<li>The carbene complexes lie below the separate molecules in electronic energy but above them in G at 1 mol/L, " +
      "because binding costs entropy; the barriers are therefore measured from the separate molecules, as in the PDF.</li></ul>" +
      probs + "<h2>Still to come</h2><ul>" + pend + "</ul>" +
      "<p class='note'>Built " + esc(M.built) + " by tools/build_site.py from the calculation outputs in calculations/. " +
      "3D views use 3Dmol.js (BSD-3-Clause).</p>";
  }

  // ------------------------------------------------------------------ start
  routeList();
  speciesTable();
  pathList();
  energies();
  about();
  var start = (location.hash || "").replace("#", "");
  setTab(TABS.indexOf(start) >= 0 ? start : "routes");
  var rt;
  window.addEventListener("resize", function () {
    clearTimeout(rt);
    rt = setTimeout(function () { drawDiagram(); drawPathPlot(); Object.keys(views).forEach(function (k) { views[k].refresh(); }); }, 150);
  });
  if (window.matchMedia) {
    var mq = window.matchMedia("(prefers-color-scheme: dark)");
    var onScheme = function () { Object.keys(views).forEach(function (k) { views[k].refresh(); }); drawDiagram(); };
    if (mq.addEventListener) mq.addEventListener("change", onScheme); else if (mq.addListener) mq.addListener(onScheme);
  }
})();
