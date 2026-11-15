/* FBW tetracyclone trapping explorer: energy diagrams, structure table and path player.
   Data come from data/data.js (window.SB_DATA) and data/paths.js (window.SB_PATHS); 3D views use 3Dmol.js. */
(function () {
  "use strict";

  var D = window.SB_DATA;
  var PATHS = window.SB_PATHS || [];
  if (!D) { document.body.innerHTML = "<p style='padding:20px'>data/data.js did not load.</p>"; return; }
  var SP = {};
  D.species.forEach(function (s) { SP[s.id] = s; });
  var ST = D.states;
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
  function kindName(k) { return k === "ts" ? "TS" : (k === "complex" ? "complex" : "minimum"); }
  function pill(k) { return "<span class='pill " + k + "'>" + kindName(k) + "</span>"; }
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
    if (this.v) { this.v.setBackgroundColor(cssVar("--mol-bg") || "white"); this.v.resize(); this.v.render(); }
  };
  var views = {};
  function view(id) { if (!views[id]) views[id] = new MolView(document.getElementById(id)); return views[id]; }

  // ------------------------------------------------------------------ tabs
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
  function drawDiagram() {
    var el = $("#diagram");
    if (!el || !el.offsetParent) return;
    $("#routetitle").textContent = curRoute.label;
    $("#routedesc").textContent = curRoute.desc;
    var W = Math.max(520, el.clientWidth || 800), H = 430;
    var m = { l: 60, r: 14, t: 40, b: 34 };
    var ids = curRoute.states;
    function val(id) { var v = ST[id][mode]; return (v === null || v === undefined) ? null : v; }
    var present = ids.filter(function (id) { return val(id) !== null; });
    var vals = present.map(val);
    var lo = Math.min.apply(null, vals), hi = Math.max.apply(null, vals);
    var span = (hi - lo) || 20;
    lo -= span * 0.10; hi += span * 0.16;
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
                      (mode === "G" ? "qh-G" : "Electronic E") + " (kJ/mol)"));
    s.appendChild(svg("text", { x: W - m.r, y: H - 8, "text-anchor": "end", "class": "axislab" }, "reaction co-ordinate →"));
    // links, bridging over states without a value in this mode
    var links = curRoute.links.map(function (l) { return [l[0], l[1], false]; });
    ids.filter(function (id) { return val(id) === null; }).forEach(function (mid) {
      var ins = links.filter(function (l) { return l[1] === mid; }).map(function (l) { return l[0]; });
      var outs = links.filter(function (l) { return l[0] === mid; }).map(function (l) { return l[1]; });
      links = links.filter(function (l) { return l[0] !== mid && l[1] !== mid; });
      ins.forEach(function (a) { outs.forEach(function (b) { links.push([a, b, true]); }); });
    });
    links.forEach(function (l) {
      var ia = ids.indexOf(l[0]), ib = ids.indexOf(l[1]);
      if (ia < 0 || ib < 0 || val(l[0]) === null || val(l[1]) === null) return;
      var left = ia < ib, xa = X(ia) + (left ? bw / 2 : -bw / 2), xb = X(ib) + (left ? -bw / 2 : bw / 2);
      s.appendChild(svg("line", { x1: xa, y1: Y(val(l[0])), x2: xb, y2: Y(val(l[1])), "class": l[2] ? "link bridge" : "link" }));
    });
    // levels
    var maxChars = Math.max(8, Math.floor(slot / 6.6));
    ids.forEach(function (id, i) {
      var v = val(id);
      if (v === null) return;
      var st = ST[id];
      var g = svg("g", { "class": "lvl" + (id === selState ? " sel" : ""), "data-id": id, tabindex: 0, role: "button",
                         "aria-label": st.label + ", " + signed(v) + " kJ/mol" });
      g.appendChild(svg("rect", { x: X(i) - slot / 2 + 2, y: Y(v) - 46, width: slot - 4, height: 66, "class": "hit" }));
      g.appendChild(svg("line", { x1: X(i) - bw / 2, x2: X(i) + bw / 2, y1: Y(v), y2: Y(v), "class": "lvlbar k-" + st.kind }));
      var lines = wrapLabel(st.short, maxChars);
      lines.forEach(function (ln, k) {
        g.appendChild(svg("text", { x: X(i), y: Y(v) - 8 - (lines.length - 1 - k) * 13, "text-anchor": "middle", "class": "lab" }, ln));
      });
      g.appendChild(svg("text", { x: X(i), y: Y(v) + 16, "text-anchor": "middle", "class": "val k-" + st.kind }, signed(v)));
      g.addEventListener("click", function () { selectState(id); });
      g.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); selectState(id); } });
      s.appendChild(g);
    });
    el.innerHTML = "";
    el.appendChild(s);
    var missing = ids.filter(function (id) { return val(id) === null; }).map(function (id) { return ST[id].label; });
    $("#diagramnote").textContent = (missing.length ? "Not shown in qh-G (no frequencies yet): " + missing.join("; ") + ". " : "") +
      D.meta.reference + ". " + (mode === "G" ? "qh-G at 313.15 K, 1 mol/L." : "M06-2X electronic energies (no ZPE).");
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
  function speciesBlock(sp) {
    var atoms = speciesAtoms(sp.id);
    var kt = keyText(atoms, sp.keys);
    return "<div class='member'><b>" + esc(sp.label) + "</b>" + pill(sp.kind) +
      "<div class='kv'>" +
      (sp.imag && sp.imag.length ? "<div>Imaginary frequency</div><div>" + sp.imag.map(function (f) { return f.toFixed(1); }).join(", ") + " cm⁻¹</div>" : "") +
      (kt ? "<div>" + (sp.kind === "ts" ? "Forming/breaking bonds" : "Key distances") + "</div><div class='mono'>" + esc(kt) + "</div>" : "") +
      "<div>E</div><div class='mono'>" + sp.E_Eh.toFixed(6) + " Eh</div>" +
      "<div>qh-G</div><div class='mono'>" + (sp.G_Eh === null ? "– (no frequencies)" : sp.G_Eh.toFixed(6) + " Eh") + "</div>" +
      "<div>Atoms</div><div>" + sp.natoms + "</div>" +
      "<div>Source</div><div>" + esc(sp.source) + "</div>" +
      "<div>File</div><div><a href='" + esc(sp.file) + "' download>" + esc(sp.file) + "</a></div>" +
      "</div></div>";
  }
  function routeInfo() {
    var box = $("#routeinfo");
    if (!selState) { box.innerHTML = ""; return; }
    var st = ST[selState];
    box.innerHTML = "<h2>" + esc(st.label) + pill(st.kind) + "</h2>" +
      "<div class='kv'><div>qh-G</div><div>" + signed(st.G) + (st.G === null ? "" : " kJ/mol") + "</div>" +
      "<div>Electronic E</div><div>" + signed(st.E) + " kJ/mol</div></div>" +
      st.members.map(function (mid) { return speciesBlock(SP[mid]); }).join("");
  }

  // ------------------------------------------------------------------ structures table
  var selSpecies = null;
  function speciesTable() {
    var t = $("#sptable");
    var rows = D.species.map(function (s) {
      var stl = s.states.map(function (k) { return ST[k].short + " (" + signed(ST[k].G) + ")"; }).join("; ");
      return "<tr data-id='" + s.id + "' tabindex='0'" + (s.id === selSpecies ? " class='on'" : "") + ">" +
        "<td>" + esc(s.label) + pill(s.kind) + "<br><small>" + esc(s.source) + "</small></td>" +
        "<td class='num'>" + s.natoms + "</td>" +
        "<td class='num mono'>" + s.E_Eh.toFixed(6) + "</td>" +
        "<td class='num mono'>" + (s.G_Eh === null ? "–" : s.G_Eh.toFixed(6)) + "</td>" +
        "<td class='num'>" + (s.imag.length ? s.imag.map(function (f) { return f.toFixed(1); }).join(", ") : "–") + "</td>" +
        "<td><small>" + esc(stl) + "</small></td>" +
        "<td><a href='" + esc(s.file) + "' download>xyz</a></td></tr>";
    }).join("");
    t.innerHTML = "<thead><tr><th>Structure</th><th>Atoms</th><th>E (Eh)</th><th>qh-G (Eh)</th>" +
      "<th>Imag. freq. (cm⁻¹)</th><th>In diagram states (qh-G, kJ/mol)</th><th></th></tr></thead><tbody>" + rows + "</tbody>";
    $$("tbody tr", t).forEach(function (tr) {
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
    $("#spinfo").innerHTML = "<h2>" + esc(sp.label) + pill(sp.kind) + "</h2>" + speciesBlock(sp).replace("<div class='member'>", "<div>");
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
      "<div>Energy</div><div>" + (e === null ? "–" : signed(e) + " kJ/mol") + "</div>" +
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

  // ------------------------------------------------------------------ about
  function about() {
    var M = D.meta;
    $("#about").innerHTML =
      "<h2>What this shows</h2>" +
      "<p>The [2.2.2] carbene (the 6→7 system) can be trapped by tetracyclone in two ways. It can add directly across a " +
      "C=C bond of tetracyclone ([1+2], on either face, giving the S,S,Z and S,S,E cyclopropanes). Or it can first undergo " +
      "the Fritsch–Buttenberg–Wiechell (FBW) rearrangement to the strained [3.2.2] alkyne, which tetracyclone traps by a " +
      "Diels–Alder ([4+2]) reaction; the adduct then loses CO. This is section B of the study's ORCA work summary.</p>" +
      "<h2>Level of theory</h2><p>" + esc(M.level) + ". " + esc(M.thermo) + ". " + esc(M.reference) + ".</p>" +
      "<h2>Where the structures come from</h2><ul>" +
      "<li><b>Minima of the original study</b>: the carbene, tetracyclone, alkyne, the two cyclopropanes, the [4+2] adduct, " +
      "the decarbonylated product and CO.</li>" +
      "<li><b>The four trapping TSs</b> (marked NOT FOUND in the original summary) were located by OptTS with an analytic " +
      "Hessian, each with one imaginary frequency along the expected bond. Relaxed scans peak within 0.04–1.4 kJ/mol of " +
      "each TS, and NEB-TS runs from the pre-reaction complexes converged onto the same S,S,Z and S,S,E TSs.</li>" +
      "<li><b>Pre-reaction complexes</b> were optimised for the NEB runs. The two carbene complexes have frequencies; the " +
      "alkyne·tetracyclone and product·CO complexes were stopped before frequencies (the product·CO complex just short " +
      "of full convergence), so they have electronic energies only.</li>" +
      "<li><b>The FBW TS</b> is the original study's; an NEB-TS run reproduced it.</li></ul>" +
      "<h2>Reading the views</h2><ul>" +
      "<li>Dashed orange lines mark the bonds that form or break in each step, labelled with their length.</li>" +
      "<li>Atom numbers are ORCA's, counting from 0 (C15 is the carbene carbon in the trapping steps).</li>" +
      "<li>States made of separate molecules (for example carbene + tetracyclone) show them side by side.</li>" +
      "<li>The carbene complexes lie below the separate molecules in electronic energy but above them in qh-G at 1 mol/L, " +
      "because binding costs entropy; the barriers are therefore measured from the separate molecules, as in the PDF.</li></ul>" +
      "<h2>Still running</h2><p>Fresh IRCs from all five TSs, each followed by optimisation of both ends, are running and " +
      "will be added to the Paths tab when they finish.</p>" +
      "<p class='note'>Built " + esc(M.built) + " from the calculation outputs. 3D views use 3Dmol.js (BSD-3-Clause).</p>";
  }

  // ------------------------------------------------------------------ start
  routeList();
  speciesTable();
  pathList();
  about();
  var start = (location.hash || "").replace("#", "");
  setTab(["routes", "species", "paths", "about"].indexOf(start) >= 0 ? start : "routes");
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
