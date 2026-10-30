/*
 * The SDP/SDPN calculator as a component: used inline on the /sdp page and in the
 * pop-up on the search page. Needs sdpcalc.js (the conversions) and sdpcalc.css.
 *
 *   const calc = mountSdpCalc(element, {prefix: "dlg", onChange: n => ...});
 *   calc.set(5461);   calc.sdpn();   calc.focus();
 *
 * Founder bases are the state; every field is a view of them and can be edited. A field
 * being typed in is left as typed; it snaps to its standard form when you leave it.
 */
"use strict";

function mountSdpCalc(root, opts = {}) {
  const P = opts.prefix || "sdp";
  const STRAINS = ["A/J", "C57BL/6J", "129S1", "NOD", "NZO", "CAST", "PWK", "WSB"];
  // display order: the bases everything is calculated from, then each encoding with its
  // groups (SDP before SDPN, as in the explanation); alleles follow last
  const FIELDS = [
    ["seq",     "Bases A–H",   "one per founder, e.g. ACCCCCCC",          t => parseSequence(t),       b => b.join("")],
    ["sdp",     "SDP",         "0–255",                                   (t, b) => parseSdp(t, b),     b => String(sdpFromBases(b)), "numeric"],
    ["split",   "SDP groups",  "reference : different",                   (t, b) => parseSplit(t, b),   b => sdpToPartition(sdpFromBases(b))],
    ["sdpn",    "SDPN",        "0–65535",                                 t => parseSdpn(t),           b => String(basesToSdpn(b)), "numeric"],
    ["pattern", "SDPN groups", "founders sharing a base, e.g. A:BC:DEFGH", (t, b) => parsePattern(t, b), b => sdpnToDisplay(basesToSdpn(b))],
  ];
  let bases = sdpnToBases(5461);

  root.classList.add("sdpcalc");
  root.innerHTML =
    `<div class="founders" role="group" aria-label="Base carried by each founder"></div>
     <div class="fields">${FIELDS.map(([k, label, hint, , , mode]) =>
       `<label for="${P}-${k}">${label}</label>
        <input id="${P}-${k}" data-k="${k}" autocomplete="off" spellcheck="false"${mode ? ` inputmode="${mode}"` : ""}>
        <span class="hint">${hint}</span>`).join("")}
       <label>Alleles</label><output class="alleles"></output><span class="hint">reference | others</span>
     </div>
     <div class="msg" role="status" aria-live="polite"></div>
     <div class="examples"><span class="hint">Examples:</span>
       <button type="button" data-sdpn="5461" title="A GWAS peak SNP; only A/J differs">rs29358506</button>
       <button type="button" data-sdpn="21843" title="Three bases: PWK and WSB differ from the reference and from each other">rs237124876</button>
       <button type="button" data-sdpn="21845" title="All 8 founders carry C">no difference</button>
     </div>`;
  const q = s => root.querySelector(s);

  SDP_LETTERS.split("").forEach((l, i) => {
    const col = document.createElement("div");
    col.className = "fcol";
    col.innerHTML = `<span class="let">${l}</span><span class="str" title="${STRAINS[i]}">${STRAINS[i]}</span>` +
                    `<span class="ref">${i === SDP_REF ? "reference" : ""}</span>`;
    for (const b of SDP_BASES) {
      const btn = document.createElement("button");
      btn.type = "button"; btn.textContent = b; btn.className = `b-${b}`;
      btn.dataset.founder = i; btn.dataset.base = b;
      btn.setAttribute("aria-label", `${STRAINS[i]} carries ${b}`);
      btn.onclick = () => { const nb = bases.slice(); nb[i] = b; update(nb); };
      col.appendChild(btn);
    }
    q(".founders").appendChild(col);
  });

  function msg(text, kind) { const m = q(".msg"); m.textContent = text || ""; m.className = "msg" + (text ? " " + kind : ""); }

  // redraw from `bases`; `skip` is the field being typed in
  function render(skip) {
    root.querySelectorAll(".fcol button").forEach(btn =>
      btn.setAttribute("aria-pressed", String(bases[+btn.dataset.founder] === btn.dataset.base)));
    for (const [k, , , , show] of FIELDS) {
      if (k === skip) continue;
      const el = q(`[data-k="${k}"]`);
      el.value = show(bases); el.removeAttribute("aria-invalid");
    }
    q(".alleles").textContent = allelesString(bases);
    if (opts.onChange) opts.onChange(basesToSdpn(bases));
  }
  function update(nb, skip, note) { bases = nb; render(skip); msg(note, "note"); }

  for (const [k, , , parse] of FIELDS) {
    const el = q(`[data-k="${k}"]`);
    el.addEventListener("input", () => {
      const r = parse(el.value, bases);
      if (r.error) { el.setAttribute("aria-invalid", "true"); msg(r.error, "err"); return; }
      el.removeAttribute("aria-invalid");
      update(r.bases, k, r.note);
    });
    el.addEventListener("change", () => { render(); if (q(".msg").classList.contains("err")) msg(""); });
  }
  root.querySelectorAll(".examples button").forEach(b => b.onclick = () => update(sdpnToBases(+b.dataset.sdpn)));

  render();
  return {
    set(n) { update(sdpnToBases(n)); },
    // apply a parse result from outside (e.g. a URL); shows its error or note
    apply(r) { if (r.error) { render(); msg(r.error, "err"); } else update(r.bases, null, r.note); },
    sdpn() { return basesToSdpn(bases); },
    bases() { return bases.slice(); },
    focus() { q(".fields input").focus(); },   // the first field
  };
}
