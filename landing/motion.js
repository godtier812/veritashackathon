(function () {
  var root = document.documentElement;
  root.classList.add("js");
  var items = document.querySelectorAll("[data-anim]");
  var calm = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  var hit = document.querySelector(".lt.hit");
  if (hit && !calm) {
    var hb = hit.getBoundingClientRect();
    var hx = hb.left + hb.width / 2, hy = hb.top + hb.height / 2;
    document.querySelectorAll(".lt:not(.hit)").forEach(function (el) {
      var b = el.getBoundingClientRect();
      var d = Math.hypot(b.left + b.width / 2 - hx, b.top + b.height / 2 - hy);
      el.style.setProperty("--rd", Math.round(2880 + d * 1.1) + "ms");
    });
  }
  function tick(el, delay, span) {
    var to = parseInt(el.dataset.to, 10);
    if (isNaN(to)) return;
    var t0 = performance.now() + delay;
    el.textContent = "0";
    function step(t) {
      var p = Math.max(0, Math.min(1, (t - t0) / span));
      el.textContent = Math.round(to * p).toLocaleString("en-US");
      if (p < 1) requestAnimationFrame(step);
    }
    requestAnimationFrame(step);
  }
  function counts(field) {
    field.querySelectorAll(".tick").forEach(function (el) {
      var from = parseInt(el.dataset.from || "0", 10), of = parseInt(el.dataset.of, 10);
      tick(el, from * 1.6, Math.max(300, of * 1.6 + 300));
    });
  }
  if (calm || !("IntersectionObserver" in window)) {
    items.forEach(function (el) { el.classList.add("in"); });
    return;
  }
  var io = new IntersectionObserver(function (entries) {
    entries.forEach(function (en) {
      if (en.isIntersecting) {
        en.target.classList.add("in");
        io.unobserve(en.target);
        if (en.target.classList.contains("dots")) counts(en.target.closest(".field"));
      }
    });
  }, { threshold: 0.25, rootMargin: "0px 0px -8% 0px" });
  items.forEach(function (el) { io.observe(el); });
})();
