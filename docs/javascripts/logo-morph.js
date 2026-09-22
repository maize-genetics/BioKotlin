/*
 * Crossfades the header logo from biokotlin_minimal*.svg to biokotlin*.svg on
 * scroll. Uses *_light.svg variants when the slate (dark) palette is active.
 */
(function () {
  "use strict";

  var SCROLL_END = 96;
  var REDUCED_MOTION_THRESHOLD = 48;
  var raf = 0;
  var themeObserver = null;

  function prefersReducedMotion() {
    return window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  }

  function isDarkScheme() {
    var scheme =
      document.body.getAttribute("data-md-color-scheme") ||
      document.documentElement.getAttribute("data-md-color-scheme") ||
      "";
    return scheme === "slate";
  }

  function logoSrc(img, dark) {
    if (dark) {
      return img.getAttribute("data-logo-dark") || img.getAttribute("src");
    }
    return img.getAttribute("data-logo-default") || img.getAttribute("src");
  }

  function syncLogoSources() {
    var dark = isDarkScheme();
    var images = document.querySelectorAll(
      ".md-logo__image[data-logo-dark][data-logo-default]"
    );

    for (var i = 0; i < images.length; i++) {
      var next = logoSrc(images[i], dark);
      if (next && images[i].getAttribute("src") !== next) {
        images[i].setAttribute("src", next);
      }
    }
  }

  function scrollProgress() {
    var y = window.scrollY || document.documentElement.scrollTop || 0;
    if (prefersReducedMotion()) {
      return y > REDUCED_MOTION_THRESHOLD ? 1 : 0;
    }
    return Math.min(1, Math.max(0, y / SCROLL_END));
  }

  function measure(container) {
    var minimal = container.querySelector(".md-logo__image--minimal");
    var full = container.querySelector(".md-logo__image--full");
    if (!minimal || !full || !minimal.naturalWidth || !full.naturalWidth) {
      return null;
    }
    var height = minimal.getBoundingClientRect().height;
    if (!height) {
      return null;
    }
    return {
      minW: (minimal.naturalWidth / minimal.naturalHeight) * height,
      fullW: (full.naturalWidth / full.naturalHeight) * height,
    };
  }

  function applyWidth(container, progress, sizes) {
    if (!sizes) {
      return;
    }
    var width = sizes.minW + (sizes.fullW - sizes.minW) * progress;
    container.style.width = width.toFixed(2) + "px";
  }

  function applyProgress(progress) {
    document.documentElement.style.setProperty(
      "--biokotlin-logo-morph",
      String(progress)
    );

    var containers = document.querySelectorAll(".md-logo__morph");
    for (var i = 0; i < containers.length; i++) {
      var sizes = measure(containers[i]);
      applyWidth(containers[i], progress, sizes);
    }
  }

  function update() {
    raf = 0;
    applyProgress(scrollProgress());
  }

  function scheduleUpdate() {
    if (raf) {
      return;
    }
    raf = window.requestAnimationFrame(update);
  }

  function whenImagesReady(callback) {
    var images = document.querySelectorAll(".md-logo__image");
    var pending = 0;

    function done() {
      pending -= 1;
      if (pending <= 0) {
        callback();
      }
    }

    for (var i = 0; i < images.length; i++) {
      if (images[i].complete) {
        continue;
      }
      pending += 1;
      images[i].addEventListener("load", done, { once: true });
      images[i].addEventListener("error", done, { once: true });
    }

    if (pending === 0) {
      callback();
    }
  }

  function watchTheme() {
    if (themeObserver || !document.body) {
      return;
    }
    themeObserver = new MutationObserver(function () {
      syncLogoSources();
      whenImagesReady(scheduleUpdate);
    });
    themeObserver.observe(document.body, {
      attributes: true,
      attributeFilter: ["data-md-color-scheme"],
    });
  }

  function bind() {
    if (!document.querySelector(".md-logo__morph")) {
      return;
    }

    syncLogoSources();
    whenImagesReady(function () {
      applyProgress(scrollProgress());
    });
    scheduleUpdate();
  }

  function init() {
    watchTheme();
    document.addEventListener("scroll", scheduleUpdate, { passive: true, capture: true });
    window.addEventListener("resize", scheduleUpdate);
    window.addEventListener("load", scheduleUpdate);

    var motion = window.matchMedia("(prefers-reduced-motion: reduce)");
    if (motion.addEventListener) {
      motion.addEventListener("change", scheduleUpdate);
    } else if (motion.addListener) {
      motion.addListener(scheduleUpdate);
    }

    bind();
  }

  /* `document$` fires on load and after instant navigation (see playground.js). */
  if (window.document$ && typeof window.document$.subscribe === "function") {
    window.document$.subscribe(bind);
    init();
  } else if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
