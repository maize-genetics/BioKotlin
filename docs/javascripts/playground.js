/*
 * TODO:
 * Turns code fences marked `{.runnable}` into Kotlin Playground editors.
 *
 * Running BioKotlin needs a JVM, so the snippets are compiled and executed by
 * a Kotlin compiler server that has BioKotlin on its classpath. Its URL comes
 * from `extra.playground_server` in zensical.toml, which overrides/main.html
 * writes into the page. While that is empty the snippets stay as ordinary
 * highlighted code, minus the `//sampleStart` and `//sampleEnd` markers that
 * only mean something to Playground.
 *
 * See website/playground-server/README.md for how to stand a server up.
 */
(function () {
  "use strict";

  var BLOCK_SELECTOR = ".runnable";
  var CODE_SELECTOR = ".runnable code";
  var MARKER = /^\s*\/\/\s*sample(Start|End)\s*$/;

  function server() {
    return (window.BIOKOTLIN_PLAYGROUND_SERVER || "").trim();
  }

  /* Pygments wraps every line in its own span, so dropping a line is just a
   * matter of removing the span whose text is a marker. */
  function stripMarkers(block) {
    var lines = block.querySelectorAll('code span[id^="__span-"]');
    for (var i = 0; i < lines.length; i++) {
      if (MARKER.test(lines[i].textContent)) {
        lines[i].remove();
      }
    }
  }

  function renderStatic() {
    var blocks = document.querySelectorAll(BLOCK_SELECTOR);
    for (var i = 0; i < blocks.length; i++) {
      if (blocks[i].dataset.playgroundState) {
        continue;
      }
      stripMarkers(blocks[i]);
      blocks[i].dataset.playgroundState = "static";
    }
  }

  function renderInteractive(url) {
    /* Playground replaces each matched element, so anything already converted
     * no longer matches CODE_SELECTOR and cannot be picked up twice. */
    var code = document.querySelectorAll(CODE_SELECTOR);
    if (!code.length) {
      return;
    }
    for (var i = 0; i < code.length; i++) {
      /* The "Open in Playground" link points at play.kotlinlang.org, which has
       * no BioKotlin on its classpath, so every snippet here would fail there. */
      code[i].setAttribute("data-crosslink", "disabled");
    }
    window.KotlinPlayground(CODE_SELECTOR, {
      server: url,
      callback: function (target) {
        target.dataset.playgroundState = "interactive";
      },
    });
  }

  function enhance() {
    var url = server();
    if (url && typeof window.KotlinPlayground === "function") {
      renderInteractive(url);
    } else {
      renderStatic();
    }
  }

  /* `document$` fires on load and again after every instant navigation, so
   * snippets on pages loaded without a full reload get enhanced too. */
  if (window.document$ && typeof window.document$.subscribe === "function") {
    window.document$.subscribe(enhance);
  } else {
    document.addEventListener("DOMContentLoaded", enhance);
  }
})();
