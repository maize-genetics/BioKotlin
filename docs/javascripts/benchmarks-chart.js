/*
 * Interactive throughput chart for docs/benchmarks.md.
 *
 * Expects Chart.js (UMD) to be loaded first and a canvas with id
 * benchmarks-throughput-chart. Does nothing on other pages.
 */
(function () {
  "use strict";

  var CANVAS_ID = "benchmarks-throughput-chart";
  var chartInstance = null;

  var DATA = {
    labels: ["Palindrome", "Translate", "Complement"],
    biokotlin: [143.0, 607.0, 1650.0],
    biopython: [1.6, 7.8, 625.2],
    speedups: ["89×", "76×", "2.6×"],
  };

  var COLORS = {
    biokotlin: "#536eb7",
    biopython: "#e69138",
    speedup: "#00c853",
  };

  function colorScheme() {
    var scheme =
      document.body.getAttribute("data-md-color-scheme") ||
      document.documentElement.getAttribute("data-md-color-scheme") ||
      "";
    var dark = scheme === "slate";
    return {
      text: dark ? "#e8eaed" : "#1a1a1a",
      muted: dark ? "#9aa0a6" : "#5f6368",
      grid: dark ? "rgba(255, 255, 255, 0.12)" : "rgba(0, 0, 0, 0.08)",
      tooltipBg: dark ? "#2d3139" : "#ffffff",
      tooltipBorder: dark ? "#536eb7" : "#d0d7de",
    };
  }

  function formatThroughput(value) {
    if (value >= 1000) {
      return value.toFixed(0) + " Mbp/s";
    }
    if (value >= 100) {
      return value.toFixed(1) + " Mbp/s";
    }
    return value.toFixed(1) + " Mbp/s";
  }

  function speedupPlugin(colors) {
    return {
      id: "benchmarkSpeedupLabels",
      afterDatasetsDraw: function (chart) {
        var ctx = chart.ctx;
        var xScale = chart.scales.x;
        var yScale = chart.scales.y;
        if (!xScale || !yScale) {
          return;
        }

        ctx.save();
        ctx.fillStyle = COLORS.speedup;
        ctx.font = "600 16px system-ui, sans-serif";
        ctx.textAlign = "center";
        ctx.textBaseline = "bottom";

        for (var i = 0; i < DATA.labels.length; i++) {
          var x = xScale.getPixelForValue(i);
          var topY = yScale.getPixelForValue(
            Math.max(DATA.biokotlin[i], DATA.biopython[i])
          );
          ctx.fillText(DATA.speedups[i], x, topY - 10);
        }
        ctx.restore();
      },
    };
  }

  function buildConfig(colors) {
    return {
      type: "bar",
      data: {
        labels: DATA.labels,
        datasets: [
          {
            label: "BioKotlin",
            data: DATA.biokotlin,
            backgroundColor: COLORS.biokotlin,
            borderColor: COLORS.biokotlin,
            borderWidth: 1,
            borderRadius: 4,
          },
          {
            label: "BioPython",
            data: DATA.biopython,
            backgroundColor: COLORS.biopython,
            borderColor: COLORS.biopython,
            borderWidth: 1,
            borderRadius: 4,
          },
        ],
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        interaction: { mode: "index", intersect: false },
        plugins: {
          legend: {
            position: "top",
            labels: {
              color: colors.text,
              usePointStyle: true,
              pointStyle: "rectRounded",
            },
          },
          tooltip: {
            backgroundColor: colors.tooltipBg,
            titleColor: colors.text,
            bodyColor: colors.muted,
            borderColor: colors.tooltipBorder,
            borderWidth: 1,
            callbacks: {
              label: function (context) {
                var label = context.dataset.label || "";
                return label + ": " + formatThroughput(context.parsed.y);
              },
            },
          },
        },
        scales: {
          x: {
            ticks: { color: colors.text },
            grid: { display: false },
          },
          y: {
            type: "logarithmic",
            title: {
              display: true,
              text: "Millions of bases per second (log₁₀ scale)",
              color: colors.text,
              font: { size: 13, weight: "500" },
            },
            ticks: {
              color: colors.muted,
              callback: function (value) {
                var n = Number(value);
                if (n === 1 || n === 10 || n === 100 || n === 1000 || n === 10000) {
                  return n;
                }
                return "";
              },
            },
            grid: { color: colors.grid },
          },
        },
      },
      plugins: [speedupPlugin(colors)],
    };
  }

  function refreshTheme() {
    if (!chartInstance || typeof Chart === "undefined") {
      return;
    }
    var colors = colorScheme();
    chartInstance.options.plugins.legend.labels.color = colors.text;
    chartInstance.options.scales.x.ticks.color = colors.text;
    chartInstance.options.scales.y.title.color = colors.text;
    chartInstance.options.scales.y.ticks.color = colors.muted;
    chartInstance.options.scales.y.grid.color = colors.grid;
    chartInstance.options.plugins.tooltip.backgroundColor = colors.tooltipBg;
    chartInstance.options.plugins.tooltip.titleColor = colors.text;
    chartInstance.options.plugins.tooltip.bodyColor = colors.muted;
    chartInstance.options.plugins.tooltip.borderColor = colors.tooltipBorder;
    chartInstance.update("none");
  }

  function init() {
    var canvas = document.getElementById(CANVAS_ID);
    if (!canvas || typeof Chart === "undefined") {
      return;
    }

    var colors = colorScheme();
    chartInstance = new Chart(canvas.getContext("2d"), buildConfig(colors));

    var observer = new MutationObserver(refreshTheme);
    observer.observe(document.body, {
      attributes: true,
      attributeFilter: ["data-md-color-scheme"],
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
