window.MathJax = {
  tex: {
    //packages: ["bm", "amsmath", "amsfonts", "amsthm"],
    tags: "ams",
    inlineMath: [
      ["$", "$"],
      ["\\(", "\\)"],
    ],
    macros: {
      bold: ['\\boldsymbol{#1}',1] ,     // this macro has one parameter
      ddx: ['\\frac{d #2}{d #1}', 2, 'x'], // this macro has an optional parameter that defaults to 'x'
      ppx: ['\\frac{\\partial #2}{\\partial #1}', 2, 'x'], 
    },
    environments: {
      eqn: ['\\begin{equation} \\begin{aligned}', '\\end{aligned}\\end{equation}'],
      eqn_n: ['\\begin{equation*} \\begin{aligned}', '\\end{aligned}\\end{equation*}'],
    }
  },
  options: {
    renderActions: {
      addCss: [
        200,
        function (doc) {
          const style = document.createElement("style");
          style.innerHTML = `
          .mjx-container {
            color: inherit;
          }
        `;
          document.head.appendChild(style);
        },
        "",
      ],
    },
  },
  chtml: {
    scale: 1,                      // global scaling factor for all expressions
    minScale: .25,                 // smallest scaling factor to use
    matchFontHeight: false,        // true to match ex-height of surrounding font
    mtextInheritFont: false,       // true to make mtext elements use surrounding font
    merrorInheritFont: true,       // true to make merror text use surrounding font
    mathmlSpacing: false,          // true for MathML spacing rules, false for TeX rules
    exFactor: .5,                  // default size of ex in em units
    displayAlign: 'center',        // default for indentalign when set to 'auto'
    displayIndent: '0',            // default for indentshift when set to 'auto'
    adaptiveCSS: true              // true means only produce CSS that is used in the processed equations
  },
};
