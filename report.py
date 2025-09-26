import os, html
import pint, math

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

class HtmlReport:
    def __init__(self, path="report.html", title="Calculation report"):
        self.path = path
        self._opened = False
        self._open(title)

    def _open(self, title):
        head = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>{html.escape(title)}</title>
<style>
 body {{ font-family: system-ui, Segoe UI, Roboto, Arial, sans-serif; line-height: 1.45; margin: 24px; }}
 h1 {{ font-size: 1.4rem; margin: 0 0 1rem; }}
 .item {{ margin: 0.35rem 0; }}
 .label {{ font-weight: 600; }}
 code {{ background: #f6f8fa; padding: 0.1rem 0.3rem; border-radius: 4px; }}
 hr {{ border: none; border-top: 1px solid #e5e5e5; margin: 1rem 0; }}
 .formula {{ margin: 0.5rem 0; }}
</style>
<!-- MathJax -->
<script id="MathJax-script" async
 src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
</head>
<body>
<h1>{html.escape(title)}</h1>
"""
        with open(self.path, "w", encoding="utf-8") as f:
            f.write(head)
        self._opened = True

    def write_html(self, raw_html):
        if not self._opened:
            self._open("Calculation report")
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(raw_html)

    def show_text(self, text):
        self.write_html(f"<div class='item'>{html.escape(text)}</div>\n")

    def show_kv(self, label, value_str):
        self.write_html(
            f"<div class='item'><span class='label'>{html.escape(label)}</span>: "
            f"{value_str}</div>\n"
        )

    def show_formula(self, latex_block):
        # support either inline or block latex; wrap in \[ \] if not already
        s = latex_block.strip()
        if not (s.startswith(r"\[") or s.startswith(r"$$")):
            s = r"\[" + s + r"\]"
        self.write_html(f"<div class='formula'>{s}</div>\n")

    def section(self, title):
        self.write_html(f"<hr/><h1>{html.escape(title)}</h1>\n")

    def close(self):
        if self._opened:
            with open(self.path, "a", encoding="utf-8") as f:
                f.write("\n</body>\n</html>")
            self._opened = False

# create the report once
REPORT = HtmlReport(path="report.html", title="Chimney Calculations")

# ================ helpers replacing display()/Markdown/Math =================
def _format_quantity(q):
    """Return HTML-safe string for a pint.Quantity or plain value."""
    if isinstance(q, ureg.Quantity):
        q = q.to_base_units()
        # special-case temperature to °C when it's K
        if q.units == ureg.kelvin:
            q = q.to("degC")
        mag = float(q.magnitude)
        if abs(mag) < 1e-2 or abs(mag) >= 1e4:
            mag_str = f"{mag:.2e}"
        else:
            mag_str = f"{mag:.3f}"
        return f"{mag_str} <code>{html.escape(str(q.units))}</code>"
    else:
        # plain numbers/strings
        if isinstance(q, float):
            return f"{q:.6g}"
        return html.escape(str(q))

def show_value(value, label=None):
    """Append a value to the HTML report."""
    label = label or "Value"
    REPORT.show_kv(label, _format_quantity(value))

def show_latex(latex_code):
    """Append a LaTeX formula (rendered by MathJax) to the HTML report."""
    REPORT.show_formula(latex_code)

def show_text(text):
    REPORT.show_text(text)

def end_report():
    REPORT.close()
    print(f"HTML report written to: {os.path.abspath(REPORT.path)}")