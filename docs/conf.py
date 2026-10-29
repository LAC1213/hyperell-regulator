"""Sphinx configuration for the package documentation."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

project = "hyperell-regulator"
author = "Lambert A'Campo"
copyright = "2026, Lambert A'Campo"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
]
root_doc = "index"
exclude_patterns = ["_build", "README.md"]
default_role = "math"  # Sage docstrings use single backticks for mathematics.
autodoc_member_order = "bysource"
autodoc_typehints = "description"
add_module_names = False

mathjax3_config = {
    "tex": {
        "macros": {
            "QQ": r"\mathbb{Q}",
            "ZZ": r"\mathbb{Z}",
            "RR": r"\mathbb{R}",
            "CC": r"\mathbb{C}",
            "PP": r"\mathbb{P}",
            "CH": r"\operatorname{CH}",
            "reg": r"\operatorname{reg}",
        },
    },
}

html_theme = "alabaster"
html_title = "hyperell-regulator documentation"
html_theme_options = {
    "description": "Periods, regulators and standard L-functions in SageMath",
    "github_user": "LAC1213",
    "github_repo": "hyperell-regulator",
    "github_button": False,
    "fixed_sidebar": True,
    "page_width": "1180px",
    "sidebar_width": "270px",
}
html_sidebars = {"**": ["about.html", "navigation.html", "searchbox.html"]}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
