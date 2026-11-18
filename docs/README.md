# Documentation

Build in a SageMath environment with the package installed:

```sh
sage -pip install -r docs/requirements.txt
sage -python -m sphinx -W --keep-going -b html docs docs/_build/html
```

Preview at <http://127.0.0.1:8000>:

```sh
sage -python -m http.server 8000 --bind 127.0.0.1 --directory docs/_build/html
```

API pages read the package docstrings. Generate pages for additional modules with:

```sh
sage -python -m sphinx.ext.apidoc --separate -o docs/api hyperell_regulator
```

The `Documentation` workflow checks builds on pull requests and deploys
relevant pushes to `main`. Builds fail on Sphinx warnings.
Set **Settings → Pages → Build and deployment → Source → GitHub Actions**
to enable publishing at <https://lac1213.github.io/hyperell-regulator/>.
