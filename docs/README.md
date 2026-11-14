# Building the documentation

From the repository root, in a SageMath environment with the package's
dependencies installed:

```sh
sage -pip install -r docs/requirements.txt
sage -python -m sphinx -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html`, or serve the result locally:

```sh
sage -python -m http.server 8000 --bind 127.0.0.1 --directory docs/_build/html
```

Then visit <http://127.0.0.1:8000>.

API pages were generated with:

```sh
sage -python -m sphinx.ext.apidoc --separate -o docs/api hyperell_regulator
```

The package overview pages deliberately leave members to their defining
modules, to avoid duplicate API entries. Sphinx reads current docstrings on
every build; only newly added modules need new reference pages.

## GitHub Pages

The `Documentation` workflow in `.github/workflows/docs.yml` builds the HTML
in a SageMath environment. Warnings fail the build. Pull requests build the
documentation without publishing; relevant pushes to `main` build and deploy
to GitHub Pages. The workflow can also be run manually on `main`.

Enable this once in the repository's **Settings → Pages → Build and deployment
→ Source → GitHub Actions**. No `gh-pages` branch or personal access token is
needed. Deployment uses the `github-pages` environment and the workflow's
GitHub token.

The published documentation is at
<https://lac1213.github.io/hyperell-regulator/>.
