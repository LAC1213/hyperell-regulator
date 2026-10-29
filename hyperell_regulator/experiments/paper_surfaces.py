"""Add independently validated display meshes to the saved paper explorer.

    sage -python -m hyperell_regulator.experiments.paper_surfaces --jobs 4

No regulator or line integral is recomputed or changed.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", type=Path, default=Path("docs/_extra/explorer"))
    parser.add_argument("--density", type=int, default=20)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--worker", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.worker:
        from hyperell_regulator.regulator.surfaces import surface_mesh
        source = args.site/args.worker
        data = json.loads(source.read_text())
        path = args.site/"meshes"/source.name
        module = Path(__file__).resolve().parents[1]/"regulator"/"surfaces.py"
        input_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        source_hash = hashlib.sha256(module.read_bytes()).hexdigest()
        if path.exists():
            cached = json.loads(path.read_text())
            if (cached.get("integral_data_sha256") == input_hash
                    and cached.get("surface_source_sha256") == source_hash
                    and cached.get("density") == args.density):
                print(source.stem, "validated cached mesh", flush=True)
                return
        output = surface_mesh(data, args.density)
        output["integral_data_sha256"] = input_hash
        output["surface_source_sha256"] = source_hash
        output["density"] = args.density
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(output, separators=(",", ":"), allow_nan=False)+"\n")
        print(source.stem, output["validation"], flush=True)
        return
    catalog_path = args.site/"catalog.json"
    catalog = json.loads(catalog_path.read_text())
    functions = [f for c in catalog["curves"] for f in c["functions"]]
    def run(f):
        command = [sys.executable, "-m", "hyperell_regulator.experiments.paper_surfaces",
                   "--site", str(args.site), "--density", str(args.density), "--worker", f["file"]]
        log = Path("/tmp")/("hyperell-surface-"+Path(f["file"]).stem+".log")
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
        with log.open("w") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                    env=env, timeout=1200, check=False)
        if result.returncode:
            raise RuntimeError(log.read_text())
        print(log.read_text().strip(), flush=True)
        return f
    errors = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run, f) for f in functions]
        for future in as_completed(futures):
            try:
                f = future.result()
                f["mesh"] = "meshes/"+Path(f["file"]).name
            except Exception as exc:
                errors.append(str(exc)); print(str(exc), file=sys.stderr, flush=True)
    if errors:
        raise SystemExit("%d surface exports need attention" % len(errors))
    catalog_path.write_text(json.dumps(catalog, indent=2)+"\n")


if __name__ == "__main__":
    main()
