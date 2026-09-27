"""Run bounded regressions under Sage without requiring pytest.

Usage: sage -python -m tests.run_regressions [--skip TEST_NAME_FRAGMENT]
"""
import argparse
import importlib
import time

MODULES = [
    'test_smoke', 'test_diagnostics', 'test_walk_coordinates',
    'test_endpoint_ramification', 'test_weierstrass_endpoint',
    'test_quadrature_refinement', 'test_finite_common_pole',
    'test_numerical_failures', 'test_regulator_assembly', 'test_infinity_ray',
    'test_beilinson_float64', 'test_pinch_continuation',
    'test_389_603_continuation', 'test_beilinson_report',
    'test_arbitrary_precision', 'test_family_n', 'test_latex_report',
]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip', action='append', default=[],
                        help='skip function names containing this fragment')
    args = parser.parse_args(argv)
    passed, skipped = [], []
    start = time.monotonic()
    for name in MODULES:
        module = importlib.import_module('tests.' + name)
        for attr in sorted(dir(module)):
            if not attr.startswith('test_'):
                continue
            label = name + '.' + attr
            if any(fragment in attr for fragment in args.skip):
                skipped.append(label)
                print('SKIP', label, flush=True)
                continue
            print('RUN', label, flush=True)
            getattr(module, attr)()
            passed.append(label)
            print('PASS', label, flush=True)
    print('TOTAL', len(passed), 'passed;', len(skipped), 'skipped;',
          round(time.monotonic()-start, 2), 'seconds', flush=True)


if __name__ == '__main__':
    main()
