"""Execute the comparison notebook in order and export local review artefacts.

Uses the capsule's existing NumPy/SciPy/Matplotlib environment; no Jupyter
runtime is required because the portable notebook contains ordinary Python.
"""
import argparse
import base64
import contextlib
import io
import json
import os
from pathlib import Path
import platform
import re
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import hyp2f1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    notebook = json.loads((root / "comparison_experiments.ipynb").read_text())
    original = json.loads((root / "provenance/javi-comparison-2026-09-09/comparison_experiments.original.ipynb").read_text())
    matplotlib.rcParams.update({"font.serif": ["Computer Modern Roman"],
                               "text.latex.preamble": r"\usepackage{amsmath}"})
    namespace = {"__name__": "__main__"}
    printed = {}
    figures = {6: "zeff_length_time", 9: "zeff_velocity_height", 13: "cattaneo_velocity_height"}
    exponents = []
    previous = Path.cwd()
    try:
        os.chdir(root)
        count = 0
        for index, cell in enumerate(notebook["cells"]):
            if cell["cell_type"] != "code" or not "".join(cell["source"]).strip():
                continue
            count += 1
            stream = io.StringIO()
            with contextlib.redirect_stdout(stream), warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="FigureCanvasAgg is non-interactive.*")
                exec(compile("".join(cell["source"]), f"comparison_experiments.ipynb:cell{index}", "exec"), namespace)
            cell["execution_count"] = count
            cell["outputs"] = []
            text = stream.getvalue()
            if text:
                printed[str(index)] = text
                cell["outputs"].append({"output_type": "stream", "name": "stdout", "text": text.splitlines(keepends=True)})
            if index in (4, 7, 11):
                beta, nu, alpha = (float(namespace[k]) for k in ("beta", "nu", "alpha"))
                z = (1 + np.cos(beta)) / 2
                independent_nu = brentq(lambda v: hyp2f1(v+1, -v, 1, z), 0, 1)
                assert abs(nu-independent_nu) < 1e-10
                exponents.append({"beta_deg": float(np.rad2deg(beta)), "nu": nu, "alpha": alpha,
                                  "independent_nu": independent_nu})
            if index in figures:
                figure = plt.gcf()
                figure.savefig(output / (figures[index]+".pdf"), bbox_inches="tight", dpi=300)
                figure.savefig(output / (figures[index]+".png"), bbox_inches="tight", dpi=160)
                encoded = base64.b64encode((output / (figures[index]+".png")).read_bytes()).decode()
                cell["outputs"].append({"output_type": "display_data", "data": {"image/png": encoded}, "metadata": {}})
        source_prints = {str(i): "".join("".join(o.get("text", [])) for o in c.get("outputs", []) if o.get("output_type")=="stream")
                         for i,c in enumerate(original["cells"]) if any(o.get("output_type")=="stream" for o in c.get("outputs", []))}
        matches = {key: printed.get(key)==value for key,value in source_prints.items()}
        # NumPy versions can wrap or format arrays differently. Compare every
        # reported number at the precision retained in the supplied output.
        number_pattern = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
        numeric_matches = {}
        for key, value in source_prints.items():
            expected = np.array([float(v) for v in re.findall(number_pattern, value)])
            actual = np.array([float(v) for v in re.findall(number_pattern, printed[key])])
            numeric_matches[key] = bool(expected.shape == actual.shape and np.allclose(expected,actual,rtol=1e-8,atol=1e-8))
        assert numeric_matches and all(numeric_matches.values()), numeric_matches
        datasets = {"zeff_length_time": np.column_stack((namespace["t"],namespace["lc"])),
                    "zeff_velocity_height": namespace["data_Fig4Zeff"],
                    "cattaneo_velocity_height": namespace["data_Cattaneo"]}
        for values in datasets.values():
            assert np.isfinite(values).all() and (values>0).all()
        summary = {"versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
                   "exponents": exponents, "stored_stdout_exact_matches": matches,
                   "stored_numeric_output_matches": numeric_matches,
                   "data_rows": {key:len(value) for key,value in datasets.items()},
                   "cattaneo": {"Oh_c": float(namespace["Ohc"]), "V_c_m_s": float(namespace["Vc"]),
                                "We_min": float(np.min(namespace["We"])), "We_max": float(np.max(namespace["We"]))},
                   "scope": "Arithmetic and supplied-output reproduction; manual parameters retained, no refitting or independent validation of published equations."}
        (output/"comparison_experiments.executed.ipynb").write_text(json.dumps(notebook,indent=1,ensure_ascii=False)+"\n")
        (output/"light-check.json").write_text(json.dumps(summary,indent=2)+"\n")
        print(json.dumps(summary,indent=2))
    finally:
        os.chdir(previous)
        plt.close("all")


if __name__ == "__main__":
    main()
