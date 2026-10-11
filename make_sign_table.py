"""
make_sign_table.py -- LaTeX table of the sign data (genus4_sign_rules_q8.json, with q = 4 for the types absent at q = 8)
inserted into paper/b28_torsor.tex at the marker %%SIGNTABLE%% (or replacing the previously generated table).
Usage: python make_sign_table.py
"""
import json
import os
import re

ORDER = ["7 1", "8", "5 2 1", "4 3 1", "6 1 1", "5 3", "5 1 1 1", "4 2 1 1", "3 2 2 1", "6 2", "3 3 2",
         "4 2 2/4 1 1 1 1", "4 4", "3 3 1 1", "3 2 1 1 1", "2 2 2 1 1", "2 2 2 2/2 2 1 1 1 1", "3 1 1 1 1 1", "2 1 1 1 1 1 1",
         "1 1 1 1 1 1 1 1"]


def tex_val(s):
    s = s.replace("i", "i")
    return "$" + s.replace("+", "+").replace("-", "-") + "$"


def load(q):
    fn = f"genus4_sign_rules_q{q}.json"
    return json.load(open(fn, encoding="utf-8")) if os.path.exists(fn) else {}


def main():
    data = {q: load(q) for q in (2, 4, 8)}
    rows = []
    seen = set()
    for q in (8, 4, 2):
        for cls in sorted(data[q], key=lambda c: (c.startswith("non"), ORDER.index(c) if c in ORDER else 99)):
            if cls in seen:
                continue
            seen.add(cls)
            dd = data[q][cls]
            first = True
            for d in sorted(dd, key=int):
                a = dd[d]
                nz = ", ".join(f"{k}" + (f"($\\times${v})" if len(a['nz']) > 1 else "") for k, v in sorted(a["nz"].items(), key=lambda kv: int(kv[0])))
                mag = int(round(a["mag"]))
                W = ", ".join(f"{tex_val(k)}" + (f"($\\times${v})" if len(a['W']) > 1 else "") for k, v in sorted(a["W"].items()))
                name = cls.replace(" ", "\\,") if first else ""
                if first and cls.startswith("non"):
                    name = "non-octad " + cls[cls.index("("):]
                rows.append(f"{name} & {q if first else ''} & {a['n'] if first else ''} & {d} & {a['f']} & {a['r']} & {nz} & {mag} & {W}\\\\")
                first = False
    table = r"""\begin{table}[p]
\centering\scriptsize\renewcommand{\arraystretch}{0.88}
\caption{Gauss sums $s_{2d}(\lambda)$ by Frobenius type, $d\mid n$, all parameters $(\lambda_5,\lambda_3)\in\F_q^2$ of the type and all $\lambda_1\in\F_q$ ($q=8$; a type absent at $q=8$ is shown at the largest $q\le8$ where it occurs). $f_d=\dim K_d$, $r_d=\dim\rad(\omega|_{K_d})$; ``$\#\lambda_1$'' is the number of $\lambda_1\in\F_q$ with $s_{2d}\ne0$ (with multiplicity over the parameters when it varies); $|s_{2d}|/q^d=2^{(f_d+r_d)/2}$ when nonzero; $W_d=\sum_{\lambda_1\in\F_q}s_{2d}/q^{d+1}$ with its multiplicity over the parameters of the type.}\label{tab:signs}
\begin{tabular}{llrrrrlrl}
\toprule
type & $q$ & \#par. & $d$ & $f_d$ & $r_d$ & $\#\lambda_1$ & $|s_{2d}|/q^d$ & $W_d$\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}"""
    tex = open("paper/b28_torsor.tex", encoding="utf-8").read()
    if "%%SIGNTABLE%%" in tex:
        tex = tex.replace("%%SIGNTABLE%%", "%%SIGNTABLE-BEGIN%%\n" + table + "\n%%SIGNTABLE-END%%")
    else:
        tex = re.sub(r"%%SIGNTABLE-BEGIN%%.*?%%SIGNTABLE-END%%", lambda m: "%%SIGNTABLE-BEGIN%%\n" + table + "\n%%SIGNTABLE-END%%", tex, flags=re.S)
    open("paper/b28_torsor.tex", "w", encoding="utf-8").write(tex)
    print(f"{len(rows)} rows written")


if __name__ == "__main__":
    main()
