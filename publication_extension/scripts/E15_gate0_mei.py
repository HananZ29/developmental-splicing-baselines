"""E15 Gate 0 at the minimum effect of interest (decision D3-B, 2026-11-06).

Writes publication_extension/results/E15_devstage/gate0_mei.json.

Provenance: this file is a verbatim transcription of the inline code that produced gate0_mei.json
on 2026-11-06 (session transcript, 19:31 UTC). It was saved as a script on 2026-11-07 for
reproducibility and has not been re-executed since; gate0_mei.json was not regenerated.

MEI = 0.10 * sqrt(2 * (1 - r)), where r is the Pearson correlation of rank-normal-transformed F1 and
A1 predictions (gate0_power.json). The fetal-favourable scenario is calibrated so that rho_F1 = 0.20
and delta_rho(F1 - A1) = MEI; the rate of declaring F1 better is the power at the MEI.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import E15_gate0_power as g  # noqa: E402

R = str(Path(__file__).resolve().parents[1] / "results" / "E15_devstage") + "/"


def main():
    r = json.load(open(R + 'gate0_power.json'))['pearson_rank_normal_F1_A1']
    mei = round(0.10 * np.sqrt(2 * (1 - r)), 4)
    df = pd.read_csv(R + 'predictions.tsv', sep='\t').rename(columns={'component_P': 'component'})
    X = {B: g.rank_normal(df[f'pred_{B}'].to_numpy()) for B in ('F1', 'A1', 'G', 'U')}
    comps = pd.factorize(df.component)[0]
    rng = np.random.default_rng(g.SEED + 1)
    c = g.calibrate(X, 0.20, mei, rng)
    res = g.simulate(X, comps, {'F1': c[0], 'A1': c[1]}, rng=rng)
    out = {'decision': 'D3-B (user, 2026-11-06)', 'r_rank_normal_F1_A1': r, 'MEI': mei,
           'definition': '0.10*sqrt(2*(1-r))',
           'scenario': 'fetal-favourable, rho_F1=0.20, delta_rho=MEI',
           'mean_delta_rho': res['F1-A1']['mean_delta_rho'],
           'power_F1-A1': res['F1-A1']['rate_declare_F1_better'], 'mc_se': res['F1-A1']['mc_se'],
           'typeI_max_from_gate0_power.json': 0.018,
           'verdict': 'PASS' if res['F1-A1']['rate_declare_F1_better'] >= 0.8 else 'FAIL',
           'R_SIMS': g.R_SIMS, 'B_BOOT': g.B_BOOT, 'seed': g.SEED + 1}
    open(R + 'gate0_mei.json', 'w').write(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
