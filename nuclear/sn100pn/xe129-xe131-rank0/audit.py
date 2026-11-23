#!/usr/bin/env python3
"""Read-only AFH/SN100PN monopole response comparison; see README.md."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import tarfile
from urllib.request import urlopen

import numpy as np
from scipy.special import eval_genlaguerre, gamma, genlaguerre, hyp1f1, roots_genlaguerre, spherical_jn

A_VALUES = (129, 131)
Z = 54
TWO_J_PLUS_ONE = {129: 2, 131: 4}
AUTHORS_SHA = '1395f0c161905341fc89a2978148b1ce7c9ee1dd'
BASELINE_SHA = '5204f1a78b32b6ff5b423cc5d278d781e88add5e'
AFH_DIR = 'nuclear/afh/data/native_v6_2026-11-22/native'
RAW = 'https://raw.githubusercontent.com'
ARXIV_SOURCE = 'https://arxiv.org/src/2311.15764v1'
EXPECTED_SHA = {
    129: ('26ef95e96f79bcd883cc3374e26aac62fdb80be23439e1604a96703468e8dd91',
          '80684a3f1201bc312c1c60af3b5b263e1480ab5a05874743a738cce9b86dc532'),
    131: ('2dd36c7aa74a88f8fc979b2423ac2970a605300d0eca973d94cb139bda438bd0',
          '585ab0f483be0d8131566cefc4cb420c6db7e07bebbcd2b887da568cda2f4ffb'),
}
TEX_SHA = '00d6b232355a140d70b64fae5ef01c5a58203866b3137d9bae6a61e8ec21cc38'  # canonical member, not transport metadata
OCC_KEYS = ((4, 7), (4, 5), (4, 3), (4, 1), (5, 11))
OCC_REF = {
    129: [(2.41, 7.53), (.88, 6.00), (.33, 2.06), (.19, 1.26), (.19, 8.14)],
    131: [(2.80, 7.71), (.62, 5.60), (.22, 2.39), (.09, 1.33), (.26, 9.95)],
}
CHANNELS = ('M', 'PhiPP', 'MPhiPP')


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def load(directory, name, url, expected=None):
    path = directory / name
    if not path.exists():
        with urlopen(url, timeout=40) as handle:
            raw = handle.read()
        if expected is not None and sha256(raw) != expected:
            raise ValueError(f'SHA-256 mismatch for {url}')
        path.write_bytes(raw)
    else:
        raw = path.read_bytes()
        if expected is not None and sha256(raw) != expected:
            raise ValueError(f'SHA-256 mismatch for {path}')
    return path, raw


def occupations(raw, a):
    blocks = {}
    for line in raw.decode('ascii').splitlines():
        if line.startswith('THE'):
            rank, twice_tau = map(int, re.fullmatch(
                r'THE ONE-BODY DENSITY MATRIX FOR 2JO = (\d+) TO =, (\d+)', line).groups())
            if rank == 0:
                if twice_tau not in (0, 2) or twice_tau in blocks:
                    raise ValueError('invalid monopole OBDME header')
                blocks[twice_tau] = {}
        elif rank == 0:
            mark, na, ja, nb, jb, density = line.split()
            key = int(na), int(ja)
            if mark != '.' or key != (int(nb), int(jb)) or key in blocks[twice_tau]:
                raise ValueError(f'invalid monopole OBDME row: {line}')
            blocks[twice_tau][key] = float(density)
    if len(blocks[0]) != 16 or tuple(sorted(blocks[2])) != tuple(sorted(OCC_KEYS)):
        raise ValueError('unexpected core/valence content')
    t = (a - 2 * Z) / 2
    norm = TWO_J_PLUS_ONE[a] * (2 * t + 1)
    occ = {}
    for (n_major, twice_j), value in blocks[0].items():
        total = value * np.sqrt(2 * (twice_j + 1) / norm)
        diff = blocks[2].get((n_major, twice_j), 0) * np.sqrt(
            6 * (twice_j + 1) / norm * t / (t + 1))
        occ[n_major, twice_j] = np.array([(total - diff) / 2, (total + diff) / 2])
    return occ


def orbit(n_major, twice_j):
    candidates = [(twice_j - 1) // 2, (twice_j + 1) // 2]
    ell, = (l for l in candidates if 0 <= l <= n_major and (n_major - l) % 2 == 0)
    j = twice_j / 2
    return (n_major - ell) // 2, ell, j * (j + 1) - ell * (ell + 1) - 3/4


def radial(y, n, ell, method, nodes=0):
    alpha = ell + 1 / 2
    norm = float(gamma(n + 1) / gamma(n + alpha + 1))
    if method == 'laguerre':
        x, w = roots_genlaguerre(nodes, alpha)
        z = 2 * np.sqrt(y[:, None] * x)
        weights = norm * w * eval_genlaguerre(n, alpha, x)**2
        m = spherical_jn(0, z)
        phi = np.divide(spherical_jn(1, z), z, out=np.full_like(z, 1 / 3), where=z != 0)
        return np.stack([(m * weights).sum(axis=1), (phi * weights).sum(axis=1)])
    if method == 'analytic':
        # j0(z) = 0F1(3/2;-z^2/4), j1(z)/z = 0F1(5/2;-z^2/4)/3.
        # Integrate monomials of the Laguerre square using the Laplace transform to 1F1.
        square = np.polynomial.polynomial.polypow(genlaguerre(n, alpha).c[::-1], 2)
        p = alpha + 1 + np.arange(square.size)
        c = norm * square * gamma(p)
        return np.stack([sum(ci * hyp1f1(pi, b, -y) for ci, pi in zip(c, p)) / pref
                         for b, pref in ((1.5, 1), (2.5, 3))])
    raise ValueError(method)


def amplitudes(y, occ, method, nodes=0):
    amp = np.zeros((len(y), 4))
    for (n_major, twice_j), pn in occ.items():
        n, ell, sigma_l = orbit(n_major, twice_j)
        m, phi = radial(y, n, ell, method, nodes)
        amp[:, :2] += m[:, None] * pn
        amp[:, 2:] -= sigma_l * phi[:, None] * pn
    return amp


def printed_response(tex, a, y):
    # Read rounded Appendix II.G polynomials from original source; not production input.
    section = tex.split(r'\label{nuclear responses}', 1)[1].split(r'\subsection*{$^{' + str(a) + r'}$Xe}' + '\n')
    fields = section[1].split(r'\subsection*', 1)[0]
    symbols = {'M': 0, r"\Phi''": 2, r"M \Phi''": -1}
    pairs = {'p': 0, 'n': 1}
    matrix = np.zeros((y.size, 4, 4))
    seen = set()
    for label, p, n, terms in re.findall(
            r"\$F_\{([^{}]+)\}\^\{\(([pn]),([pn])\)\}\s*=\s*e\^\{-2\s*y\}\s*\(([^)]+)\)\s*\$",
            fields, flags=re.S):
        if label not in symbols:
            continue
        value = terms.replace(r'\times', '*').replace('10^{', '10**{')
        value = re.sub(r'\{([+-]?\d+)\}', r'\1', value)
        value = re.sub(r'\s+', '', value).replace('^', '**').replace('y', '*y')
        if re.fullmatch(r'[0-9.+*y-]+', value) is None:
            raise ValueError(f'unsupported response polynomial syntax for {a}: {terms}')
        def term_value(term):
            coeff, separator, power = term.partition('*y')
            factors = re.fullmatch(r'([+-]?[0-9.]+)(?:\*10\*\*([+-]?[0-9]+))?', coeff)
            if factors is None:
                raise ValueError(term)
            v = float(factors[1]) * 10**int(factors[2] or 0)
            return v if not separator else v * y**int(power.removeprefix('**') or 1)
        pieces = re.split(r'(?<!\*)(?=[+-])', value)
        r = np.exp(-2*y)*sum(term_value(t) for t in pieces if t)
        if symbols[label] == -1:
            i, j = pairs[p], 2 + pairs[n]
        else:
            i, j = symbols[label] + pairs[p], symbols[label] + pairs[n]
        matrix[:, i, j] = matrix[:, j, i] = r
        seen.add((i, j))
    if len(seen) != 10:
        raise ValueError(f'missing printed Appendix II.G responses for {a}: {seen}')
    return matrix


def iso_response(normalized):
    m = normalized[:, :2, :2].sum(axis=(1, 2))
    p = normalized[:, 2:, 2:].sum(axis=(1, 2))
    i = normalized[:, :2, 2:].sum(axis=(1, 2))
    return np.column_stack((m, p, i))


def run(inputs, output):
    inputs.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    archive, _ = load(inputs, '2311.15764v1', ARXIV_SOURCE)
    with tarfile.open(archive, 'r:gz') as tar:
        tex_raw = tar.extractfile('tek/bodysupplementary.tex').read()
    if sha256(tex_raw) != TEX_SHA:
        raise ValueError('Supplementary TeX member SHA-256 mismatch')
    tex = tex_raw.decode()
    summary = {'scope': 'provisional SN100PN rank 0 only, no cross-section or package promotion',
               'authors_sha': AUTHORS_SHA, 'afh_sha': BASELINE_SHA,
               'supplement_tex_sha256': TEX_SHA, 'nodes': [32, 80], 'isotopes': {}}
    for a in A_VALUES:
        obdme_name = f'{a}Xe_OBDME_SN100PN.dat'
        afh_name = f'afh_v6_xe{a}_eight_response_native.csv'
        _, raw = load(inputs, obdme_name, f'{RAW}/Raghda-AK/DM-Form-Factor/{AUTHORS_SHA}/{obdme_name}', EXPECTED_SHA[a][0])
        _, base_raw = load(inputs, afh_name,
                           f'{RAW}/lm846802297-bit/radiative-dark-matter-nuclear-currents/{BASELINE_SHA}/{AFH_DIR}/{afh_name}',
                           EXPECTED_SHA[a][1])
        base = np.genfromtxt(base_raw.decode().splitlines(), names=True, delimiter=',')
        q, y = base['q_GeV'], base['y']
        assert (len(q), q[0], q[-1]) == (2401, 0, 0.6)
        assert np.max(abs(q - np.arange(2401)/4000)) < 2e-15
        assert np.max(abs(y - (q * 5.0677/2)**2 * (41.467 / (45*a**(-1/3)-25*a**(-2/3))))) < 2e-13
        occ = occupations(raw, a)
        nr_error = np.array(list(occ.values())).sum(axis=0) - [Z, a-Z]
        occ_error = np.array([occ[k] for k in OCC_KEYS]) - OCC_REF[a]
        assert np.max(abs(nr_error)) < 0.04  # published input deficit, not quadrature error
        assert np.max(abs(occ_error)) < 0.011  # supplement independently printed to 2 decimals
        for (n, two_j), pn in occ.items():
            assert np.min(pn) > -1e-6 and np.max(pn) < two_j + 1 + 2e-6
        amp = amplitudes(y, occ, 'laguerre', 80)
        quad_delta = np.max(np.abs(amp - amplitudes(y, occ, 'laguerre', 32)), axis=0)
        analytic_delta = np.max(np.abs(amp - amplitudes(y, occ, 'analytic')), axis=0)
        assert max(quad_delta) < 3e-11 and max(analytic_delta) < 3e-11
        assert np.allclose(radial(y, 0, 0, 'analytic')[0], np.exp(-y), rtol=2e-13, atol=1e-14)
        assert np.allclose(amp[0, :2], np.array(list(occ.values())).sum(axis=0), rtol=0, atol=3e-12)
        assert np.allclose(amp[0, 2:], -sum(orbit(*k)[-1]*v for k, v in occ.items())/3, rtol=0, atol=3e-12)
        exact = amp[:, :, None] * amp[:, None, :]
        printed = printed_response(tex, a, y)
        eig_exact = np.linalg.eigvalsh(exact)
        eig_printed = np.linalg.eigvalsh(printed)
        assert eig_exact.min() > -5e-12
        assert np.max(abs(printed[0] - exact[0]) / np.maximum(abs(exact[0]), 1)) < 0.005
        # There is no A^2/4-only native normalization: the raw AFH W also retains (2Ji+1)/4pi.
        native_to_physical = 16*np.pi / TWO_J_PLUS_ONE[a]
        afh_equal = native_to_physical*np.column_stack([base[f'W_{k}00_AFH'] for k in CHANNELS])
        assert abs(afh_equal[0, 0]/a**2 - 1) < 2e-4
        sn_equal = iso_response(exact)
        p_equal = iso_response(printed)
        err = np.abs(p_equal-sn_equal)
        sn_zero = int(np.flatnonzero(np.signbit(amp[:-1, :2].sum(axis=1)) != np.signbit(amp[1:, :2].sum(axis=1)))[0])
        idx = [int(round(x*4000)) for x in [0, .05, .10, .15, .20, .30, .40, .60]]
        info = {
            'input_sha256': dict(zip((obdme_name, afh_name), EXPECTED_SHA[a])),
            'counts_p_n': amp[0, :2].tolist(), 'count_error_p_n': nr_error.tolist(),
            'max_abs_printed_valence_occupation_error': float(np.max(abs(occ_error))),
            'max_abs_amplitude_delta_32_vs_80': quad_delta.tolist(),
            'max_abs_amplitude_delta_analytic_vs_80': analytic_delta.tolist(),
            'min_reconstructed_gram_eigenvalue': float(eig_exact.min()),
            'AFH_q0_native_W_M00': float(base['W_M00_AFH'][0]),
            'AFH_q0_physical_F_M_minus_A2': float(afh_equal[0, 0] - a**2),
            'printed_min_gram_eigenvalue': float(eig_printed.min()),
            'printed_negative_min_eigenvalue_nodes_at_1e-8': int((eig_printed[:, 0] < -1e-8).sum()),
            'printed_negative_diagonal_nodes_at_1e-8': int((np.min(np.diagonal(printed, axis1=1, axis2=2), axis=1) < -1e-8).sum()),
            'max_abs_printed_minus_rebuilt_equalpn_over_AFH_q0_charge': (err.max(axis=0)/afh_equal[0, 0]).tolist(),
            'first_zero_of_equal_pn_M0_amplitude_bracket_GeV': q[[sn_zero,sn_zero+1]].tolist(),
            'landmarks': [{'q_GeV': q[k], 'SN100PN_rank0': sn_equal[k].tolist(),
                           'AFH_BD_all_allowed_multipoles': afh_equal[k].tolist(),
                           'SN100PN_printed_all_allowed_multipoles': p_equal[k].tolist()} for k in idx],
        }
        summary['isotopes'][a] = info
        upper = [(i, j) for i in range(4) for j in range(i, 4)]
        tags = ('Mp','Mn','PhiPPp','PhiPPn')
        with (output/f'xe{a}_rank0.csv').open('w', newline='') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerow(['q_GeV','y']+[f'SN100PN_rank0_amp_{k}' for k in tags]
                            +[f'SN100PN_rank0_F_{tags[i]}_{tags[j]}' for i, j in upper]
                            +[f'printed_F_{tags[i]}_{tags[j]}' for i, j in upper]
                            +[f'AFH_BD_F_{c}_eqpn_full' for c in CHANNELS]+['printed_min_eig'])
            writer.writerows([q[k], y[k], *amp[k], *(exact[k,i,j] for i,j in upper),
                              *(printed[k,i,j] for i,j in upper), *afh_equal[k], eig_printed[k,0]]
                             for k in range(q.size))
    (output/'summary.json').write_text(json.dumps(summary, indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', required=True, type=Path, help='cache for pinned source files')
    parser.add_argument('--output', required=True, type=Path)
    options = parser.parse_args()
    run(options.inputs, options.output)
