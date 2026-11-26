import platform
import unittest

import numpy as np

METRIC = np.asarray([1.0, -1.0, -1.0, -1.0])


def bilinear(L, A, B, JA):
    return np.einsum(
        "mn,mfi,nfi->",
        L,
        METRIC[:, None, None] * A,
        (METRIC[:, None, None] * B).conjugate(),
    ) / (2.0 * JA + 1.0)


def contract_afh_continuity_projection(L, current, JA, alpha):
    full = np.asarray(current["full"], dtype=complex)
    if L.shape != (4, 4) or full.ndim != 3 or full.shape[0] != 4:
        raise ValueError("invalid Lorentz tensor/current dimensions")
    rho, jz = full[0], full[3]
    d = jz - alpha * rho
    long_actual = np.zeros_like(full)
    long_actual[0], long_actual[3] = rho, jz
    long_projected = np.zeros_like(full)
    long_projected[0], long_projected[3] = rho, alpha * rho
    transverse = np.zeros_like(full)
    transverse[1:3] = full[1:3]
    k0z_full = bilinear(L, long_actual, long_actual, JA)
    k0z_lc = bilinear(L, long_projected, long_projected, JA)
    k_transverse = bilinear(L, transverse, transverse, JA)
    k_lt_projected = 2.0 * bilinear(L, long_projected, transverse, JA).real
    k_lt_actual = 2.0 * bilinear(L, long_actual, transverse, JA).real
    k_full = bilinear(L, full, full, JA)
    residual_terms = (
        -(L[0, 3] - alpha * L[3, 3]) * rho * d.conjugate()
        -(L[3, 0] - alpha * L[3, 3]) * d * rho.conjugate()
        + L[3, 3] * np.abs(d) ** 2
    )
    residual_formula = float(np.real(np.sum(residual_terms)) / (2.0 * JA + 1.0))
    residual_direct = float(np.real(k0z_full - k0z_lc))
    formula_relative_residual = abs(residual_direct - residual_formula) / max(
        abs(residual_direct), abs(residual_formula), 1.0
    )
    k_afh6 = float(np.real(k0z_lc + k_transverse))
    k_projected_lt = float(k_lt_projected)
    k_nonconserved_lt = float(k_lt_actual - k_lt_projected)
    k_nonafh_residual = float(k_full.real - k_afh6)
    reconstructed = k_afh6 + k_projected_lt + residual_direct + k_nonconserved_lt
    decomposition_relative_residual = abs(float(k_full.real) - reconstructed) / max(
        abs(float(k_full.real)), abs(reconstructed), 1.0
    )
    return {
        "K_full_raw": float(k_full.real),
        "K_LC_longitudinal_projection_raw": float(k0z_lc.real),
        "K_transverse_xy_projection_raw": float(k_transverse.real),
        "K_AFH6_projection_raw": k_afh6,
        "K_projected_LT_cross_raw": k_projected_lt,
        "K_longitudinal_0z_nonconserved_residual_raw": residual_direct,
        "K_longitudinal_transverse_nonconserved_residual_raw": k_nonconserved_lt,
        "K_nonAFH_projection_residual_raw": k_nonafh_residual,
        "nonconserved_formula_relative_residual": float(formula_relative_residual),
        "projection_decomposition_relative_residual": float(decomposition_relative_residual),
    }


class ContinuityProjectionTests(unittest.TestCase):
    def test_conserved_current_reduces_to_lc(self):
        rng = np.random.default_rng(20261125)
        for JA in (0.0, 0.5, 1.5):
            n = int(round(2 * JA)) + 1
            x = rng.normal(size=(4, 7)) + 1j * rng.normal(size=(4, 7))
            L = x @ x.conjugate().T
            current = rng.normal(size=(4, n, n)) + 1j * rng.normal(size=(4, n, n))
            alpha = 0.0047
            current[3] = alpha * current[0]
            got = contract_afh_continuity_projection(L, {"full": current}, JA, alpha)
            self.assertLess(abs(got["K_longitudinal_0z_nonconserved_residual_raw"]), 1e-12)
            self.assertLess(got["nonconserved_formula_relative_residual"], 1e-12)
            self.assertLess(got["projection_decomposition_relative_residual"], 1e-12)
            lc = np.real(L[0, 0] - alpha * (L[0, 3] + L[3, 0]) + alpha * alpha * L[3, 3])
            expected = lc * np.sum(np.abs(current[0]) ** 2) / (2.0 * JA + 1.0)
            self.assertLess(abs(got["K_LC_longitudinal_projection_raw"] - expected), 1e-11 * max(abs(expected), 1.0))

    def test_nonconserved_residual_matches_independent_formula(self):
        rng = np.random.default_rng(20261126)
        for JA in (0.0, 0.5, 1.5):
            n = int(round(2 * JA)) + 1
            x = rng.normal(size=(4, 9)) + 1j * rng.normal(size=(4, 9))
            L = x @ x.conjugate().T
            current = rng.normal(size=(4, n, n)) + 1j * rng.normal(size=(4, n, n))
            alpha = 0.011
            rho, jz = current[0], current[3]
            d = jz - alpha * rho
            lower = np.stack((rho, -jz))
            block = L[np.ix_([0, 3], [0, 3])]
            k_full_0z = np.real(np.einsum("ab,afi,bfi->", block, lower, lower.conjugate())) / (2.0 * JA + 1.0)
            lc = np.real(L[0, 0] - alpha * (L[0, 3] + L[3, 0]) + alpha * alpha * L[3, 3])
            k_lc = lc * np.sum(np.abs(rho) ** 2) / (2.0 * JA + 1.0)
            expected = (
                -(L[0, 3] - alpha * L[3, 3]) * rho * d.conjugate()
                -(L[3, 0] - alpha * L[3, 3]) * d * rho.conjugate()
                + L[3, 3] * np.abs(d) ** 2
            )
            expected = np.real(np.sum(expected)) / (2.0 * JA + 1.0)
            got = contract_afh_continuity_projection(L, {"full": current}, JA, alpha)
            self.assertLess(abs(got["K_longitudinal_0z_nonconserved_residual_raw"] - (k_full_0z-k_lc)), 1e-11 * max(abs(k_full_0z), abs(k_lc), 1.0))
            self.assertLess(abs(got["K_longitudinal_0z_nonconserved_residual_raw"] - expected), 1e-11 * max(abs(expected), 1.0))
            self.assertLess(got["nonconserved_formula_relative_residual"], 1e-12)
            self.assertLess(got["projection_decomposition_relative_residual"], 1e-12)

    def test_full_kernel_reconstructs_from_projection_and_residuals(self):
        rng = np.random.default_rng(20261127)
        for JA in (0.0, 0.5, 1.5):
            n = int(round(2 * JA)) + 1
            x = rng.normal(size=(4, 8)) + 1j * rng.normal(size=(4, 8))
            L = x @ x.conjugate().T
            current = rng.normal(size=(4, n, n)) + 1j * rng.normal(size=(4, n, n))
            got = contract_afh_continuity_projection(L, {"full": current}, JA, 0.006)
            reconstructed = (
                got["K_AFH6_projection_raw"]
                + got["K_projected_LT_cross_raw"]
                + got["K_longitudinal_0z_nonconserved_residual_raw"]
                + got["K_longitudinal_transverse_nonconserved_residual_raw"]
            )
            self.assertLess(abs(got["K_full_raw"] - reconstructed), 1e-12 * max(abs(got["K_full_raw"]), 1.0))
            self.assertLess(got["projection_decomposition_relative_residual"], 1e-12)

    def test_20000_random_hermitian_psd_tensor_stress(self):
        rng = np.random.default_rng(20261125)
        max_formula = 0.0
        max_decomp = 0.0
        max_conserved = 0.0
        for _ in range(20000):
            x = rng.normal(size=(4, 7)) + 1j * rng.normal(size=(4, 7))
            L = x @ x.conjugate().T
            alpha = float(rng.uniform(0.0, 0.02))
            J = rng.normal(size=(4, 4, 3)) + 1j * rng.normal(size=(4, 4, 3))
            JA = float(rng.choice([0.0, 0.5, 1.5]))
            got = contract_afh_continuity_projection(L, {"full": J}, JA, alpha)
            max_formula = max(max_formula, got["nonconserved_formula_relative_residual"])
            max_decomp = max(max_decomp, got["projection_decomposition_relative_residual"])
            Jc = J.copy()
            Jc[3] = alpha * Jc[0]
            cons = contract_afh_continuity_projection(L, {"full": Jc}, JA, alpha)
            max_conserved = max(max_conserved, abs(cons["K_longitudinal_0z_nonconserved_residual_raw"]) / max(abs(cons["K_full_raw"]), 1.0))
        self.assertLess(max_formula, 1e-12)
        self.assertLess(max_decomp, 1e-12)
        self.assertLess(max_conserved, 1e-12)
        print({"samples": 20000, "max_formula_relative_residual": max_formula, "max_decomposition_relative_residual": max_decomp, "max_conserved_current_residual": max_conserved})


if __name__ == "__main__":
    unittest.main(verbosity=2)
