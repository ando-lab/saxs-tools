"""Estimate the Earth Mover's Distance between model and data P(r) functions

Compute the regularized indirect Fourier transform (IFT) of the residual between the model and data scattering profiles,
and then compute the Wasserstein metric (Earth Mover's Distance) between the two P(r) functions.

The input is a .fit file containing q, I(q), error(I(q)), model(I(q)) columns, and the maximum dimension to use in the
IFT.

Optionally, the user can specify the smoothness regularization parameter, the number of points in real space, and the
 q-range to consider.

The output is the Wasserstein metric value, printed to standard output.
"""

# 2025-11-18, Steve, added optional --qmin --qmax

# TODO: consider adding option to output P(r) function, iqreg profile, etc.

import argparse

import numpy as np
import scipy.sparse as sp
from scipy.linalg import lstsq


def ift(q, intensity, intensity_error, dmax, smoothness, npoints):
    # x-axis
    r = np.linspace(0, dmax, npoints)
    dr = r[1] - r[0]

    # compute transform
    F = 4 * np.pi * dr * np.sinc(np.outer(q, r) / np.pi)
    F[:, 0] *= 0.5
    F[:, 1] *= 0.5

    # compute regularizing operator (smoothness)
    n = np.arange(0, npoints - 2)
    o = np.ones(npoints - 2)
    L = (
        sp.csr_matrix((-0.5 * o, (n, n)), shape=(npoints - 2, npoints))
        + sp.csr_matrix((o, (n, n + 1)), shape=(npoints - 2, npoints))
        + sp.csr_matrix((-0.5 * o, (n, n + 2)), shape=(npoints - 2, npoints))
    )

    # calculate the least-squares problem
    w = 1 / intensity_error
    A = sp.diags(w, 0) @ F
    b = sp.diags(w, 0) @ intensity
    AA = A.T @ A
    Ab = A.T @ b
    H = L.T @ L

    lambda_reg = smoothness * np.trace(AA) / H.trace()

    # solve the least squares problem
    pofr, residuals, rank, s = lstsq(AA + lambda_reg * H, Ab)

    iqreg = F @ pofr

    return r, pofr, iqreg


def wasserstein(r, pofr):
    dr = r[1] - r[0]
    cs = 4 * np.pi * dr * np.cumsum(pofr)
    ws = dr * np.sum(np.abs(cs))
    return ws


def load_fit(fit_file, normalize=True):
    values = np.loadtxt(fit_file, skiprows=1)
    q = values[:, 0]
    intensity = values[:, 1]
    intensity_error = values[:, 2]
    model = values[:, 3]

    if normalize is True:
        # compute I(0) of the model by Guinier extrapolation
        npts = 5
        y = np.log(model[:npts])
        x = q[:npts] ** 2
        coeffs = np.linalg.pinv(np.stack((x, np.ones_like(x))).T) @ y
        i0 = np.exp(coeffs[1])
        intensity /= i0
        intensity_error /= i0
        model /= i0

    return q, intensity, intensity_error, model


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fit_file", type=str, help="The .fit file to process")
    parser.add_argument("dmax", type=float, help="The maximum dimension for IFT")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("-s", "--smoothness", type=float, default=1.0, help="Smoothness for IFT")
    parser.add_argument("-n", "--numpoints", type=int, default=51, help="Number of real space points for IFT")
    parser.add_argument("--qmin", type=float, default=0.0, help="Minimum q-value")
    parser.add_argument("--qmax", type=float, default=float("inf"), help="Maximum q-value")
    args = parser.parse_args()

    if args.verbose:
        print(f"Loading model and data profiles from: {args.fit_file}")
    q, intensity, intensity_error, model = load_fit(args.fit_file)
    if args.verbose:
        print(f"  In file:   {len(q)} points from {q[0]:.6f} to {q[-1]:.6f}")
    idxmin = np.searchsorted(q, args.qmin, side="left")
    idxmax = np.searchsorted(q, args.qmax, side="right")

    q = q[idxmin:idxmax]
    intensity = intensity[idxmin:idxmax]
    intensity_error = intensity_error[idxmin:idxmax]
    model = model[idxmin:idxmax]

    if args.verbose:
        print(f"  Truncated: {len(q)} points from {q[0]:.6f} to {q[-1]:.6f}")

    residual = intensity - model
    residual_error = intensity_error

    x2_fit = (((residual) / residual_error) ** 2).mean()

    if args.verbose:
        print(f"  Calculating IFT (dmax={args.dmax}, numpoints={args.numpoints}, smoothness={args.smoothness})")

    r, pofr, iqreg = ift(q, residual, residual_error, args.dmax, args.smoothness, args.numpoints)
    x2 = (((iqreg - residual) / residual_error) ** 2).mean()

    emd = wasserstein(r, pofr)
    if args.verbose:
        print(f"  Reduced chi-squared of the model-data fit: {x2_fit:.4f}")
        print(f"  Reduced chi-squared of the regularized residual: {x2:.4f}")
        print(f"  Earth Mover's Distance (EMD): {emd:.4f} Å")
    else:
        print(f"emd (Å): {emd:.4f}, x2: {x2_fit:.4f}, x2_ift: {x2:.4f}")


if __name__ == "__main__":
    main()
