# saxs-tools

Simple command-line tools for SAXS data processing and analysis

## Tool Chest

### `saxs.emd`

Compute the Earth Mover's Distance (EMD, or Wasserstein-1 metric) between model and data P(r) functions.

#### Usage

```bash
# --help or -h flags to print usage info
saxs.emd --help
```

```
positional arguments:
  fit_file              The .fit file to process
  dmax                  The maximum dimension for IFT

options:
  -h, --help            show this help message and exit
  -v, --verbose         Enable verbose output (default: False)
  -s, --smoothness SMOOTHNESS
                        Smoothness for IFT (default: 1.0)
  -n, --numpoints NUMPOINTS
                        Number of real space points for IFT (default: 51)
  --qmin QMIN           Minimum q-value (default: 0.0)
  --qmax QMAX           Maximum q-value (default: inf)
```

#### Example 1: simple output

```bash
# download example dataset from the SASBDB
curl -O https://www.sasbdb.org/media/fitting_files/SASDPQ4_fit2.fit 

# Estimate the EMD
saxs.emd SASDPQ4_fit2.fit 92
```

```
emd (Å): 0.1737, x2: 40.6738, x2_ift: 2.6933
```

#### Example 2: verbose output, custom q-range

```bash
# --verbose or -v flag prints more information
# --qmin and --qmax to trim the q-range
saxs.emd SASDPQ4_fit2.fit 92 -v --qmin 0.02 --qmax 0.25
```

```
Loading model and data profiles from: SASDPQ4_fit2.fit
  In file:   498 points from 0.006000 to 1.000000
  Truncated: 116 points from 0.020000 to 0.250000
  Calculating IFT (dmax=92.0, numpoints=51, smoothness=1.0)
  Reduced chi-squared of the model-data fit: 134.2172
  Reduced chi-squared of the regularized residual: 4.4268
  Earth Mover's Distance (EMD): 0.1799 Å
```