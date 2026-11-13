# Recursive Subsample Testing for Co-Explosivity Breakdowns

**Taejun Kim · Job Market Paper · 2026**

This repository contains Python code from the working research script
accompanying *Recursive Subsample Testing for Co-Explosivity Breakdowns*.
The empirical example uses the Nasdaq Computer Index and Microsoft and
includes the author's functions `sBIC`, `sPSY_df`, `sPSY`, `weight_one`,
`recursive_residual`, and `recursive_residual2`.

> **Research status:** This is a minimally packaged version of the code
> supplied by the author. It is **not** a verified replication package for
> the paper. Numerical methods and research parameters have been preserved;
> the original input dataset was not provided.

## Files

- `analysis.py` — original research functions and empirical analysis, set
  up to run on any machine without the author's former Windows `os.chdir`.
- `requirements.txt` — Python dependencies.
- `data/README.md` — required input filename and CSV-column positions.
- `.gitignore` — excludes local data and generated files by default.

## Run in Spyder

1. Download and extract this repository.
2. In Spyder, open `analysis.py` (File → Open).
3. Install dependencies into the same Python environment used by Spyder:

   ```python
   %pip install -r requirements.txt
   ```

   If the IPython Console is not running in the project folder, instead
   install `numpy`, `pandas`, and `matplotlib` through Anaconda Prompt.
4. Add your dataset as **`data/index_vs_msft.csv`**. Confirm its layout using
   `data/README.md`.
5. Press **F5** to execute the analysis. Because the file locates the CSV
   relative to its own directory, no `os.chdir(...)` is needed.

## Run from a terminal

```bash
python -m pip install -r requirements.txt
python analysis.py
```

## Research implementation notes

The original script specifies a sign-based BIC procedure, an `sPSY`
statistic, recursive-residual calculations, bandwidth/lag settings, and
critical values. These are retained from the supplied script; no
independent econometric verification has been performed. The computation
can be intensive on longer time series.

Two non-methodological execution updates were made: (1) the empirical
analysis runs inside `main()` to prevent execution on import, and (2) the
plot handles the case where no rejection date is found. After filtering,
the date index is reset to keep positional plotting consistent.

## Citation

Kim, Taejun (2026). *Recursive Subsample Testing for Co-Explosivity
Breakdowns.* Job Market Paper (working draft).

## Related links

- [Academic website](https://teejuun123-dev.github.io/foreignworker/)
- [GitHub profile](https://github.com/teejuun123-dev)

## Data sharing and license

The `index_vs_msft.csv` dataset is not included. No software license has
been selected; if you intend to permit reuse, add a `LICENSE` file of your
choosing. Check with coauthors or supervisors before publishing code or
data that are not yet ready for public release.
