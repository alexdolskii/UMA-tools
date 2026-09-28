# Functional assay: stitching

Combine nine single-channel ND2 Z-stacks per well into one TIFF stack for
subsequent counting of green-labelled cancer cells. The full Z-stack is
preserved; this program does not count cells.

## Install and run

In your existing UMA environment, from the repository root:

```bash
conda activate uma_tools_new
python -m pip install --no-deps ./functional_assay
uma_stitching -i input_paths.json
```

Use your actual environment name. UMA-tools must already be installed in
that environment. No additional scientific dependencies are needed; follow
the main README's [Java/Fiji check](../README.md#verify-java-and-fiji-after-installation)
on a new computer. `uma_stitching --help` and `--version` do not start Fiji.

The default overlap is **32.8%**. Override it with `--overlap 30`.
The JSON uses the existing `folder_paths` list:

```json
{"folder_paths": ["/absolute/path/to/original_images"]}
```

Each folder may contain several wells. Each well requires exactly nine
unique frames `0000`–`0008`, named like
`sample__WellB02_PointB02_0000_ChannelFITC_Seq0000.nd2`.
Each file must contain one channel. The fixed 3 × 3 grid is:

| Row | Left | Center | Right |
|---|---|---|---|
| Top | 0000 | 0001 | 0002 |
| Middle | 0005 | 0004 | 0003 |
| Bottom | 0006 | 0007 | 0008 |

The original Linear Blending settings and **Sharpen on all slices** are
retained. Output: `Stitched_Results/WellB02_stitched.tif` for each well.
**A new run deletes and replaces the previous `Stitched_Results`** when a
folder contains a valid well. Original ND2 files are retained.

Hidden/`._` files are ignored; explicitly selected `._…json` files are
rejected. Incomplete or duplicate frame sets are reported and skipped.
Fiji startup and worker shutdown use the existing UMA runtime. The Java
heap cap remains 16 GiB, as in the supplied script.

Tests: `python -m unittest discover -s functional_assay/tests -v`.
