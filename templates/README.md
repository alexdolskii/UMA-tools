# Shared FIA / UMA plate scheme

`FIA_96_well_Template.xlsx` is the empty shared **Plate Map / Channels /
PLOT LAYOUT** template, copied without changes from
[FIA-tools tech_dev, b547772](https://github.com/alexdolskii/FIA-tools/blob/b547772ae3bf735fb77d606632010d29b69af74e/templates/FIA_96_well_Template.xlsx).

Fill and save one `FIA_<experiment>.xlsx` beside the original images in each
JSON folder before running the five main UMA commands:

```bash
uma_alignment -i input_paths.json -a 15
uma_thickness -i input_paths.json
area_analysis -i input_paths.json -t 2000
uma_collect_results -i input_paths.json
uma_report -i input_paths.json --fn-threshold 20
```

Use `Fiber: Fibronectin` or another `Fiber: marker` at its actual C1–C4
position in each well. Different source channels are analyzed and plotted
separately. The same completed workbook can be used by current FIA tech_dev.
[UMA format, routing and migration rules](../README.md#shared-experimental-scheme).

This format applies to UMA's main `code` workflow (0.2.29+).
Functional assays keep the root `UMA_96_well_plate_template.xlsx` and their
existing discovery/order rules in this release.
