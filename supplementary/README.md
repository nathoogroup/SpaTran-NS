# Supplementary Material Package

This directory contains a Genome Biology-oriented supplementary package for the
SpaTran-NS manuscript.

## Main files

- `Supplementary_Information.tex` and `Supplementary_Information.pdf`: source and
  compiled Supplementary Information document.
- `additional_files/`: submission-ready Additional files 1-10.
- `additional_file_descriptions_for_manuscript.tex`: LaTeX block to paste into
  the main manuscript so each uploaded additional file has a caption/description.
- `availability_and_code_declarations_draft.tex`: draft declarations text for
  Genome Biology's required data/code availability sections.
- `Genome_Biology_submission_checklist.md`: remaining items to verify before
  submission.

## Regenerating

Run from the repository root:

```bash
Rscript supplementary/scripts/create_supplementary_materials.R
cd supplementary
latexmk -pdf -interaction=nonstopmode -halt-on-error Supplementary_Information.tex
cp Supplementary_Information.pdf additional_files/Additional_file_1_Supplementary_Information.pdf
latexmk -c Supplementary_Information.tex
```

The generator script reads current files in `results/`, `simulation/output/`,
and `data/spatial_datasets/dataset_manifest.csv`, then refreshes the CSV
additional files and generated LaTeX table fragments.
