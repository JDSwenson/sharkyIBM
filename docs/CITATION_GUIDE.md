# Citation and References in sharkyIBM

This document explains how citations are structured and used in the sharkyIBM package.

## File Structure

The package contains three files related to citations:

### 1. **CITATION** (root directory)
- R package citation metadata file
- Read by `citation("sharkyIBM")` command
- Contains package citation info plus key methodological references
- Users can access with: `citation("sharkyIBM")`
- Can export to BibTeX: `print(citation("sharkyIBM"), bibtex=TRUE)`

### 2. **inst/REFERENCES.bib**
- BibTeX bibliography database with full reference details
- Used as backup/reference for building citations
- Available to users who need to build a full bibliography
- Can be copied and used in Quarto, R Markdown, or LaTeX documents

### 3. **roxygen2 @references tags** (in R/*.R files)
- Each function documentation includes `@references` section
- These appear in the help files (e.g., `?create.stable.pop`)
- Generated from roxygen2 comments in the source code

## How to Use

### For package users:

**View package-level citations:**
```r
citation("sharkyIBM")
```

**View function-specific citations:**
```r
?create.stable.pop
?simulate.pop
?sample.pop
```
(Scroll to the "References" section at the bottom)

**Export to BibTeX:**
```r
print(citation("sharkyIBM"), bibtex=TRUE)
```

**Include in Quarto/R Markdown:**
```yaml
bibliography: /path/to/sharkyIBM/inst/REFERENCES.bib
```

### For package developers:

**To add/modify citations:**

1. Edit `CITATION` file for package-level references
2. Edit `inst/REFERENCES.bib` for BibTeX entries
3. Add/modify `@references` in roxygen2 comments in R/*.R files
4. Regenerate documentation: `devtools::document()`

**To update roxygen2 documentation:**
```r
devtools::document()
```
This regenerates the .Rd files in `man/` with updated references.

## Key References in sharkyIBM

### Population Dynamics Framework
- **Caswell, H. (2001).** Matrix Population Models. Foundation for Leslie matrix approach.
- **Hoyle & Maunder (2004).** Bayesian stock assessment. Density dependence framework.

### Density Dependence Mechanism
- **Pella & Tomlinson (1973).** Generalized stock production model. The PT compensation mechanism.

### Marine Mammal Life History
- **Barlow & Boveng (1991).** Age-specific mortality for marine mammals. Siler mortality model.

### Breeding Cycle / Half-Sibling Dynamics
- **Jacobson (2026, in prep).** Markov breeding cycles in cetacean dynamics. Coupling between calf survival and maternal breeding.

## Accessing Bibliographies

### Full package bibliography (BibTeX):
```bash
cat inst/REFERENCES.bib
```

### Individual function help:
```r
help(create.stable.pop)
help(simulate.pop)
help(sample.pop)
```

### Programmatic access:
```r
# Extract help as text
tools::Rd2txt(utils:::.getHelpFile(help(create.stable.pop)))
```

## Maintenance

When updating the package with new methods or significant changes:

1. Update `CITATION` with new references
2. Add entries to `inst/REFERENCES.bib` in BibTeX format
3. Add/modify `@references` roxygen2 tags in the relevant R files
4. Run `devtools::document()` to regenerate .Rd files
5. Commit changes to git

---

**Last updated:** September 2, 2026
