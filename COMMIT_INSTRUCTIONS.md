# Commit instructions

This package is designed to be copied over a local clone of PANOPTICON and reviewed before committing.

## 1. Create a branch

```bash
cd /path/to/PANOPTICON
git switch main
git pull
git switch -c ecology-v0.2
```

## 2. Copy package contents

From the directory containing this package:

```bash
rsync -av --exclude COMMIT_INSTRUCTIONS.md PANOPTICON_v0.2_commit_package/ /path/to/PANOPTICON/
```

The package renames the generic stratification script to:

```text
scripts/pangenome_stratification.py
```

Keep the old `pangenome_host_association.py` temporarily as a compatibility wrapper or delete it only after checking existing workflows.

## 3. Copy the repository logo

The README expects:

```text
panopticon_logo_variant2_partial.svg
```

Keep the existing repository logo in place.

## 4. Validate

```bash
conda activate panopticon
python -m py_compile scripts/*.py
pytest -q
```

Check CLI help:

```bash
python scripts/pangenome_io.py --help
python scripts/pangenome_plotting.py --help
python scripts/pangenome_stratification.py --help
python scripts/pangenome_ecology.py --help
```

## 5. Review changed files

```bash
git status
git diff --stat
git diff
```

## 6. Commit

```bash
git add README.md environment.yml CHANGELOG.md CONTRIBUTING.md \
  CITATION.cff LICENSE_CHOICE.md scripts docs examples tests .github

git commit -m "Add population-aware ecology and allele-sharing workflow"
git push -u origin ecology-v0.2
```

Open a pull request into `main`.

## 7. Wiki

The file:

```text
docs/Staphylococcus-aureus-host-association.md
```

can also replace the existing GitHub wiki page. GitHub wikis are separate Git repositories, so either paste the Markdown through the wiki editor or clone the wiki repository and replace the page there.
