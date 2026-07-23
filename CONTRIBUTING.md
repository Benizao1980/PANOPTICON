# Contributing

Please open an issue before substantial changes.

Include with bug reports:

- PANOPTICON commit
- complete command
- full traceback or scheduler error
- upstream tool and version
- available input filenames and sizes
- number of genomes and gene families
- operating system or HPC scheduler

Before submitting changes:

```bash
python -m py_compile scripts/*.py
pytest -q
```

New input formats should be implemented through the shared pangenome loader rather than parsed independently by downstream scripts.
