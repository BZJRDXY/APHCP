# APHCP

APHCP builds a reusable resource of experimental protein–ligand structures. The first stage catalogs Holo observations, their source revisions, extracted systems, and quality evidence. Apo–Holo matching and task-specific views are subsequent stages.

This repository contains the project source, workflows, migrations, tests, and public configuration templates. Production data, run records, and local deployment settings are managed separately.

- **Independent review:** [review guide](docs/public/review.md) · [verification snapshot](docs/public/status.md).
- `src/` contains the control interface and local Dashboard.
- `scripts/` contains metadata, structure, annotation and registration tools.
- `workflow/`, `migrations/`, and `tests/` contain workflow definitions, the SQLite schema and verification code.

For the local management tools, use the dependencies in `requirements-local.txt` and install the project with `python -m pip install -e .`. Server deployment requires site-specific settings; use `config/deployment.example.yaml` as the template. This public snapshot does not include the production registry or its run receipts.

Project code changes follow the [GitHub synchronization rule](docs/governance/代码同步规则.md).
