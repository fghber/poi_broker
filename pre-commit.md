# Pre-commit Setup Guide

## Overview
Pre-commit is a framework for managing Git hooks. It runs checks automatically before you commit code, ensuring consistent style, type safety, and security.

**Current Status**: Create a `.pre-commit-config.yaml` configured with Ruff, Mypy, and Bandit hooks. Supporting config files for Bandit and tool configurations for Ruff/Mypy in `pyproject.toml` need to be added as well.

## Installation
1. Install the tool:
   ```bash
   pip install pre-commit
   ```
   (or add `pre-commit` to `requirements-dev.txt` / `pyproject.toml`)

2. Install hooks into Git:
   ```bash
   pre-commit install
   ```

## Pre-Commit Configuration (`.pre-commit-config.yaml`)

At the repo root create a `.pre-commit-config.yaml` that contains this configuration:
```yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.5.0   # use latest Ruff release
    hooks:
      - id: ruff
        args: [--fix]   # auto-fix lint issues

  - repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.10.0   # use latest mypy release
    hooks:
      - id: mypy
        additional_dependencies: [types-requests]  # add stubs as needed

  - repo: https://github.com/PyCQA/bandit
    rev: 1.7.9
    hooks:
      - id: bandit
        args: ["-c", "bandit.yaml", "-r", "."]  # recursive scan
```

### 1. `bandit.yaml` — Required by current pre-commit config

The Bandit hook references `-c bandit.yaml`. Create it at repo root:
```yaml
# bandit.yaml
exclude_dirs:
  - tests
  - .venv
  - htmlcov
  - poi_broker_frontend.egg-info
skips:
  - B101  # assert_used
  - B601  # paramiko_calls (not used)
```

### 2. Ruff configuration in `pyproject.toml` — Recommended

Add to `pyproject.toml` for consistent linting:

```toml
[tool.ruff]
line-length = 100
target-version = "py312"
select = [
    "E",   # pycodestyle errors
    "W",   # pycodestyle warnings
    "F",   # pyflakes
    "I",   # isort
    "N",   # pep8-naming
    "UP",  # pyupgrade
    "B",   # flake8-bugbear
    "C4",  # flake8-comprehensions
    "T20", # flake8-print
]
ignore = [
    "E501",  # line too long (handled by formatter)
    "B008",  # function calls in default args (common in Flask)
]
per-file-ignores = {
    "tests/*" = ["S101", "S106", "S311"],  # test-specific allows
}

[tool.ruff.format]
quote-style = "double"
indent-style = "space"
skip-magic-trailing-comma = false
```

### 3. Mypy configuration in `pyproject.toml` — Recommended

Add to `pyproject.toml` for consistent type checking:
```toml
[tool.mypy]
python_version = "3.12"
warn_return_any = true
warn_unused_configs = true
disallow_untyped_defs = false
disallow_incomplete_defs = false
check_untyped_defs = true
no_implicit_optional = true
strict_optional = true
show_error_codes = true
pretty = true

[[tool.mypy.overrides]]
module = [
    "poi_broker.*",
]
ignore_missing_imports = true

[[tool.mypy.overrides]]
module = "tests.*"
ignore_missing_imports = true
```

## Usage

- Run hooks on staged files (default on commit):
  ```bash
  git commit -m "message"
  ```
- Run hooks manually on all files:
  ```bash
  pre-commit run --all-files
  ```
- Update hooks to latest versions:
  ```bash
  pre-commit autoupdate
  ```

## Benefits
- **Ruff** → linting & auto-fixes (fast, replaces flake8/isort/black)
- **Mypy** → static type checking
- **Bandit** → security scanning

## Adding Custom Hooks
You can add tools like pytest or additional linters:
```yaml
- repo: https://github.com/pytest-dev/pytest
  rev: 8.2.0
  hooks:
    - id: pytest
      args: ["-x", "-q"]
```

---

## 📝 Pre-commit Onboarding Checklist

Follow these steps to set up and use `pre-commit` in this project.

---

### 1. Install Pre-commit
- [ ] Ensure you have Python 3.12+ installed (project requirement).
- [ ] Install `pre-commit`:
  ```bash
  pip install pre-commit
  ```
- [ ] (Optional) Add `pre-commit` to your dev dependencies:
  - `requirements-dev.txt` (already has pytest stack)
  - or `pyproject.toml` under `[project.optional-dependencies] dev`

---

### 2. Configure Git Hooks
- [ ] Verify `.pre-commit-config.yaml` exists at the **repo root** (already present).
- [ ] **Create missing `bandit.yaml`** at repo root (see config above).
- [ ] **Add Ruff and Mypy config to `pyproject.toml`** (see config above).
- [ ] Install hooks into Git:
  ```bash
  pre-commit install
  ```
- [ ] Confirm installation:
  ```bash
  pre-commit run --all-files
  ```

---

### 3. Daily Workflow
- [ ] Stage changes with `git add`.
- [ ] Commit normally:
  ```bash
  git commit -m "Your message"
  ```
  → Hooks will run automatically before the commit is finalized.
- [ ] If a hook fails, fix issues and re-stage files.

---

### 4. Updating Hooks
- [ ] Periodically update hook versions:
  ```bash
  pre-commit autoupdate
  ```
- [ ] Re-run checks:
  ```bash
  pre-commit run --all-files
  ```

---

### 5. Common Hooks in This Project
- **Ruff** → Linting & auto-fixes (`--fix`)
- **Mypy** → Static type checking
- **Bandit** → Security scanning

---

### 6. Adding New Hooks
- [ ] Edit `.pre-commit-config.yaml` and add new repos/hooks.
- Example (pytest):
  ```yaml
  - repo: https://github.com/pytest-dev/pytest
    rev: 8.2.0
    hooks:
      - id: pytest
        args: ["-x", "-q"]
  ```
- [ ] Run:
  ```bash
  pre-commit install
  pre-commit run --all-files
  ```

---

## ✅ Final Check
- [ ] `bandit.yaml` exists at repo root.
- [ ] Ruff and Mypy configs added to `pyproject.toml`.
- [ ] You can commit without errors.
- [ ] Hooks run automatically on every commit.
- [ ] Team members share consistent linting, typing, and security checks.
