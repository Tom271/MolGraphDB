# python-project-template

A template for projects in Python

## Getting Started

If you are setting up the project for the first time, you will need to install the following:

1. [uv](https://docs.astral.sh/uv/getting-started/installation/), a Python project and package manager built in Rust.
2. [python] (https://apps.microsoft.com/detail/9ncvdn91xzqp?launch=true&mode=full&hl=en-us&gl=gb&ocid=bingwebsearch) will be used to install pipx
3. [pipx](https://pypi.org/project/pipx/), will be used to install pre-commit
4. [pre-commit](https://pre-commit.com/), will be used to manage pre-commit hooks. ()

What Needs to Be Installed Again for Each New Project

1. Pre-commit hooks (pre-commit install).
2. Virtual environment (uv sync and uv venv) and then activating the virtual environment
3. Project dependencies (uv add <package>).


Recommend:

1. [ruff VSCode extension](https://marketplace.visualstudio.com/items?itemName=charliermarsh.ruff). Add the following to your settings.json (Press Ctrl+Shift+P, type settings, click Preference: Open User Settins (JSON)):

```json
"[python]": {
        "editor.defaultFormatter": "charliermarsh.ruff",
        "editor.codeActionsOnSave": {
            "source.fixAll": "explicit",
            "source.organizeImports": "explicit"
        },
        "editor.tabSize": 4
    },
```

2. [Mypy type checker extension](https://marketplace.visualstudio.com/items?itemName=matangover.mypy) for VSCode.

I recommend enabling Format on Save in your settings for all file types.

## Create a repo based on this template

1. If you've got the [GitHub CLI](https://cli.github.com/) installed, you can run `gh repo create` and select 'Create repository from template'.
2. Alternatively, follow instructions [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template).
3. Run `python -m pip install --user pipx` to install pipx
4. After installing pipx ,Add the directory to the path
5. Run `pipx install pre-commit` to install pre-commit globally.
6. Run `pre-commit install` to install pre-commit hooks
7. Run `uv venv` to create a virtual environment
8. Run `.\.venv\Scripts\activate` (on Windows. on MacOS/Linux run `source .venv/bin/activate`) to activate the environment.
9. Run `uv sync` to install the dependencies listed in the `uv.lock` file.
10. To add requirements, run `uv add my-requirement` e.g. for NumPy, it would be `uv add numpy`.
11. To open a python shell, run `uv run python`.
12. To run tests with a coverage report, run `uv run pytest tests --cov=src`. To just run the tests then run `uvx pytest`.
13. Anything else, read `uv --help`.
14. Update [`pyproject.toml`](/pyproject.toml) with the name, version and description of your new project/package.
15. When you've finished making changes, commit them as normal using git. The first commit will take a little while while it initialises. If something doesn't pass a pre-commit test, you will be notified and required to edit the commit until it conforms to the standard. Note that **no commit will occur** until the checks pass.

## Includes

- uv for package and project management
- pre-commit to run pre-commit hooks
- pytest for testing.
- ruff for formatting and linting
- mypy for type checking


