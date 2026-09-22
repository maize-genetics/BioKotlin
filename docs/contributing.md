# Contributing

Please join our effort. We welcome new contributors to assist with the
development of BioKotlin.

## Where to start

- [Discussions](https://github.com/maize-genetics/BioKotlin/discussions)
  - ask questions and propose ideas.
- [Issues](https://github.com/maize-genetics/BioKotlin/issues) - report
  bugs and pick up work.
- [Trello board](https://trello.com/b/QVu8r335/biokotlin) - what the team
  is working on.
- New to pull requests? Atlassian has a good
  [walkthrough](https://www.atlassian.com/git/tutorials/making-a-pull-request).

Code examples are especially welcome. If you have a notebook that shows off
something BioKotlin does well, open a discussion or a pull request - the
[tutorials](tutorials/index.md) on this site are generated directly from
notebooks in the repository.

## Working on the library

```shell
git clone https://github.com/maize-genetics/BioKotlin.git
cd BioKotlin
./gradlew build
```

Tests run with `./gradlew test` and use [Kotest](https://kotest.io) and JUnit
5. Opening a pull request against `master` runs the suite in CI.

## Working on this website

The site lives in the same repository and is built with
[Zensical](https://zensical.org). Pages are Markdown under `docs/`, and the
tutorials are generated from the notebooks under `notebooks/`.

```shell
python -m venv .venv && source .venv/bin/activate
pip install -r website/requirements.txt

python website/convert_notebooks.py --no-execute   # generate tutorial pages
zensical serve                                     # preview at localhost:8000
```

`website/README.md` has the full details, including how to re-execute the
notebooks and how the runnable code snippets are wired up.
