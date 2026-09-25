# The BioKotlin website

[biokotlin.org](https://www.biokotlin.org/) is built from this repository with
[Zensical](https://zensical.org/) and published to GitHub Pages by
[`.github/workflows/docs.yml`](../.github/workflows/docs.yml).

Nothing in this directory is site content. The pages live in [`docs/`](../docs);
this directory holds the tooling that turns them, the notebooks and the Dokka
output into a site.




## Building locally

```sh
python3 -m venv .venv && source .venv/bin/activate
pip install -r website/requirements.txt

python website/convert_notebooks.py   # writes docs/tutorials/*.md
zensical serve                        # default ---> http://127.0.0.1:8000
```

`zensical serve` watches `docs/` and `zensical.toml` and reloads the browser. It
does not watch the notebooks; rerun the converter after editing one.

For a one-off build into `build/site/`:

```sh
zensical build --clean
```

Both the venv and the generated tutorial pages are gitignored.

### Skipping the notebooks

Executing the notebooks needs a JDK and pulls BioKotlin and its dependencies
from Maven Central, which takes a few minutes on a cold cache. When you are
editing prose rather than tutorials, convert from the saved outputs instead:

```sh
python website/convert_notebooks.py --no-execute
```

Most notebooks are checked in without outputs, so those pages will be code-only.
You can also name the tutorials you care about:

```sh
python website/convert_notebooks.py basic-sequence sequence-io
```

### The API reference

The `/api/` section of the site is Dokka output, which is not part of the
Zensical build:

```sh
./gradlew dokkaHtml
mkdir -p build/site/api && cp -r build/dokka/html/. build/site/api/
```

Without that copy, `docs/api-reference.md` links to `../api/` resolve to a 404
in a local build. The workflow does the same two steps before uploading.



## Adding a page

Write the Markdown in `docs/`, then add it to `nav` in `zensical.toml` - Zensical
does not pick up files that are not in the nav. Pages support admonitions,
content tabs, footnotes and the rest of the
[Python Markdown Extensions](https://facelessuser.github.io/pymdown-extensions/)
enabled under `[project.markdown_extensions]`.

### Runnable snippets

Tag a Kotlin fence with `{.runnable}` and wrap the interesting lines in sample
markers:

````markdown
```kotlin {.runnable}
import biokotlin.seq.*

fun main() {
//sampleStart
val dna = NucSeq("GCAGAT")
println(dna.reverse_complement())   // ATCTGC
//sampleEnd
}
```
````



## Adding a tutorial

Tutorials come from real notebooks so that the published pages cannot drift from
the published library. Add a `Tutorial(...)` entry to `NOTEBOOKS` in
`convert_notebooks.py` and a matching `nav` entry in `zensical.toml`; the slug
in both must match.

Notebooks must load the library with `%use biokotlin`, which resolves the
release from Maven Central, rather than a `@file:DependsOn` on a local build
output. Clear the outputs before committing - CI regenerates them - and prefer
`.print()` over bare DataFrame expressions, since the HTML widget output does
not survive conversion to Markdown.


