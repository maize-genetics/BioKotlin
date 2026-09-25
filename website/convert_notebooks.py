#!/usr/bin/env python3


"""Turn the tutorial notebooks into Markdown pages for the website.

Each notebook in NOTEBOOKS is executed with the Kotlin Jupyter kernel and
written to docs/tutorials/. The notebooks load BioKotlin with `%use biokotlin`,
which resolves the released artifact from Maven Central, so the pages show what
the published library actually does rather than a hand-written transcript. This
idea is heavily inspired by how `pkgdown` works for R vignettes to web page
production

Usage:
    python website/convert_notebooks.py               # execute, then convert
    python website/convert_notebooks.py --no-execute  # use saved outputs
    python website/convert_notebooks.py basic-sequence

The generated pages are not checked in; docs/tutorials/index.md is.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbconvert import MarkdownExporter


ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "docs" / "tutorials"
ASSET_DIR = OUT_DIR / "assets"
REPO_URL = "https://github.com/maize-genetics/BioKotlin/blob/master"



KERNEL = "kotlin"
CELL_TIMEOUT = 900 # add timeout for potentially frozen states


@dataclass(frozen=True)
class Tutorial:
    """One notebook and the page it becomes."""

    slug: str
    title: str
    notebook: str
    summary: str


# Order matches the Tutorials section of zensical.toml.
NOTEBOOKS: tuple[Tutorial, ...] = (
    Tutorial(
        "basic-sequence",
        "Basic sequence",
        "documentation_resources/raw_tutorials/basic_sequence.ipynb",
        "Creating nucleotide and protein sequences and the operations that come with them.",
    ),
    Tutorial(
        "sequence-operations",
        "Sequence operations",
        "notebooks/TutorialSeq.ipynb",
        "Complementing, transcribing, translating, and searching, including at scale.",
    ),
    Tutorial(
        "nucleotides-and-residues",
        "Nucleotides and residues",
        "notebooks/NucleotidesDNARNA.ipynb",
        "The NUC enum: IUPAC codes, molecular weights, and complements.",
    ),
    Tutorial(
        "amino-acids-and-proteins",
        "Amino acids and proteins",
        "notebooks/AminoAcidProtein.ipynb",
        "Building peptides and reading residue properties from the AminoAcid enum.",
    ),
    Tutorial(
        "sequence-io",
        "Sequence IO",
        "notebooks/SeqIO.ipynb",
        "Streaming FASTA and FASTQ records with NucSeqIO.",
    ),
    Tutorial(
        "genomic-ranges",
        "Genomic ranges",
        "notebooks/RangesTutorial.ipynb",
        "Positions, intervals, flanking, intersection, and reading BED files.",
    ),
    Tutorial(
        "feature-tree",
        "Feature tree (GFF)",
        "notebooks/FeatureTree.ipynb",
        "Parsing GFF3 into an immutable gene/transcript/exon tree, and mutating it.",
    ),
    Tutorial(
        "genomic-features",
        "Genomic features",
        "notebooks/Biokotlin_GenomicFeatures_Tutorial.ipynb",
        "A tabular, DataFrame-backed view of a GFF annotation.",
    ),
    Tutorial(
        "maf-processing",
        "MAF processing",
        "notebooks/MAFProcessingTutorial.ipynb",
        "Coverage and identity from MAF alignments, exported as BED and wiggle.",
    ),
)

# Dokka shorthand for API links, as in `[complement](.seq.NucSeq.complement)`.
# Only the Dokka build resolves these; on the website they would be dead links,
# so the link text becomes inline code instead.
DOKKA_LINK = re.compile(r"\[([^\]]+)\]\(\.[A-Za-z][\w.]*\)")

MARKDOWN_IMAGE = re.compile(r"!\[([^\]]*)\]\((?!https?:)([^)]+)\)")





def quiet_kernel_logging() -> None:
    """Turn down the kernel's logging before it is launched.

    See website/kernel-logback.xml. Respect the variable if the caller already
    set it, so this stays overridable.
    """
    variable = "KOTLIN_JUPYTER_JAVA_OPTS_EXTRA"
    if variable in os.environ:
        return
    config = Path(__file__).resolve().parent / "kernel-logback.xml"
    os.environ[variable] = f"-Dlogback.configurationFile={config}"


def strip_widget_outputs(nb) -> None:
    """Drop notebook outputs that a static page cannot render.

    `%use biokotlin` loads the Kotlin DataFrame integration, which renders
    tables as an HTML fragment that fetches a script from a CDN and fills in an
    empty <table> with JavaScript. Outside a live kernel that leaves a blank
    table, and the loader emits a full <html> document that would corrupt the
    page it is spliced into. The notebooks call .print() so every table also has
    a plain text form; here we keep that and discard the HTML.

    stderr is dropped for the same reason it is hidden in a rendered notebook:
    it carries the libraries' progress logging, not results. Genuine failures
    arrive as `error` outputs and are left in place so they are impossible to
    miss.
    """
    for cell in nb.cells:
        if cell.get("cell_type") != "code":
            continue
        kept = []
        for output in cell.get("outputs", []):
            if output.get("output_type") in ("display_data", "execute_result"):
                data = output.get("data", {})
                data.pop("text/html", None)
                if not data:
                    continue
            elif output.get("name") == "stderr":
                continue
            kept.append(output)
        cell["outputs"] = kept


def plain_text(output) -> str | None:
    """The text of an output that carries nothing but text, else None."""
    if output.get("output_type") == "stream":
        text = output.get("text", "")
    elif output.get("output_type") in ("display_data", "execute_result"):
        data = output.get("data", {})
        if set(data) - {"text/plain"}:
            return None
        text = data.get("text/plain", "")
    else:
        return None
    return "".join(text) if isinstance(text, list) else text


def merge_text_outputs(nb) -> None:
    """Fold a cell's consecutive text outputs into one.

    nbconvert indents each output as its own block, and adjacent indented
    blocks separated by blank lines are one Markdown code block with the blank
    lines inside it. A cell that prints and then returns a value therefore
    renders as a tall box with a gap in the middle, where the notebook shows an
    uninterrupted transcript. Trailing blank lines go at the same time, since
    they would otherwise pad the bottom of the box.
    """
    for cell in nb.cells:
        if cell.get("cell_type") != "code":
            continue

        merged: list = []
        run: list[str] = []

        def flush() -> None:
            text = "".join(run).rstrip()
            run.clear()
            if text:
                merged.append(
                    nbformat.v4.new_output("stream", name="stdout", text=text + "\n")
                )

        for output in cell.get("outputs", []):
            text = plain_text(output)
            if text is None:
                flush()
                merged.append(output)
                continue
            if run and not run[-1].endswith("\n"):
                run.append("\n")
            run.append(text)

        flush()
        cell["outputs"] = merged


def copy_images(nb, notebook_path: Path) -> None:
    """Copy images a notebook references off disk into the site assets."""
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    for cell in nb.cells:
        if cell.get("cell_type") != "markdown":
            continue

        def replace(match: re.Match) -> str:
            alt, src = match.group(1), match.group(2)
            source = (notebook_path.parent / src).resolve()
            if not source.is_file():
                print(f"  warning: missing image {src}", file=sys.stderr)
                return match.group(0)
            shutil.copyfile(source, ASSET_DIR / source.name)
            return f"![{alt}](assets/{source.name})"

        cell.source = MARKDOWN_IMAGE.sub(replace, cell.source)


def rewrite_links(nb) -> None:
    for cell in nb.cells:
        if cell.get("cell_type") == "markdown":
            cell.source = DOKKA_LINK.sub(r"`\1`", cell.source)


def convert(tutorial: Tutorial, execute: bool) -> None:
    notebook_path = ROOT / tutorial.notebook
    print(f"{tutorial.notebook} -> docs/tutorials/{tutorial.slug}.md")

    nb = nbformat.read(notebook_path, as_version=4)
    if execute:
        # Relative paths in the notebooks resolve against their own directory.
        NotebookClient(
            nb,
            timeout=CELL_TIMEOUT,
            kernel_name=KERNEL,
            resources={"metadata": {"path": str(notebook_path.parent)}},
        ).execute()

    strip_widget_outputs(nb)
    merge_text_outputs(nb)
    copy_images(nb, notebook_path)
    rewrite_links(nb)

    body, _ = MarkdownExporter().from_notebook_node(nb)

    # Quoted because summaries contain colons, which YAML would otherwise read
    # as a nested mapping.
    header = (
        "---\n"
        f'title: "{tutorial.title}"\n'
        f'description: "{tutorial.summary}"\n'
        "---\n\n"
        f"<!-- Generated from {tutorial.notebook} by website/convert_notebooks.py.\n"
        "     Edit the notebook, not this file. -->\n\n"
        f"# {tutorial.title}\n\n"
        f"[Open this tutorial as a notebook]({REPO_URL}/{tutorial.notebook})\n\n"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / f"{tutorial.slug}.md").write_text(header + body.strip() + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "slugs",
        nargs="*",
        help="only convert these tutorials (default: all)",
    )
    parser.add_argument(
        "--no-execute",
        action="store_true",
        help="convert the notebooks as saved, without running them",
    )
    args = parser.parse_args()

    selected = NOTEBOOKS
    if args.slugs:
        selected = tuple(t for t in NOTEBOOKS if t.slug in args.slugs)
        unknown = set(args.slugs) - {t.slug for t in selected}
        if unknown:
            parser.error(f"unknown tutorial(s): {', '.join(sorted(unknown))}")

    if not args.no_execute:
        quiet_kernel_logging()

    for tutorial in selected:
        convert(tutorial, execute=not args.no_execute)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
