# Getting started

BioKotlin is published to Maven Central as `org.biokotlin:biokotlin`. Pick
whichever of the following fits how you work.

!!! info "You need Java 21 or newer"

    BioKotlin is compiled to Java 21 bytecode, so a JDK 21+ runtime is
    required. Check yours with `java -version`.

## JupyterLab with the Kotlin kernel

This is the quickest way to try BioKotlin interactively, and it is how the
[tutorials](tutorials/index.md) in this site are written.

1.  Install Java 21 or newer.

2.  Install the Kotlin kernel for Jupyter, either with conda:

    ```shell
    conda install -c jetbrains kotlin-jupyter-kernel
    ```

    or with pip:

    ```shell
    pip install kotlin-jupyter-kernel
    ```

3.  Install [JupyterLab](https://jupyter.org/install), start it, and create a
    notebook with the **Kotlin** kernel.

4.  Load BioKotlin in the first cell:

    ```kotlin
    %use biokotlin
    ```

    The `%use` line magic resolves BioKotlin from Maven Central and sets up its
    default imports. To pin a version, write `%use biokotlin(1.0.0)`.

## Kotlin script

Add a dependency annotation at the top of a `.main.kts` file:

```kotlin
@file:DependsOn("org.biokotlin:biokotlin:1.0.0")

import biokotlin.seq.*

val dna = NucSeq("GCAGAT")
println(dna.reverse_complement())
```

Run it with `kotlin my-script.main.kts`.

## Gradle or Maven project

=== "Gradle (Kotlin DSL)"

    ```kotlin title="build.gradle.kts"
    repositories {
        mavenCentral()
    }

    dependencies {
        implementation("org.biokotlin:biokotlin:1.0.0")
    }

    kotlin {
        jvmToolchain(21)
    }
    ```

=== "Gradle (Groovy DSL)"

    ```groovy title="build.gradle"
    repositories {
        mavenCentral()
    }

    dependencies {
        implementation 'org.biokotlin:biokotlin:1.0.0'
    }
    ```

=== "Maven"

    ```xml title="pom.xml"
    <dependency>
        <groupId>org.biokotlin</groupId>
        <artifactId>biokotlin</artifactId>
        <version>1.0.0</version>
    </dependency>
    ```

All released versions are listed on
[Maven Central](https://central.sonatype.com/artifact/org.biokotlin/biokotlin).

## Command line tools

The command line interface lives in a separate project,
[biokotlin-tools](https://github.com/maize-genetics/biokotlin-tools). Download
an installable tar file from its
[releases page](https://github.com/maize-genetics/biokotlin-tools/releases).

## Building from source

```shell
git clone https://github.com/maize-genetics/BioKotlin.git
cd BioKotlin
./gradlew build
```

`./gradlew shadowJar` produces a fat jar under `build/libs/` if you want to
load a development build into a notebook with `@file:DependsOn`.

## Next steps

- Work through the [tutorials](tutorials/index.md).
- Look up types and functions in the [API reference](api-reference.md).
- Coming from Python? See the [BioPython comparison](biopython.md).
