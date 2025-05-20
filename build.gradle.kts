import org.jetbrains.kotlin.gradle.tasks.KotlinCompile
import java.io.ByteArrayOutputStream
import java.nio.file.Files
import java.nio.file.Paths
import kotlin.io.path.isRegularFile
import kotlin.jvm.optionals.getOrNull

// This is used to get the version from the git tag
// The version is expected to be in the format X.Y.Z
// JReleaser will use this version to create the release
fun getVersionName(): String {
    val stdout = ByteArrayOutputStream()
    exec {
        commandLine = listOf("git", "describe", "--tags", "--abbrev=0")
        standardOutput = stdout
    }
    val versionStr = stdout.toString().trim()
    val parts = versionStr.removePrefix("v").removePrefix("V").split('.')
    val normalizedStr = when (parts.size) {
        0 -> throw IllegalArgumentException("Version string is empty")
        1 -> "${parts[0]}.0.0"
        2 -> "${parts[0]}.${parts[1]}.0"
        else -> versionStr
    }
    return normalizedStr
}

group = "org.biokotlin"
version = getVersionName()

/*
This build script is need to use the early access
 */
buildscript {
    val kotlinVersion by extra("1.9.24")

    repositories {
        mavenCentral()
        gradlePluginPortal()
    }

    dependencies {
        classpath("org.jetbrains.kotlin:kotlin-gradle-plugin:$kotlinVersion")
        classpath(kotlin("serialization", version = kotlinVersion))
        classpath("org.jetbrains.dokka:dokka-gradle-plugin:1.9.20")
    }
}

plugins {
    val kotlinVersion = "1.9.24"
    java
    kotlin("jvm") version kotlinVersion
    kotlin("plugin.serialization") version kotlinVersion
    // Shadow allows for the creation of fat jars (all dependencies)
    id("com.github.johnrengelman.shadow") version "8.1.1"

    //This is used for code generation for DataFrame Schema, however, it does not work
    //https://github.com/Kotlin/dataframe/tree/eb9ec4fb90f906f6a98e69b9c5a0369009d34bbb/plugins/gradle/codegen
    //id("org.jetbrains.kotlinx.dataframe") version "1.0-SNAPSHOT"

    application
    id("org.jetbrains.dokka") version "1.9.20"
    `java-library`
    `maven-publish`
    signing
    id("org.jreleaser") version "1.18.0"
}

apply {
    plugin("kotlinx-serialization")
    plugin("org.jetbrains.dokka")
    plugin("org.jreleaser")
}

repositories {
    mavenCentral()
    gradlePluginPortal()
    maven("https://maven.imagej.net/content/groups/public/")
    maven("https://jitpack.io")
    maven("https://dl.bintray.com/kotlin/kotlin-eap")
    maven("https://kotlin.bintray.com/kotlinx")
}

dependencies {
    val kotlinVersion = rootProject.extra["kotlinVersion"]

    implementation("org.apache.logging.log4j:log4j-core:2.23.1")
    implementation("org.apache.logging.log4j:log4j-api:2.23.1")

    implementation("org.jetbrains.kotlin:kotlin-reflect:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-script-runtime:${kotlinVersion}")
    implementation("org.jetbrains.kotlinx:kotlinx-coroutines-core:1.8.1")
    implementation("org.jetbrains.kotlinx:kotlinx-serialization-json:1.2.2")

    implementation("org.nield:kotlin-statistics:1.2.1")
    implementation("org.jetbrains.kotlinx:dataframe:0.8.0-rc-7")

    // Biology possible dependencies
    // Support fasta, bam, sam, vcf, bcf support
    implementation("com.github.samtools:htsjdk:4.0.1")

    // implementation("org.jetbrains.bio:bioinf-commons:0.0.9")

    // TASSEL provides strength in marker analysis and mapping
    // This version of tassel bring along a wide range of kotlin stdlib dependencies.  See last paragraph
    // https://kotlinlang.org/docs/reference/evolution/compatibility-modes.html
    // implementation("net.maizegenetics:tassel:5.2.60")

    // Wide range of sequence tools in Java - API is dated
    // implementation("org.biojava:biojava:5.3.0")
    // implementation("org.biojava:biojava-genome:5.3.0")

    implementation("org.graalvm.sdk:graal-sdk:21.2.0")
    implementation("org.apache.commons:commons-csv:1.8")

    implementation("com.google.guava:guava:33.1.0-jre")
    implementation("org.apache.tinkerpop:gremlin-core:3.5.1")


    implementation(group = "ch.qos.logback", name = "logback-classic", version = "1.2.6")
    implementation("it.unimi.dsi:fastutil:8.5.12")
    implementation("org.lz4:lz4-java:1.8.0")

    testImplementation("org.junit.jupiter:junit-jupiter:5.8.0")

    val kotestVersion = "5.6.2"
    listOf("runner-junit5", "assertions-core", "property", "framework-datatest").forEach {
        testImplementation("io.kotest:kotest-$it-jvm:$kotestVersion")
    }

}

// This is used for code generation for DataFrame Schema, however, it does not work
// https://github.com/Kotlin/dataframe/tree/eb9ec4fb90f906f6a98e69b9c5a0369009d34bbb/plugins/gradle/codegen
// kotlin.sourceSets.getByName("main").kotlin.srcDir("build/generated/ksp/main/kotlin/")

java {
    withJavadocJar()
    withSourcesJar()
}

tasks.withType<KotlinCompile>().configureEach {
    kotlinOptions.jvmTarget = "17"
}

tasks {
    println("Source directories: ${sourceSets["main"].allSource.srcDirs}")
}

application {
    mainClass.set("biokotlin.cli.BiokotlinKt")

    // Set name of generated scripts in bin/
    applicationName = "biokotlin"
}

/**
 * Generates HTML files based on Javadoc-style comments. Supports automatic insertion of Jupyter notebook tutorials,
 * (see [tutorialInjector] for details). Supports insertion of images (see [imageInjector] for details).
 */
val dokkaHtml by tasks.getting(org.jetbrains.dokka.gradle.DokkaTask::class) {
    dokkaSourceSets {
        configureEach {
            includes.from("src/main/kotlin/biokotlin/packages.md")
        }
    }
    doLast {
        tutorialInjector()
        imageInjector()
    }
}

/**
 * Replaces occurences of ```!tutorial tutorial_name``` within the Javadocs with am embedded tutorial.
 *
 * ```!tutorial tutorial_name``` must be on its own line, surrounded by a blank line on either side.
 *
 * Do not include the file extension in tutorial_name. The tutorial must be present in documentation_resources/raw_tutorials
 * as a Jupyter notebook file.
 */
fun tutorialInjector() {
    //Convert raw notebooks into HTML files
    val raw = File("${rootProject.projectDir}/documentation_resources/raw_tutorials")
    val html = File("${rootProject.projectDir}/build/dokka/html")
    val tutorials = File("${rootProject.projectDir}/build/dokka/html/tutorials")
    for (notebook in raw.listFiles()!!) {
        if (notebook.name.endsWith(".ipynb")) {
            tutorials.mkdirs()
            try {
                ProcessBuilder(
                    "jupyter",
                    "nbconvert",
                    notebook.absolutePath,
                    "--stdout",
                    "--to",
                    "html",
                )
                    .redirectOutput(File("build/dokka/html/tutorials/${notebook.name.substringBefore(".")}.html"))
                    .start()
                    .waitFor()
            } catch (e: java.io.IOException) {
                println("Warning: Jupyter must be installed through conda for tutorials to be injected.")
                e.printStackTrace()
                return
            }

        }
    }

    val optionalPackageList = Files.walk(Paths.get("./build/dokka/"))
        .filter { it.isRegularFile() }
        .filter { it.toString().endsWith("package-list") }
        .findFirst()

    val packageList = optionalPackageList.getOrNull()?.let {
        it.toFile().readText()
    } ?: ""

    //Replaces shorthand links with accurate ones and inserts links into code
    val linkFinder = Regex("<a href=\"\\..*?>")
    for (file in tutorials.listFiles()!!) {
        if (!file.isFile) continue
        //Matches <a href=". until >
        val text = file.readText()
        val fixedLinks = linkFinder.replace(text) { a ->
            a.value.replace(Regex("\".*\"")) {
                "\"../${
                    parseLink(
                        it.value.drop(1).dropLast(1),
                        packageList
                    )
                }\" target=\"_blank\""
            }
        }

        file.writeText(fixedLinks)
    }

    //Replaces any <p></p> blocks in the dokka files that contain !tutorial with their tutorial
    recursivelyInjectTutorials(html, -1)
}

/**
 * Parses a link from a package-list file generated by Dokka. Returns a relative link.
 *
 * The path should be formatted as .package_name.class_name.member_name
 *
 * If a class/interface and a function have the same name and are in the same package,
 * the class will be returned by default. To return the function, end the path with a single
 * exclamation mark.
 *
 * Example call for the Seq interface:
 * ```kotlin
 * parseLink(".seq.Seq")
 * ```
 *
 * Example call for the Seq function
 * ```kotlin
 * parseLink(".seq.Seq!")
 * ```
 *
 * @param path the pathname of the kotlin member
 * @param text the text of a package-list file
 */
fun parseLink(path: String, text: String): String {
    //\$dokka\.location:biokotlin[\./]*seq[\./]*Seq[\./]*(#|(PointingToDeclaration)).*?$

    var isFunction = false

    // Parsing to see if it ends with exclamation and is therefore a function
    val correctedPath = if (path.endsWith("!")) {
        isFunction = true
        path.substring(0, path.length - 1)
    } else {
        path
    }

    val delimiters = "[./]*"
    val regexString = "\\\$dokka\\.location:biokotlin" +
            correctedPath.replace(".", delimiters) + "/*(#|(PointingToDeclaration)).*"

    val regex = Regex(regexString)
    val matches = regex.findAll(text)

    val functions = matches.asSequence().filter { it.value.contains("#") }
    val nonFunctions = matches.asSequence().filter { !it.value.contains("#") }

    return if (isFunction && functions.count() >= 1) {
        linkFromMatch(functions.first())
    } else if (!isFunction && nonFunctions.count() >= 1) {
        linkFromMatch(nonFunctions.first())
    } else if (matches.count() >= 1) {
        linkFromMatch(matches.first())
    } else {
        println("Warning: Could not find link for $path. Empty link placed instead")
        ""
    }
}

/**
 * Parses link from match result
 */
fun linkFromMatch(match: MatchResult): String {
    return match.value.substringAfter('\u001F')
}

/**
 * Recursively injects the tutorials into iframes in the dokka files.
 */
fun recursivelyInjectTutorials(file: File, depth: Int) {
    if (file.isDirectory) file.listFiles()!!.forEach { recursivelyInjectTutorials(it, depth + 1) }
    val tutorialIdentifier = Regex("<p [^<]*?!tutorial .*?<\\/p>")

    if (file.name.endsWith(".html")) {
        val text = file.readText()
        val injected = tutorialIdentifier.replace(text) {
            val name = it.value.substringAfter(">").substringAfter(" ").substringBefore("<")
            "<div class=\"iframeDiv\" style=\"width:100%; padding-bottom:56.25%; position:relative;\">\n" +
                    "  <iframe src=\"${"../".repeat(depth)}tutorials/$name.html\" style=\"position:absolute; top:0px; left:0px; \n" +
                    "  width:100%; height:100%; border: none; overflow: hidden;\"></iframe>\n" +
                    "</div>"
        }
        file.writeText(injected)
    }
}

/**
 * Copies all files from the documentation_resources/images directory to the build/dokka/html directory.
 * Supports nested folders within the images directory. Treat the image directory as the root in your
 * src tag. Example, to access an image in the images directory, simply use the filename "myImage.svg", but to
 * access one in a subdirectory, use "mySubdirectory/myImage.svg".
 *
 * Then, recursively descends the html folder, finds all files with a .html extensions and replaces
 * their <img> tags with the proper src attribute, taking into account the depth of the
 * folder. Ignores img src that are on the internet. If an image lacks an alt tag, it will be
 * inserted based on the content of a text file with the same name as the image (in the same
 * directory).
 *
 *
 * Prefer using svg files for images.
 *
 * To allow for flexibility, all styling for the image must be included in the style tag of the
 * image; the injector will not attempt to inject styling.
 */
fun imageInjector() {
    val images = File("${rootProject.projectDir}/documentation_resources/images")
    val html = File("${rootProject.projectDir}/build/dokka/html")
    val documentationImages = File("${rootProject.projectDir}/build/dokka/html/documentation_images")
    documentationImages.mkdirs()
    images.copyRecursively(documentationImages, overwrite = true)
    recursivelyInjectImages(html, -1)
}

fun recursivelyInjectImages(file: File, depth: Int) {
    if (file.isDirectory) file.listFiles()!!.forEach { recursivelyInjectImages(it, depth + 1) }
    if (file.name.endsWith(".html")) {
        val text = file.readText()
        val injected = Regex("<img [^<]*?>").replace(text) {
            val src = it.value.substringAfter("src=\"").substringBefore("\"")
            if (src.startsWith("http")) {
                it.value
            } else {
                val correctedSource = it.value.replace(src, "${"../".repeat(depth)}documentation_images/$src")
                if (it.value.contains("alt")) {
                    correctedSource
                } else {
                    val alt =
                        File("${rootProject.projectDir}/documentation_resources/images/${src.substringBefore(".")}.txt").readText()
                    correctedSource.replace(">", "alt=\"$alt\">")
                }
            }
        }
        file.writeText(injected)
    }
}

val dokkaJar by tasks.creating(Jar::class) {
    dependsOn(dokkaHtml)
    mustRunAfter(dokkaHtml)
    group = JavaBasePlugin.DOCUMENTATION_GROUP
    description = "BioKotlin: ${project.version}"
    archiveClassifier.set("javadoc")
    from(dokkaHtml.outputDirectory)
}

tasks.test {
    useJUnitPlatform()
    testLogging {
        events("passed", "skipped", "failed")
    }
}

tasks.jar {
    from(sourceSets.main.get().output)
}


publishing {
    publications {

        create<MavenPublication>("maven") {
            artifactId = "biokotlin"
            description = "org.biokotlin:biokotlin:$version"

            from(components["java"])
            artifact(dokkaJar)

            versionMapping {
                usage("java-api") {
                    fromResolutionOf("runtimeClasspath")
                }
                usage("java-runtime") {
                    fromResolutionResult()
                }
            }

            pom {
                name.set("BioKotlin")
                description.set("BioKotlin aims to be a high-performance bioinformatics library that brings the power and speed of compiled programming languages to scripting and big data environments.")
                url.set("http://www.biokotlin.org/")
                licenses {
                    license {
                        name.set("The Apache License, Version 2.0")
                        url.set("http://www.apache.org/licenses/LICENSE-2.0.txt")
                    }
                }
                developers {
                    developer {
                        name.set("Ed Buckler")
                        email.set("esb33@cornell.edu")
                    }
                    developer {
                        name.set("Vaishnavi Gupta")
                        email.set("vg222@cornell.edu")
                    }
                    developer {
                        id.set("tmc46")
                        name.set("Terry Casstevens")
                        email.set("tmc46@cornell.edu")
                    }
                    developer {
                        name.set("Zack Miller")
                        email.set("zrm22@cornell.edu")
                    }
                    developer {
                        name.set("Lynn Johnson")
                        email.set("lcj34@cornell.edu")
                    }
                    developer {
                        name.set("Brandon Monier")
                        email.set("bm646@cornell.edu")
                    }
                    developer {
                        name.set("Peter Bradbury")
                        email.set("pjb39@cornell.edu")
                    }
                    developer {
                        name.set("Jeffrey Morse")
                        email.set("jbm249@cornell.edu")
                    }
                    developer {
                        name.set("Ana Berthel")
                        email.set("ahb232@cornell.edu")
                    }
                    developer {
                        name.set("Mathavan Alagu Ganesan")
                        email.set("gam283@cornell.edu")
                    }
                }
                scm {
                    connection.set("scm:git:git://github.com/maize-genetics/BioKotlin.git")
                    developerConnection.set("scm:git:ssh://github.com/maize-genetics/BioKotlin.git")
                    url.set("https://github.com/maize-genetics/BioKotlin")
                }
            }
        }
    }
}

signing {
    useInMemoryPgpKeys(System.getenv("GPG_SIGNING_KEY"), System.getenv("GPG_SIGNING_PASSWORD"))
    sign(publishing.publications["maven"])
}

tasks.javadoc {
    dependsOn("dokkaJavadoc")
    if (JavaVersion.current().isJava9Compatible) {
        (options as StandardJavadocDocletOptions).addBooleanOption("html5", true)
    }
}

jreleaser {
    signing {
        setActive("ALWAYS")
        armored.set(true)
    }
    deploy {
        maven {
            // Portal Publisher API via Central Publishing Portal
            mavenCentral {
                setActive("ALWAYS")
                uri("https://central.sonatype.com/api/v1/publisher")
            }
        }
    }
}

tasks.publish {
    dependsOn(dokkaJar)
    mustRunAfter(dokkaJar)
}
