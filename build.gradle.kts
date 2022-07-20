import org.gradle.api.JavaVersion.VERSION_11
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

// Note Kotlin version needs to be updated in both the buildscript and plugins.
// Dependencies will follow the buildscript

group = "org.biokotlin"
version = "0.04"

/*
This build script is need to use the early access
 */
buildscript {
    val kotlinVersion by extra ("1.5.31")

    repositories {
        mavenCentral()
        gradlePluginPortal()
        maven("https://dl.bintray.com/kotlin/kotlin-eap")
    }

    dependencies {
        classpath("org.jetbrains.kotlin:kotlin-gradle-plugin:$kotlinVersion")
        classpath(kotlin("serialization", version = kotlinVersion))
        classpath("org.jetbrains.dokka:dokka-gradle-plugin:1.6.21")
    }
}


plugins {
    val kotlinVersion = "1.5.31"
    java
    kotlin("jvm") version kotlinVersion
    kotlin("plugin.serialization") version kotlinVersion
    // Shadow allows for the creation of fat jars (all dependencies)
    id("com.github.johnrengelman.shadow") version "5.2.0"

    //This is used for code generation for DataFrame Schema, however, it does not work
    //https://github.com/Kotlin/dataframe/tree/eb9ec4fb90f906f6a98e69b9c5a0369009d34bbb/plugins/gradle/codegen
    //id("org.jetbrains.kotlinx.dataframe") version "1.0-SNAPSHOT"

    id("org.jetbrains.dokka") version "1.6.21"
    `java-library`
    `maven-publish`
    signing
}
apply {
    plugin("kotlinx-serialization")
    plugin("org.jetbrains.dokka")
}


repositories {
    mavenCentral()
    maven("https://maven.imagej.net/content/groups/public/")
    maven("https://jitpack.io")
    maven("https://dl.bintray.com/kotlin/kotlin-eap")
    maven("https://kotlin.bintray.com/kotlinx")
    maven("https://oss.sonatype.org/content/repositories/snapshots/")
}

dependencies {
    val kotlinVersion = rootProject.extra["kotlinVersion"]

    implementation("org.jetbrains.kotlin:kotlin-reflect:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-script-runtime:${kotlinVersion}")
    implementation("org.jetbrains.kotlinx:kotlinx-coroutines-core:1.5.2")
    implementation("org.jetbrains.kotlinx:kotlinx-serialization-json:1.2.2")

    implementation("org.nield:kotlin-statistics:1.2.1")
    implementation("de.mpicbg.scicomp:krangl:0.13")
    implementation("org.jetbrains.kotlinx:dataframe:0.8.0-rc-7")

    // Biology possible dependencies
    // Support fasta, bam, sam, vcf, bcf support
    implementation("com.github.samtools:htsjdk:2.24.1")

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
    implementation("io.ktor:ktor-client-core:1.6.3")
    implementation("io.ktor:ktor-client-cio:1.6.3")

    implementation("com.google.guava:guava:30.1.1-jre")
    implementation("org.apache.tinkerpop:gremlin-core:3.5.1")
    implementation("org.jgrapht:jgrapht-core:1.5.1")

    testImplementation("org.junit.jupiter:junit-jupiter:5.8.0")

    val kotestVersion = "4.2.6"
    listOf("runner-junit5", "assertions-core", "property").forEach {
        testImplementation("io.kotest:kotest-$it-jvm:$kotestVersion")
    }
    //consider adding Kotlintest
}
//This is used for code generation for DataFrame Schema, however, it does not work
//https://github.com/Kotlin/dataframe/tree/eb9ec4fb90f906f6a98e69b9c5a0369009d34bbb/plugins/gradle/codegen
//kotlin.sourceSets.getByName("main").kotlin.srcDir("build/generated/ksp/main/kotlin/")

java {
    sourceCompatibility = VERSION_11
    targetCompatibility = VERSION_11
    withSourcesJar()
}

tasks.withType<KotlinCompile>().configureEach {
    kotlinOptions.jvmTarget = "11"
}

tasks {
    println("Source directories: ${sourceSets["main"].allSource.srcDirs}")
}

val dokkaHtml by tasks.getting(org.jetbrains.dokka.gradle.DokkaTask::class) {
    dokkaSourceSets {
        configureEach {
            includes.from("src/main/kotlin/biokotlin/packages.md")
        }
    }
}

val tutorialInjector by tasks.registering {
    doLast {
        //Convert raw notebooks into HTML files
        val raw = File("${rootProject.projectDir}/tutorials/raw")
        val html = File("${rootProject.projectDir}/build/dokka/html")
        for (notebook in raw.listFiles()!!) {
            if (notebook.name.endsWith(".ipynb")) {
                //TODO look at how to manage this executable properly
                File("${rootProject.projectDir}/build/dokka/html").mkdirs()
                ProcessBuilder(
                    "/home/jeff/miniconda3/bin/jupyter",
                    "nbconvert",
                    notebook.absolutePath,
                    "--stdout",
                    "--to",
                    "html",
                )
                    .redirectOutput(File("build/dokka/html/${notebook.name.substringBefore(".")}.html"))
                    .start()
                    .waitFor()
            }
        }

        val packageList = File("build/dokka/html/biokotlin/package-list").readText()

        //Replaces shorthand links with accurate ones and inserts links into code
        val linkFinder = Regex("<a href=\"\\..*?>")
        for (file in html.listFiles()!!) {
            if (!file.isFile) continue
            //Matches <a href=". until >
            val text = file.readText()
            val fixedLinks = linkFinder.replace(text) { a ->
                a.value.replace(Regex("\".*\"")) { "\"${parseLink(it.value.drop(1).dropLast(1), packageList)}\" target=\"_blank\"" }
            }

            file.writeText(fixedLinks)
        }

        //Replaces any <p></p> blocks in the dokka files that contain !tutorial with their tutorial
        recursivelyInjectTutorials(html, -1)
    }
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
    println("==============")
    matches.forEach { println(it.value) }

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

fun recursivelyInjectTutorials(file: File, depth: Int) {
    if (file.isDirectory) file.listFiles()!!.forEach { recursivelyInjectTutorials(it, depth + 1) }
    val tutorialIdentifier = Regex("<p [^<]*?!tutorial .*?<\\/p>")

    if (file.name.endsWith(".html")) {
        val text = file.readText()
        val injected = tutorialIdentifier.replace(text) {
            val name = it.value.substringAfter(">").substringAfter(" ").substringBefore("<")
            "<div class=\"iframeDiv\" style=\"width:100%; padding-bottom:56.25%; position:relative;\">\n" +
                    "  <iframe src=\"${"../".repeat(depth)}$name.html\" style=\"position:absolute; top:0px; left:0px; \n" +
                    "  width:100%; height:100%; border: none; overflow: hidden;\"></iframe>\n" +
                    "</div>"
        }
        file.writeText(injected)
    }
}

val dokkaJar by tasks.creating(Jar::class) {
    dependsOn(dokkaHtml)
    doLast {
        group = JavaBasePlugin.DOCUMENTATION_GROUP
        description = "BioKotlin: ${property("version")}"
        archiveClassifier.set("javadoc")
        from(dokkaHtml.outputDirectory)
    }
    finalizedBy(tutorialInjector)
}

tasks.test {
    useJUnitPlatform()
    testLogging {
        events("passed", "skipped", "failed")
    }
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

            repositories {
                maven {
                    val releasesRepoUrl = "https://oss.sonatype.org/service/local/staging/deploy/maven2"
                    val snapshotsRepoUrl = "https://oss.sonatype.org/content/repositories/snapshots"
                    url = uri(if (version.toString().endsWith("SNAPSHOT")) snapshotsRepoUrl else releasesRepoUrl)
                    credentials {
                        try {
                            username = property("ossrhUsername") as String?
                            password = property("ossrhPassword") as String?
                        } catch (e: Exception) {
                            println("Unable to get username and password for nexus maven central release.")
                        }
                    }
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
                }
                scm {
                    connection.set("scm:git:git://bitbucket.org:bucklerlab/biokotlin.git")
                    developerConnection.set("scm:git:ssh://bitbucket.org:bucklerlab/biokotlin.git")
                    url.set("https://bitbucket.org/bucklerlab/biokotlin/src")
                }
            }
        }
    }
}

signing {
    useGpgCmd()
    sign(publishing.publications["maven"])
}

tasks.javadoc {
    dependsOn("dokkaJavadoc")
    if (JavaVersion.current().isJava9Compatible) {
        (options as StandardJavadocDocletOptions).addBooleanOption("html5", true)
    }
}


tasks.jar {
    configurations["compileClasspath"].forEach { file: File ->
        from(zipTree(file.absoluteFile))
    }
    duplicatesStrategy = DuplicatesStrategy.INCLUDE
}