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
    implementation("com.github.holgerbrandl:krangl:0.18")
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

    //Multidimensional matrices
    implementation("org.jetbrains.kotlinx:multik-api:0.1.1")
    implementation("org.jetbrains.kotlinx:multik-jvm:0.1.1")

    testImplementation("org.junit.jupiter:junit-jupiter:5.8.0")

    val kotestVersion = "4.2.6"
    listOf("runner-junit5", "assertions-core", "property").forEach {
        testImplementation("io.kotest:kotest-$it-jvm:$kotestVersion")
    }
    //consider adding Kotlintest
    //mockk
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

val dokkaJar by tasks.creating(Jar::class) {
    dependsOn(dokkaHtml)
    group = JavaBasePlugin.DOCUMENTATION_GROUP
    description = "BioKotlin: ${property("version")}"
    archiveClassifier.set("javadoc")
    from(dokkaHtml.outputDirectory)
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
