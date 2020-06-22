import org.gradle.api.JavaVersion.VERSION_11
import org.jetbrains.dokka.gradle.DokkaTask
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

//Note Kotlin version nears to be updated in both the buildscript and plugins.
//Dependencies will follow the buildscript

/*
This build script is need to use the early access
 */
buildscript {
    val kotlinVersion by extra ("1.3.70")
   // val kotlinVersion by extra ("1.4-M1")

    repositories {
        mavenCentral()
        gradlePluginPortal()
        maven ( "https://dl.bintray.com/kotlin/kotlin-eap" )
    }

    dependencies {
        classpath ("org.jetbrains.kotlin:kotlin-gradle-plugin:$kotlinVersion")
    }
}


plugins {
    val kotlinVersion = "1.3.70"
    //val kotlinVersion = "1.4-M1"
    java
    kotlin("jvm") version kotlinVersion
    //Shadow allows for the creation of fat jars (all dependencies)
    id("com.github.johnrengelman.shadow") version "5.2.0"

    id("org.jetbrains.dokka") version "0.10.1"
}

repositories {
    mavenCentral()
    jcenter()
    maven("https://maven.imagej.net/content/groups/public/")
    maven("https://jitpack.io")
    maven ("https://dl.bintray.com/kotlin/kotlin-eap")
    maven ("https://kotlin.bintray.com/kotlinx")
    maven("https://oss.sonatype.org/content/repositories/snapshots/")
}

dependencies {
    val kotlinVersion = rootProject.extra["kotlinVersion"]

    implementation("org.jetbrains.kotlin:kotlin-reflect:${kotlinVersion}")
    implementation("org.jetbrains.kotlin:kotlin-script-runtime:${kotlinVersion}")


    implementation( "org.nield:kotlin-statistics:1.2.1")
    implementation("de.mpicbg.scicomp:krangl:0.11")

    //Biology possible dependencies
    //Support fasta, bam, sam, vcf, bcf support
    implementation("com.github.samtools:htsjdk:2.21.3")

    //    implementation("org.jetbrains.bio:bioinf-commons:0.0.9")

    //TASSEL provides strength in marker analysis and mapping
    //This version of tassel bring along a wide range of kotlin stdlib dependencies.  See last paragraph
    //https://kotlinlang.org/docs/reference/evolution/compatibility-modes.html
//    implementation("net.maizegenetics:tassel:5.2.60")

    //Wide range of sequence tools in Java - API is dated
//    implementation("org.biojava:biojava:5.3.0")
//    implementation("org.biojava:biojava-genome:5.3.0")

    implementation("org.graalvm.sdk:graal-sdk:20.0.0")
    implementation("org.apache.commons:commons-csv:1.8")
    implementation("com.github.jkcclemens:khttp:0.1.0")


    implementation("com.google.guava:guava:29.0-jre")

    testImplementation("org.junit.jupiter:junit-jupiter:5.6.2")

    val kotestVersion = "4.1.0.293-SNAPSHOT"
    listOf("runner-junit5", "assertions-core", "runner-console", "property").forEach {
        testImplementation("io.kotest:kotest-$it-jvm:$kotestVersion")
    }
    //consider adding Kotlintest
    //mockk
}

java {
    sourceCompatibility = VERSION_11
    targetCompatibility = VERSION_11
}

//kotlin {
//    sourceSets{
//        main. += "src/main/kotlin/"
//        testsrcDirs += src/test/kotlin/'
//
//    }
//}

tasks.withType<KotlinCompile>().configureEach {
    kotlinOptions.jvmTarget = "11"
}

tasks {
    println("Ed say ${sourceSets["main"].allSource.srcDirs}")
    val dokka by getting(DokkaTask::class) {
        outputFormat = "html"
        outputDirectory = "$buildDir/dokka"
        configuration {
            includes = listOf("/Users/edbuckler/Code/biokotlin/src/main/kotlin/biokotlin/packages.md")
        }
    }
}
tasks.test {
    useJUnitPlatform()
    testLogging {
        events("passed", "skipped", "failed")
    }
}

//kotlin {
//}
//compileKotlin {
//    kotlinOptions.jvmTarget = "11"
//}
//compileTestKotlin {
//    kotlinOptions.jvmTarget = "11"
//}
//
//sourceSets {
//    main.java.srcDirs += 'src/main/kotlin/'
//    test.java.srcDirs += 'src/test/kotlin/'
//}
//
