import org.gradle.api.JavaVersion.VERSION_11
import org.jetbrains.dokka.gradle.DokkaTask

plugins {
    java
    kotlin("jvm") version "1.3.72"
    //Shadow allows for the creation of fat jars (all dependencies)
    id("com.github.johnrengelman.shadow") version "5.2.0"

    id("org.jetbrains.dokka") version "0.10.1"
}

repositories {
    mavenCentral()
    jcenter()
    maven("http://maven.imagej.net/content/groups/public/")
    maven("https://jitpack.io")
}

dependencies {
    implementation(kotlin("stdlib"))
    testImplementation("org.junit.jupiter:junit-jupiter:5.6.2")
//    implementation("org.jetbrains.bio:bioinf-commons:0.0.9")
    implementation("com.github.samtools:htsjdk:2.21.3")
    implementation("de.mpicbg.scicomp:krangl:0.11")
    implementation("net.maizegenetics:tassel:5.2.60")
    implementation("org.jetbrains.kotlin:kotlin-script-runtime:1.3.71")
//    implementation("org.biojava:biojava:5.3.0")
//    implementation("org.biojava:biojava-genome:5.3.0")
    implementation("org.graalvm.sdk:graal-sdk:20.0.0")
    implementation("org.apache.commons:commons-csv:1.8")
//    implementation("com.github.jkcclemens:khttp:0.1.0")

}

java {
    sourceCompatibility = VERSION_11
    targetCompatibility = VERSION_11
}

kotlin {
    sourceSets{

    }
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
