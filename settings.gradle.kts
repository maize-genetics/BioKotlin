pluginManagement {
    repositories {
        mavenCentral()
        gradlePluginPortal()
        maven ("https://dl.bintray.com/kotlin/kotlin-eap")
        //jitpack required for DataFrame Gradle integration
        //This is used for code generation for DataFrame Schema, however, it does not work
        //https://github.com/Kotlin/dataframe/tree/eb9ec4fb90f906f6a98e69b9c5a0369009d34bbb/plugins/gradle/codegen
        maven(url="https://jitpack.io")
    }
}