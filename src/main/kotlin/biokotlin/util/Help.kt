package biokotlin.util

fun Any.help():String {
    ///Logic to go to Kotlin help StdLib, perhaps a few key imports Guava, Krangl, etc.
    //This is over
    //By reflection we should be able to get lots of packages and paths.
    return """
            To learn more about ${this::class} go to:
            https://javadoc.io/static/com.google.guava/guava/29.0-jre/com/google/common/collect/Range.html
        """.trimIndent()
}

fun Any.javadoc():String {
    ///Logic to go to Kotlin help StdLib, perhaps a few key imports Guava, Krangl, etc.
    //This is over
    //By reflection we should be able to get lots of packages and paths.
    return TODO()
}

fun Any.example():String {
    //Provide example code for using this class
    return TODO()
}