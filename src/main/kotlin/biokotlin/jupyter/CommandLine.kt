package biokotlin.jupyter

import java.io.File
import java.lang.ProcessBuilder.*
import java.util.concurrent.*
import java.io.IOException

/**
Runs a string command line
e.g. "ls".runCommand(File("tmp/wkdir"))
 */
fun String.runCommand(workingDir: File): String? {
    try {
        val parts = this.split("\\s".toRegex())
        val proc = ProcessBuilder(*parts.toTypedArray())
            .directory(workingDir)
            .redirectOutput(ProcessBuilder.Redirect.PIPE)
            .redirectError(ProcessBuilder.Redirect.PIPE)
            .start()
        proc.waitFor(60, TimeUnit.MINUTES)
        return proc.inputStream.bufferedReader().readText()
    } catch(e: IOException) {
        e.printStackTrace()
        return null
    }
}