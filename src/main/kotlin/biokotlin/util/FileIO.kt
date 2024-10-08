@file:JvmName("FileIO")

package biokotlin.util

import org.apache.logging.log4j.LogManager
import java.io.*
import java.net.URL
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

private val myLogger = LogManager.getLogger("biokotlin.util.FileIO")

fun bufferedReader(filename: String): BufferedReader {
    return try {
        if (filename.startsWith("http")) {
            if (filename.endsWith(".gz")) {
                GZIPInputStream(URL(filename).openStream()).bufferedReader()
            } else {
                URL(filename).openStream().bufferedReader()
            }
        } else if (filename.endsWith(".gz")) {
            GZIPInputStream(File(filename).inputStream()).bufferedReader()
        } else {
            File(filename).bufferedReader()
        }
    } catch (e: Exception) {
        throw IllegalStateException("FileIO: getBufferedReader: problem getting reader: " + e.message)
    }
}

fun bufferedWriter(filename: String, append: Boolean = false): BufferedWriter {
    return try {
        if (filename.endsWith(".gz")) {
            BufferedWriter(OutputStreamWriter(GZIPOutputStream(FileOutputStream(File(filename), append))))
        } else {
            BufferedWriter(OutputStreamWriter(FileOutputStream(File(filename), append)))
        }
    } catch (e: java.lang.Exception) {
        myLogger.error("bufferedWriter: problem getting writer for: $filename")
        myLogger.error("bufferedWriter: ${e.message}")
        throw IllegalArgumentException("bufferedWriter: problem getting writer for: $filename  error: ${e.message}")
    }
}