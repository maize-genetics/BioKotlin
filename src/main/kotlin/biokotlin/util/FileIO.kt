@file:JvmName("FileIO")

package biokotlin.util

import java.io.BufferedReader
import java.io.File
import java.net.URL
import java.util.zip.GZIPInputStream

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