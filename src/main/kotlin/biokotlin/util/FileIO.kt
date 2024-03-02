@file:JvmName("FileIO")

package biokotlin.util

import java.io.*
import java.net.URL
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

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
        throw IllegalStateException("FileIO: bufferedReader: problem getting reader: " + e.message)
    }
}

fun bufferedWriter(filename: String, append: Boolean = false): BufferedWriter {
    return try {
        if (filename.endsWith(".gz")) {
            BufferedWriter(OutputStreamWriter(GZIPOutputStream(FileOutputStream(filename, append))))
        } else {
            BufferedWriter(OutputStreamWriter(FileOutputStream(filename, append)))
        }
    } catch (e: Exception) {
        throw IllegalStateException("FileIO: bufferedWriter: problem getting writer: " + e.message)
    }
}