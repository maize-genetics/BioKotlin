@file:JvmName("BEDtoGFF")

package biokotlin.seqIO

import biokotlin.util.*

fun bedToGff3(bedfile: String, gffFile: String) {
    val writer = bufferedWriter(gffFile, false)

    /*
    0: chrom
    1: chromStart (starts at 0)
    2: chromEnd (EXCLUSIVE)
    3: name
    5: score
    6: strand
    */
    bufferedReader(bedfile).forEachLine { line ->
        val split = line.split("\t")
        //Ignore comments and headers
        if (!(line.startsWith("#") || line.startsWith("browser") || line.startsWith("track"))) {
            /*
            0: seqid (equivalent to chrom)
            1: source (leave as .)
            2: type (comes from name)
            3: start (1-BASE OFFSET)
            4: end (1-BASE OFFSET)
            5: score
            7: strand
            8: phase (leave as .)
            9: attributes (leave as .)
            */
            writer.write(split[0] + "\t")
            writer.write(".\t")
            val truncatedName = split[3].substringBeforeLast(".")
            if (truncatedName.endsWith("_C")) {
                writer.write("exon\t")
            } else if (truncatedName.endsWith("_I")) {
                writer.write("intron\t")
            } else {
                writer.write("gene\t")
            }
            writer.write(split[1] + 1 + "\t")
            writer.write(split[1] + "\t")
            writer.write(split[5] + "\t")
            writer.write(split[6] + "\t")
            writer.write(".")
            writer.write(".")
        }
    }
}