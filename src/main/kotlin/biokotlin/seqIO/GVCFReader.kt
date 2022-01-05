package biokotlin.seqIO

import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import biokotlin.seq.SeqRecord
import com.google.common.collect.ImmutableMap
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import java.io.File

/**
 * This reads a GVCF formatted file. The file is required to only have one sample.
 * The file can only be read one time. Either as an iterator (hasNext(), next()),
 * repeatedly calling read(), or using readAll().
 */
class GVCFReader(val gvcfFile: String, val referenceFile: String) : SequenceIterator {

    private val resultChannel = Channel<Deferred<Pair<String, SeqRecord>>>(5)
    private val resultIterator = resultChannel.iterator()

    private var reference: Map<String, SeqRecord> = reader(referenceFile).readAll()
    private var taxon: String

    init {

        VCFFileReader(File(gvcfFile), false).use { reader ->

            val header = reader.fileHeader

            val altHeaderLines = header.idHeaderLines
                .filterIsInstance<VCFAltHeaderLine>()
                .map { Pair(it.id, it.toString().substringAfter("Description=\"").substringBefore("\"")) }
                .toMap()

            val samples = header.sampleNamesInOrder

            require(samples.size == 1) { println("Filename: $gvcfFile: Number of samples should be 1 but is: ${samples.size}") }

            taxon = samples[0]

            CoroutineScope(Dispatchers.IO).launch {
                readFile(reader)
            }

        }

    }

    /**
     * Returns next SeqRecord or null if finished reading.
     */
    override fun read(): SeqRecord? {
        return runBlocking {
            resultChannel.receiveCatching().getOrNull()?.await()?.second
        }
    }

    /**
     * Reads a GVCF file (.gvcf) and returns a map
     * of Chromosomes (String) to Sequence (SeqRecord)
     */
    override fun readAll(): Map<String, SeqRecord> {
        return runBlocking { combineResults() }
    }

    override fun hasNext(): Boolean {
        return runBlocking {
            resultIterator.hasNext()
        }
    }

    override fun next(): SeqRecord {
        return runBlocking {
            resultIterator.next().await().second
        }
    }

    private suspend fun readFile(reader: VCFFileReader) =
        withContext(Dispatchers.IO) {

            var chromosomeChannel = Channel<VariantContext>()

            var currentChromosome = ""

            reader.forEach { context ->

                val chromosome = context.contig
                if (currentChromosome == "") {
                    currentChromosome = chromosome
                    resultChannel.send(CoroutineScope(Dispatchers.Default).async {
                        sequenceFromContexts(currentChromosome, chromosomeChannel)
                    })
                }

                if (chromosome != currentChromosome) {

                    chromosomeChannel.close()

                    currentChromosome = chromosome
                    // Start new Channel for next chromosome
                    chromosomeChannel = Channel()

                    resultChannel.send(CoroutineScope(Dispatchers.Default).async {
                        sequenceFromContexts(currentChromosome, chromosomeChannel)
                    })

                }
                chromosomeChannel.send(context)

            }

            chromosomeChannel.close()

            resultChannel.close()

        }

    // This method constructs the sequence for a single contig (i.e. chromosome)
    private suspend fun sequenceFromContexts(
        chromosome: String,
        channel: Channel<VariantContext>
    ): Pair<String, SeqRecord> {

        // Holds resulting sequence
        val result = StringBuilder()

        for (context in channel) {

            val contig = context.contig
            val start = context.start
            val end = context.end
            val genotypes = context.genotypes[0]
            val alleles = genotypes.alleles
            val refAllele = context.reference.baseString

            val allele = alleles.first()

            if (allele.isReference) {

                val referenceSeq = reference[contig]
                if (referenceSeq == null) {
                    throw IllegalStateException("GVCFReader: sequenceFromContexts: contig: $contig doesn't exist in reference.")
                } else {
                    // subtracting 1 from start and end to convert 1-based to 0-based
                    val refStr = (referenceSeq as NucSeq)[start - 1, end - 1].toString()
                    check(refStr.startsWith(refAllele))
                    { "GVCFReader: sequenceFromContexts: reference allele: $refAllele for contig: $contig start: $start doesn't match start sequence from reference file: $refStr" }
                    result.append(refStr)
                }

            } else {
                val alleleStr = allele.baseString
                if (alleleStr != "*") result.append(allele.baseString)
            }

        }

        return Pair(chromosome, NucSeqRecord(NucSeq(result.toString()), chromosome))

    }

    private suspend fun combineResults(): Map<String, SeqRecord> {

        val result = ImmutableMap.builder<String, SeqRecord>()
        for (deferred in resultChannel) {
            val chrSeq = deferred.await()
            println("GVCFReader: combineResults: adding chromosome ${chrSeq.first} to results")
            result.put(chrSeq.first, chrSeq.second)
        }
        return result.build()

    }

}

fun main() {
    val gvcfFile =
        "/Users/tmc46/hackathon/20210916_tribe/su1_gvcfs_v2/AN20TSCR000225_CKDL200169516-1a-AK6517-AK6687_HLTF5DSXY_L4_srt_Su1.gvcf"
    val reference = "/Users/tmc46/hackathon/20210916_tribe/Sbicolor_454_v3_numeric.0.1.fa"
    val result = GVCFReader(gvcfFile, reference)
    result
        .forEach {
            println("${it.id}: ${it.toString().substring(0, 100)}")
        }
}