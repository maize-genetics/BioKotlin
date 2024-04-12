package biokotlin.cli

import com.github.ajalt.clikt.core.CliktCommand
import org.apache.logging.log4j.LogManager

class MergeGVCFs : CliktCommand(help = "Merge GVCF files into Single VCF file") {

    private val myLogger = LogManager.getLogger(MergeGVCFs::class.java)

    override fun run() {
        TODO("Not yet implemented")
    }

}