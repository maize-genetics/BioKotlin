package biokotlin.cli

import biokotlin.util.setupDebugLogging
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.subcommands

/**
 * This class is the main class for the BioKotlin command line interface.
 * It is a subclass of CliktCommand and is used to provide users a BioKotlin command line interface.
 *
 */
class Biokotlin : CliktCommand() {

    init {
        setupDebugLogging()
    }

    override fun run() = Unit
}

fun main(args: Array<String>) = Biokotlin()
    .subcommands(
        MafToGvcfConverter(),
        ValidateGVCFs(),
        MergeGVCFs(),
        ValidateVCFs(),
        MutateProteins()
    )
    .main(args)