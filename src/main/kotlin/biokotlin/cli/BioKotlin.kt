package biokotlin.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.subcommands

/**
 * This class is the main class for the BioKotlin command line interface.
 * It is a subclass of CliktCommand and is used to provide users a BioKotlin command line interface.
 *
 */
class BioKotlin : CliktCommand() {
    override fun run() = Unit
}

fun main(args: Array<String>) = BioKotlin()
    .subcommands(MafToGvcfConverter(), ValidateGVCFs())
    .main(args)