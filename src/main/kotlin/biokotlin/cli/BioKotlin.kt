package biokotlin.cli

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.core.subcommands

class BioKotlin : CliktCommand() {
    override fun run() = Unit
}

fun main(args: Array<String>) = BioKotlin()
    .subcommands(MafToGvcfConverter())
    .main(args)