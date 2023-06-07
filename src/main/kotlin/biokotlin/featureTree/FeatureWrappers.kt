//Per Kotlin style convention, libraries should have redundant visibility modifiers
@file:Suppress("RedundantVisibilityModifier")

package biokotlin.featureTree

/**
 * TODO Interfaces
 * 8. gffForNewSpecies()?
 * 9. parent aliases
 *
 * TODO Safe mutability
 * 1. Ensure attributes handled safely
 * 1a. Do not allow modifications of Parent attribute
 * 1b. Propagate modifications of ID attribute
 * 1c. Do not allow rep exposure of underlying map
 * 2. Ensure start/end invariants maintained
 * 3. Ensure that CDS cannot have its phase be unspecified
 * 4. Rich error messages
 *
 * TODO I/O
 * 1. Parse from file
 * 2. Connect to FASTA
 * 3. Use FASTA to output BioKotlin Seq data
 * 4. Handle escape sequences correctly
 */
public sealed interface Feature {
    public val seqid: String
    public val source: String
    public val type: FeatureType
    public val start: UInt
    public val end: UInt
    public val score: Double?
    public val strand: Strand
    public val phase: Phase
    public val parent: Parent
    public val genome: Genome
    public val id: String?

    /**
     * The length of the feature, equivalent to `end - start + 1`
     */
    public val length: UInt
        get() = end - start + 1u

    /**
     * The [String] representation of this [Feature] as it appears as a row in a GFF file.
     */
    public fun asRow(): String {
        val scoreString = score?.toString() ?: "."
        val phaseString = phase.gffName
        val strandString = strand.gffName

        val attributesString = StringBuilder()

        allAttributes().forEach { (tag, values) ->
            attributesString.append("$tag = ")
            attributesString.append(values.fold("") { acc, elem ->
                acc + elem
            }.trimEnd(','))
        }

        return "$seqid\t$source\t${type.gffName}\t$start\t$end\t$scoreString\t$strandString\t$phaseString\t${attributesString}\n"
    }

    public fun attribute(tag: String): List<String>
    public fun allAttributes(): Map<String, List<String>>
}

public sealed interface MutableFeature : Feature {
    public override var seqid: String
    public override var source: String
    public override var start: UInt
    public override var end: UInt
    public override var score: Double?
    public override var strand: Strand
    public override var phase: Phase
    public override val parent: MutableParent
    public override val genome: MutableGenome
    public override var id: String?

    /**
     * Deletes this [MutableFeature] from its [MutableGenome]. Subsequent attempts to read information from this
     * will result in [DeletedAccessException].
     */
    public fun delete()
    public fun addAttribute(tag: String, value: String)
    public fun setAttribute(tag: String, value: String)
    public fun overwriteAttributes(tag: String, values: Iterable<String>)
    public fun clearAttribute(tag: String)

}

public sealed interface GenomeChild : Feature {
    public override val parent: Genome

    public fun copyTo(genome: MutableGenome): MutableGenomeChild
}

public sealed interface MutableGenomeChild : GenomeChild, MutableFeature {
    public override val parent: MutableGenome
}

public sealed interface AssemblyUnit : Feature, GenomeChild {
    public fun genes(): Sequence<Gene> {
        TODO()
    }

    public override fun copyTo(genome: MutableGenome): MutableAssemblyUnit
}

public sealed interface MutableAssemblyUnit : AssemblyUnit, MutableGenomeChild, MutableFeature {
    public override val parent: MutableGenome
    public override fun genes(): Sequence<MutableGene> {
        TODO()
    }
}

public sealed interface Chromosome : AssemblyUnit {
    public override fun copyTo(genome: MutableGenome): MutableChromosome
}

public sealed interface MutableChromosome : MutableAssemblyUnit, Chromosome
public sealed interface Contig : AssemblyUnit {
    public override fun copyTo(genome: MutableGenome): MutableContig
}

public sealed interface MutableContig : MutableAssemblyUnit, Contig
public sealed interface Scaffold : AssemblyUnit {
    public override fun copyTo(genome: MutableGenome): MutableScaffold
}

public sealed interface MutableScaffold : MutableAssemblyUnit, Scaffold
public sealed interface Gene : Feature, GenomeChild, Parent, TranscriptAncestor, TranscriptChildAncestor {
    public override val parent: Genome
    public override val children: List<Transcript>
    public val transcripts
        get() = children

    public override fun copyTo(genome: MutableGenome): MutableGene
}

public sealed interface MutableGene : MutableParent, MutableGenomeChild, Gene, MutableTranscriptAncestor, MutableTranscriptChildAncestor {
    public override val parent: MutableGenome
    public override val children: List<MutableTranscript>

    public override val transcripts
        get() = children

    /**
     * Inserts transcript with defined data into this gene's children.
     * @return the inserted transcript.
     */
    public fun insertTranscript(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): MutableTranscript
}

public sealed interface Transcript : Feature, Parent, TranscriptChildAncestor {
    public override val parent: Gene
    public override val children: List<TranscriptChild>

    public val leaders
        get() = filterChildren<Leader>()
    public val exons
        get() = filterChildren<Exon>()
    public val codingSequences
        get() = filterChildren<CodingSequence>()
    public val terminators
        get() = filterChildren<Terminator>()

    public fun copyTo(gene: MutableGene): MutableTranscript
}

public sealed interface MutableTranscript : MutableFeature, MutableParent, Transcript, MutableTranscriptChildAncestor {
    public override val parent: MutableGene
    public override val children: List<MutableTranscriptChild>

    public override val leaders
        get() = filterChildren<MutableLeader>()
    public override val exons
        get() = filterChildren<MutableExon>()
    public override val codingSequences
        get() = filterChildren<MutableCodingSequence>()
    public override val terminators
        get() = filterChildren<MutableTerminator>()

    /**
     * Inserts a leader into this transcript's children.
     * @return the inserted leader
     */
    public fun insertLeader(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): Leader

    /**
     * Inserts a new exon into this transcript's children.
     * @return the inserted exon
     */
    public fun insertExon(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): Exon

    /**
     * Inserts a new coding sequence into this transcript's children.
     * @return the inserted coding sequence
     */
    public fun insertCodingSequence(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): CodingSequence

    /**
     * Inserts a new terminator into this transcript's children.
     * @return the inserted terminator
     */
    public fun insertTerminator(
        seqid: String,
        source: String,
        start: UInt,
        end: UInt,
        score: Double?,
        strand: Strand,
        phase: Phase,
        attributes: Attributes
    ): Terminator
}

public sealed interface TranscriptChild : Feature {
    public override val parent: Transcript
    public val transcript: Transcript
        get() = parent
}

public sealed interface MutableTranscriptChild : MutableFeature, TranscriptChild {
    public override val parent: MutableTranscript

    public override val transcript: MutableTranscript
        get() = parent
}

public sealed interface Leader : TranscriptChild
public sealed interface MutableLeader : MutableTranscriptChild, Leader
public sealed interface Exon : TranscriptChild
public sealed interface MutableExon : MutableTranscriptChild, Exon
public sealed interface CodingSequence : TranscriptChild
public sealed interface MutableCodingSequence : MutableTranscriptChild, CodingSequence
public sealed interface Terminator : TranscriptChild
public sealed interface MutableTerminator : MutableTranscriptChild, Terminator
