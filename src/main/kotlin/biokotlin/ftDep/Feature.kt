////Per Kotlin style convention, libraries should have redundant visibility modifiers
//@file:Suppress("RedundantVisibilityModifier")
//
//package biokotlin.ftDep
//
//import biokotlin.featureTree.*
//
//public sealed interface Feature {
//    public val seqid: String
//    public val source: String
//    public val type: FeatureType
//    public val start: Int
//        get() = ranges.minimum()
//    public val end: Int
//        get() = ranges.maximum()
//    public val range: IntRange
//        get() = ranges[0]
//    public val ranges: List<IntRange>
//    public val score: Double?
//    public val strand: Strand
//    /**
//     * `phases[0]`
//     */
//    public val phase: Phase
//        get() = phases[0]
//    public val phases: List<Phase>
//    public val parent: Parent
//    public val genome: Genome
//    public val id: String?
//
//    /**
//     * The length of the feature, equivalent to `end - start + 1`
//     */
//    public val length: Int
//        get() = end - start + 1
//
//    /**
//     * Discontinuous annotations are represented in [biokotlin.featureTree] as a single [Feature] object
//     * with a [multiplicity] equal to the number of discontinuous segments of the annotation. [ranges] represents
//     * all the start-end ranges of this discontinuous feature while [phases] represents all the phases.
//     */
//    public val multiplicity: Int
//        get() = ranges.size
//
//    public val discontinuous: Boolean
//        get() = multiplicity > 1
//
//    /**
//     * The [String] representation of this [Feature] as it appears as a row in a GFF file.
//     */
//    public fun asRow(): String
//
//    public fun attribute(tag: String): List<String>
//    public fun allAttributes(): Map<String, List<String>>
//}
//
//public sealed interface MutableFeature : Feature {
//    public override var seqid: String
//    public override var source: String
//    public override var start: Int
//    public override var end: Int
//    public override var range: IntRange
//    public override var score: Double?
//    public override var strand: Strand
//    public override var phase: Phase
//    public override val parent: MutableParent
//    public override val genome: MutableGenome
//    public override var id: String?
//
//
//    /**
//     * Deletes this [MutableFeature] from its [MutableGenome]. Subsequent attempts to read information from this
//     * will result in [DeletedAccessException].
//     */
//    public fun delete()
//
//    /**
//     * Associates [value] with [tag] in the attributes of this feature.
//     * @throws IllegalArgumentException if [tag] is "Parent". The "Parent" attribute is determined by the
//     * actual topology of the tree and cannot be directly modified. Hint: see [copyTo].
//     * @throws IllegalArgumentException if [tag] is "ID". Hint: modify the [id] property instead.
//     */
//    public fun addAttribute(tag: String, value: String)
//
//    /**
//     * Associates all elements of [values] with [tag] in the attributes of this feature.
//     * @throws IllegalArgumentException if [tag] is "Parent". The "Parent" attribute is determined by the
//     * actual topology of the tree and cannot be directly modified. Hint: see [copyTo].
//     * @throws IllegalArgumentException if [tag] is "ID". Hint: modify the [id] property instead.
//     */
//    public fun addAttributes(tag: String, values: Iterable<String>)
//
//    /**
//     * Associates [tag] with [value] in the attributes of this feature. Removes all other associations with [tag]
//     * that may have existed prior.
//     * @throws IllegalArgumentException if [tag] is "Parent". The "Parent" attribute is determined by the
//     * actual topology of the tree and cannot be directly modified. Hint: see [copyTo].
//     * @throws IllegalArgumentException if [tag] is "ID". Hint: modify the [id] property instead.
//     */
//    public fun setAttribute(tag: String, value: String)
//
//    /**
//     * Associates [tag] with the elements of [values] in the attributes of this feature. Removes all other associations
//     * with [tag] that may have existed prior.
//     */
//    public fun overwriteAttributes(tag: String, values: Iterable<String>)
//
//    /**
//     * Clears all associations with [tag] in the attributes of this feature.
//     */
//    public fun clearAttribute(tag: String)
//    public fun addDiscontinuity(range: IntRange, phase: Phase)
//    public fun setDiscontinuity(index: Int, range: IntRange, phase: Phase)
//    public fun overwriteDiscontinuities(ranges: Iterable<IntRange>, phases: Iterable<Phase>)
//
//}
//
//public sealed interface GenomeChild : Feature {
//    public override val parent: Genome
//
//    public fun copyTo(genome: MutableGenome): MutableGenomeChild
//}
//
//public sealed interface MutableGenomeChild : GenomeChild, MutableFeature {
//    public override val parent: MutableGenome
//}
//
//public sealed interface AssemblyUnit : Feature, GenomeChild {
//    public fun genes(): List<Gene> = genome.genes
//
//    public override fun copyTo(genome: MutableGenome): MutableAssemblyUnit
//}
//
//public sealed interface MutableAssemblyUnit : AssemblyUnit, MutableGenomeChild, MutableFeature {
//    public override val parent: MutableGenome
//    public override fun genes(): List<MutableGene> = genome.genes
//}
//
//public sealed interface Chromosome : AssemblyUnit {
//    public override fun copyTo(genome: MutableGenome): MutableChromosome
//}
//
//public sealed interface MutableChromosome : MutableAssemblyUnit, Chromosome
//public sealed interface Contig : AssemblyUnit {
//    public override fun copyTo(genome: MutableGenome): MutableContig
//}
//
//public sealed interface MutableContig : MutableAssemblyUnit, Contig
//public sealed interface Scaffold : AssemblyUnit {
//    public override fun copyTo(genome: MutableGenome): MutableScaffold
//}
//
//public sealed interface MutableScaffold : MutableAssemblyUnit, Scaffold
//public sealed interface Gene : Feature, GenomeChild, Parent, TranscriptAncestor, TranscriptChildAncestor {
//    public override val parent: Genome
//    public override val children: List<Transcript>
//    public val transcripts
//        get() = children
//
//    public override fun copyTo(genome: MutableGenome): MutableGene
//}
//
//public sealed interface MutableGene : MutableParent, MutableGenomeChild, Gene, MutableTranscriptAncestor, MutableTranscriptChildAncestor {
//    public override val parent: MutableGenome
//    public override val children: List<MutableTranscript>
//
//    public override val transcripts
//        get() = children
//
//    /**
//     * Inserts transcript with defined data into this gene's children.
//     * @return the inserted transcript.
//     */
//    public fun insertTranscript(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): MutableTranscript
//}
//
//public sealed interface Transcript : Feature, Parent, TranscriptChildAncestor {
//    public override val parent: Gene
//    public override val children: List<TranscriptChild>
//
//    public val leaders
//        get() = filterChildren<Leader>()
//    public val exons
//        get() = filterChildren<Exon>()
//    public val codingSequences
//        get() = filterChildren<CodingSequence>()
//    public val terminators
//        get() = filterChildren<Terminator>()
//
//    public fun copyTo(gene: MutableGene): MutableTranscript
//}
//
//public sealed interface MutableTranscript : MutableFeature, MutableParent, Transcript, MutableTranscriptChildAncestor {
//    public override val parent: MutableGene
//    public override val children: List<MutableTranscriptChild>
//
//    public override val leaders
//        get() = filterChildren<MutableLeader>()
//    public override val exons
//        get() = filterChildren<MutableExon>()
//    public override val codingSequences
//        get() = filterChildren<MutableCodingSequence>()
//    public override val terminators
//        get() = filterChildren<MutableTerminator>()
//
//    /**
//     * Inserts a leader into this transcript's children.
//     * @return the inserted leader
//     */
//    public fun insertLeader(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): Leader
//
//    /**
//     * Inserts a new exon into this transcript's children.
//     * @return the inserted exon
//     */
//    public fun insertExon(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): Exon
//
//    /**
//     * Inserts a new coding sequence into this transcript's children.
//     * @return the inserted coding sequence
//     */
//    public fun insertCodingSequence(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): CodingSequence
//
//    /**
//     * Inserts a new terminator into this transcript's children.
//     * @return the inserted terminator
//     */
//    public fun insertTerminator(
//        seqid: String,
//        source: String,
//        ranges: Iterable<IntRange>,
//        score: Double?,
//        strand: Strand,
//        phases: Iterable<Phase>,
//        attributes: Attributes
//    ): Terminator
//}
//
//public sealed interface TranscriptChild : Feature {
//    public override val parent: Transcript
//    public val transcript: Transcript
//        get() = parent
//}
//
//public sealed interface MutableTranscriptChild : MutableFeature, TranscriptChild {
//    public override val parent: MutableTranscript
//
//    public override val transcript: MutableTranscript
//        get() = parent
//}
//
//public sealed interface Leader : TranscriptChild
//public sealed interface MutableLeader : MutableTranscriptChild, Leader
//public sealed interface Exon : TranscriptChild
//public sealed interface MutableExon : MutableTranscriptChild, Exon
//public sealed interface CodingSequence : TranscriptChild
//public sealed interface MutableCodingSequence : MutableTranscriptChild, CodingSequence
//public sealed interface Terminator : TranscriptChild
//public sealed interface MutableTerminator : MutableTranscriptChild, Terminator