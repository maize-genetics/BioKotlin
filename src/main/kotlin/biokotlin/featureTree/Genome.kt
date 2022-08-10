package biokotlin.featureTree

import java.io.File
import biokotlin.featureTree.FeatureType.*

/**
 * TODO import documentation
 */
public sealed interface Genome: Ancestor {
    public override val children: List<GenomeChild>

    /**
     * The assembly units (contigs, chromosomes, or scaffolds) that are direct children of the genome.
     * TODO include diagram.
     */
    public val assemblyUnits: List<AssemblyUnit>

    /**
     * The chromosomes that are direct children of the genome.
     * TODO include diagram.
     */
    public val chromosomes: List<AssemblyUnit>

    /**
     * The scaffolds that are direct children of the genome.
     * TODO include diagram.
     */
    public val scaffolds: List<AssemblyUnit>

    /**
     * The contigs that are direct children of the genome.
     * TODO include diagram.
     */
    public val contigs: List<AssemblyUnit>

    /**
     * Returns a deeply immutable clone of this [Genome].
     */
    public fun immutable(): Genome

    /**
     * Returns a mutable clone of this [Genome].
     */
    public fun mutable(): MutableGenome

    /**
     * Provides factory methods for [Genome] objects, including parsing them from GFF files and creating them
     * from lists.
     */
    public companion object {
        /**
         * Parses the GFF file with path [gff] and returns an immutable [Genome] representation of it.
         */
        public fun fromGFF(gff: String): Genome {
            val registry = mutableMapOf<String, MutableFeature>() //id to feature
            val orphans = mutableListOf<MutableFeature>() //allows for out-of-order features to find parents
            val genome = MutableGenomeImpl() //will be returned in its immutable form

            File(gff).useLines { lines ->
                for (line in lines) {
                    if (line.startsWith("#")) continue //ignores comments

                    val split = line.split("\t") //split into columns
                    val feature = try {
                        val attributes = split[8].split(";") //divides attributes
                        val attributeMap = mutableMapOf<String, String>()
                        for (attribute in attributes) {
                            val tagValue = attribute.split("=")
                            if (tagValue.size < 2) continue //this prevents exceptions from trailing semicolons
                            attributeMap[tagValue[0]] = tagValue[1]
                        }
                        val score = split[5].toDoubleOrNull()
                        val strand = if (split[5][0] == '.') null else split[5][0]
                        val phase = split[7].toIntOrNull()

                        when (val type = FeatureType.convert(split[2])) {
                            CHROMOSOME, SCAFFOLD, CONTIG -> MutableAssemblyUnitImpl(
                                MutableFeatureImpl(
                                    split[0], //seqid
                                    split[1], //source
                                    type,
                                    split[3].toInt(), //start
                                    split[4].toInt(), //end
                                    score,
                                    strand,
                                    phase,
                                    attributeMap
                                )
                            )
                            GENE -> MutableGeneImpl(
                                MutableAncestralFeature(
                                    emptyList<MutableTranscript>(),
                                    split[0],
                                    split[1],
                                    type,
                                    split[3].toInt(),
                                    split[4].toInt(),
                                    score,
                                    strand,
                                    phase,
                                    attributeMap
                                )
                            )
                            TRANSCRIPT -> MutableTranscriptImpl(
                                MutableAncestralFeature(
                                    emptyList<MutableTranscriptChild>(),
                                    split[0],
                                    split[1],
                                    type,
                                    split[3].toInt(),
                                    split[4].toInt(),
                                    score,
                                    strand,
                                    phase,
                                    attributeMap
                                )
                            )
                            LEADER, EXON, CODING_SEQUENCE, TERMINATOR -> MutableTranscriptChildImpl(
                                MutableFeatureImpl(
                                    split[0],
                                    split[1],
                                    type,
                                    split[3].toInt(),
                                    split[4].toInt(),
                                    score,
                                    strand,
                                    phase,
                                    attributeMap
                                )
                            )
                        }
                    } catch (e: IllegalArgumentException) {
                        //this occurs when a type cannot be converted properly
                        print(e.message)
                        println(" Therefore, the feature represented by the following line has been skipped:")
                        println(line)
                        continue
                    } catch (ex: IndexOutOfBoundsException) {
                        println("Warning: The following row in your follow does not contain 9 tab-delimited columns, so it has been skipped:")
                        println(line)
                        continue
                        //TODO document this exception handling
                    } catch (nf: NumberFormatException) {
                        println("Start/end must be valid integers.")
                        throw nf
                    }

                    //Enroll features with an ID into the registry
                    if (feature.attributes.contains("ID")) {
                        registry[feature.attributes["ID"]!!] = feature
                    }

                    //Match a feature to its parent or put it in the orphan list
                    //or throw an exception if the stated parent/child relationship is invalid
                    if (feature.attributes.contains("Parent")) {
                        val parentID = feature.attributes["Parent"]

                        if (registry.contains(parentID)) {
                            val parent = registry[parentID]!!
                            //since the parent is supposed to be an ancestor and feature, it should only be
                            //a gene or transcript
                            try {
                                addChildOrThrow(feature, parent as MutableAncestor)
                            } catch (e: ClassCastException) {
                                //if the cast failed then parent cannot support children
                                throw IllegalChild(feature, parent)
                            }
                        } else {
                            orphans.add(feature)
                            //TODO warning message
                        }
                    } else { //not containing a Parent attribute means it's a direct child of the genome
                        addChildOrThrow(feature, genome)
                    }
                }
                for (orphan in orphans) {
                    val parentID = orphan.attributes["Parent"]
                    if (registry.contains(parentID)) {
                        val parent = registry[parentID]!!
                        try {
                            addChildOrThrow(orphan, parent as MutableAncestor)
                        } catch (_: ClassCastException) {
                            throw IllegalChild(orphan, parent)
                        }
                    } else println(
                        "Warning: A feature's parent (Parent=${orphan.attributes["Parent"]}) could not be located, so " +
                                "this feature and its descendants will not be added to the tree."
                    )
                }
                return genome.immutable()
            }
        }

        /**
         * Moves [child] to [parent].
         * Throws [IllegalParentChild] if the parent/child relationship is invalid.
         */
        private fun addChildOrThrow(child: MutableFeature, parent: MutableAncestor) {
            //This is exhaustive! The compiler just doesn't know it.
            when (parent) {
                is MutableGenome ->
                    if (child is MutableGenomeChild) child.moveTo(parent)
                    else throw IllegalParentChild(child, parent)
                is MutableGene ->
                    if (child is MutableTranscript) child.moveTo(parent)
                    else throw IllegalParentChild(child, parent)
                is MutableTranscript ->
                    if (child is MutableTranscriptChild) child.moveTo(parent)
                    else throw IllegalParentChild(child, parent)
            }
        }

        /**
         * Creates a [Genome] from a list of features. The [Genome] will contain all features in [list] as well
         * as their ancestors. If [includeDescendants] is ```true```, then all descendants of features within [list]
         * will be included as well. Overlaps in the features in [list] or their ancestors or descendants
         * will be resolved so that there are no duplicates in the returned [Genome].
         * All features within the returned [Genome] will be deep copies of the originals.
         */
        public fun fromList(list: List<Feature>, includeDescendants: Boolean = true): Genome {
            TODO("Not yet implemented")
        }

    }
}

/**
 * A mutable representation of [Genome]. Prefer using an immutable [Genome] unless mutability is actually needed.
 * TODO more documentation of the mutability framework, how to use it effectively, etc.
 */
public sealed interface MutableGenome: Genome, MutableAncestor {
    public override val children: List<MutableGenomeChild>
    public override val assemblyUnits: List<MutableAssemblyUnit>
    public override val chromosomes: List<MutableAssemblyUnit>
    public override val scaffolds: List<MutableAssemblyUnit>
    public override val contigs: List<MutableAssemblyUnit>

    /**
     * Inserts and returns new, mutable chromosome as a direct child of the [MutableGenome].
     * Any value with the key "Parent" within [attributes] will not be included in the attributes property
     * of the chromosome.
     */
    public fun insertChromosome(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Char?,
        phase: Int?,
        attributes: Map<String, String>
    ): MutableAssemblyUnit

    /**
     * Inserts and returns new, mutable scaffold as a direct child of the [MutableGenome].
     * Any value with the key "Parent" within [attributes] will not be included in the attributes property
     * of the scaffold.
     */
    public fun insertScaffold(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Char?,
        phase: Int?,
        attributes: Map<String, String>
    ): MutableAssemblyUnit

    /**
     * Inserts and returns new, mutable scaffold as a direct child of the [MutableGenome].
     * Any value with the key "Parent" within [attributes] will not be included in the attributes property
     * of the scaffold.
     */
    public fun insertContig(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Char?,
        phase: Int?,
        attributes: Map<String, String>
    ): MutableAssemblyUnit

    /**
     * Removes [child] from this genome's children, if it exists. Returns true iff the child was removed.
     */
    public fun removeChild(child: MutableGenomeChild): Boolean
}


//TODO mutable genome blank constructor, mutable genome children constructor,
//TODO proper interaction with interfaac

internal open class GenomeImpl protected constructor (
    protected open val delegate: AncestorImpl
) : Ancestor by delegate,
    Genome {
    internal constructor(): this(AncestorImpl(emptyList()))

    /* INVARIANTS:
    1. The children of this genome are all GenomeChildren
     */

    public override val children: List<GenomeChild> get() = delegate.children.map { it as GenomeChild }
    public override val assemblyUnits: List<AssemblyUnit> get() = children.filterIsInstance<AssemblyUnit>()
    public override val chromosomes: List<AssemblyUnit> get() = assemblyUnits.filter { it.type == CHROMOSOME }
    public override val scaffolds: List<AssemblyUnit> get() = assemblyUnits.filter { it.type == SCAFFOLD }
    public override val contigs: List<AssemblyUnit> get() = assemblyUnits.filter { it.type == CONTIG }

    public override fun mutable(): MutableGenome {
        TODO()
    }
    public override fun immutable(): Genome {
        TODO()
    }

    //TODO special genome behaviors

    /**
     * Creates a mutable clone of the genome and returns the cloned MutableFeature that corresponds to [feature]
     */
    internal fun mutable(feature: Feature): MutableFeature {
        TODO()
    }

    /**
     * Creates a deeply immutable clone of the genome and returns the cloned Feature that corresponds to [feature]
     */
    internal fun immutable(feature: Feature): Feature {
        TODO()
    }

    /**
     * Returns a deeply immutable [Genome] that contains all the features of this genome and of [other].
     */
    public fun appendedWith(other: Genome): Genome {
        TODO()
    }

    /**
     * Returns a deeply immutable [Genome] that contains all the features of this genome and of [other].
     */
    public operator fun plus(other: Genome): Genome = appendedWith(other)

    /**
     * Represents the [Genome] as a GFF file.
     */
    public override fun toString(): String = descendantsToString()

}
internal class MutableGenomeImpl private constructor (
    override val delegate: MutableAncestorImpl
) : GenomeImpl(delegate),
    MutableGenome {
    internal constructor(): this(MutableAncestorImpl(emptyList()))

    /* INVARIANTS:
    1. The children of this genome are all MutableGenomeChildren
     */
    public override val children: List<MutableGenomeChild> get() = delegate.children.map { it as MutableGenomeChild }
    public override val assemblyUnits: List<MutableAssemblyUnit> get() = children.filterIsInstance<MutableAssemblyUnit>()
    public override val chromosomes: List<MutableAssemblyUnit> get() = assemblyUnits.filter { it.type == CHROMOSOME }
    public override val scaffolds: List<MutableAssemblyUnit> get() = assemblyUnits.filter { it.type == SCAFFOLD }
    public override val contigs: List<MutableAssemblyUnit> get() = assemblyUnits.filter { it.type == CONTIG }

    /**
     * Used for the implementations of inserts and moveTo.
     */
    internal fun addChild(child: MutableGenomeChild): Unit = delegate.addChild(child)
    public override fun insertChromosome(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Char?,
        phase: Int?,
        attributes: Map<String, String>
    ): MutableAssemblyUnit {
        TODO("Not yet implemented")
    }

    public override fun insertScaffold(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Char?,
        phase: Int?,
        attributes: Map<String, String>
    ): MutableAssemblyUnit {
        TODO("Not yet implemented")
    }

    public override fun insertContig(
        seqid: String,
        source: String,
        start: Int,
        end: Int,
        score: Double?,
        strand: Char?,
        phase: Int?,
        attributes: Map<String, String>
    ): MutableAssemblyUnit {
        TODO("Not yet implemented")
    }

    public override fun flatten(): List<MutableFeature>  = super<MutableGenome>.flatten()

    public override fun removeChild(child: MutableGenomeChild): Boolean = delegate.removeChild(child)
}

/*
Cloning algorithm:
immutable/mutable is called on a feature. the call is passed through its delegate and propagated up to the genome.
the genome receives the call and a reference to the original receiver object. it clones downwards and when it finds
the original receiver, it stores a reference to its clone to be returned at the end. after cloning downwards, parents
are injected. this allows for maintining immutability.
 */