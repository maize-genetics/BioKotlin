package biokotlin.featureTree

import biokotlin.featureTree.FeatureType.*

/*
TODO internal soloMutable(), internal soloImmutable()

TODO library formatting pass
 */

/**
 * Represents a chromosomes, scaffold, or contig. Provides the [genes] method as a convenience accessor.
 * TODO import documentation.
 */
public sealed interface AssemblyUnit: GenomeChild {
    /**
     * Returns a list of [Gene] within the same [Genome] as this [AssemblyUnit] that share its seqid, in order.
     * While these genes are not considered children of the assembly unit (because standard GFF formatting does
     * not consider them as such), they do have a biological connection to the assembly unit, so this function
     * enable easy access of them.
     */
    public fun genes(): List<Gene>
}

/**
 * Represents a mutable chromosome, scaffold, or contig.
 * TODO import documentation
 */
public sealed interface MutableAssemblyUnit: AssemblyUnit, MutableGenomeChild {
    public override fun genes(): List<MutableGene>
}

/**
 * Provides implementation for [AssemblyUnit]. Internal to present a clean & safe public API.
 */
internal open class AssemblyUnitImpl protected constructor(
        protected open val delegate: FeatureImpl
) : Feature by delegate,
        AssemblyUnit {
    /**
     * INVARIANTS:
     * 1. type is CHROMOSOME, CONTIG, OR SCAFFOLD
     * 2. parent is a Genome
     */
    protected open fun invariants(): Boolean {
        if (type != CHROMOSOME && type != CONTIG && type != SCAFFOLD)
            throw IllegalStateException("An AssemblyUnit may only be of type CHROMOSOME, CONTIG, or SCAFFOLD.")
        parent //simply calling the property will ensure that it fails with a ClassCastException if the invariant is broken
        return true
    }

    internal constructor(
            seqid: String,
            source: String,
            type: FeatureType,
            start: Int,
            end: Int,
            score: Double?,
            strand: Char?,
            phase: Int?,
            attributes: Map<String, String>
    ): this(FeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes))

    public override val parent: Genome get() = delegate.parent as Genome
    public override val genome: Genome get() = parent

    init {
        assert(this.invariants())
    }

    public override fun copyTo(newParent: MutableGenome): MutableAssemblyUnit = delegate.copyTo(newParent) as MutableAssemblyUnit

    public override fun genes(): List<Gene> {
        TODO("Not yet implemented")
    }

    public override fun toString(): String = asRow()

    /**
     * Returns an immutable clone of the receiver ***without propagating the cloning to descendants or ancestors***.
     * This is to be used in the implementation of making an immutable clone of a whole genome.
     */
    internal fun soloImmutable(): AssemblyUnit {
        return AssemblyUnitImpl(
                seqid = seqid,
                source = source,
                type = type,
                start = start,
                end = end,
                score = score,
                strand = strand,
                phase = phase,
                attributes = attributes
        )
    }

    /**
     * Returns an immutable clone of the receiver ***without propagating the cloning to descendants or ancestors***.
     * This is to be used in the implementation of making an immutable clone of a whole genome.
     */
    internal fun soloMutable(): MutableAssemblyUnit {
        return MutableAssemblyUnitImpl(
                seqid = seqid,
                source = source,
                type = type,
                start = start,
                end = end,
                score = score,
                strand = strand,
                phase = phase,
                attributes = attributes
        )
    }
}

// The delegate hiding the supertype overriding is deliberate to reduce boilerplate in each delegator.
@Suppress("DELEGATED_MEMBER_HIDES_SUPERTYPE_OVERRIDE")
/**
 * Provides implementation for [MutableAssemblyUnit]. Internal to present clean & safe public API.
 */
internal class MutableAssemblyUnitImpl private constructor(
        override val delegate: MutableFeatureImpl
) : AssemblyUnitImpl(delegate),
        MutableFeature by delegate,
        MutableAssemblyUnit {

    /**
     * INVARIANTS (in addition to those of super):
     * 1. parent is a MutableGenome
     */
    protected override fun invariants(): Boolean {
        super.invariants()
        parent //simply calling the property will ensure that it fails with a ClassCastException if the invariant is broken
        return true
    }

    internal constructor(
            seqid: String,
            source: String,
            type: FeatureType,
            start: Int,
            end: Int,
            score: Double?,
            strand: Char?,
            phase: Int?,
            attributes: Map<String, String>
    ): this(MutableFeatureImpl(seqid, source, type, start, end, score, strand, phase, attributes))

    public override val parent: MutableGenome = delegate.parent as MutableGenome
    public override val genome: MutableGenome = parent

    //TODO safe type mutation

    public override fun moveTo(newParent: MutableGenome): MutableAssemblyUnit {
        TODO("Not yet implemented")
    }
    public override fun copyTo(newParent: MutableGenome): MutableAssemblyUnit = delegate.copyTo(newParent) as MutableAssemblyUnit
    public override fun clearSafeAttributes(): Unit = delegate.clearSafeAttributes()
    public override fun addAttribute(tag: String, value: String): Unit = delegate.addAttribute(tag, value)
    public override fun removeAttribute(tag: String): String? = delegate.removeAttribute(tag)
    public override fun genes(): List<MutableGene> {
        TODO()
    }
}