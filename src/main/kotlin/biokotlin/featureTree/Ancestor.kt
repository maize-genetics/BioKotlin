package biokotlin.featureTree

/**
 * Any class which can have instances of [Feature] as its children. Provides the useful [flatten]
 * function, which allows the [biokotlin.featureTree] framework to easily work with the Kotlin collections
 * framework. Classes using this interface are highlighted in yellow:
 *
 * <img src="feature_tree/ancestors.svg" style="display:block; margin-left:auto; margin-right:auto;">
 */
sealed interface Ancestor {
    /**
     * The immediate children of this [Ancestor]. Prefer using more specific functions to access the children
     * of an [Ancestor], such as [Gene.transcript], [Transcript.children], etc. as these provide you with access
     * to more functionality.
     */
    fun children(): List<Feature>

    /**
     * Returns the descendants of this ancestor in depth-first order (where children on the same
     * level are sorted by [Feature.compareTo]).
     *
     * <img src="feature_tree/flatten.svg" style="display:block; margin-left:auto; margin-right:auto;" width="50%" alt="A diagram of a tree being collapsed into a single list in depth-first order by the flatten method.">
     *
     *
     * This function provides a convenient link to the Kotlin collections framework. For example,
     * use in conjunction with [kotlin.sequences.associateBy], [kotlin.sequences.groupBy] to allow
     * for access to elements of your feature tree based on specific properties:
     * ```kotlin
     * val genome = Genome.fromGFF("file.gff")
     * // Let's say that your GFF file has a custom attribute called myID and you want to access
     * // your features by this attribute
     * val associatedByMyID = genome.flatten().associateBy { it.attributes["myID"] }
     * // Access the feature with 00001 in the myID attribute
     * val myFeature = associatedByMyID["00001"]
     * ```
     *
     */
    fun flatten(): List<Feature> {
        val list = mutableListOf<Feature>()

        for (child in children()) {
            list.add(child)
            if (child is Ancestor) list.addAll(child.flatten())
        }

        return list
    }

    /**
     * Returns all descendant features whose start and end are within [start] and [end], inclusive, 1-based.
     *
     * All returned features must have a seqid equivalent to [seqid]. If left null, then the function
     * will use the seqid of the receiving [Ancestor].
     *
     * If [direct] is set to true, then only direct children will be returned, instead of all descendants.
     *
     * The resulting list will be in depth-first order, where features on the same
     * level are sorted left-to-right by [Feature.compareTo].
     *
     * @throws [IllegalArgumentException] if [seqid] is null and the receiving [Ancestor] is not a [Feature]
     * or if start is greater than end.
     * @throws [IndexOutOfBoundsException] if start or end are non-positive.
     */
    fun within(start: Int, end: Int, seqid: String? = null, direct: Boolean = false): List<Feature> {
        val correctedSeqid = seqid
            ?: if (this is Feature) this.seqid
            else throw IllegalArgumentException("Instances of ${this::class} cannot receive calls of within without specifying the seqid.")

        if (start > end) throw IllegalArgumentException("Start must be less than or equal to end for function within.")

        if (start <= 0) throw IndexOutOfBoundsException("Start and end ranges must be positive integers. Indices start at 1.")

        return if (direct) {
            children().filter { it.start >= start && it.end <= end }
        } else {
            flatten().filter { it.start >= start && it.end <= end && it.seqid == correctedSeqid }
        }
    }

    /**
     * Returns a representation of this [Ancestor] and its descendants in the
     * [dot](https://en.wikipedia.org/wiki/DOT_(graph_description_language)) language. Use with
     * [graphviz](https://en.wikipedia.org/wiki/Graphviz) to produce an image output.
     *
     * Example image output for a gene after being processed through graphviz:
     * <img src="feature_tree/visualize.svg" style="display:block; margin-left:auto; margin-right:auto;" width="70%" alt="The graphical representation of a feature tree.">
     */
    fun visualize(): String {
        val sb = StringBuilder()
        sb.append("digraph {\n")
        sb.append("rank = source\n")
        sb.append("ordering = out\n")
        sb.append("node[shape = box style = filled colorscheme = set312]\n")

        //Adds top-level label
        if (this !is Feature) {
            sb.append("\"${hashCode()}\" [label=GENOME color = gray]\n")

        } else {
            sb.append("\"${hashCode()}\" ")
            sb.append("[label=\"${type().name}\\n${start}-${end}\\n${attributes["ID"] ?: ""}\" ")
            sb.append("color = ${type().ordinal + 1}]\n")
        }

        //Adds top-level children
        for (child in children()) {
            sb.append("\"${hashCode()}\" -> ${child.hashCode()}")
        }

        for (descendant in flatten()) {
            //Label
            sb.append("\"${descendant.hashCode()}\" ")
            sb.append("[label=\"${descendant.type().name}\\n${descendant.start}-${descendant.end}\\n${descendant.attributes["ID"] ?: ""}\" ")
            sb.append("color = ${descendant.type().ordinal + 1}]\n")

            //Line to parent
            val parent = when (descendant) {
                is GenomeChild -> descendant.genome()
                is Transcript -> descendant.gene()
                is TranscriptChild -> descendant.transcript()
                else -> throw Exception() //unusual, the IDE doesn't yell at me but the compiler does. The above
                //should be exhaustive because Feature is sealed, but the compiler seems to not recognize this?
            }
            sb.append("\"${descendant.hashCode()}\" -> \"${parent.hashCode()}\"\n")

            //Line to children
            if (descendant is Ancestor) {
                for (child in descendant.children()) {
                    sb.append("\"${descendant.hashCode()}\" -> \"${child.hashCode()}\"\n")
                }
            }
        }

        sb.append("}")

        return sb.toString()
    }

    /**
     * Represents this [Ancestor] and all its descendants as they would appear in a GFF file.
     *
     * For example, calling this on a [Gene] may return
     * ```plaintext
    scaf_100	NAM	gene	10647	11026	.	+	.	ID=Zm00001eb437070;biotype=protein_coding;logic_name=cshl_gene;
    scaf_100	NAM	mRNA	10647	11026	.	+	.	ID=Zm00001eb437070_T001;Parent=Zm00001eb437070;biotype=protein_coding;transcript_id=Zm00001eb437070_T001;canonical_transcript=1;
    scaf_100	NAM	exon	10647	11026	.	+	.	Parent=Zm00001eb437070_T001;Name=Zm00001eb437070_T001.exon.1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001eb437070_T001.exon.1;rank=1;
    scaf_100	NAM	five_prime_UTR	10647	10845	.	+	.	Parent=Zm00001eb437070_T001;
    scaf_100	NAM	CDS	10846	11025	.	+	0	ID=Zm00001eb437070_P001;Parent=Zm00001eb437070_T001;protein_id=Zm00001eb437070_P001;
    scaf_100	NAM	three_prime_UTR	11026	11026	.	+	.	Parent=Zm00001eb437070_T001;
     * ```
     * while calling it on a [Genome] would return the entire content of a GFF file.
     * @see [Genome.toString]
     */
    fun descendantsToString(): String {
        val sb = StringBuilder()
        if (this is Feature) sb.append(this.toString())
        for (descendant in flatten()) sb.append(descendant.toString())
        return sb.toString()
    }
}