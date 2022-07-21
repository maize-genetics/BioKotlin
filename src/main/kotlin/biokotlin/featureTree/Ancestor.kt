package biokotlin.featureTree

/**
 * Any class which can have instances of [Feature] as its children. Provides the useful [flatten]
 * function, which allows the [biokotlin.featureTree] framework to easily work with the Kotlin collections
 * framework.
 */
sealed interface Ancestor {
    /**
     * The immediate children of this [Ancestor]. Prefer using the children property of the
     * specific instance of [Ancestor] instead as this will provide them as a specific subtype of
     * [Feature], allowing access to more functionality.
     */
    fun children(): List<Feature>

    /**
     * Returns the descendants of this ancestor in depth-first order (where children on the same
     * level are sorted left-right by [FeatureComparator]. If the receiver is a [Feature], then
     * it will be included in the list, otherwise it will be omitted.
     * This function provides a convenient link to the Kotlin collections framework.
     *
     * Use in conjunction with [kotlin.sequences.associateBy], [kotlin.sequences.groupBy] to allow
     * for access to elements of your feature tree based on specific properties.
     *
     * TODO include example of this
     */
    fun flatten(): List<Feature> {
        val list = mutableListOf<Feature>()
        if (this is Feature) list.add(this)

        for (child in children()) {
            if (child is Ancestor) list.addAll(child.flatten())
            else list.add(child)
        }

        return list
    }

    /**
     * Returns all descendant features whose start and end are within [start] and [end], inclusive, 1-based.
     *
     * All returned features must have a seqid equivalent to [seqid]. If left null, then the function
     * will use the seqid of the receiving [Ancestor]. Throws [UnsupportedOperationException] if a null
     * parameter is used on an [Ancestor] that is not a [Feature].
     *
     * If [direct] is set to true, then only direct children will be returned, instead of all descendants.
     *
     * The resulting list will be in depth-first left-to-right order, where features on the same
     * level are sorted left-to-right by [FeatureComparator].
     */
    fun within(start: Int, end: Int, seqid: String? = null, direct: Boolean = false): List<Feature> {
        val correctedSeqid = if (seqid == null) {
            if (this is Feature) this.seqid
            else throw UnsupportedOperationException("Instances of ${this::class} cannot receive calls of within without specifying the seqid.")
        } else {
            seqid
        }

        return if (direct) {
            children().filter { it.start >= start && it.end <= end }
        } else {
            flatten().filter { it.start >= start && it.end <= end && it.seqid == correctedSeqid }
        }
    }

    /**
     * Returns a representation of this [Ancestor] and its descendants in the
     * [dot](https://en.wikipedia.org/wiki/DOT_(graph_description_language)) language. Use with
     * [graphviz](https://en.wikipedia.org/wiki/Graphviz) to produce an image output. Primarily useful
     * for debugging.
     *
     * Be wary of visualizing large trees as this will make the resulting image difficult to parse
     * visually and slow for graphviz to produce.
     */
    fun visualize(): String {
        val sb = StringBuilder()
        sb.append("digraph {\n")
        sb.append("rank = source\n")
        sb.append("ordering = out\n")
        sb.append("node[shape = box style = filled colorscheme = set312]\n")

        //Adds genome label
        if (this !is Feature) {
            sb.append("\"${hashCode()}\" [label=GENOME color = gray]\n")
            for (child in children()) {
                sb.append("\"${hashCode()}\" -> ${child.hashCode()}")
            }
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
     * If this [Ancestor] is an instance of [Feature], its representation as a row in a GFF
     * file will be the first line of this string. Otherwise, it is omitted. The rest of the string
     * is the representation of all descendants of this [Ancestor] in depth-first left-to-right order
     * where the features on the same level are ordered left-to-right by [FeatureComparator] as rows
     * in a GFF file.
     */
    fun descendantsToString(): String {
        val sb = StringBuilder()
        for (descendant in flatten()) sb.append(descendant.toString())
        return sb.toString()
    }
}