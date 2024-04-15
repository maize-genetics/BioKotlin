package biokotlin.genome

/**
 * Data class to represent a position on a contig.
 * Position is 1-based
 */
data class Position(val contig: String, val position: Int) : Comparable<Position> {

    init {
        require(position >= 1) { "Position must be greater than or equal to 1. Position: $position" }
    }

    override fun compareTo(other: Position): Int {
        return if (this.contig == other.contig) {
            position - other.position
        } else {
            try {
                contig.toInt() - other.contig.toInt()
            } catch (e: NumberFormatException) {
                // If we can't convert contigs to an int, then compare the strings
                contig.compareTo(other.contig)
            }
        }
    }

    override fun toString(): String {
        return "$contig:$position"
    }

}

data class PositionRange(val contig: String, val start: Int, val end: Int) : Comparable<PositionRange> {

    constructor(start: Position, end: Position) : this(start.contig, start.position, end.position) {
        require(start.contig == end.contig) { "Start and end positions must be on the same contig." }
    }

    init {
        require(start >= 1) { "Start position must be greater than or equal to 1. Start: $start" }
        require(end >= 1) { "End position must be greater than or equal to 1. End: $end" }
        require(start <= end) { "Start position must be less than or equal to end position. Start: $start End: $end" }
    }

    val startPosition by lazy { Position(contig, start) }

    val endPosition by lazy { Position(contig, end) }

    fun length(): Int {
        return end - start + 1
    }

    override fun compareTo(other: PositionRange): Int {
        return if (this.contig == other.contig) {
            start - other.start
        } else {
            try {
                contig.toInt() - other.contig.toInt()
            } catch (e: NumberFormatException) {
                // If we can't convert contigs to an int, then compare the strings
                contig.compareTo(other.contig)
            }
        }
    }

    /**
     * Compare this position range to a position.
     *
     * @return the value 0 if the position is within the range on the same contig;
     * a value less than 0 if this position range is less than the position;
     * and a value greater than 0 if this position range is greater than the position.
     *
     */
    fun compareTo(position: Position): Int {
        return if (this.contig == position.contig) {
            when (position.position) {
                in start..end -> 0
                else -> start - position.position
            }
        } else {
            try {
                contig.toInt() - position.contig.toInt()
            } catch (e: NumberFormatException) {
                // If we can't convert contigs to an int, then compare the strings
                contig.compareTo(position.contig)
            }
        }
    }

    fun contains(position: Position): Boolean {
        return contig == position.contig && start <= position.position && position.position <= end
    }

    override fun toString(): String {
        return "${contig}:${start}-${end}"
    }

}