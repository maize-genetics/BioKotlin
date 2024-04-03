package biokotlin.genome

data class Position(val contig: String, val position: Int) : Comparable<Position> {

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
        require(start <= end) { "Start position must be less than or equal to end position." }
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

    override fun toString(): String {
        return "${contig}:${start}-${end}"
    }

}