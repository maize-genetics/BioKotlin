package biokotlin.genome

data class Position(val contig: String, val position: Int) : Comparable<Position> {

    override fun compareTo(other: Position): Int {
        return if (this.contig == other.contig) {
            this.position - other.position
        } else {
            try {
                this.contig.toInt() - other.contig.toInt()
            } catch (e: NumberFormatException) {
                // If we can't convert contigs to an int, then compare the strings
                this.contig.compareTo(other.contig)
            }
        }
    }

    override fun toString(): String {
        return "$contig:$position"
    }

}

data class PositionRange(val start: Position, val end: Position) : Comparable<PositionRange> {

    init {
        require(start.contig == end.contig) { "Start and end positions must be on the same contig." }
        require(start.position <= end.position) { "Start position must be less than or equal to end position." }
    }

    val contig = start.contig

    override fun compareTo(other: PositionRange): Int {
        return start.compareTo(other.start)
    }

    override fun toString(): String {
        return "${start.contig}:${start.position}-${end.position}"
    }

}