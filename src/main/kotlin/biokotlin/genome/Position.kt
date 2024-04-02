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
