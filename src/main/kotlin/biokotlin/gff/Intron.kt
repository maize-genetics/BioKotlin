package biokotlin.gff

class Intron(
    seqid: String,
    source: String,
    start: Int,
    end: Int
) : Feature(seqid, source, start, end) {

    override fun type(): FeatureType = FeatureType.Intron

}