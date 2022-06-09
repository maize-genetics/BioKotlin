package biokotlin.gff

class Intron(
    seqid: String,
    start: Int,
    end: Int
) : Feature(seqid, start, end) {

    override fun type(): FeatureType = FeatureType.Intron

}