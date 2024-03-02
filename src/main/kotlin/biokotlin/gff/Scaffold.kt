package biokotlin.gff

class Scaffold(
    seqid: String,
    source: String,
    start: Int,
    end: Int,
    attributes: Map<String, String> = emptyMap(),
    children: List<Feature> = emptyList()
) : Feature(seqid, source, start, end, attributes = attributes, children = children) {

    override fun type(): FeatureType = FeatureType.Scaffold

}