package biokotlin.genome

import biokotlin.seq.NucSeq
import biokotlin.seq.NucSeqRecord
import com.google.common.collect.HashMultimap
import com.google.common.collect.Range
import com.google.common.collect.SortedSetMultimap
import com.google.common.collect.TreeMultimap
import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe
import java.util.*


class GenomicRangesTest: StringSpec({
    "Test GP flank range "{

        val chr = Chromosome("1")
        val gr = chr[48..68]
        val gr2 = gr.enlarge(20)

        var flank5 = gr.flank(5, Direction.BOTH, 100)
        flank5.size shouldBe 2
        //flank5.contains(SeqRange.closed(Int1(43u), Int1(47u))) shouldBe true
        flank5.contains(chr[43..47]) shouldBe true

        // run again using other signature
        var flank5left = gr.flankLeft(5, 100)
        //flank5.contains(SeqRange.closed(Int1(43u), Int1(47u))) shouldBe true
        flank5left.contains(chr[43..47]) shouldBe true

        var flank5right = gr.flankRight(5)
        flank5right.contains(chr[69..73]) shouldBe true

        // Test right side is cropped
        var flank40 = gr.flank(40, Direction.BOTH, 100)
        flank40.size shouldBe 2
        println("\nRanges for flank40 of Range 48,68")
        for (range in flank40) {
            println(range.toString())
        }
        flank40.contains(chr[8..47]) shouldBe true

        // Test left side gets cropped
        var range2 = chr[5..58]
        flank40 = range2.flank(40, Direction.BOTH, 100)
        flank40.size shouldBe 2
        println("\nRanges for flank40 of Range 5,58")
        for (range in flank40) {
            println(range.toString())
        }
        flank40.contains(chr[1..4]) shouldBe true
        flank40.contains(chr[59..98]) shouldBe true

    }

    "Test flank on Set<SeqRange> " {
        var ranges: MutableSet<GenomePosRange> = mutableSetOf()

        val chr = Chromosome("1")
        ranges.add(chr[5..25])
        ranges.add(chr[48..65])
        ranges.add(chr[75..87])
        ranges.add(chr[170..175])
        ranges.add(chr[180..210])
        ranges.add(chr[300..320])

        var flankedRangeSet = flankGenomePosRangeSet(5, Direction.BOTH, 500, ranges)
        flankedRangeSet.size shouldBe 12

    }

    "Test shift range " {
        val chr = Chromosome("1")
        var range1 = chr[48..68]
        var shift5 = range1.shift(5)
        shift5.range.lowerEndpoint().site shouldBe 53
        shift5.range.upperEndpoint().site shouldBe 73

        var shift50 = range1.shift(50, 100) // 100 is max size for range
        shift50.range.lowerEndpoint().site shouldBe 98
        shift50.range.upperEndpoint().site shouldBe 100

        var shiftNeg10 = range1.shift(-10, 100)
        shiftNeg10.range.lowerEndpoint().site shouldBe 38
        shiftNeg10.range.upperEndpoint().site shouldBe 58

        // Shift left beyond 1
        var shiftNeg50 = range1.shift(-50, 100)
        shiftNeg50.range.lowerEndpoint().site shouldBe 1
        shiftNeg50.range.upperEndpoint().site shouldBe 18

    }

    "Test Shift on Set<SeqRange> " {
        var ranges: MutableSet<GenomePosRange> = mutableSetOf()

        val chr = Chromosome("1")
        ranges.add(chr[5..25])
        ranges.add(chr[48..65])
        ranges.add(chr[75..87])
        ranges.add(chr[170..175])
        ranges.add(chr[180..210])
        ranges.add(chr[300..320])

        // Test positive shift - shift right
        var shiftedRangeSet = shiftGenomePosRangeSet(5, ranges, 500)

        shiftedRangeSet.size shouldBe 6

        println("\nRanges for ShiftedRangeSet ,shift by 5")
        for (range in shiftedRangeSet) {
            println(range.toString())
        }
        shiftedRangeSet.contains(chr[10..30]) shouldBe true
        shiftedRangeSet.contains(chr[53..70]) shouldBe true
        shiftedRangeSet.contains(chr[80..92]) shouldBe true
        shiftedRangeSet.contains(chr[175..180]) shouldBe true
        shiftedRangeSet.contains(chr[185..215]) shouldBe true
        shiftedRangeSet.contains(chr[305..325]) shouldBe true

        // Test negative shift (shift left)
        shiftedRangeSet = shiftGenomePosRangeSet(-5, ranges, 500)
        shiftedRangeSet.contains(chr[1..20]) shouldBe true
        shiftedRangeSet.contains(chr[43..60]) shouldBe true
        shiftedRangeSet.contains(chr[70..82]) shouldBe true
        shiftedRangeSet.contains(chr[165..170]) shouldBe true
        shiftedRangeSet.contains(chr[175..205]) shouldBe true
        shiftedRangeSet.contains(chr[295..315]) shouldBe true
    }

    "Test Merge Ranges " {
        var ranges: MutableSet<GenomePosRange> = mutableSetOf()

        val chr = Chromosome("1")
        ranges.add(chr[48..65])
        ranges.add(chr[5..25])
        ranges.add(chr[27..30])
        ranges.add(chr[180..210])
        ranges.add(chr[170..175])
        ranges.add(chr[300..320])
        ranges.add(chr[69..75])

        println("Unsorted , ummerged ranges:  ")
        for (range in ranges) {
            println(range.toString())
        }

        var mergedRangeSet = mergeGenomePosRangeSet(5, ranges)

        println("\nSORTED , MERGED ranges:  ")
        for (range in mergedRangeSet) {
            println(range.toString())
        }

        mergedRangeSet.size shouldBe 4

        println("\nSorted, merged ranges:")
        for (range in mergedRangeSet) {
            println(range.toString())
        }

        mergedRangeSet.contains(chr[5..30]) shouldBe true
        mergedRangeSet.contains(chr[48..75]) shouldBe true
        mergedRangeSet.contains(chr[170..210]) shouldBe true

    }

    "Test Merge Ranges Multiple Chromosomes" {
        var ranges: MutableSet<GenomePosRange> = mutableSetOf()

        // Sorting now fixed.
        val chr1 = Chromosome("1")
        val chr2 = Chromosome("2")
        val chr3 = Chromosome("3")
        val chr4 = Chromosome("4")

        ranges.add(chr1[48..65])
        ranges.add(chr1[5..25])
        ranges.add(chr2[27..30])
        ranges.add(chr3[180..210])
        ranges.add(chr4[300..320])
        ranges.add(chr3[170..175])
        ranges.add(chr4[69..75])


        println("Unsorted , ummerged ranges:  ")
        for (range in ranges) {
            println(range.toString())
        }

        var mergedRangeSet = mergeGenomePosRangeSet(5, ranges)

        mergedRangeSet.size shouldBe 6

        println("\nSorted, merged ranges:")
        for (range in mergedRangeSet) {
            println(range.toString())
        }

        // 5..25 and 27..30 should NOT merge as they are on different chromosoems
        mergedRangeSet.contains(chr1[5..30]) shouldBe false
        mergedRangeSet.contains(chr1[48..75]) shouldBe false
        mergedRangeSet.contains(chr3[170..210]) shouldBe true

    }

    "Test Chromosome and GenomePosRange Sort" {
        var chromList = mutableListOf<Chromosome>()
        var chr1 = Chromosome("1")
        var chr2 = Chromosome("2")
        var chr3 = Chromosome("3")
        var chr4 = Chromosome("4")

        chromList.add(chr2)
        chromList.add(chr1)
        chromList.add(chr4)
        chromList.add(chr3)

        println("\nUnsorted chrom list:")
        for (chrom in chromList) {
            println(chrom.toString())
        }

        var sortedChromList = chromList.sorted()
        println("\nSorted chrom list - it works:")
        for (chrom in sortedChromList) {
            println(chrom.toString())
        }

        var ranges: MutableSet<GenomePosRange> = mutableSetOf()

        ranges.add(chr1[48..65])
        ranges.add(chr1[5..25])
        ranges.add(chr2[27..30])
        ranges.add(chr3[180..210])
        ranges.add(chr4[300..320])
        ranges.add(chr3[170..175])
        ranges.add(chr4[69..75])

        // Try sorting them.  It should merge the 3 chr3 entries,
        // Nothing else should merge
        println("\nUnsorted ranges Set:")
        for (range in ranges) {
            println(range.toString())
        }

        // var sortedRanges = ranges.toSortedSet()
        // override fun compareTo(other: GenomePosRange): Int = compareValuesBy(this,other,
        //            {this.range.lowerEndpoint().chromosome.name},{it.range.lowerEndpoint().site})

        var sortedRanges = ranges.sorted()
        // THis one fails - it uses the comparator in GenomePosRange
        println("\nSorted ranges S uses class comparator - fails:")
        for (range in sortedRanges) {
            println(range.toString())
        }

        // Add this sort to mergeGenomePosRangeSet()
        val sortedSet2 = ranges.sortedWith(compareBy({ it.range.lowerEndpoint().chromosome.name }, { it.range.lowerEndpoint().site }))
        println("\nSorted ranges Set 2 !! using explicit compareBy - works!:")
        for (range in sortedSet2) {
            println(range.toString())
        }

        // Test sorting the list:
        // THis works
        var rangeList: MutableList<GenomePosRange> = mutableListOf()
        rangeList.add(chr1[48..65])
        rangeList.add(chr1[5..25])
        rangeList.add(chr2[27..30])
        rangeList.add(chr3[180..210])
        rangeList.add(chr4[300..320])
        rangeList.add(chr3[170..175])
        rangeList.add(chr4[69..75])

        println("\nUnsorted ranges LIST:")
        for (range in rangeList) {
            println(range.toString())
        }

        val sortedList = rangeList.sorted()
        println("\nSorted ranges LIST defaulting to class comparator - fails")
        for (range in sortedList) {
            println(range.toString())
        }

        val sortedList2 = rangeList.sortedWith(compareBy({ it.range.lowerEndpoint().chromosome.name }, { it.range.lowerEndpoint().site }))
        println("\nSorted ranges LIST using explicit compareBy:")
        for (range in sortedList2) {
            println(range.toString())
        }
    }

    "Test GenomePosition sorting " {
        var ranges: MutableSet<GenomePosition> = mutableSetOf()
        val chr1 = Chromosome("1")
        val chr2 = Chromosome("2")
        val gp1 = GenomePosition(chr1, 10)

        ranges.add(GenomePosition(chr2, 10))
        ranges.add(GenomePosition(chr1, 10))
        ranges.add(GenomePosition(chr1, 5))
        ranges.add(GenomePosition(chr2, 5))
        ranges.add(GenomePosition(chr1, 7))

        println("\nUnsorted GenomePosition Range")
        for (range in ranges) {
            println(range.toString())
        }

        println("\nSorted GenomePosition Ranges")
        val sortedRanges = ranges.sorted()
        for (range in sortedRanges) {
            println(range.toString())
        }
    }

    "Test navigatble set" {

        // THis works - but must have comparator defined
        // If switch to SeqPosition, then need that comparator
        val irComparator = IRComparator()
        val ns: NavigableSet<IntRange> = TreeSet(irComparator)

        ns.add(43..58)
        ns.add(7..12)
        ns.add(82..145)
        ns.add(57..78)

        // get reverse order
        var reverseNs = ns.descendingSet()

        // Print the normal and reverse views
        println("Normal order: " + ns);
        println("Reverse order: " + reverseNs);

        // Try now with Range of SeqPosition
        //LCJ _ these need to be changed to use SeqRecord
//        var chr1 = SeqContig("1")
//        var chr2 = SeqContig("2")
//        val SeqRange = Range.closed(SeqPosition(chr1,43),SeqPosition(chr1,58))
//        val sComparator = SeqPositionRangeComparator()
//        val sR: NavigableSet<Range<SeqPosition>> = TreeSet(sComparator)
//
//        sR.add(SeqRange)
//        sR.add(Range.closed(SeqPosition(chr2,7),SeqPosition(chr2,12)))
//        sR.add(Range.closed(SeqPosition(chr1,82),SeqPosition(chr1,156)))
//        sR.add(Range.closed(SeqPosition(chr1,73),SeqPosition(chr1,78)))
//
//        var reversesR = sR.descendingSet()
//        // Print the normal and reverse views
//        println("Normal order SeqPosition Range: " + sR);
//        println("Reverse order SeqPosition Range: " + reversesR);
//

        // Try with Range<Int>
        val gRange = Range.closed(43, 58)
        val rComparator = RComparator()
        val nsR: NavigableSet<Range<Int>> = TreeSet(rComparator)

        nsR.add(gRange)
        nsR.add(Range.closed(7, 12))
        nsR.add(Range.closed(82, 156))
        nsR.add(Range.closed(73, 78))

        var reverseNsR = nsR.descendingSet()
        // Print the normal and reverse views
        println("Normal order nsR: " + ns);
        println("Reverse order nsR: " + reverseNs);

    }
    "Test NS of SeqPosRange" {

        val sR: NavigableSet<SeqPosRange> = TreeSet()
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        var spr1 = SeqPosRange(record1, 43..58)
        var spr2 = SeqPosRange(record2, 7..12)
        var spr3 = SeqPosRange(record1, 82..156)
        var spr4 = SeqPosRange(record1, 73..78)

        sR.add(spr1)
        sR.add(spr2)
        sR.add(spr3)
        sR.add(spr4)

        var reversesR = sR.descendingSet()
        // Print the normal and reverse views
        println("Normal order SeqPosition Range: " + sR);
        println("\nReverse order SeqPosition Range: " + reversesR);
    }
    "Test NS of short seq len SeqPosRange" {

        val sR: NavigableSet<SeqPosRange> = TreeSet()
        val dnaString = "ACGTGGTGAATATATATGCGCTCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGCGTATCAGTCAGTCATGCATGCATGTGTGTACACACATGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        var spr1 = SeqPosRange(record1, 43..50)
        var spr2 = SeqPosRange(record2, 7..12)
        var spr3 = SeqPosRange(record1, 32..45)
        var spr4 = SeqPosRange(record1, 22..26)

        sR.add(spr1)
        sR.add(spr2)
        sR.add(spr3)
        sR.add(spr4)

        var reversesR = sR.descendingSet()
        // Print the normal and reverse views
        println("Normal order SeqPosition Range: " + sR);
        println("\nReverse order SeqPosition Range: " + reversesR);
    }
    "Test SeqPosRange in a Map " {
        val sprMap = TreeMap<SeqPosRange, String>()
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        var spr1 = SeqPosRange(record1, 43..58)
        var spr2 = SeqPosRange(record2, 7..12)
        var spr3 = SeqPosRange(record1, 82..156)
        var spr4 = SeqPosRange(record1, 73..78)

        sprMap.put(spr1, "forward")
        sprMap.put(spr2, "forward")
        sprMap.put(spr3, "reverse")
        sprMap.put(spr4, "reverse")

        // I don't get what we want to test here.  Show me a use case
    }
    "Test SeqPosRange in a multi Map " {
        // TreeSets will fail if you give it anything that doesn't have a comparator
        //val sprMap: SortedSetMultimap<String, Int> = TreeMultimap.create()
       // val sprMap: HashMultimap<NucSeqRecord, IntRange> = HashMultimap.create()
        val sprMap = TreeMap<SeqPosRange, String>()
        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        var spr1 = SeqPosRange(record1, 43..58)
        var spr2 = SeqPosRange(record2, 7..12)
        var spr3 = SeqPosRange(record1, 82..156)
        var spr4 = SeqPosRange(record1, 73..78)

        sprMap.put(spr1, "forward")
        sprMap.put(spr2, "forward")
        sprMap.put(spr3, "reverse")
        sprMap.put(spr4, "reverse")

        // I don't get what we want to test here.  Show me a use case
    }

    "Test SequenceRecord in a Map " {
        val myTreeMultimap: SortedSetMultimap<String, String> = TreeMultimap.create()
        // Tree maps fail if no comparator  - neither NucSeqRecord nor IntRange has comparator
        //val sprMap: SortedSetMultimap<NucSeqRecord, IntRange> = TreeMultimap.create()
        val sprMap: HashMultimap<NucSeqRecord, IntRange> = HashMultimap.create()

        val dnaString = "ACGTGGTGAATATATATGCGCGCGTGCGTGGATCAGTCAGTCATGCATGCATGTGTGTACACACATGTGATCGTAGCTAGCTAGCTGACTGACTAGCTGAC"
        val dnaString2 = "ACGTGGTGAATATATATGCGCGCGTGCGTGGACGTACGTACGTACGTATCAGTCAGCTGAC"
        val record1 = NucSeqRecord(NucSeq(dnaString), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))
        val record2 = NucSeqRecord(NucSeq(dnaString2), "Sequence 1", description = "The first sequence",
                annotations = mapOf("key1" to "value1"))

        sprMap.put(record1, 43..58)
        sprMap.put(record2, 7..12)

    }

})

