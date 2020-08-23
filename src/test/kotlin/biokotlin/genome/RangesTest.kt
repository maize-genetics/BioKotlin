package biokotlin.genome

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class RangesTest: StringSpec ({
    "Test GP flank range "{

        val chr = Chromosome("1")
        val gr = chr[48..68]
        val gr2= gr.enlarge(20)

        var flank5 = gr.flank(5, Direction.BOTH, 100)
        flank5.size shouldBe 2
        //flank5.contains(SeqRange.closed(Int1(43u), Int1(47u))) shouldBe true
        flank5.contains(chr[43..47]) shouldBe true

        // run again using other signature
        var flank5left = gr.flankLeft(5,100)
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
        var ranges : MutableSet<GenomePosRange> = mutableSetOf()

        val chr = Chromosome("1")
        ranges.add(chr[5..25])
        ranges.add(chr[48..65])
        ranges.add(chr[75..87])
        ranges.add(chr[170..175])
        ranges.add(chr[180..210])
        ranges.add(chr[300..320])

        var flankedRangeSet = flankGenomePosRangeSet ( 5,  Direction.BOTH, 500, ranges)
        flankedRangeSet.size shouldBe 12

    }

    "Test shift range " {
        val chr = Chromosome("1")
        var range1 = chr[48..68]
        var shift5 = range1.shift(5  )
        shift5.range.lowerEndpoint().site shouldBe 53
        shift5.range.upperEndpoint().site shouldBe 73

        var shift50 = range1.shift(50,100) // 100 is max size for range
        shift50.range.lowerEndpoint().site shouldBe 98
        shift50.range.upperEndpoint().site shouldBe 100

        var shiftNeg10 = range1.shift(-10,100)
        shiftNeg10.range.lowerEndpoint().site shouldBe 38
        shiftNeg10.range.upperEndpoint().site shouldBe 58

        // Shift left beyond 1
        var shiftNeg50 = range1.shift(-50,100)
        shiftNeg50.range.lowerEndpoint().site shouldBe 1
        shiftNeg50.range.upperEndpoint().site shouldBe 18

    }

    "Test Shift on Set<SeqRange> " {
        var ranges : MutableSet<GenomePosRange> = mutableSetOf()

        val chr = Chromosome("1")
        ranges.add(chr[5..25])
        ranges.add(chr[48..65])
        ranges.add(chr[75..87])
        ranges.add(chr[170..175])
        ranges.add(chr[180..210])
        ranges.add(chr[300..320])

        // Test positive shift - shift right
        var shiftedRangeSet = shiftGenomePosRangeSet ( 5,    ranges,500)

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
        shiftedRangeSet = shiftGenomePosRangeSet ( -5,   ranges, 500)
        shiftedRangeSet.contains(chr[1..20]) shouldBe true
        shiftedRangeSet.contains(chr[43..60]) shouldBe true
        shiftedRangeSet.contains(chr[70..82]) shouldBe true
        shiftedRangeSet.contains(chr[165..170]) shouldBe true
        shiftedRangeSet.contains(chr[175..205]) shouldBe true
        shiftedRangeSet.contains(chr[295..315]) shouldBe true
    }

    "Test Merge Ranges " {
        var ranges : MutableSet<GenomePosRange> = mutableSetOf()

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

        var mergedRangeSet = mergeGenomePosRangeSet ( 5, ranges)

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
        var ranges : MutableSet<GenomePosRange> = mutableSetOf()

        // THis is not sorting the way I would expect.
        // ANd the test only passes because the ranges to merge
        // are next to each other.  If I have chrom 2 interspersed
        // differently I expect this will fail.  Sorting needs work
        val chr = Chromosome("1")
        val chr2 = Chromosome("2")
        ranges.add(chr[48..65])
        ranges.add(chr[5..25])
        ranges.add(chr2[27..30])
        ranges.add(chr[180..210])
        ranges.add(chr[170..175])
        ranges.add(chr[300..320])
        ranges.add(chr[69..75])

        println("Unsorted , ummerged ranges:  ")
        for (range in ranges) {
            println(range.toString())
        }

        var mergedRangeSet = mergeGenomePosRangeSet ( 5, ranges)

        println("\nSORTED , MERGED ranges:  ")
        for (range in mergedRangeSet) {
            println(range.toString())
        }

        mergedRangeSet.size shouldBe 5

        println("\nSorted, merged ranges:")
        for (range in mergedRangeSet) {
            println(range.toString())
        }

        // 5..25 and 27..30 should NOT merge as they are on different chromosoems
        mergedRangeSet.contains(chr[5..30]) shouldBe false
        mergedRangeSet.contains(chr[48..75]) shouldBe true
        mergedRangeSet.contains(chr[170..210]) shouldBe true

    }

})