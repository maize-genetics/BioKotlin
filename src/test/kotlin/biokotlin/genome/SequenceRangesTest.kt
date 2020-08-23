package biokotlin.genome

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class SequenceRangesTest: StringSpec({
    "Test flank range "{

        var range1 = (48..68).toSeqRange()
        var flank5 = range1.flank(5, Direction.BOTH, 100)
        flank5.size shouldBe 2
        //flank5.contains(SeqRange.closed(Int1(43u), Int1(47u))) shouldBe true
        flank5.contains((43..47).toSeqRange()) shouldBe true

        // run again using other signature
        var flank5left = range1.flankLeft(5,100)
        //flank5.contains(SeqRange.closed(Int1(43u), Int1(47u))) shouldBe true
        flank5left.contains((43..47).toSeqRange()) shouldBe true

        var flank5right = range1.flankRight(5)
        flank5right.contains((69..73).toSeqRange())

        // Test right side is cropped
        var flank40 = range1.flank(40, Direction.BOTH, 100)
        flank40.size shouldBe 2
        println("\nRanges for flank40 of Range 48,68")
        for (range in flank40) {
            println(range.toString())
        }
        flank40.contains((8..47).toSeqRange()) shouldBe true



        // Test left side gets cropped
        var range2 = (5..58).toSeqRange()
        flank40 = range2.flank(40, Direction.BOTH, 100)
        flank40.size shouldBe 2
        println("\nRanges for flank40 of Range 5,58")
        for (range in flank40) {
            println(range.toString())
        }
        flank40.contains((1..4).toSeqRange()) shouldBe true
        flank40.contains((59..98).toSeqRange()) shouldBe true

    }

    "Test flank on Set<SeqRange> " {
        var ranges : MutableSet<SeqRange> = mutableSetOf()

        ranges.add((5..25).toSeqRange())
        ranges.add((48..65).toSeqRange())
        ranges.add((75..87).toSeqRange())
        ranges.add((170..175).toSeqRange())
        ranges.add((180..210).toSeqRange())
        ranges.add((300..320).toSeqRange())

        var flankedRangeSet = flankSeqRangeSet ( 5,  Direction.BOTH, 500, ranges)
        flankedRangeSet.size shouldBe 12

    }

    "Test shift range " {
        var range1 = (48..68).toSeqRange()
        var shift5 = range1.shift(5  )
        shift5.range.lowerEndpoint() shouldBe 53
        shift5.range.upperEndpoint() shouldBe 73

        var shift50 = range1.shift(50,100) // 100 is max size for range
        shift50.range.lowerEndpoint() shouldBe 98
        shift50.range.upperEndpoint() shouldBe 100

        var shiftNeg10 = range1.shift(-10,100)
        shiftNeg10.range.lowerEndpoint() shouldBe 38
        shiftNeg10.range.upperEndpoint() shouldBe 58

        // Shift left beyond 1
        var shiftNeg50 = range1.shift(-50,100)
        shiftNeg50.range.lowerEndpoint() shouldBe 1
        shiftNeg50.range.upperEndpoint() shouldBe 18

    }

    "Test Shift on Set<SeqRange> " {
        var ranges : MutableSet<SeqRange> = mutableSetOf()

        ranges.add((5..25).toSeqRange())
        ranges.add((48..65).toSeqRange())
        ranges.add((75..87).toSeqRange())
        ranges.add((170..175).toSeqRange())
        ranges.add((180..210).toSeqRange())
        ranges.add((300..320).toSeqRange())

        // Test positive shift - shift right
        var shiftedRangeSet = shiftSeqRangeSet ( 5,    ranges,500)

        shiftedRangeSet.size shouldBe 6

        println("\nRanges for ShiftedRangeSet ,shift by 5")
        for (range in shiftedRangeSet) {
            println(range.toString())
        }
        shiftedRangeSet.contains((10..30).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((53..70).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((80..92).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((175..180).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((185..215).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((305..325).toSeqRange()) shouldBe true

        // Test negative shift (shift left)
        shiftedRangeSet = shiftSeqRangeSet ( -5,   ranges, 500)
        shiftedRangeSet.contains((1..20).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((43..60).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((70..82).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((165..170).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((175..205).toSeqRange()) shouldBe true
        shiftedRangeSet.contains((295..315).toSeqRange()) shouldBe true
    }

    "Test Merge Ranges " {
        var ranges : MutableSet<SeqRange> = mutableSetOf()

        ranges.add((48..65).toSeqRange())
        ranges.add((5..25).toSeqRange())
        ranges.add((27..30).toSeqRange())
        ranges.add((180..210).toSeqRange())
        ranges.add((170..175).toSeqRange())
        ranges.add((300..320).toSeqRange())
        ranges.add((69..75).toSeqRange())

        println("Unsorted , ummerged ranges:  ")
        for (range in ranges) {
            println(range.toString())
        }

        var mergedRangeSet = mergeSeqRangeSet ( 5, ranges)

        println("\nSORTED , MERGED ranges:  ")
        for (range in mergedRangeSet) {
            println(range.toString())
        }

        mergedRangeSet.size shouldBe 4

        println("\nSorted, merged ranges:")
        for (range in mergedRangeSet) {
            println(range.toString())
        }

        mergedRangeSet.contains((5..30).toSeqRange()) shouldBe true
        mergedRangeSet.contains((48..75).toSeqRange()) shouldBe true
        mergedRangeSet.contains((170..210).toSeqRange()) shouldBe true

    }
})