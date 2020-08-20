package biokotlin.genome

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class RangesTest: StringSpec ({
    "Test flank range "{

        var range1 = (48..68).toSRange()
        var flank5 = range1.flank(5, Direction.BOTH, 100)
        flank5.size shouldBe 2
        //flank5.contains(SRange.closed(Int1(43u), Int1(47u))) shouldBe true
        flank5.contains((43..47).toSRange()) shouldBe true

        // run again using other signature
        var flank5left = range1.flankLeft(5,100)
        //flank5.contains(SRange.closed(Int1(43u), Int1(47u))) shouldBe true
        flank5left.contains((43..47).toSRange()) shouldBe true

        var flank5right = range1.flankRight(5)
        flank5right.contains((69..73).toSRange())

        // Test right side is cropped
        var flank40 = range1.flank(40, Direction.BOTH, 100)
        flank40.size shouldBe 2
        println("\nRanges for flank40 of Range 48,68")
        for (range in flank40) {
            println(range.toString())
        }
        flank40.contains((8..47).toSRange()) shouldBe true



        // Test left side gets cropped
        var range2 = (5..58).toSRange()
        flank40 = range2.flank(40, Direction.BOTH, 100)
        flank40.size shouldBe 2
        println("\nRanges for flank40 of Range 5,58")
        for (range in flank40) {
            println(range.toString())
        }
        flank40.contains((1..4).toSRange()) shouldBe true
        flank40.contains((59..98).toSRange()) shouldBe true

    }

    "Test flank on Set<SRange> " {
        var ranges : MutableSet<SRange> = mutableSetOf()

        ranges.add((5..25).toSRange())
        ranges.add((48..65).toSRange())
        ranges.add((75..87).toSRange())
        ranges.add((170..175).toSRange())
        ranges.add((180..210).toSRange())
        ranges.add((300..320).toSRange())

        var flankedRangeSet = flankSRangeSet ( 5,  Direction.BOTH, 500, ranges)
        flankedRangeSet.size shouldBe 12

    }

    "Test shift range " {
        var range1 = (48..68).toSRange()
        var shift5 = range1.shift(5,  100)
        shift5.lowerEndpoint().site shouldBe 53u
        shift5.upperEndpoint().site shouldBe 73u

        var shift50 = range1.shift(50,100) // 100 is max size for range
        shift50.lowerEndpoint().site shouldBe 98u
        shift50.upperEndpoint().site shouldBe 100u

        var shiftNeg10 = range1.shift(-10,100)
        shiftNeg10.lowerEndpoint().site shouldBe 38u
        shiftNeg10.upperEndpoint().site shouldBe 58u

        // Shift left beyond 1
        var shiftNeg50 = range1.shift(-50,100)
        shiftNeg50.lowerEndpoint().site shouldBe 1u
        shiftNeg50.upperEndpoint().site shouldBe 18u

    }

    "Test Shift on Set<SRange> " {
        var ranges : MutableSet<SRange> = mutableSetOf()

        ranges.add((5..25).toSRange())
        ranges.add((48..65).toSRange())
        ranges.add((75..87).toSRange())
        ranges.add((170..175).toSRange())
        ranges.add((180..210).toSRange())
        ranges.add((300..320).toSRange())

        // Test positive shift - shift right
        var shiftedRangeSet = shiftSRangeSet ( 5,   500, ranges)

        shiftedRangeSet.size shouldBe 6

        println("\nRanges for ShiftedRangeSet ,shift by 5")
        for (range in shiftedRangeSet) {
            println(range.toString())
        }
        shiftedRangeSet.contains((10..30).toSRange()) shouldBe true
        shiftedRangeSet.contains((53..70).toSRange()) shouldBe true
        shiftedRangeSet.contains((80..92).toSRange()) shouldBe true
        shiftedRangeSet.contains((175..180).toSRange()) shouldBe true
        shiftedRangeSet.contains((185..215).toSRange()) shouldBe true
        shiftedRangeSet.contains((305..325).toSRange()) shouldBe true

        // Test negative shift (shift left)
        shiftedRangeSet = shiftSRangeSet ( -5,   500, ranges)
        shiftedRangeSet.contains((1..20).toSRange()) shouldBe true
        shiftedRangeSet.contains((43..60).toSRange()) shouldBe true
        shiftedRangeSet.contains((70..82).toSRange()) shouldBe true
        shiftedRangeSet.contains((165..170).toSRange()) shouldBe true
        shiftedRangeSet.contains((175..205).toSRange()) shouldBe true
        shiftedRangeSet.contains((295..315).toSRange()) shouldBe true
    }

    "Test Merge Ranges " {
        var ranges : MutableSet<SRange> = mutableSetOf()

        ranges.add((48..65).toSRange())
        ranges.add((5..25).toSRange())
        ranges.add((27..30).toSRange())
        ranges.add((180..210).toSRange())
        ranges.add((170..175).toSRange())
        ranges.add((300..320).toSRange())
        ranges.add((69..75).toSRange())

        println("Unsorted , ummerged ranges:  ")
        for (range in ranges) {
            println(range.toString())
        }

        var mergedRangeSet = mergeSRangeSet ( 5, ranges)
        mergedRangeSet.size shouldBe 4

        println("\nSorted, merged ranges:")
        for (range in mergedRangeSet) {
            println(range.toString())
        }

        mergedRangeSet.contains((5..30).toSRange()) shouldBe true
        mergedRangeSet.contains((48..75).toSRange()) shouldBe true
        mergedRangeSet.contains((170..210).toSRange()) shouldBe true

    }

})