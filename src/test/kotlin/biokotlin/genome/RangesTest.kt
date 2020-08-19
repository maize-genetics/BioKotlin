package biokotlin.genome

import io.kotest.core.spec.style.StringSpec
import io.kotest.matchers.shouldBe

class RangesTest: StringSpec ({
    "Test flank range "{

        var range1 = SRange.closed(Int1(48u), Int1(68u))
        var flank5 = range1.flank(5, "BOTH", 100)
        flank5.size shouldBe 2
        flank5.contains(SRange.closed(Int1(43u), Int1(73u))) shouldBe true

        // Test right side is cropped
        var flank40 = range1.flank(40, "BOTH", 100)
        flank40.size shouldBe 2
        flank40.contains(SRange.closed(Int1(8u), Int1(39u))) shouldBe true
        println("\nRanges for flank40 of Range 48,68")
        for (range in flank40) {
            println(range.toString())
        }

        // Test left side gets cropped
        var range2 = SRange.closed(Int1(5u), Int1(58u))
        flank40 = range2.flank(40, "BOTH", 100)
        flank40.size shouldBe 2
        flank40.contains(SRange.closed(Int1(1u), Int1(4u))) shouldBe true
        println("\nRanges for flank40 of Range 5,58")
        for (range in flank40) {
            println(range.toString())
        }
    }

    "Test flank on Set<SRange> " {
        var ranges : MutableSet<SRange> = mutableSetOf()

        ranges.add(SRange.closed(Int1(5u), Int1(25u)))
        ranges.add(SRange.closed(Int1(48u), Int1(65u)))
        ranges.add(SRange.closed(Int1(75u), Int1(87u)))
        ranges.add(SRange.closed(Int1(170u), Int1(175u)))
        ranges.add(SRange.closed(Int1(180u), Int1(210u)))
        ranges.add(SRange.closed(Int1(300u), Int1(320u)))

        var flankedRangeSet = flankSRangeSet ( 5,  "BOTH", 500, ranges)

    }

    "Test shift range " {
        var range1 = SRange.closed(Int1(48u), Int1(68u))
        var shift5 = range1.shift(5,  100)
        shift5.lowerEndpoint().site shouldBe 53
        shift5.upperEndpoint().site shouldBe 73

        var shift50 = range1.shift(50,100) // 100 is max size for range
        shift50.lowerEndpoint().site shouldBe 98
        shift50.upperEndpoint().site shouldBe 100

        var shiftNeg10 = range1.shift(-10,100)
        shiftNeg10.lowerEndpoint().site shouldBe 38
        shiftNeg10.upperEndpoint().site shouldBe 58

        // Shift left beyond 1
        var shiftNeg50 = range1.shift(-50,100)
        shiftNeg50.lowerEndpoint().site shouldBe 1
        shiftNeg50.upperEndpoint().site shouldBe 18

    }

    "Test Merge Ranges " {
        var ranges : MutableSet<SRange> = mutableSetOf()

        ranges.add(SRange.closed(Int1(48u), Int1(65u)))
        ranges.add(SRange.closed(Int1(5u), Int1(25u)))
        ranges.add(SRange.closed(Int1(27u), Int1(30u)))
        ranges.add(SRange.closed(Int1(180u), Int1(210u)))
        ranges.add(SRange.closed(Int1(170u), Int1(175u)))
        ranges.add(SRange.closed(Int1(300u), Int1(320u)))
        ranges.add(SRange.closed(Int1(69u), Int1(75u)))
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

        mergedRangeSet.contains(SRange.closed(Int1(5u), Int1(30u))) shouldBe true
        mergedRangeSet.contains(SRange.closed(Int1(48u), Int1(75u))) shouldBe true
        mergedRangeSet.contains(SRange.closed(Int1(170u), Int1(210u))) shouldBe true
    }

})