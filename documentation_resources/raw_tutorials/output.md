Sequences are series of characters that can represent RNA, DNA or chains of amino acids. Within BioKotlin, they are represented by the [seq](.seq) package. Below, we demonstrate the basics of creating and manipulating sequences.

Start by importing the BioKotlin module, and the [seq](.seq) package.


```kotlin
@file:DependsOn("org.biokotlin:biokotlin:0.05")
import biokotlin.seq.*
```

Next, create a nucleotide sequence and output its [complement](.seq.NucSeq.complement). DNA is inferred (instead of RNA) because of the "T" within the sequence.


```kotlin
val dna = NucSeq("GCAGAT")
dna.complement()
```




    CGTCTA



Ensure that a sequence is a specific type when there is no thymine or uracil for inference with the [NUC](.seq.NUC) enum.


```kotlin
val rna = NucSeq("AGCG", NUC.RNA)
rna.complement() //Note that this outputs uracil because the sequence is RNA!
```




    UCGC



Now, [translate](.seq.NucSeq.translate) it to a protein sequence (will print both the codon table ID, defaulting to 1, and the translated sequence).


```kotlin
dna.translate()
```

    1





    AD



Now, [transcribe](.seq.NucSeq.transcribe) it to and RNA sequence


```kotlin
dna.transcribe()
```




    GCAGAU



Create a [protein sequence](.seq.ProteinSeq), count its number of glysines, and print it in a [String template](https://kotlinlang.org/docs/basic-types.html#string-templates). Here, the [AminoAcid](.seq.AminoAcid) enum is used, which stores useful properties of the amino acids, such as their three letter encoding.


```kotlin
val protein = ProteinSeq("GCAGAT")
val gCount = protein.count(AminoAcid.G)
println("${AminoAcid.G.name3letter} count is $gCount")
```

    Gly count is 2


The + and * operators can be used to concatenate two sequences or the same sequence with itself multiple times, respectively.


```kotlin
val otherProtein = ProteinSeq("ARSQRS")
val concatenatedProtein = protein + otherProtein
println("A concatenated protein: $concatenatedProtein")
println("A repeated DNA sequence: ${dna * 3}")
```

    A concatenated protein: GCAGATARSQRS
    A repeated DNA sequence: GCAGATGCAGATGCAGAT


Finally, let's find the mass of our concatenated protein in daltons. The [brackets](.seq.ProteinSeq.get) return an [AminoAcid](.seq.AminoAcid) when a single integer is used as the index. Each [AminoAcid](.seq.AminoAcid) stores its mass in the ```weight``` property.


```kotlin
var mass = 0.0
for (index in 0 until concatenatedProtein.size()) {
    mass += concatenatedProtein[index].weight
}
println("The of $concatenatedProtein is $mass daltons")
```

    The of GCAGATARSQRS is 1362.4218999999998 daltons

