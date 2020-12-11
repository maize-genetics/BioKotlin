package biokotlin.testData

class BedFileTestPaths {
    companion object {
        private val userHome = System.getProperty("user.home")
        val outputDir = "$userHome/temp/biokotlinTest/"
        val dataDir = outputDir + "data/"

        val bedFile1 = dataDir + "testBedFile.txt"
        val bedFileEnds = dataDir + "testBedFileEnds.txt"
    }
}