package biokotlin.cli

import biokotlin.util.setupDebugLogging
import org.junit.jupiter.api.extension.BeforeAllCallback
import org.junit.jupiter.api.extension.ExtensionContext
import java.io.File

class TestExtension : BeforeAllCallback {

    companion object {

        // the next line converts Windows \ to linux / in the user home path
        private val userHome = System.getProperty("user.home").replace('\\', '/')
        val tempDir = "$userHome/temp/biokotlinTests/tempDir/"
        val testOutputGVCFDir = "${tempDir}outputGVCFDir/"
        val testOutputHVCFDir = "${tempDir}outputHVCFDir/"
        val testLineAGvcfOutput = "${testOutputGVCFDir}LineA.g.vcf/"
        val testLineAGvcfDiffRefOutput = "${testOutputGVCFDir}LineA-DiffRef.g.vcf/"

        const val gvcfInputDir = "data/test/gvcf/"

        const val refFasta = "data/test/ref/Ref.fa"
        const val lineAGvcf = "${gvcfInputDir}LineA.g.vcf"
        const val lineAGvcfDiffRef = "${gvcfInputDir}LineA-DiffRef.g.vcf"

    }

    override fun beforeAll(context: ExtensionContext) {

        // Setup the test document environments

        setupDebugLogging()

        File(tempDir).mkdirs()
        File(testOutputGVCFDir).mkdirs()
        File(testOutputHVCFDir).mkdirs()

    }
}