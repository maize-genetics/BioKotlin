package biokotlin.integration

import ch.qos.logback.classic.Level
import ch.qos.logback.classic.Logger
import org.slf4j.LoggerFactory

fun setUniProtLogging(level: Level = Level.INFO) {

    val loggers = listOf(
            "org.apache.http",
            "groovyx.net.http"
    )

    loggers.forEach { name ->
        val logger = LoggerFactory.getLogger(name) as Logger
        logger.setLevel(level)
        logger.setAdditive(false)
    }

}