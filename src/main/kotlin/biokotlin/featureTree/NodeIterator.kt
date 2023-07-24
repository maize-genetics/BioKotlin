package biokotlin.featureTree

import java.util.*

internal class NodeIterator(root: NodeParent): Iterator<Node> {
    private val stack = Stack<Node>()
    private val graph = when(root) {
        is Graph -> root
        is Node -> root.graph
    }
    private val topologicalMutations = graph.topologicalMutations

    init {
        root.belows.forEach { stack.push(it) }
    }
    override fun hasNext(): Boolean = stack.isNotEmpty()

    override fun next(): Node {
        if (graph.topologicalMutations != topologicalMutations) throw ConcurrentModificationException()
        val popped = stack.pop()
        popped.belows.forEach { stack.add(it) }
        return popped
    }

}