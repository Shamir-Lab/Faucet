#include "ContigIterator.h" 

ContigIterator::ContigIterator(ContigGraph* theGraph){
    graph = theGraph;
    nodeIt = graph->getNodeMap()->begin();
    index = 0;
    findNextContig();
}

Contig* ContigIterator::getContig(){
    if(nodeIt != graph->getNodeMap() -> end()){
        Contig* result = nodeIt->second.contigs[index];
        increment();
        findNextContig();
        return result;
    }
}

void ContigIterator::increment(){
    if(index < 4){
        index++;
    }
    else{
        if(nodeIt != graph->getNodeMap()->end()){
            index = 0;
            nodeIt++;
        }
    }
}
    
Contig* ContigIterator::findNextContig(){    
    for( ; nodeIt != graph->getNodeMap()->end(); nodeIt++){
        ContigNode* node = &nodeIt->second;
        for(index %= 5; index < 5; index++){//if index is 5, reset to 0. Otherwise use it as is
            if(node->contigs[index]){
                Contig* contig = node->contigs[index];
                if(!contig->node1_p || !contig->node2_p){ //if one side is a sink, always return
                    return contig;
                }
                else if(contig->node1_p == contig->node2_p){ //if the contig attaches to the same node twice, print when you see lower index
                    if(index == contig->getMinIndex()){    
                        return contig;
                    }
                }
                else if(contig->getSide(node,index) == 1){ //If it attaches to two distinct nodes, print when you're on side 1
                     return contig;  
                }
            }       
        }
    }
    return NULL;
}
    
bool ContigIterator::hasNextContig(){
    return nodeIt != graph->getNodeMap()->end();
}
