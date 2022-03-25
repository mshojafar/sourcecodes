import java.util.List;

public class Flow{
  public static int[] flow(Graph graph){
    double flow = 0;
    int[] edgeTo = BFS.BFS(graph, 0);
    while(edgeTo[edgeTo.length-1] != -1){
      List<Edge> path = BFS.getPath(edgeTo, graph);
      flow = augment(flow, path);
      edgeTo = BFS.BFS(graph, 0);
    }
    System.out.println("FLOW: " + flow);
    return edgeTo;
  }

  public static double augment(double flow, List<Edge> path){
    double bottleneck = Double.POSITIVE_INFINITY;
    for(Edge edge : path){ //find the bottleneck
      if(edge.getResidualCapacity() < bottleneck) {
        bottleneck = edge.getResidualCapacity();
      }
    }
    for(Edge edge : path){ //modify the flows on all edges
      edge.incrementFlow(bottleneck);
    }
    return (flow + bottleneck);
  }
}
