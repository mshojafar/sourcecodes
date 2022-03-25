import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Graph{

  private int numNodes;
  private int numEdges;
  private Map<Integer, List<Edge>> node_to_edges;

  public Graph(int numNodes, int numEdges){
    this.numNodes = numNodes;
    this.numEdges = numEdges;
    this.node_to_edges = new HashMap<>();
    List<Edge> aux = new ArrayList<>();
  }

  public void addEdge(Edge edge){
    node_to_edges.putIfAbsent(edge.getFromNode(), new ArrayList<Edge>());
    node_to_edges.get(edge.getFromNode()).add(edge);
  }

  public String toString(){
    StringBuilder s = new StringBuilder();
    for(Integer node : node_to_edges.keySet()){
      s.append("Node: " + node);
      s.append("\n");
      for(Edge e : node_to_edges.get(node)){
        s.append(e.toString() + ", ");
      }
      s.append("\n");
    }
    return s.toString();
  }

  public int V(){
    return node_to_edges.keySet().size();
  }

  public List<Integer> adj(int V){
    List<Integer> toNodes = new ArrayList<>();
    for (Edge edge : node_to_edges.get(V)){
      if(edge.getResidualCapacity() > 0) toNodes.add(edge.getToNode());
    }
    return toNodes;
  }

  public Edge getEdge(int fromNode, int toNode){
    for(Edge edge : node_to_edges.get(fromNode)){
      if(edge.getToNode() == toNode) return edge;
    }
    return null;
  }

}
