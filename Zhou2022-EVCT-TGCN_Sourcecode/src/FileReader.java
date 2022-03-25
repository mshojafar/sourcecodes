import java.util.*;
import java.io.*;

public class FileReader{

  public static Graph readIn(String filename, boolean directed){
    try{
      Scanner sc = new Scanner(new File(filename));
      int numNodes = sc.nextInt(); //number of vertices
      for(int i = 0; i <= numNodes; i++) sc.nextLine(); //scim over the vertex names
      int numEdges = Integer.parseInt(sc.nextLine());

      Graph graph = new Graph(numNodes, numEdges);

      for(int j = 0; j < numEdges; j++){
        String[] edge = sc.nextLine().split(" ");
        double capacity = Double.parseDouble(edge[2]);
        if(capacity == -1) capacity = Double.POSITIVE_INFINITY;

        Edge forwardEdge = new Edge(Integer.parseInt(edge[0]), Integer.parseInt(edge[1]), capacity, false, directed);
        Edge backwardEdge = new Edge(Integer.parseInt(edge[1]), Integer.parseInt(edge[0]), capacity, true, directed);
        forwardEdge.addPartnerEdge(backwardEdge);
        backwardEdge.addPartnerEdge(forwardEdge);
        graph.addEdge(forwardEdge);
        graph.addEdge(backwardEdge);
      }
      return graph;
    }
    catch(IOException e) {System.out.println("File not found.");}
    return null;
    }

}
