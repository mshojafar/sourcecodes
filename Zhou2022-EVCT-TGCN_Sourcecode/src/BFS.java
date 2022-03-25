import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
//modified from Sedwich and Wayne,
public class BFS{

  public static int[] BFS(Graph G, int s) {
    boolean[] marked = new boolean[G.V()];
    int[] distTo = new int[G.V()];
    int[] edgeTo = new int[G.V()];
    for (int v = 0; v < G.V(); v++) distTo[v] = Integer.MAX_VALUE;
    for (int v = 0; v < G. V(); v++) edgeTo[v] = -1;
    validateVertex(s, marked);
    bfs(G, s, marked, distTo, edgeTo);
    return edgeTo;
  }

  private static void bfs(Graph G, int s, boolean[] marked, int[] distTo, int[] edgeTo) {
    Queue<Integer> q = new LinkedList<Integer>();
    marked[s] = true;
    distTo[s] = 0;
    q.add(s);
    while (!q.isEmpty()) {
      int v = q.poll();
      for (int w : G.adj(v)) {
        if (!marked[w]) {
          edgeTo[w] = v;
          distTo[w] = distTo[v] + 1;
          marked[w] = true;
          q.add(w);
        }
      }
    }
  }
  public static List<Edge> getPath(int[] edgeTo, Graph g){
    List<Edge> path = new ArrayList<>();
    for(int node = edgeTo.length - 1; node != 0; node = edgeTo[node]){
      path.add(g.getEdge(edgeTo[node], node));
    }
    return path;
  }

  private static void validateVertex(int v, boolean[] marked) {
    int V = marked.length;
    if (v < 0 || v >= V) throw new IllegalArgumentException("vertex " + v + " is not between 0 and " + (V-1));
  }
}
