public class Edge{

  private int fromNode;
  private int toNode;
  private double capacity;
  private Edge partnerEdge;
  private boolean isResidualEdge;
  private double flow;
  private boolean directed;

  public Edge(int fromNode, int toNode, double capacity, boolean isResidualEdge, boolean directed){
    this.fromNode = fromNode;
    this.toNode = toNode;
    this.capacity = capacity;
    this.isResidualEdge = isResidualEdge;
    this.directed = directed;

    if(directed){
      if(isResidualEdge) this.flow = capacity;
      else this.flow = 0;
    }
    else{
      flow = 0;
    }
  }

  public void incrementFlow(double flow){
    this.flow += flow;
    if(!directed && !isResidualEdge) partnerEdge.setCapacity(this.flow);
    else partnerEdge.decrementFlow(flow);
  }

  public void decrementFlow(double flow){
    this.flow -= flow;
  }

  public double getResidualCapacity(){
    return capacity - flow;
  }

  public void addPartnerEdge(Edge partnerEdge){
    this.partnerEdge = partnerEdge;
  }

  public int getFromNode(){
    return fromNode;
  }

  public int getToNode(){
    return this.toNode;
  }

  public String toString(){
    return "" + fromNode + " " + toNode + " " + capacity;
  }

  public void setCapacity(double capacity){
    this.capacity = capacity;
  }
  public boolean isResidualEdge(){
    return this.isResidualEdge;
  }
}
