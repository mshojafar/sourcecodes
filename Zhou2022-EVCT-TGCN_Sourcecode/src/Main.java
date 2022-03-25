import java.util.ArrayList;
import java.util.List;

public class Main{

  public static void main(String[] args){
//    String[] files = new String[6];
//    files[0] = "E:/test1.txt";
//    files[1] = "E:/test2.txt";
//    files[2] = "E:/test3.txt";
//    files[3] = "E:/test4.txt";
//    files[4] = "E:/test5.txt";
//    files[5] = "E:/rail.txt";

//    for(int s = 0; s < files.length; s++){
//      System.out.println("\n" + "File: " + files[s]);
//      boolean directed = true;
//      if(s==5) directed = false;
      Graph graph = FileReader.readIn("E:/rail.txt", false);
      System.out.println(graph);
      int[] edgeTo = Flow.flow(graph);
      List<Integer> ACut = new ArrayList<>();
      List<Integer> BCut = new ArrayList<>();
      ACut.add(0); //manually add source to the A cut
      for(int i = 1; i < edgeTo.length; i++){
        if(edgeTo[i] != -1) {ACut.add(i);}
        else BCut.add(i);
      }
      System.out.println(ACut);
      System.out.println(BCut);
      System.out.println("Min cut edges: ");
      for(int i : ACut){
        for(int j : BCut){
          Edge cuttingEdge = graph.getEdge(i, j);
          if(cuttingEdge != null){
              System.out.println(graph.getEdge(i, j));
          }
        }
      }
  }
}
       
      
      
//      if(directed){ //print only forward edges
//            if(cuttingEdge != null && !cuttingEdge.isResidualEdge()){
//              System.out.println(graph.getEdge(i, j));
//            }
//          }
//          else{ //print all edges
//            if(cuttingEdge != null){
//              System.out.println(graph.getEdge(i, j));
//            }
//          }

//        }
//      }
//    }
//  }
//}

