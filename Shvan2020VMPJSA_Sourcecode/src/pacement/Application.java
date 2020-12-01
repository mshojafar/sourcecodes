
package pacement;

import java.util.ArrayList;
import java.util.Comparator;

/**
 *
 * @author SHVAN
 */
public class Application {
    String type;
    ArrayList<VirtualMachine> VMList;
    int noVM;
    int weightArray[][];
    int type1;
    
    public Application(String t,int noVm,ArrayList<VirtualMachine> vmlist,int weight[][])
    {
        type=t;
        if(type.equalsIgnoreCase("Critical"))
            type1=0;
        else
            type1=1;
        
        noVM=noVm;
        weightArray=weight;
        VMList=vmlist;
    }
    
    public static Comparator<Application> TypeCompareAscending = new Comparator<Application>() {

        public int compare(Application v1, Application v2) {
            double w1 = v1.type1;
            double w2 = v2.type1;
            return Double.compare(w1, w2);
        }
    };   
}
