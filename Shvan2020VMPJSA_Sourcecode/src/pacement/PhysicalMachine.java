
package pacement;

import java.util.Comparator;

/**
 *
 * @author SHVAN
 */
public class PhysicalMachine {

    double CPU;
    double RAM;
    double PMAX;
    double PMIN;
    String Name;
    String reserve;
    int count;
    double weight;
    int EdgeSwitchNo;
    int PodSwitchNo;
    double PowerEfficincy;

    public PhysicalMachine(double CPU, double RAM, double pmax, double pmin, String Name, int edge, int pod) {
        this.CPU = CPU;
        this.RAM = RAM;
        weight = CPU / RAM;
        this.PMAX = pmax;
        this.PMIN = pmin;
        this.PowerEfficincy = CPU / PMAX;
        this.Name = Name;
        reserve = "";
        count = 0;
        EdgeSwitchNo = edge;
        PodSwitchNo = pod;
    }

    public static Comparator<PhysicalMachine> WeightCompareAscending = new Comparator<PhysicalMachine>() {

        public int compare(PhysicalMachine v1, PhysicalMachine v2) {
            double w1 = v1.RAM;
            double w2 = v2.RAM;
            return Double.compare(w1, w2);
        }
    };

    public static Comparator<PhysicalMachine> CompareCPUAscending = new Comparator<PhysicalMachine>() {

        public int compare(PhysicalMachine v1, PhysicalMachine v2) {
            double w1 = v1.CPU;
            double w2 = v2.CPU;
            return Double.compare(w1, w2);
        }
    };

    public static Comparator<PhysicalMachine> ComparePmNameAscending = new Comparator<PhysicalMachine>() {

        public int compare(PhysicalMachine v1, PhysicalMachine v2) {

            return v1.Name.compareTo(v2.Name);
        }
    };

    public static Comparator<PhysicalMachine> ComparePmName = new Comparator<PhysicalMachine>() {

        public int compare(PhysicalMachine v1, PhysicalMachine v2) {

            String pm1 = v1.Name;
            int x1 = pm1.indexOf("_");
            int pmNo1 = Integer.parseInt(pm1.substring(x1 + 1, pm1.length()));

            String pm2 = v2.Name;
            int x2 = pm2.indexOf("_");
            int pmNo2 = Integer.parseInt(pm2.substring(x2 + 1, pm2.length()));

            double w1 = pmNo1;
            double w2 = pmNo2;
            return Double.compare(w1, w2);

            //return v1.Name.compareTo(v2.Name);
        }
    };

    public static Comparator<PhysicalMachine> CompareCPUAs = new Comparator<PhysicalMachine>() {
        public int compare(PhysicalMachine v2, PhysicalMachine v1) {
            double w1 = v1.CPU;
            double w2 = v2.CPU;
            return Double.compare(w2, w1);
        }
    };

    public static Comparator<PhysicalMachine> SortEnergyAsending = new Comparator<PhysicalMachine>() {
        public int compare(PhysicalMachine v1, PhysicalMachine v2) {
            double w1 = v1.PowerEfficincy;
            double w2 = v2.PowerEfficincy;
            return Double.compare(w2, w1);
        }
    };

    public static Comparator<PhysicalMachine> SortPMEdgeSwitchEnergyAsending = new Comparator<PhysicalMachine>() {
        public int compare(PhysicalMachine v1, PhysicalMachine v2) {
            double w1=0,w2=0;
            if (v1.EdgeSwitchNo == v2.EdgeSwitchNo) {
                 w1 = v1.PowerEfficincy;
                 w2 = v2.PowerEfficincy;
            }
                return Double.compare(w2, w1);
            
        }
    };
    
    
    public static Comparator<PhysicalMachine> SortPMEdgeSwitchEnergy = new Comparator<PhysicalMachine>() {
        public int compare(PhysicalMachine v1, PhysicalMachine v2) {
            double w1=0,w2=0;
            if (v1.EdgeSwitchNo == v2.EdgeSwitchNo) {
                 w1 = (v1.CPU) / (v1.PMAX);
                 w2 = (v2.CPU) / (v2.PMAX);
            }
                return Double.compare(w2, w1);
            
        }
    };

}
