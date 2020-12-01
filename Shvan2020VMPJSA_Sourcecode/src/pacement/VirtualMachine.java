
package pacement;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 *
 * @author SHVAN
 */
public class VirtualMachine {

    double CPU;
    double RAM;
    double bandwidth;
    String Name;
    double weight;
    int res;
    String pm;
    int edge;
    int pod;
    double resource;

    public VirtualMachine(double C, double R, String N) {
        CPU = C;
        RAM = R;
        Name = N;
        weight = CPU / RAM;

    }

    public VirtualMachine() {
    }

    public VirtualMachine(double C, double R, String N, double B) {
        CPU = C;
        RAM = R;
        Name = N;
        bandwidth = B;
        weight = CPU / RAM;
        resource = CPU + RAM;

    }
    public static Comparator<VirtualMachine> ResourceAscending = new Comparator<VirtualMachine>() {

        public int compare(VirtualMachine v1, VirtualMachine v2) {
            double r1 = v1.CPU + v1.RAM;
            double r2 = v2.CPU + v2.RAM;
            return Double.compare(r2, r1);
        }
    };

    public static Comparator<VirtualMachine> WeightCompareAscending = new Comparator<VirtualMachine>() {

        public int compare(VirtualMachine v1, VirtualMachine v2) {
            double w1 = v1.weight;
            double w2 = v2.weight;
            return Double.compare(w1, w2);
        }
    };

    public static Comparator<VirtualMachine> WeightCompareDescending = new Comparator<VirtualMachine>() {

        public int compare(VirtualMachine v1, VirtualMachine v2) {
            double w1 = v1.weight;
            double w2 = v2.weight;
            return Double.compare(w2, w1);
        }
    };

    public static Comparator<VirtualMachine> CPUCompareDescending = new Comparator<VirtualMachine>() {

        public int compare(VirtualMachine v1, VirtualMachine v2) {

            double w1 = v1.CPU;
            double w2 = v2.CPU;
            return Double.compare(w2, w1);
        }
    };

    @Override
    public String toString() {
        return " CPU=" + CPU + "\t RAM=" + RAM + "\t Name=" + Name + "\t Bandwidth =" + bandwidth + "\t Weight = " + weight;
    }
}
