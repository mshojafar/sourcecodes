package pacement;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

/**
 *
 * @author SHVAN
 */
public class Computation {

    ArrayList<PhysicalMachine> PM;
    ArrayList<VirtualMachine> VM;
    ArrayList<Application> app;
    ArrayList<Application> sortedApp;
    double weight[];
    int sharearray[][];

    ArrayList<PhysicalMachine> result;

    private double findResourceWastage(ArrayList<PhysicalMachine> result, double[] CPUPM, double[] RAMPM) {
        double wastage = 0;
        double TotalWastage = 0;
        int counter = 0;
        for (int l = 0; l < result.size(); l++) {
            double LC = 0;
            double LR = 0;
            double UC = 0;
            double UR = 0;
            double NUC = 0;
            double NUR = 0;
            double L, U = 0;
            double Pw = 0;
            if (result.get(l).reserve.equals("")) {
                ;
            } else {
                counter++;
                System.out.println("CPUPM " + CPUPM[l] + " result PM " + result.get(l).CPU);
                UC = Math.abs(CPUPM[l] - result.get(l).CPU);

                NUC = UC / CPUPM[l];
                System.out.println("Normalize UC " + NUC);
                UR = Math.abs(RAMPM[l] - result.get(l).RAM);
                NUR = UR / RAMPM[l];
                System.out.println("Normalize UR " + NUR);

                LC = 1 - NUC;
                LR = 1 - NUR;
                System.out.println("LC " + LC + "\nLR " + LR);
                L = Math.abs(LC - LR);
                L = L + 0.0001;
                Pw = L / (UC + UR);
                TotalWastage = TotalWastage + (Math.abs(result.get(l).CPU - result.get(l).RAM) + 0.0001) / (UC + UR);

                System.out.println(result.get(l).Name + " R W =" + Pw);
            }
        }
        //System.out.println("number of active pm "+counter);
        return TotalWastage;
    }

    private boolean placmentNormalApplication(int[] arraysortpod, double[] sortedarrayedge, ArrayList<PhysicalMachine> result, ArrayList<VirtualMachine> vms, int r, double[] CPUPM, double[] RAMPM, int j, ArrayList<Application> ap) {
        boolean apply = false;
        for (int i = 0; i < arraysortpod.length; i++) {

            for (int k = 0; k < sortedarrayedge.length; k++) {
                double MIN = Double.MAX_VALUE;
                int index = 0;
                for (int l = 0; l < result.size(); l++) {
                    double LC = 0;
                    double LR = 0;
                    double UC = 0;
                    double UR = 0;
                    double L, U = 0;
                    double Pw = 0;

                    String s = result.get(l).reserve;

                    if (s.contains(" App ") && result.get(l).EdgeSwitchNo == sortedarrayedge[k] && result.get(l).CPU >= vms.get(r).CPU && (result.get(l).CPU) - (vms.get(r).CPU) >= 0 && (result.get(l).RAM) - (vms.get(r).RAM) >= 0 && result.get(l).RAM >= vms.get(r).RAM && result.get(l).PodSwitchNo == arraysortpod[i] && result.get(l).EdgeSwitchNo == sortedarrayedge[k] && vms.get(r).res != 5) {

                        System.out.println("CPUPM " + CPUPM[l] + " result PM " + result.get(l).CPU + " VM CPU " + vms.get(r).CPU);
                        UC = Math.abs(CPUPM[l] - result.get(l).CPU + vms.get(r).CPU);

                        UC = UC / CPUPM[l];
                        System.out.println("Normalize UC " + UC);
                        UR = Math.abs(RAMPM[l] - result.get(l).RAM + vms.get(r).RAM);
                        UR = UR / RAMPM[l];
                        System.out.println("Normalize UR " + UR);

                        LC = 1 - UC;
                        LR = 1 - UR;
                        System.out.println("LC " + LC + "\nLR " + LR);
                        L = LC - LR;
                        //System.out.println("L1 " + L);
                        if (L < 0) {
                            L = Math.abs(L);
                        }
                        //System.out.println("L2 " + L);
                        if (L == 0.0) {
                            L = 0.0001;
                        } else {
                            L = L + 0.0001;
                        }
                        //System.out.println("L3 " + L);
                        double UU = (UC + UR);
                        Pw = L / UU;
                        System.out.println("L= " + L + " U= " + UU);
                        System.out.println(result.get(l).Name + " VM " + vms.get(r).Name + " R W = " + Pw);
                        if (Pw < MIN) {
                            MIN = Pw;
                            index = l;
                        }
                    }
                }

                if (MIN != Double.MAX_VALUE) {

                    System.out.println("norrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr");
                    System.out.println(" Minimume resource Wastage " + MIN + " result " + result.get(index).Name);
                    //System.out.println("Rrsource wastage "+MIN);
                    result.get(index).CPU = result.get(index).CPU - (vms.get(r).CPU);
                    result.get(index).RAM = result.get(index).RAM - (vms.get(r).RAM);
                    result.get(index).reserve = result.get(index).reserve + "\t  App " + j + " type " + ap.get(j).type + ": vm " + r + " :" + vms.get(r).Name;
                    result.get(index).count = result.get(index).count + 1;
                    vms.get(r).res = 5;
                    //counterapp = counterapp + 1;
                    vms.get(r).pm = result.get(index).Name + "";  /// Or edge Or POD 
                    vms.get(r).edge = result.get(index).EdgeSwitchNo;
                    vms.get(r).pod = result.get(index).PodSwitchNo;
                    System.out.println("VM " + r + " placed on pm " + result.get(index).Name);
                    System.out.println("pm " + result.get(index).Name + " : " + result.get(index).reserve);
                    //withoutdublicate.add(r);
                    apply = true;

                } else {
                    ;
                }

            }

        }
        if (apply == false) {

            for (int i = 0; i < arraysortpod.length; i++) {

                for (int k = 0; k < sortedarrayedge.length; k++) {
                    double MIN = Double.MAX_VALUE;
                    int index = 0;
                    for (int l = 0; l < result.size(); l++) {
                        double LC = 0;
                        double LR = 0;
                        double UC = 0;
                        double UR = 0;
                        double L, U = 0;
                        double Pw = 0;

                        if (result.get(l).EdgeSwitchNo == sortedarrayedge[k] && result.get(l).CPU >= vms.get(r).CPU && (result.get(l).CPU) - (vms.get(r).CPU) >= 0 && (result.get(l).RAM) - (vms.get(r).RAM) >= 0 && result.get(l).RAM >= vms.get(r).RAM && result.get(l).PodSwitchNo == arraysortpod[i] && result.get(l).EdgeSwitchNo == sortedarrayedge[k] && vms.get(r).res != 5) {
                            System.out.println("ssssssssssssssssssssssssssssssssssss " + result.get(l).Name);
                            System.out.println("CPUPM " + CPUPM[l] + " result PM " + result.get(l).CPU + " VM CPU " + vms.get(r).CPU);
                            UC = Math.abs(CPUPM[l] - result.get(l).CPU + vms.get(r).CPU);

                            UC = UC / CPUPM[l];
                            System.out.println("Normalize UC " + UC);
                            UR = Math.abs(RAMPM[l] - result.get(l).RAM + vms.get(r).RAM);
                            UR = UR / RAMPM[l];
                            System.out.println("Normalize UR " + UR);

                            LC = 1 - UC;
                            LR = 1 - UR;
                            System.out.println("LC " + LC + "\nLR " + LR);
                            L = LC - LR;
                            //System.out.println("L1 " + L);
                            if (L < 0) {
                                L = Math.abs(L);
                            }
                            //System.out.println("L2 " + L);
                            if (L == 0.0) {
                                L = 0.0001;
                            } else {
                                L = L + 0.0001;
                            }
                            //System.out.println("L3 " + L);
                            double UU = (UC + UR);
                            Pw = L / UU;
                            System.out.println("L= " + L + " U= " + UU);
                            System.out.println(result.get(l).Name + " VM " + vms.get(r).Name + " R W = " + Pw);
                            if (Pw < MIN) {
                                MIN = Pw;
                                index = l;
                            }
                        }
                    }

                    if (MIN != Double.MAX_VALUE) {

                        System.out.println("norrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr");
                        System.out.println(" Minimume resource Wastage " + MIN + " result " + result.get(index).Name);
                        //System.out.println("Rrsource wastage "+MIN);
                        result.get(index).CPU = result.get(index).CPU - (vms.get(r).CPU);
                        result.get(index).RAM = result.get(index).RAM - (vms.get(r).RAM);
                        result.get(index).reserve = result.get(index).reserve + "\t  App " + j + " type " + ap.get(j).type + ": vm " + r + " :" + vms.get(r).Name;
                        result.get(index).count = result.get(index).count + 1;
                        vms.get(r).res = 5;
                        //counterapp = counterapp + 1;
                        vms.get(r).pm = result.get(index).Name + "";  /// Or edge Or POD 
                        vms.get(r).edge = result.get(index).EdgeSwitchNo;
                        vms.get(r).pod = result.get(index).PodSwitchNo;
                        System.out.println("VM " + r + " placed on pm " + result.get(index).Name);
                        System.out.println("pm " + result.get(index).Name + " : " + result.get(index).reserve);
                        //withoutdublicate.add(r);
                        apply = true;

                    } else {
                        ;
                    }

                }

            }
        }

        return apply;

    }

    public ArrayList<VirtualMachine> sortVMByResourceDecreascing(ArrayList<VirtualMachine> arraylist) {

        ArrayList<VirtualMachine> temp = new ArrayList<VirtualMachine>(arraylist.size());
        try {

            Collections.sort(arraylist, VirtualMachine.ResourceAscending);

            for (VirtualMachine str : arraylist) {
                temp.add(str);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return temp;
    }

    static class Edge {

        int source;
        int destination;
        int weight;

        public Edge(int source, int destination, int weight) {
            this.source = source;
            this.destination = destination;
            this.weight = weight;
        }

        public void updateEdge(Graph g, int source, int destination, int weight) {

            g.update(source, destination, weight);
        }
    }

    static class Graph {

        List<List<Edge>> adjList = null;

        // Constructor
        Graph(List<Edge> edges, int N) {
            vertices = N;
            adjList = new ArrayList<>(N);

            for (int i = 0; i < N; i++) {
                adjList.add(i, new ArrayList<>());
            }

            // add edges to the undirected graph
            for (Edge edge : edges) {
                adjList.get(edge.source).add(edge);
            }
        }

        int vertices;

        public void addEgde(int source, int destination, int weight) {
            Edge edge = new Edge(source, destination, weight);
            adjList.get(source).add(edge);
            //adjacencylist[source].addFirst(edge); //for directed graph
        }

        public void printGraph() {
            for (int i = 0; i < vertices; i++) {
                ArrayList<Edge> list = (ArrayList<Edge>) adjList.get(i);//adjacencylist[i];
                for (int j = 0; j < list.size(); j++) {
                    System.out.println("vertex-" + i + " is connected to "
                            + list.get(j).destination + " with weight " + list.get(j).weight);
                }
            }
        }

        public void countTotalBandwidthGraph() {
            long TotalBandwidth = 0;
            for (int i = 0; i < vertices; i++) {
                ArrayList<Edge> list = (ArrayList<Edge>) adjList.get(i);//adjacencylist[i];
                for (int j = 0; j < list.size(); j++) {
                    if (list.get(j).weight < 0) {
                        TotalBandwidth = TotalBandwidth + list.get(j).weight;
                    }
                    if (list.get(j).weight >= 0 && list.get(j).weight != 10) {

                        TotalBandwidth = TotalBandwidth - (10 - list.get(j).weight);
                    }
                    //System.out.println("vertex-" + i + " is connected to "+ list.get(j).destination + " with weight " + list.get(j).weight);
                }
            }
            System.out.println("Total used Banwidth = " + TotalBandwidth);
        }

//        public void updateEdge(Graph g, int sourc, int destinat, int weight) {
//
//            if (sourc == destinat) {
//                System.out.println("same source and destination");
//            } else {
//                System.out.println("size " + ListPath.size());
//                for (int i = 0; i < ListPath.size() - 1; i++) {
//
//                    //g.update(ListPath.get(i), ListPath.get(i + 1), (g.getWeight(i, i + 1)) - weight);
//                    g.update(ListPath.get(i), ListPath.get(i + 1), weight);
//                }
//
//            }
//
//        }
//        public void updateEdgeOurAlgorithm(Graph g, int sourc, int destinat, int weight) {
//
//            if (sourc == destinat) {
//                System.out.println("same source and destination");
//            } else {
//                System.out.println("size " + ListPath.size());
//                for (int i = 0; i < ListPath.size() - 1; i++) {
//
//                    g.updateOurAlgorithm(ListPath.get(i), ListPath.get(i + 1),weight);
//                    //g.updateOurAlgorithm(ListPath.get(i), ListPath.get(i + 1), (g.getWeight(i, i + 1)) - weight);
//                }
//
//            }
//
//        }
        public int getWeight(int s, int d) {
            for (int i = 0; i < vertices; i++) {
                ArrayList<Edge> list = (ArrayList<Edge>) adjList.get(i);
                //LinkedList<Edge> list = adjacencylist[i];
                for (int j = 0; j < list.size(); j++) {
                    if (s == i && d == list.get(j).destination) {

                        //System.out.println("vertex-" + i + " is connected to " + list.get(j).destination + " with weight " + list.get(j).weight);
                        return list.get(j).weight;
                    }
                }
            }
            //System.out.println("00000000000000000000000000000000000000000000000");
            return 0;
        }

        public boolean find(int s, int d) {
            for (int i = 0; i < vertices; i++) {
                ArrayList<Edge> list = (ArrayList<Edge>) adjList.get(i);
                //LinkedList<Edge> list = adjacencylist[i];
                for (int j = 0; j < list.size(); j++) {
                    if (s == i && d == list.get(j).destination) {

                        System.out.println("vertex-" + i + " is connected to " + list.get(j).destination + " with weight " + list.get(j).weight);
                        return true;
                    }
                }
            }
            return false;
        }

        public boolean update(int s, int d, int w) {
            for (int i = 0; i < vertices; i++) {
                ArrayList<Edge> list = (ArrayList<Edge>) adjList.get(i);
                //LinkedList<Edge> list = adjacencylist[i];
                for (int j = 0; j < list.size(); j++) {
                    if (s == i && d == list.get(j).destination) {
                        list.get(j).weight = (list.get(j).weight) - w;
                        //list.get(i).weight=w;
                        System.out.println("vertex-" + i + " is connected to " + list.get(j).destination + " with weight " + list.get(j).weight);
                        return true;
                    }
                }
            }
            return false;
        }

//        public boolean updateOurAlgorithm(int s, int d, int w) {
//            for (int i = 0; i < vertices; i++) {
//                ArrayList<Edge> list = (ArrayList<Edge>) adjList.get(i);
//                //LinkedList<Edge> list = adjacencylist[i];
//                for (int j = 0; j < list.size(); j++) {
//                    if (s == i && d == list.get(j).destination) {
//                        list.get(j).weight = (list.get(j).weight) + w;
//                        //list.get(i).weight=w;
//                        System.out.println("vertex-" + i + " is connected to " + list.get(j).destination + " with weight " + list.get(j).weight);
//                        return true;
//                    }
//                }
//            }
//            return false;
//        }
    }

    public ArrayList<VirtualMachine> sortVMByCPUDescending(ArrayList<VirtualMachine> arraylist) {

        ArrayList<VirtualMachine> temp = new ArrayList<VirtualMachine>(arraylist.size());
        try {
            Collections.sort(arraylist, VirtualMachine.CPUCompareDescending);
            int count = 0;
            for (VirtualMachine str : arraylist) {
                temp.add(str);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return temp;
    }

    public ArrayList<PhysicalMachine> sortPMByName(ArrayList<PhysicalMachine> arraylist) {

        ArrayList<PhysicalMachine> temp = new ArrayList<PhysicalMachine>(arraylist.size());
        try {

            Collections.sort(arraylist, PhysicalMachine.ComparePmName);

            int count = 0;
            for (PhysicalMachine str : arraylist) {
                temp.add(str);

            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return temp;
    }

    public ArrayList<PhysicalMachine> sortPMByCPUAsending(ArrayList<PhysicalMachine> arraylist) {

        ArrayList<PhysicalMachine> temp = new ArrayList<PhysicalMachine>(arraylist.size());
        try {

            Collections.sort(arraylist, PhysicalMachine.CompareCPUAs);

            int count = 0;
            for (PhysicalMachine str : arraylist) {
                temp.add(str);
                //System.out.println(count++ + " - " + str);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return temp;
    }

    public void printPM(ArrayList<PhysicalMachine> P) {

        for (int i = 0; i < P.size(); i++) {

            System.out.println(i + 1 + " - " + P.get(i).Name + "\t CPU = " + P.get(i).CPU + "\t RAM = " + P.get(i).RAM + "\t Pmax = " + P.get(i).PMAX + "\t Pmin = " + P.get(i).PMIN + "\t Power Efficincy = " + P.get(i).PowerEfficincy);
        }
    }

    public void printPM() {
        for (int i = 0; i < PM.size(); i++) {

            System.out.println(i + 1 + " - " + PM.get(i).Name + "\t No of Core = " + PM.get(i).CPU + "\t RAM = " + PM.get(i).RAM + "\t  By " + PM.get(i).reserve);
        }
    }

    public void printVM() {
        try {
            for (int i = 0; i < VM.size(); i++) {
                System.out.println(i + " - " + VM.get(i).Name + " \t CPU = " + VM.get(i).CPU + "\t RAM = " + VM.get(i).RAM);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void printListVM(ArrayList<VirtualMachine> vmlist) {
        try {
            for (int i = 0; i < vmlist.size(); i++) {
                System.out.println(i + " - " + vmlist.get(i).Name + " \t CPU = " + vmlist.get(i).CPU + "\t RAM = " + vmlist.get(i).RAM);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    int countdown2 = 0;

    int countdown = 0;

    Map VMmap;
    int indexmap[];

    private ArrayList<VirtualMachine> findVMWithMaxConnection(int[][] matrix, ArrayList<VirtualMachine> vms) {

        if (vms.size() == 1) {
            return vms;
        } else {

            VMmap = new HashMap();
            indexmap = new int[vms.size()];
            int indexmapcounter = 0;
            ArrayList<Integer> al = new ArrayList<Integer>();
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix.length; j++) {
                    //System.out.println(matrix[i][j]);
                    matrix[j][i] = matrix[i][j];
                }
            }

            int vm = 0;
            int max = 0;
            int index = -1;
            for (int i = 0; i < matrix.length; i++) {
                max = 0;
                for (int j = 0; j < matrix[i].length; j++) {
                    max = max + matrix[i][j];
                }
                if (max > vm) {
                    vm = max;
                    index = i;
                }

            }
            System.out.println("first index = " + index);
            ArrayList<VirtualMachine> vmlist = new ArrayList<VirtualMachine>();
            vmlist.add(vms.get(index));

            //int inindex = -1;
            int ii = -1;
            int jj = -1;
            while (vmlist.size() != vms.size()) {
                int inmax = 0;
                //System.out.println("test");
                for (int i = 0; i < vmlist.size(); i++) {
                    //System.out.println("inside for 1");
                    for (int j = 0; j < matrix[i].length; j++) {
                        //System.out.println("inside for 2");
                        String VmName = vmlist.get(i).Name;

                        int y = VmName.indexOf("_");
                        int No = Integer.parseInt(VmName.substring(y + 1, VmName.length()));

                        //System.out.println("matrix value " + matrix[No][j]);
                        if (matrix[No][j] >= inmax) {
                            inmax = matrix[No][j];
                            //inindex = j;
                            ii = No;
                            jj = j;
                        }
                    }
                }
                if (!vmlist.contains(vms.get(jj))) {
                    System.out.println("addddd " + vms.get(jj));
                    vmlist.add(vms.get(jj));
                    System.out.println("jj " + jj + " ii " + ii);
                    String vmap = " vm " + ii;
                    indexmap[indexmapcounter++] = jj;
                    VMmap.put(jj, vmap);
                    matrix[ii][jj] = matrix[jj][ii] = 0;
                } else {
                    matrix[ii][jj] = matrix[jj][ii] = 0;
                }
                for (int i = 0; i < matrix.length; i++) {
                    for (int j = 0; j < matrix[i].length; j++) {
                        System.out.print(matrix[i][j] + "\t");
                    }
                    System.out.println("");

                }
                for (int i = 0; i < vmlist.size(); i++) {
                    System.out.println("vm data " + vmlist.get(i));
                }
                for (int i = 0; i < VMmap.size(); i++) {
                    String element1 = (String) VMmap.get(indexmap[i]);
                    System.out.println(indexmap[i] + " by " + element1);
                }
                //System.out.println("VMlist size "+vmlist.size());
            }
            return vmlist;
        }
    }

    public int[][] generateWeight(int n) {
        int array[][] = new int[n][n];

        double percentage = 0.3;
        int ch[] = new int[n * n];
        for (int i = 0; i < ch.length; i++) {
            if (i <= percentage * ch.length) {
                ch[i] = 1;
            } else {
                ch[i] = 0;
            }
        }
        ch = RandomizeArray(ch);
        int chcounter = 0;
        for (int i = 0; i < array.length; i++) {
            for (int j = i; j < array[i].length; j++) {

                if (i == j) {
                    chcounter++;
                    array[i][j] = 0;
                } else {

                    //int chk[] = new int[n * n];
                    if (ch[chcounter++] == 1) {
                        array[i][j] = (int) ((Math.random() * ((100 - 1))) + 1);
                    } else if (ch[chcounter++] == 0) {
                        array[i][j] = 0;
                    }
                }
            }

        }

        return array;
    }

    public int[][] checkMatrix(int array[][]) {
        for (int i = 0; i < array.length - 1; i++) {
            int zerocounter = 0;
            for (int j = 0; j < array[i].length; j++) {
                if (i == array[i].length) {
                    continue;
                }
                if (array[i][j] == 0) {
                    zerocounter++;
                }
            }
            if (zerocounter == array[i].length) {
                Random r = new Random();
                int low = i + 1;
                int high = array[i].length;
                int indexchange = r.nextInt(high - low) + low;
                array[i][indexchange] = r.nextInt(100 - 2) + 2;
            }
        }
        return array;
    }

    public void printWeightMatrix(int array[][]) {
        System.out.println("");
        for (int i = 0; i < array.length; i++) {
            System.out.print("\tVM" + i);
        }
        System.out.println("");
        //System.out.println("\tVM0 \t VM1 \t VM2 \t VM3 \t VM4 \t VM5 \t");
        for (int i = 0; i < array.length; i++) {
            System.out.print("VM" + i);
            for (int j = 0; j < array[i].length; j++) {
                System.out.print("\t" + array[i][j]);
                //array[i][j]=(int) ((Math.random()*((10-0)+10))+0);

            }
            System.out.println("");

        }
    }

    public ArrayList<VirtualMachine> createVMApp(int n) {

        ArrayList<VirtualMachine> vmList = new ArrayList<>(n);

        for (int i = 0; i < n; i++) {

            Random random = new Random();
            int check = random.nextInt(12);//((i * Math.random()) % 16);

            if (check == 0) {
                vmList.add(new VirtualMachine(1, 1, "StandardB1s"));

            } else if (check == 1) {
                vmList.add(new VirtualMachine(1, 2, "StandardB1ms"));

            } else if (check == 2) {
                vmList.add(new VirtualMachine(2, 4, "StandardB2s"));

            } else if (check == 3) {
                vmList.add(new VirtualMachine(2, 8, "StandardB2ms"));

            } else if (check == 4) {
                vmList.add(new VirtualMachine(4, 16, "StandardB4ms"));

            } else if (check == 5) {
                vmList.add(new VirtualMachine(8, 32, "StandardB8ms"));

            } else if (check == 6) {
                vmList.add(new VirtualMachine(12, 48, "StandardB12ms"));

            } else if (check == 7) {
                vmList.add(new VirtualMachine(16, 64, "StandardB16ms"));

            } else if (check == 8) {
                vmList.add(new VirtualMachine(20, 80, "StandardB20ms"));

            } else if (check == 9) {
                vmList.add(new VirtualMachine(4, 8, "m4.large"));

            } else if (check == 10) {
                vmList.add(new VirtualMachine(8, 16, "m4.xlarge"));

            } else if (check == 11) {
                vmList.add(new VirtualMachine(16, 32, "m4.2xlarge"));

            }
        }
        return vmList;
    }

    int NormalCounter = 0;
    int CriticalCounter = 0;

    public void printList(int array[][]) {
        ArrayList<Integer>[] al = new ArrayList[array.length];

        // initializing 
        for (int i = 0; i < al.length; i++) {
            al[i] = new ArrayList<Integer>();
        }

        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[i].length; j++) {
                if (array[i][j] != 0) {
                    al[i].add(j);
                }
            }
        }
        for (int i = 0; i < al.length; i++) {
            System.out.print("VM" + i + "\t");
            for (int j = 0; j < al[i].size(); j++) {
                System.out.print("VM" + al[i].get(j) + "\t");
            }
            System.out.println("");

        }

    }

    public void printValue(int array[][]) {
        ArrayList<Integer>[] al = new ArrayList[array.length];

        // initializing 
        for (int i = 0; i < al.length; i++) {
            al[i] = new ArrayList<Integer>();
        }

        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[i].length; j++) {
                if (array[i][j] != 0) {
                    al[i].add(+array[i][j]);
                }
            }
        }

        for (int i = 0; i < al.length; i++) {
            System.out.print("VM" + i + "\t");
            for (int j = 0; j < al[i].size(); j++) {
                System.out.print(al[i].get(j) + "\t");
            }
            System.out.println("");

        }

    }

    public double[] orderIndexAscending(double[] array) {
        double[] minimumIndexes = new double[array.length];
        double[] sortedArray = array.clone();

        Arrays.sort(sortedArray);
        Double a[] = new Double[array.length];
        for (int i = 0; i < array.length; i++) {
            a[i] = array[i];
        }
        Arrays.sort(a, Collections.reverseOrder());
        for (int i = 0; i < a.length; i++) {
            sortedArray[i] = a[i];
        }
        Set<Integer> savedIndexes = new HashSet<>();

        for (int index = 0; index < array.length; index++) {
            int minIndex = 0;
            // Add the index in ascending order, we need to keep the indexes already
            for (double number : array) {
                if (number == sortedArray[index] && !savedIndexes.contains(minIndex)) {
                    savedIndexes.add(minIndex);
                    minimumIndexes[index] = minIndex;
                    break;
                }
                minIndex++;
            }
        }
        return minimumIndexes;
    }

    public void placementComputationPAVA(ArrayList<PhysicalMachine> pm, ArrayList<Application> ap, int noCriticalApp) {
        try {

            result = new ArrayList<PhysicalMachine>(pm.size());
            result = pm;

            double CPUPM[] = new double[pm.size()];
            double RAMPM[] = new double[pm.size()];
            double PowerEff[] = new double[pm.size()];
            for (int i = 0; i < CPUPM.length; i++) {
                CPUPM[i] = pm.get(i).CPU;
                RAMPM[i] = pm.get(i).RAM;
                PowerEff[i] = pm.get(i).PowerEfficincy;
            }

            ArrayList<ArrayList<VirtualMachine>> ListOfList = new ArrayList<ArrayList<VirtualMachine>>(ap.size());

            double pmsize = pm.size();
            double podsize = Math.cbrt(pmsize * 4);
            double NoOfEdgeSwitchs = (podsize / 2) * podsize;
            System.out.println(" No of Edge Switches " + NoOfEdgeSwitchs);

            double EdgeEnergy[] = new double[(int) NoOfEdgeSwitchs];

            int startEdge = result.get(0).EdgeSwitchNo;
            int endEdge = 0;

            System.out.println(" pod " + podsize);
            System.out.println("NoOfEdgeSwitchs " + NoOfEdgeSwitchs);

            double noPODPm = Math.pow((podsize / 2), 2);
            //System.out.println("NoPOrt" +NoPort);
            double noEdgePM = (podsize / 2);
            ArrayList<VirtualMachine> vms = new ArrayList<VirtualMachine>();

            List<int[][]> weightlist = new ArrayList<int[][]>();

            for (int o = 0; o < ap.size(); o++) {

                if (ap.get(o).type.equalsIgnoreCase("Critical")) {

                    System.out.println("Application " + o + " Critical");

                    result = sortPMByName(result);
                    vms = ap.get(o).VMList;
                    ListOfList.add(vms);
                    int TempShareArray[][] = ap.get(o).weightArray;
                    weightlist.add(TempShareArray);

                    sharearray = ap.get(o).weightArray;

                    boolean ch = false;

                    int counterapp = 0;

                    int notservised = 0;
                    //result = QH;
                    ArrayList<Integer> withoutdublicate = new ArrayList<Integer>();
                    ArrayList<Integer> listnotservised = new ArrayList<Integer>();
                    for (int vmi = 0; vmi < vms.size(); vmi++) {

                        ArrayList<PhysicalMachine> QH = new ArrayList<PhysicalMachine>();
                        ArrayList<PhysicalMachine> Happ = new ArrayList<PhysicalMachine>();
                        double capacity[] = new double[(int) NoOfEdgeSwitchs];
                        int arraysortpod[] = new int[capacity.length];
                        int arraysortedge[] = new int[capacity.length];
                        ArrayList<Integer> checklist = new ArrayList<Integer>();
                        ArrayList<Integer> remain = new ArrayList<Integer>();
                        int pmsortindex = 0;

                        ///// line 10 algorithm ////
                        for (int ha = 0; ha < result.size(); ha++) {
                            //System.out.println("reserrrrrrrrrrrrrrrrrrr " + result.get(ha).reserve);
                            String s = (result.get(ha).reserve);

                            if (s.contains(" App " + o)) {

                                System.out.println("zyad");
                                if (!Happ.contains(result.get(ha))) {
                                    Happ.add(result.get(ha));
                                }

                            }
                        }
                    ////////////////////////end line 10/////

                        //////line 11 until line 21///
                        if (!Happ.isEmpty()) {
                            for (int i = 0; i < Happ.size(); i++) {
                                if (!QH.contains(Happ.get(i))) {
                                    QH.add(Happ.get(i));  //line 12
                                }
                            }
                            for (int Happindex = 0; Happindex < Happ.size(); Happindex++) {
                                System.out.println("Happ.get(Happindex).EdgeSwitchNo  " + Happ.get(Happindex).EdgeSwitchNo);
                                for (int indexedge = 0; indexedge < result.size(); indexedge++) {
                                    if (result.get(indexedge).EdgeSwitchNo == Happ.get(Happindex).EdgeSwitchNo && !QH.contains(result.get(indexedge))) {
                                        QH.add(result.get(indexedge));
                                        System.out.println("Edsge " + result.get(indexedge).Name);
                                    }
                                }
                            }
                            for (int Happindex = 0; Happindex < Happ.size(); Happindex++) {
                                System.out.println("Happ.get(Happindex).PodSwitchNo  " + Happ.get(Happindex).PodSwitchNo);

                                for (int indexedge = 0; indexedge < result.size(); indexedge++) {
                                    if (result.get(indexedge).PodSwitchNo == Happ.get(Happindex).PodSwitchNo && !QH.contains(result.get(indexedge))) {
                                        QH.add(result.get(indexedge));
                                        System.out.println("poddd " + result.get(indexedge).Name);
                                    }
                                }
                            }
                        }
                        System.out.println("End ! isEmpty()");
                        //////End line 11 until line 21///
                        ////////////////////// Find capacity and sort Hgroup high to low///

                        System.out.println("Startttttttttttttttt edgeeeeeeeee " + startEdge);

                        /////  resort pms ascending order by id
                        PhysicalMachine arrayofPM[] = new PhysicalMachine[result.size()];

                        result = sortPMByName(result);
                        ///////end resort

                        int pcounter = -1;
                        for (int i = 0; i < result.size(); i++) {

                            if (i % noEdgePM == 0) {
                                pcounter++;
                                //System.out.println("iiii " + i + " " + pcounter);
                            }
                            capacity[pcounter] = capacity[pcounter] + result.get(i).CPU;//CPUPM[pmNomber-1];

                        }

//////                    sort edge by capacity high to low
                        for (int i = 0; i < capacity.length; i++) {
                            System.out.println("Sum of Capacity  of Edge " + (i + startEdge) + " = " + capacity[i]);
                        }

                        double sortedEdgeCapacity[] = orderIndexAscending(capacity);

                        for (int i = 0; i < sortedEdgeCapacity.length; i++) {
                            sortedEdgeCapacity[i] = sortedEdgeCapacity[i] + startEdge;
                            System.out.println("sorted edge " + sortedEdgeCapacity[i]);
                        }

                        for (int i = 0; i < sortedEdgeCapacity.length; i++) {
                            for (int j = 0; j < result.size(); j++) {
                                if (result.get(j).EdgeSwitchNo == (int) sortedEdgeCapacity[i] && !QH.contains(result.get(j))) {
                                    QH.add(result.get(j));
                                }
                            }
                        }
                        System.out.println("length of QHH " + QH.size());
                        for (int i = 0; i < QH.size(); i++) {
                            System.out.println("Final List " + QH.get(i).Name);
                        }
                        //////////////////////End Find capacity and sort Hgroup high to low---- Lin 22///

                        int pmcounter = 0;
                        result = QH;
                        for (int p = 0; p < QH.size(); p++) {
                            pmcounter++;

                            //for (int p = 0; p < allpmsortedPM.size(); p++) {
                            boolean apply = false;

                            if (apply == false) {
                                if (result.get(p).RAM > 0 && result.get(p).CPU > 0 && result.get(p).CPU >= (vms.get(vmi).CPU) && (result.get(p).CPU) - (vms.get(vmi).CPU) >= 0 && result.get(p).RAM >= (vms.get(vmi).RAM) && (result.get(p).RAM) - (vms.get(vmi).RAM) >= 0 && vms.get(vmi).res != 7) {
                                    withoutdublicate.add(vmi);

                                    result.get(p).CPU = result.get(p).CPU - (vms.get(vmi).CPU);
                                    result.get(p).RAM = result.get(p).RAM - (vms.get(vmi).RAM);
                                    result.get(p).reserve = result.get(p).reserve + "\t  App " + o + " type " + ap.get(o).type + ": vm " + vmi + " :" + vms.get(vmi).Name;
                                    result.get(p).count = result.get(p).count + 1;
                                    //System.out.println("VM name " + vms.get(vmi).Name);
                                    vms.get(vmi).res = 7;

                                    System.out.println("vm " + vmi + " servised " + vms.get(vmi).res);
                                    counterapp = counterapp + 1;
                                    vms.get(vmi).pm = result.get(p).Name + "";  /// Or edge Or POD

                                    vms.get(vmi).edge = result.get(p).EdgeSwitchNo;

                                    vms.get(vmi).pod = result.get(p).PodSwitchNo;
                                    apply = true;
                                    System.out.println("pm " + result.get(p).Name);
                                    continue;
                                }

                            }
                        }
                    }
                    // }

                    System.out.println("Application " + o);

                    //}
                } else if (ap.get(o).type.equalsIgnoreCase("Normal")) {
                    System.out.println("Application " + o + " Normal");
                    //pm=sortPMByCPUAsending(pm);
                    result = new ArrayList<PhysicalMachine>();
                    result = pm;
                    System.out.println("check stuts of pms");
                    printPM(result);

                    ArrayList<VirtualMachine> vmlist = ap.get(o).VMList;

                    vms = ap.get(o).VMList;
                    ListOfList.add(vms);
                    int TempShareArray[][] = ap.get(o).weightArray;
                    weightlist.add(TempShareArray);

                    //vmlist = sortVMByCPUDescending(vmlist);
                    for (int k = 0; k < vmlist.size(); k++) {
                        pm = sortPMByCPUAsending(pm);
                        System.out.println("normal app sort pm");
                        printPM(pm);
                        result = new ArrayList<PhysicalMachine>();
                        result = pm;
                        for (int i = 0; i < result.size(); i++) {

                            if (result.get(i).RAM >= 0 && vmlist.get(k).res != 7 && result.get(i).CPU >= 0 && result.get(i).CPU >= vmlist.get(k).CPU && (result.get(i).CPU - vmlist.get(k).CPU) >= 0 && result.get(i).RAM >= vmlist.get(k).RAM && (result.get(i).RAM - vmlist.get(k).RAM) >= 0) {
                                {

                                    result.get(i).CPU = result.get(i).CPU - vmlist.get(k).CPU;
                                    result.get(i).RAM = result.get(i).RAM - vmlist.get(k).RAM;
                                    vmlist.get(k).res = 7;
                                    vmlist.get(k).pm = result.get(i).Name;

                                    result.get(i).reserve = result.get(i).reserve + "\t  App " + o + " type " + ap.get(o).type + ": vm " + k + " :" + vmlist.get(k).Name;
                                    result.get(i).count = result.get(i).count++;

                                    System.out.println("vm " + k + " reserved by " + result.get(i).Name);

                                }
                            }

                        }

                    }

                }

            }

            System.out.println("--------------------------------------------------------------------------------");

            ///////////////////////////////////////////////////
            BufferedReader load = new BufferedReader(new FileReader("Datasets.txt"));

            String filename = load.readLine();
            BufferedReader reader = new BufferedReader(new FileReader(filename));

            //int novm = (int) ((Math.random() * ((8 - 2) + 8)) + 2);;
            String apptype = "";
            int noofvm = 0;

            ArrayList<VirtualMachine> vmlist = new ArrayList<VirtualMachine>();
            List<int[][]> savedweight = new ArrayList<int[][]>();
            int WeightArray[][];

            String str;
            int g = 0;
            //int WeightArray2[][];

            for (int j = 0; j < ap.size(); j++) {
                BufferedReader br3 = new BufferedReader(new FileReader(filename));
                int cc = 0;
                while ((str = br3.readLine()) != null) {

                    if (str.equalsIgnoreCase("Application " + j)) {

                        apptype = br3.readLine();

                        if (apptype.equalsIgnoreCase("Critical") || apptype.equalsIgnoreCase("Normal")) {
                            noofvm = Integer.parseInt(br3.readLine());

                            WeightArray = new int[noofvm][noofvm];
                            try {
                                BufferedReader br4 = new BufferedReader(new FileReader(filename));
                                while ((str = br4.readLine()) != null) {

                                    if (str.equalsIgnoreCase("ConnectionVMApp " + j)) {

                                        int i2 = 0;
                                        while ((!(str = br4.readLine()).equalsIgnoreCase("EndConnectionVMApp " + j))) {

                                            String[] line = str.split(",");

                                            for (int j2 = 0; j2 < line.length; j2++) {
                                                //System.out.println("test  " + Integer.parseInt(line[j2]));
                                                WeightArray[i2][j2] = Integer.parseInt(line[j2]);
                                                //System.out.print(WeightArray[i2][j2] + " t ");
                                            }
                                            i2++;
                                        }

                                        break;
                                    }
                                    //continue;
                                }
                                //System.out.println("");
                                br4.close();
                            } catch (Exception e) {
                                e.printStackTrace();
                            }
                            savedweight.add(WeightArray);

                            break;
                        }
                    }
                    //System.out.println(j);

                    cc++;
                    continue;
                }
                br3.close();
            }

            double resultbandwidth = 0;
            double Criticalbandwidth = 0;
            double NormalBandwidth = 0;
            System.out.println("saved weight size " + savedweight.size());
            for (int i = 0; i < ap.size(); i++) {

                WeightArray = savedweight.get(i);//weightlist.get(i);//sortedApp.get(i).weightArray;//savedweight.get(i);//

                System.out.println("size of list of list " + ListOfList.size());
                vmlist = ListOfList.get(i);

                System.out.println("size vmlist " + vmlist.size());
                for (int j = 0; j < vmlist.size(); j++) {
                    System.out.println("vml " + vmlist.get(j).pm);
                }
                System.out.println("weghit array length " + WeightArray.length);
                System.out.println("///////////////////update bandwidth////////////// application" + i);
                for (int k = 0; k < WeightArray.length; k++) {
                    System.out.println("WeightArray[k].length " + WeightArray[k].length);
                    for (int l = 0; l < WeightArray[k].length; l++) {
                        if (WeightArray[k][l] != 0 && k < l) {

                            String pm1 = vmlist.get(l).pm;
                            System.out.println("pm1 " + pm1);
                            int pmNo1 = 0;
                            int pmNo2 = 0;
                            if (pm1 != null) {
                                //System.out.println("pm1 " + pm1);
                                int y = pm1.indexOf("_");
                                pmNo1 = Integer.parseInt(pm1.substring(y + 1, pm1.length()));

                            }

                            String pm2 = vmlist.get(k).pm;
                            if (pm2 != null) {
                                int x = pm2.indexOf("_");
                                pmNo2 = Integer.parseInt(pm2.substring(x + 1, pm2.length()));
//                               
                            } //System.out.println("pmNo " + pmNo2);
                            if (pm1.equalsIgnoreCase(pm2)) {
                                ;
                            } else if (pm.get(pmNo1 - 1).EdgeSwitchNo == pm.get(pmNo2 - 1).EdgeSwitchNo) {
                                resultbandwidth = resultbandwidth + WeightArray[k][l] * 2;
                            } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo == pm.get(pmNo2 - 1).PodSwitchNo) {
                                resultbandwidth = resultbandwidth + WeightArray[k][l] * 4;
                            } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo != pm.get(pmNo2 - 1).PodSwitchNo) {
                                resultbandwidth = resultbandwidth + WeightArray[k][l] * 6;
                            }

                            if (ap.get(i).type.equalsIgnoreCase("Critical")) {
                                if (pm1.equalsIgnoreCase(pm2)) {
                                    ;
                                } else if (pm.get(pmNo1 - 1).EdgeSwitchNo == pm.get(pmNo2 - 1).EdgeSwitchNo) {
                                    Criticalbandwidth = Criticalbandwidth + WeightArray[k][l] * 2;
                                } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo == pm.get(pmNo2 - 1).PodSwitchNo) {
                                    Criticalbandwidth = Criticalbandwidth + WeightArray[k][l] * 4;
                                } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo != pm.get(pmNo2 - 1).PodSwitchNo) {
                                    Criticalbandwidth = Criticalbandwidth + WeightArray[k][l] * 6;
                                }
                            }
                            if (ap.get(i).type.equalsIgnoreCase("Normal")) {
                                if (pm1.equalsIgnoreCase(pm2)) {
                                    ;
                                } else if (pm.get(pmNo1 - 1).EdgeSwitchNo == pm.get(pmNo2 - 1).EdgeSwitchNo) {
                                    NormalBandwidth = NormalBandwidth + WeightArray[k][l] * 2;
                                } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo == pm.get(pmNo2 - 1).PodSwitchNo) {
                                    NormalBandwidth = NormalBandwidth + WeightArray[k][l] * 4;
                                } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo != pm.get(pmNo2 - 1).PodSwitchNo) {
                                    NormalBandwidth = NormalBandwidth + WeightArray[k][l] * 6;
                                }
                            }

                            System.out.println("costttt " + WeightArray[k][l]);

                        }

                    }

                }
                //graph.printGraph();

            }

            result = sortPMByName(result);

            System.out.println("Result of PAVA Algorith");
            for (int i = 0; i < pm.size(); i++) {
                System.out.println(result.get(i).Name + " CPU =" + result.get(i).CPU + " RAM =" + result.get(i).RAM + " VMS :" + result.get(i).reserve);
                //System.out.println("PM " + i + " RAM " + result.get(i).RAM);
            }
            double resource = findResourceWastage(result, CPUPM, RAMPM);

            System.out.println("------------------------------------Power Consumption--------------------------------------------");
            double SumOfCPUUtilazation = 0;
            double TotalPower = 0;
            double CriticalPower = 0;
            double RamWastage = 0;
            double CpuWastage = 0;
            double noactivepm = 0;
            double powerconsumption = 0;
            for (int i = 0; i < pm.size(); i++) {
                if (CPUPM[i] != result.get(i).CPU) {
                    noactivepm++;
                    double Ucpu = (CPUPM[i] - result.get(i).CPU) / CPUPM[i];

                    System.out.println("CPU " + CPUPM[i] + " Remain  " + result.get(i).CPU + " CPU U " + Ucpu);
                    CpuWastage = CpuWastage + result.get(i).CPU;
                    RamWastage = RamWastage + result.get(i).RAM;

                    System.out.println("reserve " + result.get(i).reserve);
                    double PC = result.get(i).PMIN + (result.get(i).PMAX - result.get(i).PMIN) * Ucpu;
                    if (result.get(i).reserve.contains("Critical")) {
                        CriticalPower = CriticalPower + PC;
                    }
                    TotalPower = TotalPower + PC;
                    System.out.println("CPU utilazation of " + result.get(i).Name + " = " + Ucpu + "\t Power Consumption of " + result.get(i).Name + " = " + PC);
                }

            }
            powerconsumption = TotalPower;/// noactivepm;
            System.out.println("-------------------------------------------------------------------------");
            //graph.countTotalBandwidthGraph();
            System.out.println(Criticalbandwidth);
            System.out.println(NormalBandwidth);

            System.out.println(resultbandwidth);//System.out.println("Total Bandwidth : " + resultbandwidth);

            System.out.println(String.format("%.2f", powerconsumption));//System.out.println("Total Power Consumption : " + String.format("%.2f", powerconsumption));
            //System.out.println("--------------------------Resource Wastage-------------------------------");
            System.out.println(String.format("%.4f", resource));//System.out.println("Total Resource Wastage " + String.format("%.4f", resource) );

        } catch (Exception e) {
            e.printStackTrace();
        }

        int countpm = 0;
        for (int i = 0; i < pm.size(); i++) {
            if (pm.get(i).reserve == "") {
                countpm++;
            }
        }
        ArrayList<VirtualMachine> vms;
        int vmcounter = 0;
        for (int i = 0; i < ap.size(); i++) {
            vms = ap.get(i).VMList;

            for (int j = 0; j < vms.size(); j++) {
                if (vms.get(j).res != 7) {
                    vmcounter++;
                }

            }
        }
        //System.out.println("--------------------------Active Physical Machine-------------------------------");

        System.out.println((pm.size() - countpm));//System.out.println("Number of Active physical machine : " + (pm.size() - countpm));

    }

    public ArrayList<PhysicalMachine> sortPmByEnergy(ArrayList<PhysicalMachine> arraylist) {

        ArrayList<PhysicalMachine> temp = new ArrayList<PhysicalMachine>(arraylist.size());
        try {
            Collections.sort(arraylist, PhysicalMachine.SortEnergyAsending);
            int count = 0;
            for (PhysicalMachine str : arraylist) {
                temp.add(str);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return temp;
    }

    public ArrayList<PhysicalMachine> sortPmByEnergyInEdgeSwitch(ArrayList<PhysicalMachine> arraylist) {

        ArrayList<PhysicalMachine> temp = new ArrayList<PhysicalMachine>(arraylist.size());
        try {
            Collections.sort(arraylist, PhysicalMachine.SortPMEdgeSwitchEnergyAsending);
            int count = 0;
            for (PhysicalMachine str : arraylist) {
                temp.add(str);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return temp;
    }

    public void proposed(ArrayList<PhysicalMachine> pm, ArrayList<Application> ap) {
        double Criticalbandwidth = 0;
        double NormalBandwidth = 0;
        try {
            //double wpm;

            double CPUPM[] = new double[pm.size()];
            double RAMPM[] = new double[pm.size()];
            double PowerEff[] = new double[pm.size()];
            result = pm;
            ////// power efficincy 

            ArrayList<PhysicalMachine> allpmedge = new ArrayList<PhysicalMachine>();
            double pmsize = pm.size();
            double podsize = Math.cbrt(pmsize * 4);
            double NoOfEdgeSwitchs = (podsize / 2) * podsize;
            double NoPmPerEdge = podsize / 2;
            System.out.println(" pod " + podsize);

            double avgPodEnergy[] = new double[(int) podsize];
            double sumEdgeEnergy[] = new double[(int) NoOfEdgeSwitchs];
            double noPODPm = Math.pow((podsize / 2), 2);
            int pcounter = -1;
            int pcounter2 = -1;

            int startEdge = 0;
            int endEdge = 0;

            for (int i = 0; i < CPUPM.length; i++) {
                CPUPM[i] = pm.get(i).CPU;
                RAMPM[i] = pm.get(i).RAM;
                PowerEff[i] = pm.get(i).PowerEfficincy;

                if (i == 0) {
                    startEdge = pm.get(i).EdgeSwitchNo;
                }
                if (i == (pm.size() - 1)) {
                    endEdge = pm.get(i).EdgeSwitchNo;
                }

                //System.out.println("pcounter " + pcounter);
                if (i % noPODPm == 0) {
                    pcounter++;
                    //System.out.println("pcounter   1 " + i + " " + pcounter);

                }
                avgPodEnergy[pcounter] = avgPodEnergy[pcounter] + PowerEff[i];  /// power efficinecy of each POD
                //System.out.println("no port " + NoPod);
                if (i % (podsize / 2) == 0) {
                    pcounter2++;
                    //System.out.println("pcounter   2 " + i + " " + pcounter2);

                }

                sumEdgeEnergy[pcounter2] = sumEdgeEnergy[pcounter2] + PowerEff[i];  // power efficinecy of each edge switch

            }

            /////sort pm per edge ///////////
            int pmnumber = 0;

            for (int i = 0; i < NoOfEdgeSwitchs; i++) {

                ArrayList<PhysicalMachine> pmlisttemp = new ArrayList<PhysicalMachine>();
                for (int count = 0; count < NoPmPerEdge; count++) {
                    pmlisttemp.add(result.get(pmnumber++));
                }
                ArrayList<PhysicalMachine> r = sortPmByEnergy(pmlisttemp);
                for (int j = 0; j < r.size(); j++) {
                    //System.out.println("sort pm by power efficincy "+ r.get(j).Name + " Energy "+r.get(j).PowerEfficincy);
                    allpmedge.add(r.get(j));
                    //allpmedge.clone();
                }
            }

            /////End sort pm per edge ///////////
            /////sort edge per pod ///////////
            int edgenumber = 0;
            ArrayList<PhysicalMachine> sortededge = new ArrayList<PhysicalMachine>();
            double sortedarrayedge[] = new double[(int) NoOfEdgeSwitchs];
            int countsortedge = 0;
            for (int k = 0; k < podsize; k++) {

                //ArrayList<PhysicalMachine> pmlisttemp = new ArrayList<PhysicalMachine>();
                double edgeenergytemp[] = new double[(int) NoPmPerEdge];
                for (int count = 0; count < NoPmPerEdge; count++) {
                    edgeenergytemp[count] = sumEdgeEnergy[edgenumber++];
                }

                double tempsumation[] = edgeenergytemp;
                double sortedge[] = new double[(int) NoPmPerEdge];
                for (int i = 0; i < tempsumation.length; i++) {
                    double mx = tempsumation[0];
                    int ind = 0;
                    for (int j = 0; j < tempsumation.length; j++) {
                        if (tempsumation[j] > mx) {
                            mx = tempsumation[j];
                            //startEdge=startEdge+j;
                            //sortedge[i] = startEdge+j;
                            ind = j;
                        }

                    }
                    sortedge[i] = startEdge + ind;
                    tempsumation[ind] = 0;
                    //startEdge=startEdge-ind;

                }
                startEdge = startEdge + (int) NoPmPerEdge;

                for (int i = 0; i < sortedge.length; i++) {
                    sortedarrayedge[countsortedge++] = sortedge[i];
                }
            }

            ArrayList<PhysicalMachine> sortededgeinsidesamepod = new ArrayList<PhysicalMachine>();

            /////End sort edge per pod ///////////
            /////sort edge inside same pod/////////
            for (int j = 0; j < sortedarrayedge.length; j++) {
                for (int i = 0; i < allpmedge.size(); i++) {

                    if (allpmedge.get(i).EdgeSwitchNo == sortedarrayedge[j]) {
                        sortededgeinsidesamepod.add(allpmedge.get(i));
                    }
                }

            }

            /////sort edge inside same pod/////////
            //////////////////sort energy pod///////////
            int arraysortpod[] = new int[avgPodEnergy.length];
            double tempsum[] = avgPodEnergy;

            for (int i = 0; i < tempsum.length; i++) {
                double mx = 0;
                int ind = 0;
                for (int j = 0; j < tempsum.length; j++) {
                    if (tempsum[j] > mx) {
                        mx = tempsum[j];
                        arraysortpod[i] = j + 1;
                        ind = j;
                    }

                }
                tempsum[ind] = 0;

            }

            //////////////////End sort energy pod///////////
            result = new ArrayList<PhysicalMachine>(pm.size());
            result = sortededgeinsidesamepod;
            ArrayList<VirtualMachine> vms = new ArrayList<VirtualMachine>();
            ArrayList<ArrayList<VirtualMachine>> ListOfList = new ArrayList<ArrayList<VirtualMachine>>(ap.size());
            for (int j = 0; j < ap.size(); j++) {
                int TempShareArray[][] = ap.get(j).weightArray;
                sharearray = ap.get(j).weightArray;

                boolean ch = false;

                vms = ap.get(j).VMList;
                ListOfList.add(vms);
                printListVM(vms);
                //System.out.println(" Size of VMs of application "+j+" = "+vms.size());
                int counterapp = 0;
                int notservised = 0;
                ArrayList<Integer> withoutdublicate = new ArrayList<Integer>();
                ArrayList<Integer> listnotservised = new ArrayList<Integer>();
                {
                    if (ap.get(j).type.equalsIgnoreCase("Critical")) {

                        {
                            ArrayList<VirtualMachine> ar = findVMWithMaxConnection(sharearray, vms);
                            //System.out.println("pinttttttttttttttttttttt listttttttttttttttttttttt");
                            for (int i = 0; i < ar.size(); i++) {
                                System.out.println(i + " = " + ar.get(i).Name);
                            }

                            //new idea for each vm (sort vms by most connection)
                            for (int i = 0; i < ar.size(); i++) {

                                boolean apply = false;

                                if (apply == false) {

                                    for (int y = 0; y < arraysortpod.length; y++) {
                                        int ind = 0;
                                        for (int t = 0; t < sortededgeinsidesamepod.size(); t++) {
                                            if (sortededgeinsidesamepod.get(t).PodSwitchNo == arraysortpod[y]) {
                                                ind = t;
                                                System.out.println("ind " + ind);
                                                break;
                                            }

                                        }
                                        result = sortPmByEnergyInEdgeSwitch(result);
                                        int pmcounter = 0;

                                        for (int p = ind; pmcounter < podsize; p++) {
                                            pmcounter++;

                                            if (result.get(p).RAM >= 0 && result.get(p).CPU >= 0 && result.get(p).CPU >= (ar.get(i).CPU) && (result.get(p).CPU) - (ar.get(i).CPU) >= 0 && result.get(p).RAM >= (ar.get(i).RAM) && (result.get(p).RAM) - (ar.get(i).RAM) >= 0 && ar.get(i).res != 5) {
                                                withoutdublicate.add(i);

                                                result.get(p).CPU = result.get(p).CPU - (ar.get(i).CPU);
                                                result.get(p).RAM = result.get(p).RAM - (ar.get(i).RAM);
                                                result.get(p).reserve = result.get(p).reserve + "\t  App " + j + " type " + ap.get(j).type + ": vm " + i + " :" + ar.get(i).Name;
                                                result.get(p).count = result.get(p).count + 1;
                                                //System.out.println("VM name " + vms.get(r).Name);
                                                ar.get(i).res = 5;
                                                //vms.get(c).res = 5;
                                                System.out.println("vm " + ar.get(i).Name + " servised " + ar.get(i).res);
                                                //System.out.println("vm " + c + " servised " + vms.get(c).res);
                                                counterapp = counterapp + 1;
                                                ar.get(i).pm = result.get(p).Name + "";  /// Or edge Or POD
                                                //vms.get(c).pm = result.get(p).Name + "";

                                                ar.get(i).edge = result.get(p).EdgeSwitchNo;
                                                //ar.get(i).edge = result.get(p).EdgeSwitchNo;
                                                ar.get(i).pod = result.get(p).PodSwitchNo;
                                                //vms.get(r).pod = result.get(p).PodSwitchNo;
                                                apply = true;
                                                System.out.println("pm " + result.get(p).Name + " vms : " + result.get(p).reserve);
                                                continue;
                                            }

                                        }
                                    }
                                }

                            }

                        }
                    } else if (ap.get(j).type.equalsIgnoreCase("Normal")) {
                        while ((withoutdublicate.size() + listnotservised.size()) != vms.size()) {
                            for (int i = 0; i < vms.size(); i++) {
                                if (vms.get(i).res != 5) {
                                    //System.out.println("Both of vms not servised");
                                    boolean apply = false;
                                    result = sortPMByName(result);

                                    {
                                        boolean b = placmentNormalApplication(arraysortpod, sortedarrayedge, result, vms, i, CPUPM, RAMPM, j, ap);

                                        if (b == true) {
                                            counterapp = counterapp + 1;
                                            withoutdublicate.add(i);
                                            apply = true;

                                        }
                                    }
                                }
                            }
                        }
                    }

                }

                System.out.println("Application " + j);

                //////////////22222222222222222222222222222222222///////////////////    
            }

            for (int i = 0; i < ap.size(); i++) {
                ArrayList<VirtualMachine> vv = app.get(i).VMList;
                for (int j = 0; j < vv.size(); j++) {
                    if (vv.get(j).res != 5) {
                        System.out.println("vm not placed " + vv.get(j).Name);
                    }
                }

            }
            double resultbandwidth = 0;

            System.out.println("--------------------------------------------------------------------------------");

            //ap = sortAppByType(ap);
            for (int i = 0; i < ap.size(); i++) {

                //int cost[][] = ap.get(i).weightArray;
                ///////////////////////////////////////////////////
                BufferedReader load = new BufferedReader(new FileReader("Datasets.txt"));

                String filename = load.readLine();
                BufferedReader reader = new BufferedReader(new FileReader(filename));

                //int novm = (int) ((Math.random() * ((8 - 2) + 8)) + 2);;
                String apptype = "";
                int noofvm = 0;

                ArrayList<VirtualMachine> vmlist = new ArrayList<VirtualMachine>();
                //int WeightArray[][];

                String str;
                int g = 0;
                int WeightArray[][];
                List<int[][]> savedweight = new ArrayList<int[][]>();

                for (int j = 0; j < ap.size(); j++) {
                    BufferedReader br3 = new BufferedReader(new FileReader(filename));
                    int cc = 0;
                    while ((str = br3.readLine()) != null) {
                        if (str.equalsIgnoreCase("Application " + j)) {

                            apptype = br3.readLine();

                            {
                                noofvm = Integer.parseInt(br3.readLine());

                                WeightArray = new int[noofvm][noofvm];
                                try {
                                    BufferedReader br4 = new BufferedReader(new FileReader(filename));
                                    while ((str = br4.readLine()) != null) {

                                        if (str.equalsIgnoreCase("ConnectionVMApp " + j)) {

                                            int i2 = 0;
                                            while ((!(str = br4.readLine()).equalsIgnoreCase("EndConnectionVMApp " + j))) {

                                                String[] line = str.split(",");

                                                for (int j2 = 0; j2 < line.length; j2++) {
                                                    //System.out.println("test  " + Integer.parseInt(line[j2]));
                                                    WeightArray[i2][j2] = Integer.parseInt(line[j2]);
                                                    //System.out.print(WeightArray[i2][j2] + " t ");
                                                }
                                                i2++;
                                            }

                                            break;
                                        }
                                        //continue;
                                    }
                                    //System.out.println("");
                                    br4.close();
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }
                                savedweight.add(WeightArray);

                                break;
                            }
                        }
                        //System.out.println(j);

                        cc++;
                        continue;
                    }
                    br3.close();
                }
                ////////////////////////////////////////end load weight//////////////////////////////////////////
                WeightArray = savedweight.get(i);
                vmlist = ListOfList.get(i);
                System.out.println("size vmlist " + vmlist.size());
                System.out.println("///////////////////update bandwidth////////////// application" + i);
                if (vmlist.size() == 1) {
                    ;
                } //ArrayList<VirtualMachine> vmlist = ap.get(i).VMList;
                else {
                    for (int k = 0; k < WeightArray.length; k++) {
                        //System.out.println("WeightArray[k].length " + WeightArray[k].length);
                        for (int l = 0; l < WeightArray[k].length; l++) {
                            if (WeightArray[k][l] != 0 && k < l) {

                                String pm1 = vmlist.get(l).pm;

                                int pmNo1 = 0;
                                int pmNo2 = 0;
                                if (pm1 != null) {
                                    //System.out.println("pm1 " + pm1);
                                    int y = pm1.indexOf("_");
                                    pmNo1 = Integer.parseInt(pm1.substring(y + 1, pm1.length()));

                                }

                                String pm2 = vmlist.get(k).pm;
                                if (pm2 != null) {
                                    int x = pm2.indexOf("_");
                                    pmNo2 = Integer.parseInt(pm2.substring(x + 1, pm2.length()));

                                }

                                if (pm1.equalsIgnoreCase(pm2)) {
                                    ;
                                } else if (pm1 == null && pm2 == null) {
                                    ;
                                } else if (pm.get(pmNo1 - 1).EdgeSwitchNo == pm.get(pmNo2 - 1).EdgeSwitchNo) {
                                    resultbandwidth = resultbandwidth + WeightArray[k][l] * 2;
                                } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo == pm.get(pmNo2 - 1).PodSwitchNo) {
                                    resultbandwidth = resultbandwidth + WeightArray[k][l] * 4;
                                } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo != pm.get(pmNo2 - 1).PodSwitchNo) {
                                    resultbandwidth = resultbandwidth + WeightArray[k][l] * 6;
                                }
                                if (ap.get(i).type.equalsIgnoreCase("Critical")) {
                                    if (pm1.equalsIgnoreCase(pm2)) {
                                        ;
                                    } else if (pm.get(pmNo1 - 1).EdgeSwitchNo == pm.get(pmNo2 - 1).EdgeSwitchNo) {
                                        Criticalbandwidth = Criticalbandwidth + WeightArray[k][l] * 2;
                                    } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo == pm.get(pmNo2 - 1).PodSwitchNo) {
                                        Criticalbandwidth = Criticalbandwidth + WeightArray[k][l] * 4;
                                    } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo != pm.get(pmNo2 - 1).PodSwitchNo) {
                                        Criticalbandwidth = Criticalbandwidth + WeightArray[k][l] * 6;
                                    }
                                }
                                if (ap.get(i).type.equalsIgnoreCase("Normal")) {
                                    if (pm1.equalsIgnoreCase(pm2)) {
                                        ;
                                    } else if (pm.get(pmNo1 - 1).EdgeSwitchNo == pm.get(pmNo2 - 1).EdgeSwitchNo) {
                                        NormalBandwidth = NormalBandwidth + WeightArray[k][l] * 2;
                                    } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo == pm.get(pmNo2 - 1).PodSwitchNo) {
                                        NormalBandwidth = NormalBandwidth + WeightArray[k][l] * 4;
                                    } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo != pm.get(pmNo2 - 1).PodSwitchNo) {
                                        NormalBandwidth = NormalBandwidth + WeightArray[k][l] * 6;
                                    }
                                }

                            }

                        }

                    }
                }
            }
            result = sortPMByName(result);
            System.out.println("Result of our Method");
            for (int i = 0; i < pm.size(); i++) {
                System.out.println(result.get(i).Name + " CPU =" + result.get(i).CPU + " RAM =" + result.get(i).RAM + " VMS :" + result.get(i).reserve);
                //System.out.println("PM " + i + " RAM " + result.get(i).RAM);
            }

            double resource = findResourceWastage(result, CPUPM, RAMPM);

            System.out.println("------------------------------------Power Consumption--------------------------------------------");
            double SumOfCPUUtilazation = 0;
            double TotalPower = 0;
            double CriticalPower = 0;
            double RamWastage = 0;
            double CpuWastage = 0;
            double noactivepm = 0;
            double powerconsumption = 0;
            //double resource = findResourceWastage(result, CPUPM, RAMPM);
            for (int i = 0; i < pm.size(); i++) {
                if (CPUPM[i] != result.get(i).CPU) {
                    noactivepm++;

                    double Ucpu = (CPUPM[i] - result.get(i).CPU) / CPUPM[i];
                    double PC = result.get(i).PMIN + (result.get(i).PMAX - result.get(i).PMIN) * Ucpu;
                    if (result.get(i).reserve.contains("Critical")) {
                        CriticalPower = CriticalPower + PC;
                    }
                    TotalPower = TotalPower + PC;
                    System.out.println("CPU utilazation of " + result.get(i).Name + " = " + Ucpu + "\t Power Consumption of " + result.get(i).Name + " = " + PC);
                }

            }
            powerconsumption = TotalPower;/// noactivepm;
            System.out.println("-------------------------------------------------------------------------");
            System.out.println(Criticalbandwidth);
            System.out.println(NormalBandwidth);
            System.out.println(resultbandwidth);//System.out.println("Total Bandwidth : " + resultbandwidth);
            System.out.println(String.format("%.2f", powerconsumption));//System.out.println("Total Power Consumption : " + String.format("%.2f", powerconsumption));
            //System.out.println("--------------------------Resource Wastage-------------------------------");
            System.out.println(String.format("%.4f", resource));//System.out.println("Total Resource Wasstage : " +  String.format("%.4f", resource));
        } catch (Exception e) {
            e.printStackTrace();
        }

        int countpm = 0;
        for (int i = 0; i < pm.size(); i++) {
            if (pm.get(i).reserve == "") {
                countpm++;
            }
        }
        ArrayList<VirtualMachine> vms;
        int vmcounter = 0;
        for (int i = 0; i < ap.size(); i++) {
            vms = ap.get(i).VMList;

            for (int j = 0; j < vms.size(); j++) {
                if (vms.get(j).res != 5) {
                    vmcounter++;
                }

            }
        }
        //System.out.println("-------------------------------------------------------");
        System.out.println((pm.size() - countpm));//System.out.println("Number of Active physical machine : " + (pm.size() - countpm));

    }

    public void placementComputationFirstFitDecreasing(ArrayList<PhysicalMachine> pm, ArrayList<Application> ap) {

        result = new ArrayList<PhysicalMachine>(pm.size());

        result = pm;

        ///////////////power////////////////
        double CPUPM[] = new double[pm.size()];
        double PowerEff[] = new double[pm.size()];
        //double CPUPM[] = new double[pm.size()];
        double RAMPM[] = new double[pm.size()];
        //double PowerEff[] = new double[pm.size()];
        for (int i = 0; i < CPUPM.length; i++) {
            CPUPM[i] = pm.get(i).CPU;
            RAMPM[i] = pm.get(i).RAM;
            PowerEff[i] = pm.get(i).PowerEfficincy;
        }

        double pmsize = pm.size();
        double podsize = Math.cbrt(pmsize * 4);
        System.out.println(" pod " + podsize);

        double avgPodEnergy[] = new double[(int) podsize];
        double noPODPm = Math.pow((podsize / 2), 2);
        int pcounter = -1;
        for (int i = 0; i < CPUPM.length; i++) {
            CPUPM[i] = pm.get(i).CPU;
            PowerEff[i] = pm.get(i).PowerEfficincy;

            System.out.println("pcounter " + pcounter);
            if (i % noPODPm == 0) {
                pcounter++;
                System.out.println("iiii " + i + " " + pcounter);

            }
            avgPodEnergy[pcounter] = avgPodEnergy[pcounter] + PowerEff[i];

        }

        /////power////
        for (int j = 0; j < ap.size(); j++) {

            ArrayList<VirtualMachine> vmlist = ap.get(j).VMList;

            vmlist = sortVMByCPUDescending(vmlist);
            for (int k = 0; k < vmlist.size(); k++) {
                System.out.println("sort VMlist by CPU " + vmlist.get(k).CPU);
            }
            System.out.println("------------");

            for (int k = 0; k < vmlist.size(); k++) {
                for (int i = 0; i < pm.size(); i++) {

                    if (result.get(i).RAM >= 0 && vmlist.get(k).res != 2 && result.get(i).CPU >= 0 && result.get(i).CPU >= vmlist.get(k).CPU && (result.get(i).CPU - vmlist.get(k).CPU) >= 0 && result.get(i).RAM >= vmlist.get(k).RAM && (result.get(i).RAM - vmlist.get(k).RAM) >= 0) {
                        {
                            result.get(i).CPU = result.get(i).CPU - vmlist.get(k).CPU;
                            result.get(i).RAM = result.get(i).RAM - vmlist.get(k).RAM;
                            vmlist.get(k).res = 2;
                            vmlist.get(k).pm = result.get(i).Name;
                            result.get(i).reserve = result.get(i).reserve + " " + " ";

                            result.get(i).reserve = result.get(i).reserve + "\t  App " + j + " type " + ap.get(j).type + ": vm " + k + " :" + vmlist.get(k).Name;
                            result.get(i).count = result.get(i).count++;
                            continue;
                        }
                    }

                }

            }
        }

        System.out.println("--------------------------------------------------------------------------------");
        double resultbandwidth = 0;
        double Criticalbandwidth = 0;
        double NormalBandwidth = 0;
        for (int i = 0; i < ap.size(); i++) {

            int cost[][] = ap.get(i).weightArray;
            ArrayList<VirtualMachine> vmlist = ap.get(i).VMList;
            vmlist = sortVMByCPUDescending(vmlist);
            System.out.println("vmlist size " + vmlist.size());
            System.out.println("///////////////////update bandwidth////////////// application" + i);
            for (int k = 0; k < cost.length; k++) {

                for (int l = 0; l < cost[k].length; l++) {
                    if (cost[k][l] != 0 && k < l) {
                        int pmNo1 = 0;
                        int pmNo2 = 0;
                        String pm1 = vmlist.get(k).pm;
                        //if (pm1 != null) 
                        {
                            int y = pm1.indexOf("_");
                            pmNo1 = Integer.parseInt(pm1.substring(y + 1, pm1.length()));

                        }

                        String pm2 = vmlist.get(l).pm;
                        {
                            int x = pm2.indexOf("_");
                            pmNo2 = Integer.parseInt(pm2.substring(x + 1, pm2.length()));
                        }
                        //System.out.println("pmNo " + pmNo2);
                        if (pm1.equalsIgnoreCase(pm2)) {
                            ;
                        } else if (pm.get(pmNo1 - 1).EdgeSwitchNo == pm.get(pmNo2 - 1).EdgeSwitchNo) {
                            resultbandwidth = resultbandwidth + cost[k][l] * 2;
                        } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo == pm.get(pmNo2 - 1).PodSwitchNo) {
                            resultbandwidth = resultbandwidth + cost[k][l] * 4;
                        } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo != pm.get(pmNo2 - 1).PodSwitchNo) {
                            resultbandwidth = resultbandwidth + cost[k][l] * 6;
                        }
                        if (ap.get(i).type.equalsIgnoreCase("Critical")) {
                            if (pm1.equalsIgnoreCase(pm2)) {
                                ;
                            } else if (pm.get(pmNo1 - 1).EdgeSwitchNo == pm.get(pmNo2 - 1).EdgeSwitchNo) {
                                Criticalbandwidth = Criticalbandwidth + cost[k][l] * 2;
                            } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo == pm.get(pmNo2 - 1).PodSwitchNo) {
                                Criticalbandwidth = Criticalbandwidth + cost[k][l] * 4;
                            } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo != pm.get(pmNo2 - 1).PodSwitchNo) {
                                Criticalbandwidth = Criticalbandwidth + cost[k][l] * 6;
                            }
                        }
                        if (ap.get(i).type.equalsIgnoreCase("Normal")) {
                            if (pm1.equalsIgnoreCase(pm2)) {
                                ;
                            } else if (pm.get(pmNo1 - 1).EdgeSwitchNo == pm.get(pmNo2 - 1).EdgeSwitchNo) {
                                NormalBandwidth = NormalBandwidth + cost[k][l] * 2;
                            } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo == pm.get(pmNo2 - 1).PodSwitchNo) {
                                NormalBandwidth = NormalBandwidth + cost[k][l] * 4;
                            } else if (pm.get(pmNo1 - 1).EdgeSwitchNo != pm.get(pmNo2 - 1).EdgeSwitchNo && pm.get(pmNo1 - 1).PodSwitchNo != pm.get(pmNo2 - 1).PodSwitchNo) {
                                NormalBandwidth = NormalBandwidth + cost[k][l] * 6;
                            }
                        }
                    }

                }

            }
        }
        result = sortPMByName(result);
        System.out.println("--------------------------------------------------------------------------------");
        System.out.println("Result of First Fit Decreasing Algorithm");
        for (int i = 0; i < result.size(); i++) {
            System.out.println("PM " + result.get(i).Name + " \t CPU " + result.get(i).CPU + " RAM " + result.get(i).RAM + " VMS :" + result.get(i).reserve);
            //System.out.println("PM " + i + " RAM " + result.get(i).RAM);
        }

        double resource = findResourceWastage(result, CPUPM, RAMPM);

        System.out.println("------------------------------------Power Consumption--------------------------------------------");
        double SumOfCPUUtilazation = 0;
        double TotalPower = 0;
        double CriticalPower = 0;
        double RamWastage = 0;
        double CpuWastage = 0;
        double noactivepm = 0;
        double powerconsumption = 0;
        for (int i = 0; i < pm.size(); i++) {
            if (CPUPM[i] != result.get(i).CPU) {
                noactivepm++;
                double Ucpu = (CPUPM[i] - result.get(i).CPU) / CPUPM[i];

                System.out.println("CPU " + CPUPM[i] + " Remain  " + result.get(i).CPU + " CPU U " + Ucpu);
                CpuWastage = CpuWastage + result.get(i).CPU;
                RamWastage = RamWastage + result.get(i).RAM;
                System.out.println("reserve " + result.get(i).reserve);
                double PC = result.get(i).PMIN + (result.get(i).PMAX - result.get(i).PMIN) * Ucpu;
                if (result.get(i).reserve.contains("Critical")) {
                    CriticalPower = CriticalPower + PC;
                }
                TotalPower = TotalPower + PC;
                System.out.println("CPU utilazation of " + result.get(i).Name + " = " + Ucpu + "\t Power Consumption of " + result.get(i).Name + " = " + PC);
            }

        }
        powerconsumption = TotalPower;/// noactivepm;
        System.out.println("-------------------------------------------------------------------------");
        //graph.countTotalBandwidthGraph();

        System.out.println(Criticalbandwidth);
        System.out.println(NormalBandwidth);
        System.out.println(resultbandwidth);

        //System.out.println("Power Consumption for Critical app : " + String.format("%.2f", CriticalPower));
        System.out.println(String.format("%.2f", powerconsumption));
        System.out.println(String.format("%.4f", resource));

        int countpm = 0;
        for (int i = 0; i < pm.size(); i++) {
            if (pm.get(i).reserve == "") {
                countpm++;
            }
        }
        ArrayList<VirtualMachine> vms;
        int vmcounter = 0;
        for (int i = 0; i < app.size(); i++) {
            vms = app.get(i).VMList;

            for (int j = 0; j < vms.size(); j++) {
                if (vms.get(j).res != 2) {
                    vmcounter++;
                }

            }
        }
        //System.out.println("-------------------------------------------------------");
        System.out.println((pm.size() - countpm));
    }

    public int[] RandomizeArray(int[] array) {
        Random rgen = new Random(); // Random number generator
        for (int i = 0; i < array.length; i++) {
            int randomPosition = rgen.nextInt(array.length);
            int temp = array[i];
            array[i] = array[randomPosition];
            array[randomPosition] = temp;
        }
        return array;
    }

    public void EQVMP(ArrayList<PhysicalMachine> pmss, ArrayList<Application> ap) {
        // Stores block id of the block allocated to a 
        // process 
        ArrayList<VirtualMachine> vmlist = new ArrayList<VirtualMachine>();
        ArrayList<PhysicalMachine> result = new ArrayList<PhysicalMachine>(pmss.size());
        result = pmss;

        double CPUPM[] = new double[pmss.size()];
        double PowerEff[] = new double[pmss.size()];

        //double CPUPM[] = new double[pm.size()];
        double RAMPM[] = new double[pmss.size()];
        //double PowerEff[] = new double[pm.size()];
        for (int i = 0; i < CPUPM.length; i++) {
            CPUPM[i] = pmss.get(i).CPU;
            RAMPM[i] = pmss.get(i).RAM;
            PowerEff[i] = pmss.get(i).PowerEfficincy;
        }

        double pmsize = pmss.size();
        double podsize = Math.cbrt(pmsize * 4);
        System.out.println(" pod " + podsize);

        double avgPodEnergy[] = new double[(int) podsize];
        double noPODPm = Math.pow((podsize / 2), 2);
        int pcounter = -1;
        for (int i = 0; i < CPUPM.length; i++) {
            CPUPM[i] = pmss.get(i).CPU;
            PowerEff[i] = pmss.get(i).PowerEfficincy;

            //System.out.println("pcounter " + pcounter);
            if (i % noPODPm == 0) {
                pcounter++;
                //System.out.println("iiii " + i + " " + pcounter);

            }
            avgPodEnergy[pcounter] = avgPodEnergy[pcounter] + PowerEff[i];

        }

        /////////////////power/////
        // Initially 
        for (int a = 0; a < ap.size(); a++) {
            vmlist = ap.get(a).VMList;
            int cost[][] = ap.get(a).weightArray;

            for (int i = 0; i < cost.length; i++) {
                for (int j = 0; j < cost[i].length; j++) {
                    cost[j][i] = cost[i][j];
                }
            }
            int m = pmss.size();
            int n = vmlist.size();
            int perfectarray1[];
            int perfectarray2[];
            int lengthpart1 = 0;
            int finalcost = Integer.MAX_VALUE;

            ArrayList<VirtualMachine> vmlist1 = new ArrayList<VirtualMachine>();
            ArrayList<VirtualMachine> vmlist2 = new ArrayList<VirtualMachine>();

            int array[] = new int[n];
            for (int i = 0; i < array.length; i++) {
                array[i] = i;
            }
            double nn = (double) n / 2;
            double counter = 0;
            //double counter = Math.pow(nn, n) * n * n;

            if (n < 8) {
                counter = Math.pow(nn, n) * n * n;
            } else {
                counter = Math.pow(n, 4) * n * n;
            }
            for (int fac = 0; fac < counter; fac++) {
                int minimumecost = 0;
                /////vm parttion
                int arrayre[] = RandomizeArray(array);
                int len = arrayre.length;
                int mid = len / 2;
                int part1[] = new int[mid];
                int part2[] = new int[arrayre.length - mid];
                int indexpart2 = 0;

                for (int i = 0; i < arrayre.length; i++) {
                    if (i < part1.length) {
                        part1[i] = arrayre[i];
                    } else {
                        part2[indexpart2++] = arrayre[i];
                    }
                }

                for (int i = 0; i < part1.length; i++) {
                    for (int j = 0; j < part2.length; j++) {
                        //System.out.println("cost "+part1[i] +" "+part2[j]+" = "+cost[part1[i]][part2[j]]);
                        minimumecost = minimumecost + cost[part1[i]][part2[j]];
                    }
                }
                if (minimumecost < finalcost) {
                    finalcost = minimumecost;

                    perfectarray1 = new int[part1.length];
                    perfectarray1 = part1;
                    vmlist1.clear();
                    vmlist2.clear();

                    for (int i = 0; i < perfectarray1.length; i++) {
                        vmlist1.add(vmlist.get(perfectarray1[i]));
                    }
                    //lengthpart1=part1.length;

                    perfectarray2 = new int[part2.length];
                    perfectarray2 = part2;
                    for (int i = 0; i < perfectarray2.length; i++) {
                        vmlist2.add(vmlist.get(perfectarray2[i]));
                    }
                }
            }
            for (int i = 0; i < vmlist1.size(); i++) {
                System.out.println("list 1 " + vmlist1.get(i).Name);
            }
            for (int i = 0; i < vmlist2.size(); i++) {
                System.out.println("list 2 " + vmlist2.get(i).Name);
            }
            System.out.println("application " + a);

            ////ta era best hop reduction
            System.out.println("Sortttttttttttttttttttttttt Sublist 1");

            vmlist1 = sortVMByResourceDecreascing(vmlist1);
            for (int i = 0; i < vmlist1.size(); i++) {
                System.out.println(vmlist1.get(i).Name);
            }
            System.out.println("Sortttttttttttttttttttttttt Sublist 2");
            vmlist2 = sortVMByResourceDecreascing(vmlist2);
            for (int i = 0; i < vmlist2.size(); i++) {
                System.out.println(vmlist2.get(i).Name);
            }
            for (int i = 0; i < result.size(); i++) {

                for (int t = 0; t < vmlist1.size(); t++) {
                    int perfectindex = -1;
                    for (int j = 0; j < vmlist1.size(); j++) {
                        if (result.get(i).CPU >= vmlist1.get(j).CPU && result.get(i).RAM >= vmlist1.get(j).RAM && vmlist1.get(j).res != 8) {
                            if (perfectindex == -1) {
                                perfectindex = j;
                            } else if (vmlist1.get(perfectindex).CPU < vmlist1.get(j).CPU) {
                                perfectindex = j;
                            }
                        }
                    }
                    if (perfectindex != -1) {

                        result.get(i).CPU -= vmlist1.get(perfectindex).CPU;
                        result.get(i).RAM -= vmlist1.get(perfectindex).RAM;
                        vmlist1.get(perfectindex).res = 8;
                        //result.get(i).reserve = result.get(i).reserve + "   " + vmlist1.get(perfectindex).Name;
                        result.get(i).reserve = result.get(i).reserve + "\t  App " + a + " type " + ap.get(a).type + ": vm " + perfectindex + " :" + vmlist1.get(perfectindex).Name;
                        result.get(i).count = result.get(i).count++;

                        vmlist1.get(perfectindex).pm = result.get(i).Name;
                    }
                }
            }
            for (int i = 0; i < result.size(); i++) {
                //////
                for (int t = 0; t < vmlist2.size(); t++) {
                    int perfectindex = -1;
                    for (int j = 0; j < vmlist2.size(); j++) {
                        if (result.get(i).CPU >= vmlist2.get(j).CPU && result.get(i).RAM >= vmlist2.get(j).RAM && vmlist2.get(j).res != 8) {
                            if (perfectindex == -1) {
                                perfectindex = j;
                            } else if (vmlist2.get(perfectindex).CPU < vmlist2.get(j).CPU) {
                                perfectindex = j;
                            }
                        }
                    }
                    if (perfectindex != -1) {

                        result.get(i).CPU -= vmlist2.get(perfectindex).CPU;
                        result.get(i).RAM -= vmlist2.get(perfectindex).RAM;
                        vmlist2.get(perfectindex).res = 8;
                        //result.get(i).reserve = result.get(i).reserve + "   " + vmlist2.get(perfectindex).Name;
                        result.get(i).reserve = result.get(i).reserve + "\t  App " + a + " type " + ap.get(a).type + ": vm " + perfectindex + " :" + vmlist2.get(perfectindex).Name;
                        result.get(i).count = result.get(i).count++;

                        vmlist2.get(perfectindex).pm = result.get(i).Name;
                    }
                }

            }
        }

        System.out.println("--------------------------------------------------------------------------------");
        double resultbandwidth = 0;
        double Criticalbandwidth = 0;
        double NormalBandwidth = 0;

        for (int i = 0; i < ap.size(); i++) {

            int cost[][] = ap.get(i).weightArray;
            vmlist = ap.get(i).VMList;

            System.out.println("///////////////////update bandwidth////////////// application" + i);
            for (int k = 0; k < cost.length; k++) {

                for (int l = 0; l < cost[k].length; l++) {
                    if (cost[k][l] != 0 && k < l) {

                        String pm1 = vmlist.get(k).pm;
                        //System.out.println("pm1 " + pm1);
                        int y = pm1.indexOf("_");
                        int pmNo1 = Integer.parseInt(pm1.substring(y + 1, pm1.length()));

                        String pm2 = vmlist.get(l).pm;
                        int x = pm2.indexOf("_");
                        int pmNo2 = Integer.parseInt(pm2.substring(x + 1, pm2.length()));

                        if (pm1.equalsIgnoreCase(pm2)) {
                            ;
                        } else if (result.get(pmNo1 - 1).EdgeSwitchNo == result.get(pmNo2 - 1).EdgeSwitchNo) {
                            resultbandwidth = resultbandwidth + cost[k][l] * 2;
                        } else if (result.get(pmNo1 - 1).EdgeSwitchNo != result.get(pmNo2 - 1).EdgeSwitchNo && result.get(pmNo1 - 1).PodSwitchNo == result.get(pmNo2 - 1).PodSwitchNo) {
                            resultbandwidth = resultbandwidth + cost[k][l] * 4;
                        } else if (result.get(pmNo1 - 1).EdgeSwitchNo != result.get(pmNo2 - 1).EdgeSwitchNo && result.get(pmNo1 - 1).PodSwitchNo != result.get(pmNo2 - 1).PodSwitchNo) {
                            resultbandwidth = resultbandwidth + cost[k][l] * 6;
                        }

                        if (ap.get(i).type.equalsIgnoreCase("Critical")) {
                            if (pm1.equalsIgnoreCase(pm2)) {
                                ;
                            } else if (result.get(pmNo1 - 1).EdgeSwitchNo == result.get(pmNo2 - 1).EdgeSwitchNo) {
                                Criticalbandwidth = Criticalbandwidth + cost[k][l] * 2;
                            } else if (result.get(pmNo1 - 1).EdgeSwitchNo != result.get(pmNo2 - 1).EdgeSwitchNo && result.get(pmNo1 - 1).PodSwitchNo == result.get(pmNo2 - 1).PodSwitchNo) {
                                Criticalbandwidth = Criticalbandwidth + cost[k][l] * 4;
                            } else if (result.get(pmNo1 - 1).EdgeSwitchNo != result.get(pmNo2 - 1).EdgeSwitchNo && result.get(pmNo1 - 1).PodSwitchNo != result.get(pmNo2 - 1).PodSwitchNo) {
                                Criticalbandwidth = Criticalbandwidth + cost[k][l] * 6;
                            }
                        }

                        if (ap.get(i).type.equalsIgnoreCase("Normal")) {
                            if (pm1.equalsIgnoreCase(pm2)) {
                                ;
                            } else if (result.get(pmNo1 - 1).EdgeSwitchNo == result.get(pmNo2 - 1).EdgeSwitchNo) {
                                NormalBandwidth = NormalBandwidth + cost[k][l] * 2;
                            } else if (result.get(pmNo1 - 1).EdgeSwitchNo != result.get(pmNo2 - 1).EdgeSwitchNo && result.get(pmNo1 - 1).PodSwitchNo == result.get(pmNo2 - 1).PodSwitchNo) {
                                NormalBandwidth = NormalBandwidth + cost[k][l] * 4;
                            } else if (result.get(pmNo1 - 1).EdgeSwitchNo != result.get(pmNo2 - 1).EdgeSwitchNo && result.get(pmNo1 - 1).PodSwitchNo != result.get(pmNo2 - 1).PodSwitchNo) {
                                NormalBandwidth = NormalBandwidth + cost[k][l] * 6;
                            }
                        }
                    }

                }

            }
        }
        result = sortPMByName(result);
        System.out.println(" Result of EQVMP Algorithm : ");
        for (int i = 0; i < result.size(); i++) {
            System.out.println("Name : " + result.get(i).Name + " CPU " + result.get(i).CPU + " RAM " + result.get(i).RAM + " VMS: " + result.get(i).reserve);
        }
        double resource = findResourceWastage(result, CPUPM, RAMPM);

        System.out.println("------------------------------------Power Consumption--------------------------------------------");
        double SumOfCPUUtilazation = 0;
        double TotalPower = 0;
        double CriticalPower = 0;
        double RamWastage = 0;
        double CpuWastage = 0;
        double noactivepm = 0;
        double powerconsumption = 0;
        for (int i = 0; i < result.size(); i++) {
            if (CPUPM[i] != result.get(i).CPU) {
                noactivepm++;
                double Ucpu = (CPUPM[i] - result.get(i).CPU) / CPUPM[i];

                System.out.println("CPU " + CPUPM[i] + " Remain  " + result.get(i).CPU + " CPU U " + Ucpu);
                CpuWastage = CpuWastage + result.get(i).CPU;
                RamWastage = RamWastage + result.get(i).RAM;
                System.out.println("reserve " + result.get(i).reserve);
                double PC = result.get(i).PMIN + (result.get(i).PMAX - result.get(i).PMIN) * Ucpu;
                if (result.get(i).reserve.contains("Critical")) {
                    CriticalPower = CriticalPower + PC;
                }
                TotalPower = TotalPower + PC;
                System.out.println("CPU utilazation of " + result.get(i).Name + " = " + Ucpu + "\t Power Consumption of " + result.get(i).Name + " = " + PC);
            }

        }
        powerconsumption = TotalPower;/// noactivepm;
        System.out.println("-------------------------------------------------------------------------");

        System.out.println(Criticalbandwidth);
        System.out.println(NormalBandwidth);
        System.out.println(resultbandwidth);
        System.out.println(String.format("%.2f", powerconsumption));
        System.out.println(String.format("%.4f", resource));
        //System.out.println("------------------------------------END Power Consumption--------------------------------------------");
        int countpm = 0;
        for (int i = 0; i < result.size(); i++) {
            if (result.get(i).reserve == "") {
                countpm++;
            }
        }
        int vmcounter = 0;
        for (int i = 0; i < app.size(); i++) {
            vmlist = app.get(i).VMList;

            for (int j = 0; j < vmlist.size(); j++) {
                if (vmlist.get(j).res != 8) {
                    vmcounter++;
                }

            }
        }
        System.out.println((result.size() - countpm));

    }

    static List<Integer> dist;

    static List<Integer> width;

    public static int randomNumberGenerator() {

        Random r = new Random();
        int low = 10;
        int high = 11;
        int createdRanNum = r.nextInt(high - low) + low;

        return (createdRanNum);
    }

    static double no = 0;
    static int NoPod = 0;
    static int NoPort = 0;
    static double NoCoreSwitch = 0;
    static int NoAggSwitch = 0;
    Graph graph;

    public static void main(String[] args) {
        Computation ob = new Computation();

        try {
            Scanner in = new Scanner(System.in);
            //double no = 0;
            String filename = "";
            System.out.println("Please choice one of them :");
            System.out.println("1-Load dataset");
            System.out.println("2-Create dataset and save the dataset");

            int choice_dataset = in.nextInt();
            if (choice_dataset == 1) {

                ArrayList<PhysicalMachine> arrayList3 = new ArrayList<PhysicalMachine>();
                try {

                    File yourFile = new File("Datasets.txt");
                    if (yourFile.exists()) {
                        ;
                    } else {
                        yourFile.createNewFile();
                    }

                    BufferedReader load = new BufferedReader(new FileReader(yourFile));

                    filename = load.readLine();
                    BufferedReader reader = new BufferedReader(new FileReader(filename));

                    reader.readLine();
                    int NoPort = Integer.parseInt(reader.readLine());//in.nextInt();
                    //System.out.println("portttttttttt " + NoPort);
                    reader.readLine();
                    int noapp = Integer.parseInt(reader.readLine());
                    //System.out.println("Appppppppppp " + noapp);

                    reader.close();

                    int NoPod = NoPort;
                    double MaxNoServer = (Math.pow(NoPort, 3)) / 4;
                    no = no + MaxNoServer;
                    double NoCoreSwitch = Math.pow((NoPort / 2), 2);

                    no = no + NoCoreSwitch;
                    int NoAggSwitch = (NoPort / 2) * NoPod;
                    no = no + NoAggSwitch;
                    int NoEdgeSwitch = (NoPort / 2) * NoPod;
                    no = no + NoEdgeSwitch;
                    //each aggigation connect to Noport/2 core and NoPort/2 edge
                    int NoConnAggNode = NoPort / 2 + NoPort / 2;
                    //each edge connect to Noport/2 server and NoPort/2 aggigation
                    int NoConnEdgeNode = NoPort / 2 + NoPort / 2;

                    System.out.println("no " + no);
                    //int vertices = 6;

                    List<Edge> edges = new ArrayList<Edge>();
                    int counter1 = 0;
                    //int k = (int) (no-NoCoreSwitch);
                    int value1 = (int) (no - NoCoreSwitch);
                    for (int j = (int) no; counter1 < (int) NoCoreSwitch; j--) {
                        int aggcounter = 0;
                        int k = (int) (no - NoCoreSwitch);

                        if (j % (NoPort / 2) == 0 && j != no) //2
                        {
                            k = k - 1;
                            value1 = k;
                        }
                        for (k = value1; aggcounter < NoPod; k = k - NoPort / 2) {

                            int RandomNo1 = randomNumberGenerator();
                            edges.add(new Edge(j, k, RandomNo1));
                            edges.add(new Edge(k, j, RandomNo1));

                            aggcounter++;
                            ////////////////////////////////

                        }
                        counter1++;

                    }

                    //int l = 9999;
                    int i = (int) (no - NoCoreSwitch);
                    int j = 0;

                    j = (int) (no - NoCoreSwitch - NoAggSwitch);
                    int aggrCcounter = 0;
                    int val = (int) (no - NoCoreSwitch - NoAggSwitch);
                    for (i = (int) (no - NoCoreSwitch); aggrCcounter < NoPort * 2; i--) {
                        // j = 10000;

                        int edgecounter = 0;

                        if (i % (NoPort / 2) == 0 && i != (int) (no - NoCoreSwitch)) {
                            //i = i + NoPort / 2;

                            j = val;
                            j = j - NoPort / 2;
                            val = j;
                        }

                        for (j = val; edgecounter < (NoPod / 2); j--) {

                            int RandomNo2 = randomNumberGenerator();
                            edges.add(new Edge(i, j, RandomNo2));
                            edges.add(new Edge(j, i, RandomNo2));

                            edgecounter++;

                        }
                        aggrCcounter++;

                    }
                    int servercounter = 0;

                    int server = (int) (no - NoCoreSwitch - NoAggSwitch - NoEdgeSwitch);
                    for (int m = (int) (no - NoCoreSwitch - NoAggSwitch); servercounter < NoEdgeSwitch; m--) {
                        int linkcounter = 0;
                        for (; linkcounter < NoPort / 2; server--) {

                            int RandomNo3 = randomNumberGenerator();
                            edges.add(new Edge(m, server, RandomNo3));
                            edges.add(new Edge(server, m, RandomNo3));

                            linkcounter++;
                        }
                        servercounter++;

                    }

                    ob.graph = new Graph(edges, (int) no + 1);
                    //ob.graph.printGraph();

                    System.out.println("--------------------------------------------");

                    System.out.println("Number of Core switch = " + NoCoreSwitch);
                    System.out.println("Number of Aggigation switch = " + NoAggSwitch);
                    System.out.println("Number of Edge switch = " + NoEdgeSwitch);
                    System.out.println("Number of Server  = " + MaxNoServer);
                    ArrayList<PhysicalMachine> ArrayListTemp;
                    try {
                        ob.PM = new ArrayList<PhysicalMachine>();
                        ArrayListTemp = new ArrayList<PhysicalMachine>();
                        BufferedReader br = new BufferedReader(new FileReader(filename));

                        String str;
                        int testwhile = 0;
                        while ((str = br.readLine()) != null) {
                            //System.out.println(str);
                            if (str.equalsIgnoreCase("PMList")) {
                                //System.out.println("-2222");
                                //br.readLine();
                                while (!(str = br.readLine()).equalsIgnoreCase("EndPMList")) {
                                    //System.out.println(str);
                                    String array[] = str.split(",");
                                    ArrayListTemp.add(new PhysicalMachine(Double.parseDouble(array[1]), Double.parseDouble(array[2]), Double.parseDouble(array[3]), Double.parseDouble(array[4]), array[0], Integer.parseInt(array[5]), Integer.parseInt(array[6])));
                                }
                                break;
                            }
                            //System.out.println("testwhile " + testwhile++);
                            continue;

                        }
                        ob.printPM(ArrayListTemp);
                        System.out.println("--------------------------------");
                        System.out.println(" List of VMS");
                        br.close();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    ob.app = new ArrayList<Application>(noapp);
                    ob.sortedApp = new ArrayList<Application>(noapp);
                    //int StringCount = -22;

                    for (int d = 0; d < noapp; d++) {
                        //int novm = (int) ((Math.random() * ((8 - 2) + 8)) + 2);;
                        String apptype = "";
                        int noofvm = 0;

                        ArrayList<VirtualMachine> vmlists = new ArrayList<VirtualMachine>();
                        //int WeightArray[][];

                        String str;
                        int g = 0;
                        int WeightArray[][];
                        BufferedReader br3 = new BufferedReader(new FileReader(filename));
                        while ((str = br3.readLine()) != null) {
                            if (str.equalsIgnoreCase("Application " + d)) {
                                apptype = br3.readLine();
                                noofvm = Integer.parseInt(br3.readLine());
                                break;
                            }
                            //break;
                        }

                        WeightArray = new int[noofvm][noofvm];
                        br3.close();

                        BufferedReader br2 = new BufferedReader(new FileReader(filename));
                        while ((str = br2.readLine()) != null) {
                            //System.out.println(str);
                            if (str.equalsIgnoreCase("VMApp " + d)) {
                                //System.out.println("VMApp   ppp");
                                while ((!(str = br2.readLine()).equalsIgnoreCase("EndVMApp " + d))) {

                                    String array[] = str.split(",");
                                    vmlists.add(new VirtualMachine(Double.parseDouble(array[0]), Double.parseDouble(array[1]), array[2]));

                                }
                            }

                        }
                        br2.close();
                        try {
                            BufferedReader br4 = new BufferedReader(new FileReader(filename));
                            while ((str = br4.readLine()) != null) {

                                if (str.equalsIgnoreCase("ConnectionVMApp " + d)) {

                                    int i2 = 0;
                                    while ((!(str = br4.readLine()).equalsIgnoreCase("EndConnectionVMApp " + d))) {

                                        String[] line = str.split(",");

                                        for (int j2 = 0; j2 < line.length; j2++) {
                                            //System.out.println("test  " + Integer.parseInt(line[j2]));
                                            WeightArray[i2][j2] = Integer.parseInt(line[j2]);
                                            //System.out.print(WeightArray[i2][j2] + " t ");
                                        }
                                        i2++;
                                    }

                                    break;
                                }
                                continue;
                            }
                            System.out.println("");
                            br4.close();
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                        //System.out.println("finish");

                        ob.app.add(new Application(apptype, noofvm, vmlists, WeightArray));
                        ob.printListVM(vmlists);

                        ob.printWeightMatrix(WeightArray);

                        System.out.println("---------------------List Virtual machines Connections-------------------");
                        ob.printList(WeightArray);
                        System.out.println("---------------------List of Values of Connection (Weight)-------------------");
                        ob.printValue(WeightArray);
                        System.out.println("******************************************");

                    } //End For Application

                    System.out.println("start sorttttttttttttttttttttttttttttttttttttttttttt application  read all critical app after that read normal app");
                    ////////////  start sortting application //////
                    for (int d = 0; d < noapp; d++) {
                        boolean countercheck = false;
                        //int novm = (int) ((Math.random() * ((8 - 2) + 8)) + 2);;
                        String apptype = "";
                        int noofvm = 0;

                        ArrayList<VirtualMachine> vmlists = new ArrayList<VirtualMachine>();
                        //int WeightArray[][];

                        String str;
                        int g = 0;
                        int WeightArray[][];
                        BufferedReader br3 = new BufferedReader(new FileReader(filename));
                        while ((str = br3.readLine()) != null) {
                            if (str.equalsIgnoreCase("Application " + d)) {
                                apptype = br3.readLine();
                                if (apptype.equalsIgnoreCase("Critical")) {
                                    ob.CriticalCounter++;

                                    noofvm = Integer.parseInt(br3.readLine());
                                    System.out.println("no of vm " + noofvm);
                                    countercheck = true;
                                    break;
                                }

                            }
                            //break;
                        }
                        if (countercheck == true) {

                            WeightArray = new int[noofvm][noofvm];
                            br3.close();

                            BufferedReader br2 = new BufferedReader(new FileReader(filename));
                            while ((str = br2.readLine()) != null) {
                                //System.out.println(str);
                                if (str.equalsIgnoreCase("VMApp " + d)) {
                                    //System.out.println("VMApp   ppp");
                                    while ((!(str = br2.readLine()).equalsIgnoreCase("EndVMApp " + d))) {

                                        String array[] = str.split(",");
                                        vmlists.add(new VirtualMachine(Double.parseDouble(array[0]), Double.parseDouble(array[1]), array[2]));
                                    }
                                }
                            }
                            br2.close();
                            try {
                                BufferedReader br4 = new BufferedReader(new FileReader(filename));
                                while ((str = br4.readLine()) != null) {

                                    if (str.equalsIgnoreCase("ConnectionVMApp " + d)) {

                                        int i2 = 0;
                                        while ((!(str = br4.readLine()).equalsIgnoreCase("EndConnectionVMApp " + d))) {

                                            String[] line = str.split(",");
                                            //System.out.println("line "+line[0]);

                                            for (int j2 = 0; j2 < line.length; j2++) {
                                                //System.out.println("test  " + Integer.parseInt(line[j2]));
                                                WeightArray[i2][j2] = Integer.parseInt(line[j2]);
                                                //System.out.print(WeightArray[i2][j2] + " t ");
                                            }
                                            i2++;
                                        }

                                        break;
                                    }
                                    continue;
                                }
                                System.out.println("");
                                br4.close();
                            } catch (Exception e) {
                                e.printStackTrace();
                            }
                            //System.out.println("finish");

                            ob.sortedApp.add(new Application(apptype, noofvm, vmlists, WeightArray));
                            ob.printListVM(vmlists);

                            ob.printWeightMatrix(WeightArray);

                            System.out.println("---------------------List Virtual machines Connections-------------------");
                            ob.printList(WeightArray);
                            System.out.println("---------------------List of Values of Connection (Weight)-------------------");
                            ob.printValue(WeightArray);
                            System.out.println("******************************************");
                        }

                    }

                    ///////// normal application//////
                    for (int d = 0; d < noapp; d++) {
                        boolean countercheck = false;
                        //int novm = (int) ((Math.random() * ((8 - 2) + 8)) + 2);;
                        String apptype = "";
                        int noofvm = 0;

                        ArrayList<VirtualMachine> vmlists = new ArrayList<VirtualMachine>();
                        //int WeightArray[][];

                        String str;
                        int g = 0;
                        int WeightArray[][];
                        BufferedReader br3 = new BufferedReader(new FileReader(filename));
                        while ((str = br3.readLine()) != null) {
                            if (str.equalsIgnoreCase("Application " + d)) {
                                apptype = br3.readLine();
                                if (apptype.equalsIgnoreCase("Normal")) {
                                    //System.out.println("yessssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss");
                                    System.out.println("application " + d + " Normal");
                                    noofvm = Integer.parseInt(br3.readLine());
                                    System.out.println("no of vm " + noofvm);
                                    countercheck = true;
                                    break;
                                }

                            }
                            //break;
                        }
                        if (countercheck == true) {

                            WeightArray = new int[noofvm][noofvm];
                            br3.close();

                            BufferedReader br2 = new BufferedReader(new FileReader(filename));
                            while ((str = br2.readLine()) != null) {
                                //System.out.println(str);
                                if (str.equalsIgnoreCase("VMApp " + d)) {
                                    //System.out.println("VMApp   ppp");
                                    while ((!(str = br2.readLine()).equalsIgnoreCase("EndVMApp " + d))) {

                                        String array[] = str.split(",");
                                        vmlists.add(new VirtualMachine(Double.parseDouble(array[0]), Double.parseDouble(array[1]), array[2]));

                                    }

                                }
                                //break;
                            }
                            br2.close();
                            try {
                                BufferedReader br4 = new BufferedReader(new FileReader(filename));
                                while ((str = br4.readLine()) != null) {

                                    if (str.equalsIgnoreCase("ConnectionVMApp " + d)) {

                                        int i2 = 0;
                                        while ((!(str = br4.readLine()).equalsIgnoreCase("EndConnectionVMApp " + d))) {

                                            String[] line = str.split(",");
                                            //System.out.println("line "+line[0]);

                                            for (int j2 = 0; j2 < line.length; j2++) {
                                                //System.out.println("test  " + Integer.parseInt(line[j2]));
                                                WeightArray[i2][j2] = Integer.parseInt(line[j2]);
                                                //System.out.print(WeightArray[i2][j2] + " t ");
                                            }
                                            i2++;
                                        }

                                        break;
                                    }
                                    continue;
                                }
                                System.out.println("");
                                br4.close();
                            } catch (Exception e) {
                                e.printStackTrace();
                            }
                            //System.out.println("finish");

                            ob.sortedApp.add(new Application(apptype, noofvm, vmlists, WeightArray));
                            ob.printListVM(vmlists);

                            ob.printWeightMatrix(WeightArray);

                            System.out.println("---------------------List Virtual machines Connections-------------------");
                            ob.printList(WeightArray);
                            System.out.println("---------------------List of Values of Connection (Weight)-------------------");
                            ob.printValue(WeightArray);
                            System.out.println("******************************************");
                        }

                    } //End For sortting Application

                    System.out.println("End sorttttttttttttttttttttttttttttttttttttttttttt application");

                    Scanner inp = new Scanner(System.in);
                    int choice;
                    boolean bb = false;

                    //while (true) {
                    System.out.println("Pleace chooce one of those algorithm :");
                    System.out.println("1- FFD");
                    System.out.println("2- EQVMP");
                    System.out.println("3- PAVA");
                    System.out.println("4- PROPOSED");
                    System.out.println("5- Exit");

                    choice = inp.nextInt();
                    if (choice == 1) {
                        ArrayList<PhysicalMachine> arrayList5 = new ArrayList<PhysicalMachine>();
                        try {
                            BufferedReader in1 = new BufferedReader(new FileReader(filename));

                            String str = "";

                            while ((str = in1.readLine()) != null) {
                                //System.out.println(str);
                                if (str.equalsIgnoreCase("PMList")) {

                                    while (!(str = in1.readLine()).equalsIgnoreCase("EndPMList")) {
                                        String array[] = str.split(",");
                                        arrayList5.add(new PhysicalMachine(Double.parseDouble(array[1]), Double.parseDouble(array[2]), Double.parseDouble(array[3]), Double.parseDouble(array[4]), array[0], Integer.parseInt(array[5]), Integer.parseInt(array[6])));
                                    }
                                    break;
                                }

                                continue;

                            }
                            in1.close();

                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                        ob.placementComputationFirstFitDecreasing(arrayList5, ob.app);
                    } else if (choice == 2) {
                        ArrayList<PhysicalMachine> arrayList5 = new ArrayList<PhysicalMachine>();
                        try {
                            BufferedReader in1 = new BufferedReader(new FileReader(filename));

                            String str = "";

                            while ((str = in1.readLine()) != null) {

                                if (str.equalsIgnoreCase("PMList")) {

                                    while (!((str = in1.readLine()).equalsIgnoreCase("EndPMList"))) {
                                        String array[] = str.split(",");
                                        //System.out.println("wihle  "+str);
                                        arrayList5.add(new PhysicalMachine(Double.parseDouble(array[1]), Double.parseDouble(array[2]), Double.parseDouble(array[3]), Double.parseDouble(array[4]), array[0], Integer.parseInt(array[5]), Integer.parseInt(array[6])));
                                    }
                                    //break;
                                    continue;
                                }

                                continue;

                            }
                            in1.close();

                        } catch (Exception e) {
                            e.printStackTrace();
                        }

                        //System.out.println("critical counter " + ob.CriticalCounter);
                        ob.EQVMP(arrayList5, ob.app);

                    } else if (choice == 3) {
                        ArrayList<PhysicalMachine> arrayList5 = new ArrayList<PhysicalMachine>();
                        try {
                            BufferedReader in1 = new BufferedReader(new FileReader(filename));

                            String str = "";

                            while ((str = in1.readLine()) != null) {

                                if (str.equalsIgnoreCase("PMList")) {
                                    //System.out.println("-2222");
                                    //br.readLine();
                                    while (!((str = in1.readLine()).equalsIgnoreCase("EndPMList"))) {
                                        String array[] = str.split(",");
                                        //System.out.println("wihle  "+str);
                                        arrayList5.add(new PhysicalMachine(Double.parseDouble(array[1]), Double.parseDouble(array[2]), Double.parseDouble(array[3]), Double.parseDouble(array[4]), array[0], Integer.parseInt(array[5]), Integer.parseInt(array[6])));
                                    }
                                    //break;
                                    continue;
                                }

                                continue;

                            }
                            in1.close();

                        } catch (Exception e) {
                            e.printStackTrace();
                        }

                        System.out.println("critical counter " + ob.CriticalCounter);
                        ob.placementComputationPAVA(arrayList5, ob.app, ob.CriticalCounter);

                    } else if (choice == 4) {
                        ArrayList<PhysicalMachine> arrayList5 = new ArrayList<PhysicalMachine>();
                        try {
                            BufferedReader in1 = new BufferedReader(new FileReader(filename));

                            String str = "";

                            while ((str = in1.readLine()) != null) {

                                if (str.equalsIgnoreCase("PMList")) {

                                    while (!((str = in1.readLine()).equalsIgnoreCase("EndPMList"))) {
                                        String array[] = str.split(",");
                                        //System.out.println("wihle  "+str);
                                        arrayList5.add(new PhysicalMachine(Double.parseDouble(array[1]), Double.parseDouble(array[2]), Double.parseDouble(array[3]), Double.parseDouble(array[4]), array[0], Integer.parseInt(array[5]), Integer.parseInt(array[6])));
                                    }
                                    //break;
                                    continue;
                                }

                                continue;

                            }
                            in1.close();

                        } catch (Exception e) {
                            e.printStackTrace();
                        }

                        ob.proposed(arrayList5, ob.app);

                    } else if (choice == 5) {
                        bb = true;
                        System.exit(0);
                        // break;
                    }

                    // }
                    System.out.println("--------------------------");
                } catch (Exception e) {
                    System.out.println("Please Create a Dataset First Time");
                    //e.printStackTrace();
                }

            } else if (choice_dataset == 2) {

                System.out.println("Enter No. of ports of switchs");
                int NoPort = in.nextInt();

                int NoPod = NoPort;
                double MaxNoServer = (Math.pow(NoPort, 3)) / 4;
                no = no + MaxNoServer;
                double NoCoreSwitch = Math.pow((NoPort / 2), 2);
                no = no + NoCoreSwitch;
                int NoAggSwitch = (NoPort / 2) * NoPod;
                no = no + NoAggSwitch;
                int NoEdgeSwitch = (NoPort / 2) * NoPod;
                no = no + NoEdgeSwitch;
                //each aggigation connect to Noport/2 core and NoPort/2 edge
                int NoConnAggNode = NoPort / 2 + NoPort / 2;
                //each edge connect to Noport/2 server and NoPort/2 aggigation
                int NoConnEdgeNode = NoPort / 2 + NoPort / 2;

                System.out.println("no " + no);
                //int vertices = 6;

                List<Edge> edges = new ArrayList<Edge>();
                int counter1 = 0;
                //int k = (int) (no-NoCoreSwitch);
                int value1 = (int) (no - NoCoreSwitch);
                for (int j = (int) no; counter1 < (int) NoCoreSwitch; j--) {
                    int aggcounter = 0;
                    int k = (int) (no - NoCoreSwitch);

                    if (j % (NoPort / 2) == 0 && j != no) //2
                    {
                        k = k - 1;
                        value1 = k;
                    }
                    for (k = value1; aggcounter < NoPod; k = k - NoPort / 2) {

                        int RandomNo1 = randomNumberGenerator();
                        edges.add(new Edge(j, k, RandomNo1));
                        edges.add(new Edge(k, j, RandomNo1));

                        aggcounter++;
                        ////////////////////////////////

                    }
                    counter1++;

                }

                //int l = 9999;
                int i = (int) (no - NoCoreSwitch);
                int j = 0;

                j = (int) (no - NoCoreSwitch - NoAggSwitch);
                int aggrCcounter = 0;
                int val = (int) (no - NoCoreSwitch - NoAggSwitch);
                for (i = (int) (no - NoCoreSwitch); aggrCcounter < NoPort * 2; i--) {
                    // j = 10000;

                    int edgecounter = 0;

                    if (i % (NoPort / 2) == 0 && i != (int) (no - NoCoreSwitch)) {

                        j = val;
                        j = j - NoPort / 2;
                        val = j;
                    }

                    for (j = val; edgecounter < (NoPod / 2); j--) {

                        int RandomNo2 = randomNumberGenerator();
                        edges.add(new Edge(i, j, RandomNo2));
                        edges.add(new Edge(j, i, RandomNo2));

                        edgecounter++;

                    }
                    aggrCcounter++;

                }
                int servercounter = 0;

                int server = (int) (no - NoCoreSwitch - NoAggSwitch - NoEdgeSwitch);
                for (int m = (int) (no - NoCoreSwitch - NoAggSwitch); servercounter < NoEdgeSwitch; m--) {
                    int linkcounter = 0;
                    for (; linkcounter < NoPort / 2; server--) {

                        int RandomNo3 = randomNumberGenerator();
                        edges.add(new Edge(m, server, RandomNo3));
                        edges.add(new Edge(server, m, RandomNo3));

                        linkcounter++;
                    }
                    servercounter++;

                }

                ob.graph = new Graph(edges, (int) no + 1);

                System.out.println("--------------------------------------------");

                System.out.println("Number of Core switch = " + NoCoreSwitch);
                System.out.println("Number of Aggigation switch = " + NoAggSwitch);
                System.out.println("Number of Edge switch = " + NoEdgeSwitch);
                System.out.println("Number of Server  = " + MaxNoServer);

                System.out.println("Enter number of Applications");
                int noapp = in.nextInt();
                //ob.app = ob.createApp(noapp);
                filename = "Dataset_" + NoPort + "_ports_" + noapp + "_Applications.txt";
                File file = new File(filename);

                if (!file.exists()) {
                    file.createNewFile();
                }

                File datasets = new File("Datasets.txt");

                if (!datasets.exists()) {
                    datasets.createNewFile();
                }
                PrintWriter printdatasets = new PrintWriter(datasets);
                printdatasets.write(filename);
                printdatasets.println();
                printdatasets.close();

                //ArrayList<PhysicalMachine> arrayListPM = ob.createPM((int) MaxNoServer);
                PrintWriter pw = new PrintWriter(file);
                pw.write("PortNO");
                pw.println();
                pw.write(NoPort + "");
                pw.println();
                pw.write("NumberOfApplication");
                pw.println();
                pw.write((int) noapp + "");
                pw.println();
                pw.write("PMList");
                pw.println();
                try {
                    ob.PM = new ArrayList<PhysicalMachine>();

                    System.out.println("no " + no);

                    int counter = (int) (MaxNoServer + 1);

                    int servercounter2 = 0;
                    int arraytempo[] = new int[(int) MaxNoServer];
                    int index = (int) (MaxNoServer - 1);
                    int server2 = (int) (no - NoCoreSwitch - NoAggSwitch - NoEdgeSwitch);
                    for (int m = (int) (no - NoCoreSwitch - NoAggSwitch); servercounter2 < NoEdgeSwitch; m--) {
                        int linkcounter = 0;
                        for (; linkcounter < NoPort / 2; server2--) {
                            arraytempo[index] = m;
                            index--;
                            linkcounter++;
                        }
                        servercounter2++;
                    }

                    int podcounter = 0;
                    int indexpod = 0;

                    //comment for experiment 4
                    int CPU[] = {56, 64, 128, 112, 64, 56, 128};
                    int RAM[] = {192, 128, 256, 384, 256, 192, 512};
                    double pmax[] = {347, 246, 422, 672, 250, 385, 412};
                    double pmin[] = {48.3, 99.5, 111, 121, 60.8, 48.6, 99.2};
                    //end comment experiment 4

                    Random r1 = new Random();
                    Random r2 = new Random();
                    Random rPmax = new Random();
                    Random rPmin = new Random();

                    double noofserverperpod = Math.pow((NoPort / 2), 2);

                    for (int i2 = 1; i2 <= MaxNoServer; i2++) {

                        if (indexpod % noofserverperpod == 0) {
                            podcounter++;
                        }

                        int low1 = 0;
                        int high1 = 7;
                        int cpu = r1.nextInt(high1 - low1) + low1;

                        int low2 = 0;
                        int high2 = 4;
                        int ram = r2.nextInt(high2 - low2) + low2;

                        //commnet for EXperiment 4
                        ob.PM.add(new PhysicalMachine(CPU[cpu], RAM[cpu], pmax[cpu], pmin[cpu], "PM_" + i2, arraytempo[i2 - 1], podcounter));
                        pw.write("PM_" + i2 + "," + CPU[cpu] + "," + RAM[cpu] + "," + pmax[cpu] + "," + pmin[cpu] + "," + arraytempo[i2 - 1] + "," + podcounter);
                        pw.println();

                        //end comment experiment 4
                        indexpod++;
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                }
                ob.printPM(ob.PM);

                pw.write("EndPMList");
                pw.println();

                double percentage = 0.5;
                int chk[] = new int[noapp];
                int appcounter = 0;
                int Cno = (int) Math.round(percentage * noapp);

                for (int ii = 0; ii < Cno; ii++) {
                    chk[ii] = 1;
                    appcounter++;
                }
                double Npercentage = 1 - percentage;
                int Nno = (int) Math.floor(Npercentage * noapp);

                for (int iii = appcounter; iii < Nno; iii++) {
                    chk[iii] = 0;
                    appcounter++;
                }

                chk = ob.RandomizeArray(chk);

                ob.app = new ArrayList<Application>(noapp);
                for (int d = 0; d < noapp; d++) {
                    //int novm = (int) ((Math.random() * ((8 - 2) + 8)) + 2);;
                    pw.write("Application " + d);
                    pw.println();
                    Random r = new Random();
                    int low = 1;
                    int high = 17;
                    //int novm = 16;
                    int novm = r.nextInt(high - low) + low;
                    Random r2 = new Random();
                    if (chk[d] == 1) {
                        //if (check >= 0) {
                        //novm = r.nextInt(8 - 2) + 2;
                        ob.CriticalCounter++;
                        int WeightArray[][] = ob.generateWeight(novm);
                        WeightArray = ob.checkMatrix(WeightArray);

                        ArrayList<VirtualMachine> vmlist = ob.createVMApp(novm);
                        System.out.println(d + " Critical");
                        pw.write("Critical");
                        pw.println();
                        pw.write(novm + "");
                        pw.println();
                        pw.write("VMApp " + d);
                        pw.println();
                        //for (int jjj = 0; jjj < novm; jjj++) {
                        for (int jjj = 0; jjj < vmlist.size(); jjj++) {
                            pw.write(vmlist.get(jjj).CPU + "," + vmlist.get(jjj).RAM + "," + vmlist.get(jjj).Name + "_" + jjj);
                            //System.out.println("  tt " + vmlist.get(jjj).CPU + "," + vmlist.get(jjj).RAM + "," + vmlist.get(jjj).Name + "_" + jjj);
                            pw.println();
                        }
                        pw.write("EndVMApp " + d);
                        pw.println();

                        pw.write("ConnectionVMApp " + d);
                        pw.println();
                        for (int ii = 0; ii < WeightArray.length; ii++) {
                            for (int jj = 0; jj < WeightArray[ii].length; jj++) {
                                if (jj == WeightArray[ii].length - 1) {
                                    pw.write(WeightArray[ii][jj] + "");
                                } else {
                                    pw.write(WeightArray[ii][jj] + ",");
                                }
                            }
                            pw.println();
                        }
                        pw.print("EndConnectionVMApp " + d);
                        pw.println();
                        //ob.printPM(ob.PM);

                        ob.app.add(new Application("Critical Application", novm, vmlist, WeightArray));
                        System.out.println("The " + d + "  Critical application Virtual machine's");
                        ob.printListVM(vmlist);
                        ob.printWeightMatrix(WeightArray);

                        System.out.println("---------------------List Virtual machines Connections-------------------");
                        ob.printList(WeightArray);
                        System.out.println("---------------------List of Values of Connection (Weight)-------------------");
                        ob.printValue(WeightArray);
                        System.out.println("******************************************");

                    } else {
                        ob.NormalCounter++;
                        int WeightArray[][] = ob.generateWeight(novm);
                        WeightArray = ob.checkMatrix(WeightArray);
                        ArrayList<VirtualMachine> vmlist = ob.createVMApp(novm);
                        System.out.println(d + " Normal");
                        pw.write("Normal");
                        pw.println();
                        pw.write(novm + "");
                        pw.println();
                        pw.write("VMApp " + d);
                        pw.println();

                        for (int jjj = 0; jjj < vmlist.size(); jjj++) {
                            pw.write(vmlist.get(jjj).CPU + "," + vmlist.get(jjj).RAM + "," + vmlist.get(jjj).Name + "_" + jjj);
                            pw.println();
                        }
                        pw.write("EndVMApp " + d);
                        pw.println();
                        pw.write("ConnectionVMApp " + d);
                        pw.println();

                        for (int ii = 0; ii < WeightArray.length; ii++) {
                            for (int jj = 0; jj < WeightArray[ii].length; jj++) {
                                if (jj == WeightArray[ii].length - 1) {
                                    pw.write(WeightArray[ii][jj] + "");
                                } else {
                                    pw.write(WeightArray[ii][jj] + ",");
                                }
                            }
                            pw.println();
                        }
                        pw.print("EndConnectionVMApp " + d);
                        pw.println();

                        ob.app.add(new Application("Normal Application", novm, vmlist, WeightArray));
                        System.out.println("The " + d + " Normal application Virtual machine's");
                        ob.printListVM(vmlist);
                        ob.printWeightMatrix(WeightArray);

                        System.out.println("---------------------List Virtual machines Connections-------------------");
                        ob.printList(WeightArray);
                        System.out.println("---------------------List of Values of Connection (Weight)-------------------");
                        ob.printValue(WeightArray);
                        System.out.println("******************************************");

                    }

                } //End For Application
                pw.flush();
                pw.close();
                System.out.println("Normal Application  " + ob.NormalCounter);
                System.out.println("Critical Application  " + ob.CriticalCounter);

            }

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}

class Node {

    int vertex, weight;

    public Node(int vertex, int weight) {
        this.vertex = vertex;
        this.weight = weight;
    }
};
