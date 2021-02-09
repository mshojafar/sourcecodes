Help file to run the project written in java of CloudSim toolkit, developed using Eclipse IDE.

Dr. Zhou Zhou did this implementation. 
Dr. Mohammad Shojafar helped on idea brainstorming and documentation. 
Profesors Abawajy, Alazab and Li helped in English correction and leading the team.

Contact:

Dr. Zhou Zhou
Department of mathematics and computer science, Changsha University, Changsha, China
E-mail: zhouzhou03201@126.com
or

Dr. Mohammad Shojafar
5GIC & 6GIC, Institute for Communication Systems (ICS), University of Surrey, Guildford,GU27XH, United Kingdom
E-mail: m.shojafar@surrey.ac.uk
 

Step of the running project:

* Please watch the intrusction and introduction video recorded by Dr. Zhou Zhou.

Follow these steps: 

1.	Inatall CloudSim toolkit and Eclipse IDE and set the required Machine environemnt (please see cloudsim help).
2.	Open this project in your Eclipse IDE.
3. modify some classes and  functions ( for example:  class RunnerAbstract,  PowerVmAllocationPolicyMigrationAbstract, class PowerVmAllocationPolicyMigrationAbstract,  optimizeAllocation, public PowerHost findHostForVm );
4.  Run and obtain the experimental results:  

Experiment name: 20110303_KMeansMadIQR_mmt_1.0
Number of hosts: 800
Number of VMs: 1052
Total simulation time: 86400.00 sec
Energy consumption: 73.57 kWh
Number of VM migrations: 5710
SLA: 0.00013%
SLA perf degradation due to migration: 0.00772%
SLA time per active host: 1.66%
Overall SLA violation: 0.01%
Average SLA violation: 9.50%
Number of host shutdowns: 792
Mean time before a host shutdown: 2200.33 sec
StDev time before a host shutdown: 8905.72 sec
Mean time before a VM migration: 21.65 sec
StDev time before a VM migration: 7.85 sec
Execution time - VM selection mean: 0.00057 sec
Execution time - VM selection stDev: 0.00168 sec
Execution time - host selection mean: 0.00694 sec
Execution time - host selection stDev: 0.00257 sec
Execution time - VM reallocation mean: 0.12738 sec
Execution time - VM reallocation stDev: 0.10266 sec
Execution time - total mean: 0.17585 sec
Execution time - total stDev: 0.31658 sec
