#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<conio.h>
#include<math.h>
#define max 1000
struct machine {
	int id,c,r,b,rank,loc;
	float l,h,cu,ru,bu,a,v,p_min,p_max; // pc for power consumption
};
void FFD(struct machine [],struct machine [],int,int);
void MBFD_Beloglazov(struct machine [],struct machine [],int,int);
void RAVMP(struct machine [],struct machine [],int,int);
void MinPR(struct machine [],struct machine [],int,int);
void pm_increasing_order(struct machine [],struct machine [], int);
void pm_decreasing_order(struct machine [], int);
void vm_decreasing_order(struct machine [], int);
main()
{
	srand(time(NULL));
	int i, j, code, nvm, npm, x;
	struct machine vm[max],pm[max];
	for (j=0; j<10; j++) {
	printf("\nEnter # of PMs and # of VMs (npm=nvm): \n");
	scanf("%d %d", &npm, &nvm);
	printf("\n PMs:\n");
	for (i=0; i<npm; i++) {
		x = rand()%2;  // RAVMP dataset both synthetic and real world
		if (x==0) {
			pm[i].c = 5320; 
			pm[i].r = 4096;
			pm[i].p_min = 93.7;
			pm[i].p_max = 135;
		}	
		else if (x==1) {
			pm[i].c = 3720;
			pm[i].r = 4096;
			pm[i].p_min = 86;
			pm[i].p_max = 117;
		}
	    pm[i].b = rand()%6001 + 2000;
		pm[i].id = i;
		printf("pm[%d] = %d\t%d\n", pm[i].id, pm[i].c, pm[i].r);
	}
	printf("\n VMs:\n");
	for (i=0; i<nvm; i++) {	
		x = rand()%4;  // RAVMP real world dataset
		if (x==0) {
			vm[i].c = 500;
			vm[i].r = 613;
		}	
		else if (x==1) {
			vm[i].c = 1000;
			vm[i].r = 1700;
		}
		else if (x==2) {
			vm[i].c = 2000;
			vm[i].r = 3750;
		}
		else  {
			vm[i].c = 2500;
			vm[i].r = 850;
		}
//		vm[i].c = rand()% 2001 + 500; // RAVMP dataset with a little change
//		vm[i].r = rand()% 1501 + 500;	
		vm[i].b = rand()%901 + 100;
		vm[i].id = i;
		printf("vm[%d] = %d\t%d\n", vm[i].id, vm[i].c, vm[i].r);
	}
	while (1) {	
	printf("\n1-FFD");
	printf("\n2-MBFD_Beloglazov");
	printf("\n3-RAVMP");
	printf("\n4-MinPR");
	printf("\n5-Start New Iteration\n"); 
	printf("\nEnter Your Choice: ");
	scanf("%d", &code);
	switch(code)
	{

		case 1: FFD(pm,vm,npm,nvm);
				break;
		case 2: MBFD_Beloglazov(pm,vm,npm,nvm);
				break;
		case 3: RAVMP(pm,vm,npm,nvm);
				break;
		case 4: MinPR(pm,vm,npm,nvm);
				break;
		case 5: break; 	
	}
	if (code == 5)
		break;
	}
	}
}

/////////////////////////  FFD  ////////////////////////////////
void FFD(struct machine pm[],struct machine vm[],int npm,int nvm)
{
	int i, j, idle=0, unaloc_vm=0;
	float sum_cu=0, sum_ru=0, rw=0, pc=0;
	struct machine pm_ffd[npm], pm_ffd_copy[npm], vm_ffd[nvm];
	for (i=0; i<npm; i++) 
		pm_ffd[i] = pm[i];
	for (i=0; i<nvm; i++) 
		vm_ffd[i] = vm[i];
	vm_decreasing_order(vm_ffd, nvm);
	for (i=0; i<npm; i++) 
		pm_ffd_copy[i] = pm_ffd[i];
	printf("\npm_ffd --> pm_increasing_order:\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d,%d]\n", pm_ffd[i].id, pm_ffd[i].c, pm_ffd[i].r, pm_ffd[i].b);
	printf("\nvm_ffd --> vm_decreasing_order:\n");
	for (i=0; i<nvm; i++)
		printf("[%d,%d,%d,%d]\n", vm_ffd[i].id, vm_ffd[i].c, vm_ffd[i].r, vm_ffd[i].b);
	printf("\n");	
	for (i=0; i<nvm; i++)
		for (j=0; j<npm; j++) {
			if (vm_ffd[i].c <= pm_ffd[j].c && vm_ffd[i].r <= pm_ffd[j].r && vm_ffd[i].b <= pm_ffd[j].b)
			{
				pm_ffd[j].c -= vm_ffd[i].c;
				pm_ffd[j].r -= vm_ffd[i].r;
				pm_ffd[j].b -= vm_ffd[i].b;
				vm_ffd[i].c=0;
//				printf("vm[%d] --> pm[%d]\n", vm_ffd[i].id, pm_ffd[j].id);
//				getch();
				break;
			}
		}
	printf("\npm_ffd --> final pm list :\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d,%d]\n", pm_ffd[i].id, pm_ffd[i].c, pm_ffd[i].r, pm_ffd[i].b);
	for (i=0; i<npm; i++) {
		if (pm_ffd[i].c == pm_ffd_copy[i].c)
			idle++;
		else
			rw += abs( pm_ffd[i].c-pm_ffd[i].r) / (float)((pm_ffd_copy[i].c - pm_ffd[i].c)+(pm_ffd_copy[i].r - pm_ffd[i].r)) + .0001;
		pm_ffd[i].cu = (float)(pm_ffd_copy[i].c-pm_ffd[i].c)/pm_ffd_copy[i].c;
		pm_ffd[i].ru = (float)(pm_ffd_copy[i].r-pm_ffd[i].r)/pm_ffd_copy[i].r;
		printf("pm_ffd[%d].cu=%.2f\n",pm_ffd[i].id, pm_ffd[i].cu);
		printf("pm_ffd[%d].ru=%.2f\n",pm_ffd[i].id, pm_ffd[i].ru);
	}
	for (i=0; i<npm; i++)
		if (pm_ffd[i].cu>0) {
			sum_cu += pm_ffd[i].cu;
			sum_ru += pm_ffd[i].ru;
			pc += (pm_ffd[i].p_min+(pm_ffd[i].p_max-pm_ffd[i].p_min)*pm_ffd[i].cu);
		}
	for (i=0; i<nvm; i++) 
		if (vm_ffd[i].c != 0)
			unaloc_vm++;
	printf("\n\nThe number of unallocated VMs: %d\n",unaloc_vm);
	printf("The number of active PMs: %d\n",npm-idle);
	printf("The mean CPU utilization ratio: %.2f\n",sum_cu/(npm-idle));		
	printf("The mean RAM utilization ratio: %.2f\n",sum_ru/(npm-idle));	
	printf("The resource wastage: %.2f\n",rw);
	printf("The power consumption: %.2f W\n",pc);
}	

//////////////////////////// MBFD_Beloglazov ////////////////////////////
void MBFD_Beloglazov(struct machine pm[],struct machine vm[],int npm,int nvm)
{
	int i, j, idle=0; 
	float sum_cu=0, sum_ru=0, rw=0, pc=0;
	struct machine pm_bfd[npm], pm_copy[npm], vm_bfd[nvm];
	for (i=0; i<npm; i++) 
		pm_bfd[i] = pm[i];
	for (i=0; i<npm; i++) 
		pm_copy[i] = pm_bfd[i];
	for (i=0; i<nvm; i++) 
		vm_bfd[i] = vm[i];
	vm_decreasing_order(vm_bfd, nvm);
//	for (i=0; i<npm; i++) 
//		pm_bfd_copy[i] = pm_bfd[i];
	printf("\npm_bfd:\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d]\n", pm_bfd[i].id, pm_bfd[i].c, pm_bfd[i].r);
	printf("\nvm_bfd --> vm_decreasing_order:\n");
	for (i=0; i<nvm; i++)
		printf("[%d,%d,%d]\n", vm_bfd[i].id, vm_bfd[i].c, vm_bfd[i].r);
	printf("\n");	
//////////////////////////// placement /////////////////////////////////////////////////
	int index, unaloc_vm=0;
	float minp_inc,p_inc, old_cu, new_cu; // min power increase 	
	for (i=0; i<nvm; i++) {
		minp_inc=200;
		for (j=0; j<npm; j++) {
			if (vm_bfd[i].c <= pm_bfd[j].c && vm_bfd[i].r <= pm_bfd[j].r) {
				old_cu = (pm_copy[j].c - pm_bfd[j].c)/(float)pm_copy[j].c;
				new_cu = (pm_copy[j].c - pm_bfd[j].c + vm_bfd[i].c)/(float)pm_copy[j].c;
				if (old_cu==0)
					p_inc = (new_cu - old_cu)*(pm_bfd[j].p_max - pm_bfd[j].p_min);
				else
					p_inc = (new_cu - old_cu)*(pm_bfd[j].p_max - pm_bfd[j].p_min);  
			}
			if (vm_bfd[i].c <= pm_bfd[j].c && vm_bfd[i].r <= pm_bfd[j].r && p_inc < minp_inc)	{
				index = j;
				minp_inc = p_inc;
			}
		}
		pm_bfd[index].c -= vm_bfd[i].c;
		pm_bfd[index].r -= vm_bfd[i].r;				
		vm_bfd[i].c=0;
		vm_bfd[i].loc = index;
//		printf("vm[%d] --> pm[%d]\n",vm_bfd[i].id,pm_bfd[index].id);
//		getch();
	}
///////////////////////////// end of placement ///////////////////////////////////
///////////////////////////// start of migration ////////////////////////////////
//	int index, unaloc_vm=0;
//	float minp_inc,p_inc, old_cu, new_cu; // min power increase 
/*
	int k, ans;
	printf("Do you want to have migration phase? (0 or 1):");
	scanf("%d",&ans);
	if (ans==1){
	for (i=0; i<npm; i++) {
		pm_bfd[i].cu = (float)(pm_copy[i].c-pm_bfd[i].c)/pm_copy[i].c;
		pm_bfd[i].ru = (float)(pm_copy[i].r-pm_bfd[i].r)/pm_copy[i].r;
		if (pm_bfd[i].cu>0 && pm_bfd[i].cu<=.4) 
//			if (pm_bfd[i].cu<=.4)
			for (j=0; j>nvm; j++)
				if (vm_bfd[j].loc==i) {
					minp_inc=100;
					for (k=0; k<npm; k++) {
						if (pm_bfd[k].cu>0 && vm_bfd[j].c <= pm_bfd[k].c && vm_bfd[j].r <= pm_bfd[k].r) {
							old_cu = (pm_copy[k].c - pm_bfd[k].c)/(float)pm_copy[k].c;
							new_cu = (pm_copy[k].c - pm_bfd[k].c + vm_bfd[j].c)/(float)pm_copy[k].c;	
							p_inc = (new_cu - old_cu)*(pm_bfd[k].p_max - pm_bfd[k].p_min); 	
						}
						if (vm_bfd[j].c <= pm_bfd[k].c && vm_bfd[j].r <= pm_bfd[k].r && p_inc < minp_inc) {
							index = k;
							minp_inc = p_inc;
						}
					}
					pm_bfd[index].c -= vm_bfd[j].c;
					pm_bfd[index].r -= vm_bfd[j].r;				
					vm_bfd[i].loc=index;
					pm_bfd[i].c += vm_bfd[j].c;
					pm_bfd[i].r += vm_bfd[j].r;			
				}
		}
	}
*/
////////////////////////////////// end of migration /////////////////////////
	printf("\npm_bfd --> final pm list :\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d]\n", pm_bfd[i].id, pm_bfd[i].c, pm_bfd[i].r);
	for (i=0; i<npm; i++) {
		if (pm_bfd[i].c == pm_copy[i].c)
			idle++;
		else
			rw += abs( pm_bfd[i].c-pm_bfd[i].r) / (float)((pm_copy[i].c - pm_bfd[i].c)+(pm_copy[i].r - pm_bfd[i].r)) + .0001;
		pm_bfd[i].cu = (float)(pm_copy[i].c-pm_bfd[i].c)/pm_copy[i].c;
		pm_bfd[i].ru = (float)(pm_copy[i].r-pm_bfd[i].r)/pm_copy[i].r;
		printf("pm_bfd[%d].cu=%.2f\n",pm_bfd[i].id, pm_bfd[i].cu);
		printf("pm_bfd[%d].ru=%.2f\n",pm_bfd[i].id, pm_bfd[i].ru);
	}
	for (i=0; i<npm; i++)
		if (pm_bfd[i].cu>0) {
			sum_cu += pm_bfd[i].cu;
			sum_ru += pm_bfd[i].ru;
			pc += (pm_bfd[i].p_min+(pm_bfd[i].p_max-pm_bfd[i].p_min)*pm_bfd[i].cu);	
		}
	for (i=0; i<nvm; i++) 
		if (vm_bfd[i].c != 0)
			unaloc_vm++;
	printf("\n\nThe number of unallocated VMs: %d\n",unaloc_vm);
	printf("The number of active PMs: %d\n",npm-idle);
	printf("The mean CPU utilization ratio: %.2f\n",sum_cu/(npm-idle));		
	printf("The mean RAM utilization ratio: %.2f\n",sum_ru/(npm-idle));	
	printf("The resource wastage: %.2f\n",rw);
	printf("The power consumption: %.2f W\n",pc);	
}



////////////////////// RAVMP /////////////////////////
void RAVMP(struct machine pm[],struct machine vm[],int npm,int nvm)
{
int i, j,k,r,c, idle=0, unaloc_vm=0;
	float sum_cu=0, sum_ru=0, rw=0, pc=0;
	struct machine pm_mbfd[npm], pm_copy[npm], vm_mbfd[nvm];
	for (i=0; i<npm; i++) 
		pm_copy[i] = pm[i];
	for (i=0; i<nvm; i++) 
		vm_mbfd[i] = vm[i];
	pm_decreasing_order(pm_copy, npm);
	for (i=0; i<npm; i++) 
		pm_mbfd[i] = pm_copy[i];
	printf("\nList of All PMs:\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d]\n", pm_mbfd[i].id, pm_mbfd[i].c, pm_mbfd[i].r);
	printf("\nList of All VMs:\n");
	for (i=0; i<nvm; i++)
		printf("[%d,%d,%d]\n", vm_mbfd[i].id, vm_mbfd[i].c, vm_mbfd[i].r);
	printf("\n");
	for (i=0; i<npm; i++) {
		pm_mbfd[i].cu = 0;
		pm_mbfd[i].ru = 0;
		pm_mbfd[i].l = 0;
		pm_mbfd[i].h = 0;
		pm_mbfd[i].v = 0;  // avrage utilization per vm located on a pm
		pm_mbfd[i].rank = 0; // # of VMs on each PM
	}
	
///////////////// Placement  ////////////////////////////////////////////////	
	int index;
	float max_ruf, new_cu, new_ru, new_acu, new_aru;
	for (i=0; i<npm; i++) {
	//	printf("pm[%d] is ready to allocation\n",i);
		for (k=0; k<nvm; k++) {
			max_ruf = 0;
			index = 0;
			for (j=0; j<nvm; j++)
				if (pm_mbfd[i].c >= vm_mbfd[j].c && pm_mbfd[i].r >= vm_mbfd[j].r && pm_mbfd[i].b >= vm_mbfd[j].b &&
				vm_mbfd[j].c != 0)	{
				new_cu = (pm_copy[i].c - pm_mbfd[i].c + vm_mbfd[j].c)/(float)pm_copy[i].c;
				new_ru = (pm_copy[i].r - pm_mbfd[i].r + vm_mbfd[j].r)/(float)pm_copy[i].r;
				new_acu = new_cu/(pm_mbfd[i].rank+1);
				new_aru = new_ru/(pm_mbfd[i].rank+1);
				vm_mbfd[j].v = new_cu*new_acu*(1-new_cu) + new_ru*new_aru*(1-new_ru);
			//	printf("new_cu, new_ru, new_acu, new_aru, new_v of vm[%d]: %.4f  %.4f  %.4f  %.4f  %.4f\n", vm_mbfd[j].id,new_cu,new_ru,new_acu,new_aru,pm_mbfd[j].v+1);
			//	printf("%.4f\n",new_cu*new_acu*(1-new_cu));
			//	printf("%.4f\n",new_ru*new_aru*(1-new_ru));
			//	printf("ruf of vm[%d]: %f \n", vm_mbfd[j].id, vm_mbfd[j].v);
				if (max_ruf < vm_mbfd[j].v) {
					max_ruf = vm_mbfd[j].v;
					index = j;
				}
			}
	//	printf("max_ruf=%f  , index=%d\n",max_ruf, index);
		if (max_ruf>0) {
		pm_mbfd[i].c -= vm_mbfd[index].c;
		pm_mbfd[i].r -= vm_mbfd[index].r;
		pm_mbfd[i].b -= vm_mbfd[index].b;				
		pm_mbfd[i].cu = (float)pm_mbfd[i].c/pm_copy[i].c;
		pm_mbfd[i].ru  = (float)pm_mbfd[i].r/pm_copy[i].r;
	//		printf("pm_mbfd[%d].cu=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].cu);
	//		printf("pm_mbfd[%d].ru=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].ru);
		pm_mbfd[i].rank += 1;
	//	pm_mbfd[i].l = pm_mbfd[index].cu/pm[index].c;
	//	pm_mbfd[i].h  = pm_mbfd[index].ru/pm[index].r;
		vm_mbfd[index].c=0;
	//	printf("vm[%d] was selected\n", vm_mbfd[index].id);
	//	printf("[%d,%d,%d,%d]\n", pm_mbfd[i].id, pm_mbfd[i].c, pm_mbfd[i].r, pm_mbfd[i].b);
	//	printf("--------------\n");
	//	getch();
	}
	}
}
	
///////////////////////////// end of placement ///////////////////////////////////		
	printf("\npm_mbfd --> first pm list :\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d,%d]\n", pm_copy[i].id, pm_copy[i].c, pm_copy[i].r);
	printf("\npm_mbfd --> final pm list :\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d]\n", pm_mbfd[i].id, pm_mbfd[i].c, pm_mbfd[i].r);	
	printf("\n");
	for (i=0; i<npm; i++) {
			if (pm_mbfd[i].c == pm_copy[i].c)
				idle++;
			else
				rw += abs(pm_mbfd[i].c-pm_mbfd[i].r) / (float)((pm[i].c - pm_mbfd[i].c)+(pm[i].r - pm_mbfd[i].r)) + .0001;
			pm_mbfd[i].cu = (float)(pm_copy[i].c-pm_mbfd[i].c)/pm_copy[i].c;
			pm_mbfd[i].ru = (float)(pm_copy[i].r-pm_mbfd[i].r)/pm_copy[i].r;
			printf("pm_mbfd[%d].cu=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].cu);
			printf("pm_mbfd[%d].ru=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].ru);
		//	getch();
	}
	for (i=0; i<npm; i++)
		if (pm_mbfd[i].cu>0) {
			sum_cu += pm_mbfd[i].cu;
			sum_ru += pm_mbfd[i].ru;
			pc += (pm_mbfd[i].p_min+(pm_mbfd[i].p_max-pm_mbfd[i].p_min)*pm_mbfd[i].cu);
		}
	for (i=0; i<nvm; i++) 
		if (vm_mbfd[i].c != 0)
			unaloc_vm++;
	printf("\n\nThe number of unallocated VMs: %d\n",unaloc_vm);
	printf("The number of active PMs: %d\n",npm-idle);
	printf("The mean CPU utilization ratio: %.2f\n",sum_cu/(npm-idle));		
	printf("The mean RAM utilization ratio: %.2f\n",sum_ru/(npm-idle));
	printf("The resource wastage: %.2f\n",rw);
	printf("The power consumption: %.2f W\n",pc);
}

////////////////////// MinPR /////////////////////////
void MinPR(struct machine pm[],struct machine vm[],int npm,int nvm)
{
	int i, j,k,r,c, idle=0, unaloc_vm=0;
	float sum_cu=0, sum_ru=0, rw=0, pc=0;
	struct machine pm_mbfd[npm], pm_copy[npm], vm_mbfd_c[nvm], vm_mbfd_r[nvm], vm_mbfd[nvm];
	for (i=0; i<npm; i++) 
		pm_copy[i] = pm[i];
	pm_decreasing_order(pm_copy, npm);
	for (i=0; i<npm; i++) 
		pm_mbfd[i] = pm_copy[i];
	for (i=0; i<nvm; i++) 
		vm_mbfd[i] = vm[i];
//	for (i=0; i<npm; i++) 
//		pm_mbfd_c[i] = pm[i];
//	for (i=0; i<npm; i++) 
//		pm_mbfd_r[i] = pm[i];
//	pm_decreasing_order_c(pm_mbfd_c, npm);
//	pm_decreasing_order_r(pm_mbfd_r, npm);
	//////////// calculate rank for PMs //////////////
//	for (i=0; i<npm; i++)
//		for (j=0; j<npm; j++) {
//			if (pm_mbfd_c[i].id == pm_mbfd_r[j].id) {
//				pm_mbfd_c[i].rank = i+j+2;
//				break;
//			}
//		}
//	pm_decreasing_order_prob_rank(pm_mbfd_c, npm);
//	for (i=0; i<npm; i++) {
//		pm_mbfd[i] = pm_mbfd_c[i];
//		pm[i] = pm_mbfd[i];
//	}
//	
	printf("\npm_decreasing_order:\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d]\n", pm_mbfd[i].id, pm_mbfd[i].c, pm_mbfd[i].r);
	printf("\nvms:\n");
	for (i=0; i<nvm; i++)
		printf("[%d,%d,%d]\n", vm_mbfd[i].id, vm_mbfd[i].c, vm_mbfd[i].r);
	printf("\n");
///////////////////////// Placement  ////////////////////////////
///* just improve
	int sel;
//	printf("Which one (1: just improve 2: proposed):");
//	scanf("%d",&sel);
	sel=2;
	if (sel==1){
		for (i=0; i<npm; i++) {
			pm_mbfd[i].cu = 0;
			pm_mbfd[i].ru = 0;
			pm_mbfd[i].l = 0;
			pm_mbfd[i].h = 0;
			pm_mbfd[i].v = 0;  // avrage utilization per vm located on a pm
			pm_mbfd[i].rank = 0; // # of VMs on each PM
	}	
	int index;
	float ncu, nru;
	float max_ruf, new_cu, new_ru, new_acu, new_aru;
	for (i=0; i<npm; i++) {
//		getch();
//		printf("\n\npm[%d] is ready to allocation\n",pm_mbfd[i].id);
		for (k=0; k<nvm; k++) {
			max_ruf = -2;
			index = 0;
			for (j=0; j<nvm; j++)
				if (vm_mbfd[j].c != 0 && pm_mbfd[i].c >= vm_mbfd[j].c && pm_mbfd[i].r >= vm_mbfd[j].r)	{
				new_cu = (pm_copy[i].c - pm_mbfd[i].c + vm_mbfd[j].c)/(float)pm_copy[i].c;
				new_ru = (pm_copy[i].r - pm_mbfd[i].r + vm_mbfd[j].r)/(float)pm_copy[i].r;
				new_acu = new_cu/(pm_mbfd[i].rank+1);
				new_aru = new_ru/(pm_mbfd[i].rank+1);
				ncu = pm_mbfd[k].rank+1-new_cu;
				nru = pm_mbfd[k].rank+1-new_ru;
				vm_mbfd[j].v = new_cu*ncu + new_ru*nru;
//				vm_mbfd[j].v = new_cu*new_acu*(1-new_cu) + new_ru*new_aru*(1-new_ru);
//				printf("new_cu, new_ru, new_acu, new_aru, new_v of vm[%d]: %.4f  %.4f  %.4f  %.4f  %.4f\n", vm_mbfd[j].id,new_cu,new_ru,new_acu,new_aru,pm_mbfd[j].v+1);
//				printf("%.4f\n",new_cu*new_acu*(1-new_cu));
//				printf("%.4f\n",new_ru*new_aru*(1-new_ru));
//				printf("ruf of vm[%d]: %f \n", vm_mbfd[j].id, vm_mbfd[j].v);
				
				
				vm_mbfd[j].a =  (new_cu + new_ru)/2 ; // the average of cu and ru
				vm_mbfd[j].a = (pow(vm_mbfd[j].a - new_cu,2)+pow(vm_mbfd[j].a - new_ru,2)) / 2;  // variance
//				printf("var of vm_mbfd[%d] = %.4f \n", vm_mbfd[j].id, vm_mbfd[j].a);
				if (new_cu>=.5 || new_ru>=.5) {
//					printf("penalty for variance=%.2f\n",4*vm_mbfd[j].a);
					vm_mbfd[j].v -= 8*vm_mbfd[j].a;
				}
					
				if (fabs(new_cu - new_ru)>=.5) {
//					printf("penalty for UD:%.2f\n",0.5*vm_mbfd[j].v);
					vm_mbfd[j].v -= 0.5*vm_mbfd[j].v;
				}		
				if ((new_cu + new_ru) >= 1.5) {
//					printf("reward for AD:%.2f\n",0.5*vm_mbfd[j].v);
					vm_mbfd[j].v += 0.25*vm_mbfd[j].v;
				}
				
//				getch();
				if (max_ruf < vm_mbfd[j].v) {
					max_ruf = vm_mbfd[j].v;
					index = j;
				}
			}
//		printf("max_ruf=%f  , index=%d\n",max_ruf, index);
		if (max_ruf>-2) {
			pm_mbfd[i].c -= vm_mbfd[index].c;
			pm_mbfd[i].r -= vm_mbfd[index].r;				
			pm_mbfd[i].cu = 1.0 - (float)pm_mbfd[i].c/pm_copy[i].c;
			pm_mbfd[i].ru  = 1.0 - (float)pm_mbfd[i].r/pm_copy[i].r;
//			printf("pm_mbfd[%d].cu=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].cu);
//			printf("pm_mbfd[%d].ru=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].ru);
			pm_mbfd[i].rank += 1;
		//	pm_mbfd[i].l = pm_mbfd[index].cu/pm[index].c;
		//	pm_mbfd[i].h  = pm_mbfd[index].ru/pm[index].r;
			vm_mbfd[index].c=0;
//			printf("vm[%d] was selected\n", vm_mbfd[index].id);
//			printf("[%d,%d,%d]\n", pm_mbfd[i].id, pm_mbfd[i].c, pm_mbfd[i].r);
//			printf("[%d,%.2f,%.2f]\n", pm_mbfd[i].id, pm_mbfd[i].cu, pm_mbfd[i].ru);
//			printf("--------------\n");
//			getch();
		}
	}
}
}
//*/
	
	
	
	
//	proposed
	if (sel==2){
	int index;
	float max1;
	for (i=0; i<npm; i++) {
//		printf("\n*****************************************\npm[%d,%d,%d] is ready to allocation\n",pm_mbfd[i].id,pm_mbfd[i].c,pm_mbfd[i].r);
//		getch();
		for (k=0; k<nvm; k++) {
			max1 = -2;
			index = 0;
			for (j=0; j<nvm; j++)
				if (pm_mbfd[i].c >= vm_mbfd[j].c && pm_mbfd[i].r >= vm_mbfd[j].r && vm_mbfd[j].c != 0)	{
					vm_mbfd[j].l = 1.0 - (float)(pm_mbfd[i].c-vm_mbfd[j].c)/pm_copy[i].c;
					vm_mbfd[j].h = 1.0 - (float)(pm_mbfd[i].r-vm_mbfd[j].r)/pm_copy[i].r; 
					vm_mbfd[j].a = vm_mbfd[j].l * vm_mbfd[j].h; // cu*ru
//					printf("cu and ru of vm_mbfd[%d] = %.2f,%.2f \n", vm_mbfd[j].id,vm_mbfd[j].l,vm_mbfd[j].h);
//					printf("cu * ru vm_mbfd[%d] = %.2f \n", vm_mbfd[j].id, vm_mbfd[j].a);
					vm_mbfd[j].v =  (vm_mbfd[j].l + vm_mbfd[j].h)/2 ; // the average of cu and ru
//					printf("(1-cu) * (1-ru) vm_mbfd[%d] = %.2f \n", vm_mbfd[j].id, (1-vm_mbfd[j].l)*(1-vm_mbfd[j].h));
					vm_mbfd[j].v=(pow(vm_mbfd[j].v - vm_mbfd[j].l,2)+
					pow(vm_mbfd[j].v - vm_mbfd[j].h,2)) / 2;  // variance
//					printf("var of vm_mbfd[%d] = %.4f \n", vm_mbfd[j].id, vm_mbfd[j].v);
					if (vm_mbfd[j].l>=.75 || vm_mbfd[j].h>=.75) {
						vm_mbfd[j].a -= 4*vm_mbfd[j].v;
//						printf("p variance=%.2f\n",4*vm_mbfd[j].v);
					}
					else {
						vm_mbfd[j].a -= 0*vm_mbfd[j].v;
//						printf("p variance=%.2f\n",0*vm_mbfd[j].v);
					}
					 // the effet of variance
//					vm_mbfd[j].a += (1-vm_mbfd[j].l)*(1-vm_mbfd[j].h);
					
					if (fabs(vm_mbfd[j].l - vm_mbfd[j].h)>=.5) {
//						printf("p UD:%.2f\n",0.5*vm_mbfd[j].a);
						vm_mbfd[j].a -= 0.5*vm_mbfd[j].a;
					}		
					if (vm_mbfd[j].l + vm_mbfd[j].h>=1.5) {
//						printf("r AD:%.2f\n",0.5*vm_mbfd[j].a);
						vm_mbfd[j].a += 0.25*vm_mbfd[j].a;
					}
//	//				if (fabs(vm_mbfd[j].l - vm_mbfd[j].h)>=.5 && (vm_mbfd[j].l>.75 || (vm_mbfd[j].h>.75)))
//					if (fabs(vm_mbfd[j].l - vm_mbfd[j].h)>=.5)
//							vm_mbfd[j].a -= 0.5*vm_mbfd[j].a;
//					if ( (vm_mbfd[j].l>=.75) && (vm_mbfd[j].h>=.75))
//							vm_mbfd[j].a += 0.5*vm_mbfd[j].a;
//					printf("final a vm_mbfd[%d] = %.2f \n\n", vm_mbfd[j].id, vm_mbfd[j].a);
					if (vm_mbfd[j].a > max1) {
						max1 = vm_mbfd[j].a;
						index = j;
					}
				}
			if (max1>-2) {
				pm_mbfd[i].c -= vm_mbfd[index].c;
				pm_mbfd[i].r -= vm_mbfd[index].r;
				pm_mbfd[i].b -= vm_mbfd[index].b;
				vm_mbfd[index].c = 0;
				vm_mbfd[index].l = 0;
				vm_mbfd[index].h = 0;
//				printf("vm[%d] --> pm[%d]\n", vm_mbfd[index].id, pm_mbfd[i].id);
//				printf("pm[%d,%d,%d]\n", pm_mbfd[i].id, pm_mbfd[i].c, pm_mbfd[i].r);
				pm_mbfd[i].cu = (float)(pm_copy[i].c-pm_mbfd[i].c)/pm_copy[i].c;
				pm_mbfd[i].ru = (float)(pm_copy[i].r-pm_mbfd[i].r)/pm_copy[i].r;
//				printf("pm[%d].cu=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].cu);
//				printf("pm[%d].ru=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].ru);
//				printf("--------------\n");
//				getch();
			}
		}
	}
}
//for (i=0; i<npm; i++) {
//		if ( pm_mbfd[i].c == pm_copy[i].c)
//			idle++;
//		else
//			rw += abs(pm_mbfd[i].c-pm_mbfd[i].r) / (float)((pm_copy[i].c - pm_mbfd[i].c)+(pm_copy[i].r - pm_mbfd[i].r)) + .0001;
//		pm_mbfd[i].cu = (float)(pm_copy[i].c-pm_mbfd[i].c)/pm_copy[i].c;
//		pm_mbfd[i].ru = (float)(pm_copy[i].r-pm_mbfd[i].r)/pm_copy[i].r;
//		printf("pm_mbfd[%d].cu=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].cu);
//		printf("pm_mbfd[%d].ru=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].ru);
////		if (i%50==0)
////			getch();
//	}
//	for (i=0; i<npm; i++)
//		if (pm_mbfd[i].cu>0) {
//			sum_cu += pm_mbfd[i].cu;
//			sum_ru += pm_mbfd[i].ru;
//			pc += (pm_mbfd[i].p_min+(pm_mbfd[i].p_max-pm_mbfd[i].p_min)*pm_mbfd[i].cu);
//		}
//	for (i=0; i<nvm; i++) 
//		if (vm_mbfd[i].c != 0)
//			unaloc_vm++;
//	printf("\n\nThe number of unallocated VMs: %d\n",unaloc_vm);
//	printf("The number of active PMs: %d\n",npm-idle);
//	printf("The mean CPU utilization ratio: %.2f\n",sum_cu/(npm-idle));		
//	printf("The mean RAM utilization ratio: %.2f\n",sum_ru/(npm-idle));
//	printf("The resource wastage: %.2f\n",rw);
//	printf("The power consumption: %.2f W\n",pc);
///////////////////////////// end of placement ///////////////////////////////////	
// /*
//	getch();
	idle=0, unaloc_vm=0;
	sum_cu=0, sum_ru=0, rw=0, pc=0;
	int ans;
//	printf("\nDo you want to use migration phase? (1 or 0):");
//	scanf("%d",&ans);
	ans=1;
///////////////////////////// migration ///////////////////////////////////
	if(ans==1) {

	for (i=0; i<npm; i++) {
		pm_mbfd[i].cu = (float)(pm_copy[i].c-pm_mbfd[i].c)/pm_copy[i].c;
		pm_mbfd[i].ru = (float)(pm_copy[i].r-pm_mbfd[i].r)/pm_copy[i].r;
//		printf("pm_mbfd[%d].cu=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].cu);
//		printf("pm_mbfd[%d].ru=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].ru);
	}
	for(i=0; i<npm; i++) {
		if (pm_mbfd[i].cu + pm_mbfd[i].ru>=1.5)
			pm_mbfd[i].b=0;
		else if (fabs(pm_mbfd[i].cu - pm_mbfd[i].ru)>=.5)	
			pm_mbfd[i].b=1;
		else
			pm_mbfd[i].b=2;
		}
	for(i=0; i<npm; i++) {
//		printf("The domain of pm[%d] is: %d",pm_mbfd[i].id,pm_mbfd[i].b);
//		printf("\tpm[%d,%d,%d]\n",i,pm_mbfd[i].c,pm_mbfd[i].r);
	}
	for(i=0; i<npm; i++) {
		if (pm_mbfd[i].b==1) // 
			for (j=0; j<npm; j++)
				if (pm_mbfd[j].c==3720 && pm_mbfd[j].c >= (pm_copy[i].c - pm_mbfd[i].c)) {
//					printf("pm[%d] is ready to migrate to pm[%d]\n",pm_mbfd[i].id, pm_mbfd[j].id);
//					printf("before migration: pm[%d,%d,%d]   pm[%d,%d,%d]\n",pm_mbfd[i].id,pm_mbfd[i].c,pm_mbfd[i].r,pm_mbfd[j].id,pm_mbfd[j].c,pm_mbfd[j].r);
					pm_mbfd[j].c = pm_mbfd[j].c - (pm_copy[i].c - pm_mbfd[i].c);
					pm_mbfd[j].r = pm_mbfd[j].r - (pm_copy[i].r - pm_mbfd[i].r);
					pm_mbfd[i].c = pm_copy[i].c;
					pm_mbfd[i].r = pm_copy[i].r;
//					printf("after migration: pm[%d,%d,%d]   pm[%d,%d,%d]\n\n",pm_mbfd[i].id,pm_mbfd[i].c,pm_mbfd[i].r,pm_mbfd[j].id,pm_mbfd[j].c,pm_mbfd[j].r);
					pm_mbfd[i].cu = (float)(pm_copy[i].c-pm_mbfd[i].c)/pm_copy[i].c;
					pm_mbfd[i].ru = (float)(pm_copy[i].r-pm_mbfd[i].r)/pm_copy[i].r;
					if (pm_mbfd[i].cu + pm_mbfd[i].ru>=1.5)
						pm_mbfd[i].b=0;
					else if (fabs(pm_mbfd[i].cu - pm_mbfd[i].ru)>=.5)	
						pm_mbfd[i].b=1;
					else
						pm_mbfd[i].b=2;
//					getch();
					break;	
				}
//		getch();	
	} 
	}
//	printf("\nThe domain of PMs after migration phase\n");
//	for (i=0; i<npm; i++) {
//		pm_mbfd[i].cu = (float)(pm_copy[i].c-pm_mbfd[i].c)/pm_copy[i].c;
//		pm_mbfd[i].ru = (float)(pm_copy[i].r-pm_mbfd[i].r)/pm_copy[i].r;
////		printf("pm_mbfd[%d].cu=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].cu);
////		printf("pm_mbfd[%d].ru=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].ru);
//	}
//	for(i=0; i<npm; i++) {
//		if (pm_mbfd[i].cu + pm_mbfd[i].ru>=1.5)
//			pm_mbfd[i].b=0;
//		else if (fabs(pm_mbfd[i].cu - pm_mbfd[i].ru)>=.5)	
//			pm_mbfd[i].b=1;
//		else
//			pm_mbfd[i].b=2;
//		}
//	for(i=0; i<npm; i++) {
//		printf("The domain of pm[%d] is: %d",pm_mbfd[i].id,pm_mbfd[i].b);
//		printf("\tpm[%d,%d,%d]\n",i,pm_mbfd[i].c,pm_mbfd[i].r);
//	}
//	*/
///////////////////////////// end of migration ///////////////////////////////////

	
	printf("\npm_mbfd --> first pm list :\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d,%d]\n", pm[i].id, pm[i].c, pm[i].r, pm[i].b);
	printf("\npm_mbfd --> final pm list :\n");
	for (i=0; i<npm; i++)
		printf("[%d,%d,%d,%d]\n", pm_mbfd[i].id, pm_mbfd[i].c, pm_mbfd[i].r, pm_mbfd[i].b);	
	for (i=0; i<npm; i++) {
		if ( pm_mbfd[i].c == pm_copy[i].c)
			idle++;
		else
			rw += abs(pm_mbfd[i].c-pm_mbfd[i].r) / (float)((pm_copy[i].c - pm_mbfd[i].c)+(pm_copy[i].r - pm_mbfd[i].r)) + .0001;
		pm_mbfd[i].cu = (float)(pm_copy[i].c-pm_mbfd[i].c)/pm_copy[i].c;
		pm_mbfd[i].ru = (float)(pm_copy[i].r-pm_mbfd[i].r)/pm_copy[i].r;
		printf("pm_mbfd[%d].cu=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].cu);
		printf("pm_mbfd[%d].ru=%.2f\n",pm_mbfd[i].id, pm_mbfd[i].ru);
//		if (i%50==0)
//			getch();
	}
	for (i=0; i<npm; i++)
		if (pm_mbfd[i].cu>0) {
			sum_cu += pm_mbfd[i].cu;
			sum_ru += pm_mbfd[i].ru;
			pc += (pm_mbfd[i].p_min+(pm_mbfd[i].p_max-pm_mbfd[i].p_min)*pm_mbfd[i].cu);
		}
	for (i=0; i<nvm; i++) 
		if (vm_mbfd[i].c != 0)
			unaloc_vm++;
	printf("\n\nThe number of unallocated VMs: %d\n",unaloc_vm);
	printf("The number of active PMs: %d\n",npm-idle);
	printf("The mean CPU utilization ratio: %.2f\n",sum_cu/(npm-idle));		
	printf("The mean RAM utilization ratio: %.2f\n",sum_ru/(npm-idle));
	printf("The resource wastage: %.2f\n",rw);
	printf("The power consumption: %.2f W\n",pc);
}


//////////////////////////////  Functions  //////////////////////////////
void pm_increasing_order(struct machine pm[], struct machine pm_copy[], int npm)
{
	int i , j;
	struct machine temp;
	for (i=0; i<npm-1; i++)
		pm[i].a = (pm[i].cu + pm[i].ru)/2;
	for (i=0; i<npm-1; i++)
		for (j=i+1; j<npm; j++)
			if (pm[i].b == 1 && pm[j].b ==1 && pm[i].a > pm[j].a)
			{
				temp = 	pm[i];
				pm[i] = pm[j];
				pm[j] = temp;
				
				temp = 	pm_copy[i];
				pm_copy[i] = pm_copy[j];
				pm_copy[j] = temp;
			}
//	printf("inside sorting function\n");	
//	for (i=0; i<npm; i++)
//		printf("pm[%d,%d,%d]\n",pm[i].id, pm[i].c, pm[i].r);
//	printf("\n");
}
void pm_decreasing_order(struct machine pm[], int npm)
{
	int i , j;
	struct machine temp;
	for (i=0; i<npm-1; i++)
		for (j=i+1; j<npm; j++)
			if (pm[i].c < pm[j].c)
			{
				temp = 	pm[i];
				pm[i] = pm[j];
				pm[j] = temp;
			}	
}
void vm_decreasing_order(struct machine vm[], int nvm)
{
	int i , j;
	struct machine temp;
	for (i=0; i<nvm-1; i++)
		for (j=i+1; j<nvm; j++)
			if (vm[i].c < vm[j].c)
			{
				temp = 	vm[i];
				vm[i] = vm[j];
				vm[j] = temp;
			}
}



