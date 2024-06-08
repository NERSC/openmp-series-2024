#include <iostream>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <stdio.h>
using namespace std;

void init(int* array,int length){
	srand((unsigned)time(0));
	for(int i=0;i<length;i++)
		array[i]=rand()%1234 + 1;
}

void swap(int &a, int &b){
	int temp=a;
	a=b;
	b=temp;
}


int pivot(int* array, int first, int last){
	int p = first;
	int pivotElement = array[first];
	for(int i=first+1;i<=last;i++){
		if(array[i]<=pivotElement){
			p++;
			swap(array[i],array[p]);
		}
	}
	swap(array[p],array[first]);
	return p;
}

void quicksort(int * array, int first, int last){
	int pivotElement;
	if(first<last){
		pivotElement = pivot(array,first,last);
		quicksort(array,first,pivotElement-1);
		quicksort(array,pivotElement+1,last);
	}
}

bool checkFn(int * array,int length){
	for(int i=0;i<length-1;i++){
		if(array[i]>array[i+1]){
			cout << "array[" << i << "] > array[" << i+1 << "]" << endl;
			return false;
		}
	}
	return true;
}

int main(int argc, char *argv[]){
	int length;
	int *toBeSorted;
	double t1,t2;
	if(argc < 2){
		length=6000000;
	}else{
		length=atoi(argv[1]);
	}
	toBeSorted = new int[length];

	init(toBeSorted, length);
	cout << "Sorting an array of " << length << " elements." << endl;	

	t1=omp_get_wtime();
	quicksort(toBeSorted,0,length-1);
	t2=omp_get_wtime()-t1;

	cout << "quicksort took " << t2 << " sec. to complete" << endl;
	if (!checkFn(toBeSorted, length)) {
		cout << "validation failed!" << endl;
	}
	delete [] toBeSorted;
}
