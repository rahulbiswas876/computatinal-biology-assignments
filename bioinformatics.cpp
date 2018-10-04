//============================================================================
// Name        : bioinformatics.cpp
// Author      : rahul
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iostream>
#include <time.h>
#include <fstream>
#include <math.h>
#include <string.h>
#include <bits/stdc++.h>

#define test false

using namespace std;

int bases_int(char base)
{
	if(base == 'A')
		return 0;
	else if(base == 'C')
		return 1;
	else if(base == 'G')
		return 2;
	else
		return 3;
}

/*
 * t sequences of each of lenght = l
 */
char** DNA(int l, int t)
{
	char** dna = new char*[t];
	for(int i=0;i < t; i++)
	{
		dna[i] = new char[l];

	}


}

int main() {
	//cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	int n;
	int k;
	int t;
	//motif to store start position of motif
	int motifs[t];
	int best_motifs[t];
	//initially random k-mer motifs.intially both best_motif and motifs are same
	for(int i=0; i < t; i++)
	{
		int temp = (rand() % n-k+1) + 1;
		motifs[i] = temp;
		best_motifs[i] = temp;
	}


	return 0;
}
