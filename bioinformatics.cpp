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

char int_to_base(int i)
{
	char c;
	switch(i)
	{
		case 0:
			c = 'A';
			break;
		case 1:
			c = 'C';
			break;
		case 2:
			c= 'G';
			break;
		case 3:
			c = 'T';
			break;

	}
	return c;
}

/*
 * t sequences of each of lenght = l
 */
int len = 0;
char** DNA(int t)
{
	char** dna = new char*[t];
	for(int i=0;i < t; i++)
	{
		dna[i] = new char[3000];
		char c = '#';
		int j=0;
		len = 0;
		while(c != '\n')
		{

			scanf("%c",&c);
			len++;
			dna[i][j] = c;
			j++;
		}

	}

//	dna[0] = "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA";
//	dna[1] = "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG";
//	dna[2] = "TAGTACCGAGACCGAAAGAAGTATACAGGCGT";
//	dna[3] = "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC";
//	dna[4] = "AATCCACCAGCTCCACGTGCAATGTTGGCCTA";

	return dna;


}

char** DNA2(int t)
{
	char** dna = new char*[t];
	for(int i=0;i < t; i++)
	{
		dna[i] = new char[30000];
		string s;
		cin >> s;
		int j=0;
		len = 0;
		while(s[j] != '\0')
		{
			dna[i][j] = s[j];
			len++;
			j++;
		}

	}
//	dna[0] = "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA";
//	dna[1] = "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG";
//	dna[2] = "TAGTACCGAGACCGAAAGAAGTATACAGGCGT";
//	dna[3] = "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC";
//	dna[4] = "AATCCACCAGCTCCACGTGCAATGTTGGCCTA";

	return dna;


}

int **pseudo_count(char **dna, int *motifs, int t, int k, int excluded_motif)
{
	int **p_count = new int*[4];
	for(int i=0; i<4; i++)
	{
		p_count[i] = new int[k];
	}

	for(int i=0; i<4; i++)
	{
		for(int j=0; j<k; j++)
			p_count[i][j] = 1;
	}

	//cout << "pseduo" << endl;
	for(int i=0; i< t; i++)
	{
		if(i== excluded_motif)
			continue;
		for(int j=0; j<k;j++)
		{
			char base =  dna[i][motifs[i]+j] ;
			p_count[bases_int(base)][j]++;
		}
	}

	return p_count;

}

float **profile(int **p_count, int t, int k)
{
	float **prof = new float*[4];
	for(int i=0; i<4; i++)
	{
		prof[i] = new float[k];
		for(int j=0; j<k; j++)
		{
			prof[i][j] = ((float)p_count[i][j])/((float)(t-1)+ (float)4);
		}
	}

	return prof;
}

int gibs_smapleing(float *Pr, int len)
{
	float sum = 0;
	for(int i=0; i< len; i++)
		sum += Pr[i];

	int percentage[len];
	int sum_per = 0;
	for(int i=0; i< len; i++)
	{
		percentage[i] = (Pr[i]/sum)*100;
		sum_per += percentage[i];
	}

	int distribution[100];
	int index = 0;
	for(int i=0; i< len; i++)
	{
		int times = percentage[i];
		while(times-- > 0)
		{
			distribution[index] = i;
			index++;
		}
	}

	//adjust the error at distribution[99]
	for(int i= sum_per; i < 100; i++)
	{
		distribution[i] = distribution[(rand() % sum_per)];
	}

	if(test)
	{
		cout << "prob:" << endl;
		for(int i=0; i< len;i++)
			cout << Pr[i]/sum << " ";


		cout << endl << "percentage:" << endl;
		for(int i=0; i< len;i++)
			cout << percentage[i] << " ";

		cout << endl << "dist:" << endl;
		for(int i=0; i< 100; i++)
			cout << distribution[i] << " ";
	}


	return distribution[(rand() % 100)];

}

int profile_rand_generate_k_mer(char **dna,float** profile, int excluded_motif,int t, int l, int k )
{
	//cout << "k-mers" << endl;
	//storing probability at each ith k-mer in excluced motif
	float *prob = new float[l-k+1];
	for(int i=0; i<l-k+1; i++)
	{
		//pr at ith k-mer
		float pr_i = 1;
		for(int j=0; j<k; j++)
		{
			//cout << dna[excluded_motif][j + i];
			pr_i *= profile[bases_int(dna[excluded_motif][j + i])][j];
		}
		//cout << endl;
		prob[i] = pr_i;
	}


//	for(int i=0; i<l-k+1; i++)
//	{
//		cout << prob[i] << " ";
//	}

	return gibs_smapleing(prob, l-k+1);

}
char* consensus(float **profile, int k)
{
	char *cons_str = new char[k];
	for(int i=0; i < k; i++)
	{
		int max_row_at_i_col = 0;
		for(int rows = 0; rows < 4; rows++)
		{
			if(profile[rows][i] > profile[max_row_at_i_col][i])
			{
				max_row_at_i_col = rows;
			}
		}

		cons_str[i] = int_to_base(max_row_at_i_col);
	}

	return cons_str;
}
int score(float** prof,int *motifs, char** dna, int t,int k)
{
	char *cons_str = consensus(prof,k);
	int mismatch = 0;
	for(int i = 0; i < t; i++ )
	{
		for(int cols = 0; cols < k; cols++)
		{
			char base = dna[i][motifs[i]+cols];
			if(base != cons_str[cols])
			{
				//cout << base << " ";
				mismatch++;
			}
		}

		//cout << endl;
	}

	return mismatch;
}

void print_consensus(float** profile, int k)
{
	char *str = consensus(profile,k);
	cout << str;
}
void print_prf_gen_mer(float *prob,int k)
{
	for(int i=0; i<k; i++)
	{
		cout << prob[i] << " ";
	}

}
void print_dna(char **dna, int l,int t)
{
	for(int i=0; i<t; i++)
	{
		cout << dna[i];
		cout << endl;
	}


}
void print_motifs(char **dna, int *motifs, int t, int k)
{
	for(int i=0; i< t; i++)
	{
		for(int j=0; j<k;j++)
		{
			cout << dna[i][motifs[i]+j] ;
		}
		cout << endl;
	}
}
void print_count(int **c,int k)
{
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<k; j++)
		{
			cout << c[i][j] << " ";
		}
		cout << endl;
	}
}

void print_profile(float **prof,int k)
{
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<k; j++)
		{
			cout << prof[i][j] << " ";
		}
		cout << endl;
	}
}

void test_gibs()
{
	cout << endl;
	float Pr[5] = {0.1,0.03,0.002,0.5,0.7};

	for(int i=0; i< 10; i++)
	{
		int random_motif = gibs_smapleing(Pr,5);
		cout << endl << random_motif;
	}
}

int main() {
	/*
	 * l length of each seq
	 * k-mer motif
	 * t number of seqs
	 * N number of iterations
	 */

	int k = 8;
	int t = 5;
	int N = 5000;

	cin >> k >> t >> N;//scanf("%d%d%d\n",&k,&t,&N);
	N = 5000;
	int random_starts = 20;
	srand(time(0));

	char **dna = DNA2(t);
	int l = len;
	//print_dna(dna,l,t);

	int best_motifs[t]= {6,3,0,1,6};
	int best_motifs_score = 50000;//= score(prof, best_motifs, dna, t, k);
	//initially random k-mer best_motifs.
	for(int i=0; i < t; i++)
	{
		int temp = (rand() % (l-k+1));
		best_motifs[i] = temp;
	}


	for(int i=1; i <= random_starts; i++ )
	{
		//motif to store start position of motif
		int motifs[t] = {6,3,0,1,6};
		//print_dna(dna,l,t);
		//cout << endl;
		//print_motifs(dna,motifs,t,k);

		//initially random k-mer motifs.
		for(int i=0; i < t; i++)
		{
			int temp = (rand() % (l-k+1));
			motifs[i] = temp;
		}

		for(int i=1; i<=N; i++)
		{
			int excluded_motif = (rand() % t);
			int **p_count = pseudo_count(dna, motifs,t,k, excluded_motif);
			//print_count(p_count,k);

			float **prof = profile(p_count,t,k);
			//print_profile(prof,k);

			int rand_motif_i = profile_rand_generate_k_mer(dna, prof, excluded_motif, t, l, k);
			motifs[excluded_motif] = rand_motif_i;

	//		print_profile(prof, k);
	//		cout << endl;
	//		print_consensus(prof, k);
			int motifs_score = score(prof, motifs, dna, t, k);


			if(motifs_score < best_motifs_score)
			{

				for(int i=0; i< t; i++)
				{
					best_motifs[i] = motifs[i];
					best_motifs_score = motifs_score;
				}

				//cout << best_motifs_score << " ";
			}


		}

	}
	print_motifs(dna, best_motifs, t, k);




	return 0;
}
