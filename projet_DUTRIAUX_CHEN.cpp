#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <glpk.h>

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>


using namespace std;

struct timeval start_utime, stop_utime;

void crono_start() {
	struct rusage rusage;
	getrusage(RUSAGE_SELF, &rusage);
	start_utime = rusage.ru_utime;
}

void crono_stop() {
	struct rusage rusage;
	
	getrusage(RUSAGE_SELF, &rusage);
	stop_utime = rusage.ru_utime;
}

double crono_ms() {
	return (stop_utime.tv_sec - start_utime.tv_sec) * 1000 +
    (stop_utime.tv_usec - start_utime.tv_usec) / 1000 ;
}


typedef struct {
	vector<int> vec;
	int longueur;
} tournee;

typedef struct {
	int nblieux; /* Nombre de lieux (incluant le dépôt) */
	int capacite; /* Capacité du véhicule de livraison */
	int *demande; /* Demande de chaque lieu (la case 0 est inutilisée car le dépôt n'a aucune demande à voir satisfaire) */
	int **C; /* distancier (les lignes et colonnes 0 correspondent au dépôt) */
} donnees;

/* lecture des donnees */

void lecture_data(char *file, donnees *p)
{
	int i,j;
	FILE *fin;
	
	int val;
	fin = fopen(file,"rt");
	
	/* Lecture du nombre de villes */
	
	fscanf(fin,"%d",&val);
	p->nblieux = val;

	/* Allocation mémoire pour la demande de chaque ville, et le distancier */
	
	p->demande = (int *) malloc (val * sizeof(int));
	p->C = (int **) malloc (val * sizeof(int *));
	for(i = 0;i < val;i++) p->C[i] = (int *) malloc (val * sizeof(int));
	
	/* Lecture de la capacité */
	
	fscanf(fin,"%d",&val);
	p->capacite = val;
	
	/* Lecture des demandes des clients */
	
	for(i = 1;i < p->nblieux;i++)
	{
		fscanf(fin,"%d",&val);
		p->demande[i] = val;
	}
	
	/* Lecture du distancier */

	for(i = 0; i < p->nblieux; i++)
		for(j = 0; j < p->nblieux; j++)
		{
			fscanf(fin,"%d",&val);
			p->C[i][j] = val;
		}
		
	fclose(fin);
}

/* Fonction de libération mémoire des données */

void free_data(donnees *p)
{
	int i;
	for(i = 0;i < p->nblieux;i++) free(p->C[i]);
	free(p->C);
	free(p->demande);	
}

vector<tournee>* enumeration(donnees *p, vector<tournee> *enumerate, vector<int> vec, int cursor, int cap);
string tableau_tournee(vector<tournee> v);
string tableau(vector<int> v);
bool etat_terminal(vector<int> *v, int n);
tournee longmin(donnees *p, vector<int> *v);
int calcul_longueur(donnees *p, vector<int> *v);
void crono_start();
void crono_stop();
double crono_ms();



int main(int argc, char *argv[]) {

	donnees p;
	double temps;
	vector<tournee> e;

	lecture_data(argv[1],&p);

	crono_start();

	enumeration(&p, &e, vector<int>(), 1, 0);

	crono_stop();
	temps = crono_ms()/1000,0;

	printf("temps : %f\n", temps);


	/**************
	 * GLPK BELOW *
	 **************/
	
	crono_start();
	
	// Creation du problème GLPK
	glp_prob *prob;
	prob = glp_create_prob();
	glp_set_prob_name(prob, "resolution");
	glp_set_obj_dir(prob, GLP_MIN);

	// Informations initiales
	int nombVar = e.size();
	int nombContraintes = p.nblieux - 1; 
	
	// Tableaux de contraintes
	vector<int> tab1;
	vector<int> tab2;
	vector<double> tab3;


	char nomContraintes[nombVar+1][20];
	char nomVariables[nombContraintes+1][20];
	double z;
	double x[nombVar+1];

	// Contraintes
	glp_add_rows(prob, nombContraintes);

	for (int i=1; i<=nombContraintes; i++) {
		sprintf(nomContraintes[i], "constraint n°%d", i);
		glp_set_row_name(prob, i, nomContraintes[i]);
		glp_set_row_bnds(prob, i, GLP_FX, 1.0, 0.0); // Equality constraint
	}

	// Variables
	glp_add_cols(prob, nombVar);

	for (int i=1; i<=nombVar; i++) {
		sprintf(nomVariables[i], "x%d", i);
		glp_set_col_name(prob, i, nomVariables[i]);
		glp_set_col_bnds(prob, i, GLP_DB, 0.0, 1.0); // Bounds
		glp_set_col_kind(prob, i, GLP_BV); // Binary
	}

	// Objectif
	for (int i=1; i<=nombVar; i++) {
		glp_set_obj_coef(prob, i, e[i-1].longueur);
	}

	// Matrice des contraintes
	tab1.push_back(0);
	tab2.push_back(0);
	tab3.push_back(0.0);
	for (int i=0; i<nombVar; i++) {
		for (unsigned int j=0; j<e[i].vec.size(); j++) { 
			tab1.push_back(e[i].vec[j]);
			tab2.push_back(i+1);
			tab3.push_back(1.0);
		}
	}

	glp_load_matrix(prob, tab1.size()-1, tab1.data(), tab2.data(), tab3.data()); // Use of data() because GLPK cannot deal with vectors


	// Résolution
	glp_simplex(prob, NULL);
	glp_intopt(prob, NULL);


	z = glp_mip_obj_val(prob); 
	for (int i=0; i<nombVar; i++) { 
		x[i] = glp_mip_col_val(prob, i+1);
	}

	printf("z= %lf\n",z);
	for (int i=0; i<nombVar; i++) {
		if ((int)(x[i] + 0.5) == 1)
			printf("x%d = %d, ", i, (int)(x[i] + 0.5));
	}

	crono_stop();
	temps = crono_ms()/1000,0;



	
	/************
	 * END GLPK *
	 ************/
		
	printf("Temps GLPK : %f\n", temps);

	free_data(&p);

	return 0;
}


vector<tournee>* enumeration(donnees *p, vector<tournee> *enume, vector<int> vec, int curseur, int capa) {

	if (etat_terminal(&vec, p->nblieux - 1)) {
		return enume;

	} else if (curseur > p->nblieux - 1) {
		curseur = vec.back();
		vec.pop_back();
		capa -= p->demande[curseur];

	} else if (p->demande[curseur] + capa <= p->capacite) {
		vec.push_back(curseur);
		enume->push_back(longmin(p, &vec));
		capa += p->demande[curseur];

	} return enumeration(p, enume, vec, curseur + 1, capa);
}

// Retourne le chemin le plus court
tournee longmin(donnees *p, vector<int> *v) {

	tournee t;
	int tmp = 0, min = numeric_limits<int>::max();
	
	do {
		tmp = calcul_longueur(p, v);
		if (tmp < min) {
			min = tmp;
			t.vec = *v;
			t.longueur = min;
		}
	} while (next_permutation(v->begin(),v->end()));
		
	return t;
}


int calcul_longueur(donnees *p, vector<int> *v) {
	
	int longu = p->C[0][v->at(0)];

	for (unsigned int i = 1; i < v->size(); i++) {
		longu += p->C[v->at(i-1)][v->at(i)];
	}

	longu += p->C[v->back()][0];
	
	return longu;
}

// Est etat terminal ?
bool etat_terminal(vector<int> *vec, int n) {


	if (vec->size() > 0) {
		return vec->front() == n;
	}
	return false;
}

// Retourne le tableau tournee en chaine
string tableau_tournee(vector<tournee> tabTournee) {

	string cdc;

	for (tournee & tournee : tabTournee) {
		cdc += to_string(tournee.longueur) + " ";
		cdc += tableau(tournee.vec) + "\n";
	}
	
	return cdc;
}

// Retourne un tableau tournee en chaine
string tableau(vector<int> v) {

	string cdc;

	for (int & i : v) {
		cdc += to_string(i) + " ";
	}

	return cdc;
}




