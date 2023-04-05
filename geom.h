//
// Created by David on 06/02/2023.
//

#ifndef SIM2_GEOM_H
#define SIM2_GEOM_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matriz.h"
#include <stdbool.h>
#include "jacobi.h"

#define TOL_GEOM 1E-4
#define TOL_DOUBLES 1E-4
#define MAX_ORDER 6


struct atom{
    double pos[3];
    double M;// masa atómica
    int Z; // número atómico;
    int id_atom; //Identificador del átomo
};
typedef struct atom atom;

struct mol{
    atom* arr_atoms; //array de átomos
    int N; //número de átomos
};
typedef struct mol mol;
struct elem_sim{
    double vec_dic[3];
    int ref;
    int order;
};
typedef struct elem_sim elem_sim;
void leer_geom(char* nom_fichero, mol* molecula);
double M_a (char* simbolo);
int Z (char* simbolo);
double dist (double* vecA, double* vecB);
void cons_MDIN (matrix MDIN, mol molecula);
int comp_SEA(int atom1, int atom2, matrix MDIN);
int comp_doubles (double val1, double val2, double toler);
void ad_atom_a_mol (mol* molecula, atom atomo);
int cons_sets_SEA(mol** array_SEAS, mol molecula);
void impr_mol (mol molecula);
void cons_tensor_inerc(matrix tensor, mol molecula);
void realloc_SEAS (mol** array_SEAS, size_t nuevo_tam);
void calc_cdm(mol molecula, double* r_cm);
void origen_cdm(double* cdm, mol molecula);
bool test_equiv_mol (mol molA, mol molB);
void rot_mol (mol molA, mol* molB, double* eje, double theta);
int hallar_ejes_rot_prop(elem_sim** array_elems, mol set_SEAS, mol molecula, int* n_total_elems,double* CDM_SET_REF_CDM_MOL);
struct vector {
    size_t m; //Dimensión del vector
    double* vector;
};
typedef struct vector vector;
void vector_transf (matrix T, vector A, vector A_T);
void vector_inic(vector* Vec, size_t m);
void vector_norm (vector V);
void cop_col_a_vect (matrix Mat, double* V, int col);
int hallar_planos_ref (elem_sim** array_elems, mol set_SEAS, mol molecula,int* n_total_elems);
void rest_vector (double* VecA, double* VecB, double* VecC);
void reflex_mol(mol molA, mol* molB, double* v_normal);
void normaliza_vec (double* Vec);
void impr_elems_sim(elem_sim* arr_elems, int num_total_elems_sim);
double prod_esc (double* VecA, double* VecB);
int testSIMELS (elem_sim** array_elems, mol set_SEAS, mol molecula,int* n_total_elems);
void free_mol(mol* molecula);
int alloc_mol(int N_at, mol* molecula);
void copia_mol (mol destino, mol origen);
void rest_vect (double vecA[3], double vecB[3], double vecR[3]);
int hallar_c2_perps(elem_sim** array_elems, mol set_SEAS, mol molecula, int* n_total_elems);
int hallar_inv (elem_sim** array_elems,mol molecula,int* n_total_elems);
int pos_arr_cmax(elem_sim** array_elems, int n_total_elems);
void prodvect (double* vecA, double* vecB, double* vecP);
void sum_vector (double* VecA, double* VecB, double* VecC);
#endif //SIM2_GEOM_H
