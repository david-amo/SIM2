//
// Created by David on 02/09/2021.
//

#ifndef RWMINIFOCK_JACOBI_H
#define RWMINIFOCK_JACOBI_H
#include "matriz.h"
#include <math.h>
#include <time.h>
#include <stdbool.h>


#define TOLERANCIA_JACOBI 1E-12

#define ERR_MAT_DIAG "ADVERTENCIA: LA MATRIZ DE ENTRADA YA ES DIAGONAL. SE OMITE LA EJECUCION DE LA SUBRUTINA DE DIAGONALIZACION DE JACOBI\n"
struct elemento_matriz {
    long double valor;
    int i;
    int j;
};
typedef struct elemento_matriz elemento_matriz;


struct elemento_ordenacion{
    long double valor;
    int indice_original;
};
typedef struct elemento_ordenacion elemento_ordenacion;


int eigen_jacobi(matrix matriz, double* array_eigenvalores, matrix matriz_eigenvectores, bool diagnostico);

elemento_matriz encontrar_max (matrix matriz);

int algoritmo_ordenacion (const void * a, const void * b);

#endif //RWMINIFOCK_JACOBI_H
