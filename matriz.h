//
// Created by David on 06/02/2023.
//

#ifndef SIM2_MATRIZ_H
#define SIM2_MATRIZ_H

#define ERR_ALG "ERROR DE ALGEBRA MATRICIAL\nEL PROGRAMA SE CERRARA"
#define MAX_COLS_IMPRIMIR 6

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defs.h"


struct matrix {
    int m; //Numero de filas
    int n; //Numero de columnas
    double** datos;
};
typedef struct matrix matrix;


void matriz_inicializar(matrix* matriz, int m, int n);
void matriz_imprimir(matrix matriz);
void matriz_sumar (matrix* matrizA, matrix* matrizB, matrix* matrizS);
void matriz_copiar(matrix matrizA, matrix matrizC);
void matriz_trasponer(matrix matriz, matrix matriz_traspuesta);
void matriz_liberar(matrix matriz);
void matriz_producto (matrix matrizA, matrix matrizB, matrix matrizC);
void matriz_nula(matrix matriz);
void matriz_identidad(matrix matriz);
void matriz_restar(matrix *matrizA, matrix *matrizB, matrix *matrizS);
double matriz_traza (matrix matrizA);
void num_mat_prod (double num, matrix A, matrix B);
double norma_euclidea (matrix mat);


#endif //SIM2_MATRIZ_H
