#include "jacobi.h"


#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

int eigen_jacobi(matrix matriz, double* array_eigenvalores, matrix matriz_eigenvectores, bool diagnostico) {
    double s, c, t; /*Variables trigonometricas (seno, coseno y tangente, respectivamente)*/
    double alpha;
    int p, q; /* Variables para recorrer matrices (indices de fila y columna, respectivamente)*/
    int N = matriz.m; /*Dimension de la matriz cuadrada y simetrica de entrada*/
    int r; /*Variable auxiliar para recorrer arrays*/
    int contador_iteraciones = 0;
    double *auxiliar_p;
    double *auxiliar_q;
    auxiliar_p = malloc(N * sizeof(double));
    auxiliar_q = malloc(N * sizeof(double));
    clock_t inicio, final;
    double tiempo_cpu_diag;

    matrix matriz_trabajo;
    matriz_inicializar(&matriz_trabajo,N,N);
    matriz_copiar(matriz,matriz_trabajo);

    matrix eigenvectores_trabajo;
    matriz_inicializar(&eigenvectores_trabajo,N,N);
    matriz_identidad(eigenvectores_trabajo);

    matriz_identidad(matriz_eigenvectores);
    if (diagnostico == true) {
        printf(ANSI_COLOR_GREEN"Comenzando rutina de diagonalizacion de Jacobi con tolerancia %.1e.\n"ANSI_COLOR_RESET,
               TOLERANCIA_JACOBI);
        printf("***MATRIZ DE ENTRADA***\n");
        matriz_imprimir(matriz_trabajo);
        printf("\n");
    }
    inicio = clock();
    do {
        for (p = 1; p < N; p++) {
            for (q = 0; q < p; q++) {
                if (fabsl(matriz_trabajo.datos[p][q]) > TOLERANCIA_JACOBI) {
                    /*Calculo de los parametros trigonometricos requeridos para efectuar la rotacion de Jacobi.
                     * (s,c y t son, respectivamente, el seno, el coseno y la tangente del angulo de la rotacion
                     * que anula el elemento matricial Apq) */
                    alpha = (matriz_trabajo.datos[q][q] - matriz_trabajo.datos[p][p]) / (2 * matriz_trabajo.datos[p][q]);
                    t = alpha - SIGNO(alpha) * sqrtl(powl(alpha, 2) + 1);
                    c = powl((powl(t, 2) + 1), -0.5);
                    s = t * c;

                    /*Antes de la transformacion, se copian los elementos de las columnas p y q en sendos arrays
                     * auxiliares destinados a tal efecto*/
                    for (r = 0; r < N; r++) {
                        auxiliar_p[r] = matriz_trabajo.datos[r][p];
                        auxiliar_q[r] = matriz_trabajo.datos[r][q];
                    }

                    /*Se actualizan los elementos App,Aqq,Apq,Aqp teniendo en cuenta las ecuaciones de transformacion
                     * explicitadas en la bibliografia y la anulacion de los dos ultimos como consecuencia de la rotacion*/
                    matriz_trabajo.datos[p][p] = auxiliar_p[p] + t * auxiliar_q[p];
                    matriz_trabajo.datos[p][q] = 0;
                    matriz_trabajo.datos[q][p] = 0;
                    matriz_trabajo.datos[q][q] = auxiliar_q[q] - t * auxiliar_q[p];

                    /*Finalmente, se actualizan el resto de elementos ajenos la diagonal principal*/
                    for (r = 0; r < N; r++) {
                        if (r != p && r != q) {
                            /*Transformacion de los elementos del triangulo inferior*/
                            matriz_trabajo.datos[r][p] = c * auxiliar_p[r] + s * auxiliar_q[r];
                            matriz_trabajo.datos[r][q] = c * auxiliar_q[r] - s * auxiliar_p[r];
                            /*Los elementos del triangulo superior se actualizan por simetria*/
                            matriz_trabajo.datos[p][r] = matriz_trabajo.datos[r][p];
                            matriz_trabajo.datos[q][r] = matriz_trabajo.datos[r][q];
                        }
                    }
                    // Actualizacion de los eigenvectores;
                    for (r = 0; r < N; r++) {
                        auxiliar_p[r] = eigenvectores_trabajo.datos[r][p];
                        auxiliar_q[r] = eigenvectores_trabajo.datos[r][q];
                    }
                    for (r=0;r<N;r++){
                        eigenvectores_trabajo.datos[r][p]= c * auxiliar_p[r] + s * auxiliar_q[r];
                        eigenvectores_trabajo.datos[r][q]= -s * auxiliar_p[r] + c * auxiliar_q[r];
                    }
                    contador_iteraciones++;
                }
            }
        }
    } while (encontrar_max(matriz_trabajo).valor > TOLERANCIA_JACOBI);
    final = clock();

    //Ordenacion de autovectores y autovalores
    elemento_ordenacion* array_ord_eigenvalores;
    array_ord_eigenvalores = malloc(N* sizeof(elemento_ordenacion));
    for (q = 0; q<N; q++){
        array_ord_eigenvalores[q].valor = matriz_trabajo.datos[q][q];
        array_ord_eigenvalores[q].indice_original = q;
    }

    qsort(array_ord_eigenvalores,N, sizeof(elemento_ordenacion),&algoritmo_ordenacion);


   int columna;
   for (q=0;q<N;q++){
       columna = array_ord_eigenvalores[q].indice_original;
       array_eigenvalores[q] = array_ord_eigenvalores[q].valor;
       for (p=0;p<N;p++){
           matriz_eigenvectores.datos[p][q]=eigenvectores_trabajo.datos[p][columna];
       }
   }
 

    matriz_liberar(matriz_trabajo);
    matriz_liberar(eigenvectores_trabajo);
    free (array_ord_eigenvalores);
    free(auxiliar_p);
    free(auxiliar_q);
    auxiliar_p = NULL;
    auxiliar_q = NULL;
    array_ord_eigenvalores = NULL;



    if (diagnostico == true) {
        printf(ANSI_COLOR_GREEN"eigen_jacobi: Diagonalizacion completada correctamente en %d barridos con tolerancia %.1e.\n",
               contador_iteraciones, TOLERANCIA_JACOBI);
        tiempo_cpu_diag=((double)(final-inicio))/CLOCKS_PER_SEC;
        printf("Tiempo de CPU requerido:\t%lf s\n"ANSI_COLOR_RESET,tiempo_cpu_diag);
        printf("***EIGENVALORES***\n");
        for(p=0;p<N;p++){
            printf("%lf\t",array_eigenvalores[p]);
        }
        printf("\n***AUTOVECTORES**\n");
        matriz_imprimir(matriz_eigenvectores);
    }
    return 0;
}

    int algoritmo_ordenacion (const void * a, const void * b){

        elemento_ordenacion* elemA = (elemento_ordenacion*)(a);
        elemento_ordenacion* elemB = (elemento_ordenacion*)(b);
        if(elemA->valor < elemB->valor)
            return -1;
        else if (elemA->valor > elemB->valor)
            return 1;
        else
            return 0;
}


elemento_matriz encontrar_max(matrix matriz) {
    int i, j, N;
    N = matriz.m;
    elemento_matriz elemento_max_off = {1, 0, 0};
    double max = 0;
    for (i = 1; i < N; i++) {
        for (j = 0; j < i; j++) {
            if (fabsl(matriz.datos[i][j]) > max) {
                elemento_max_off.i = i;
                elemento_max_off.j = j;
                max = fabsl(matriz.datos[i][j]);
            }
        }
    }
    elemento_max_off.valor = max;
    return elemento_max_off;
}



