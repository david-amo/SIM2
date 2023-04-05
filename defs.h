//
// Created by David on 16/09/2021.
//

#ifndef MINIFOCK_DEFS_H
#define MINIFOCK_DEFS_H

#define SIGNO(x)    ((x) >= 0.0 ? 1 : -1)

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


#define TOLERANCIA_SCF 1E-09
#define TOLERANCIA_GRAL 1E-14
// Errores genericos
#define ERR_MEM "ERROR DE ASIGNACION DINAMICA DE MEMORIA\nEL PROGRAMA SE CERRARA\n"







#endif //MINIFOCK_DEFS_H
