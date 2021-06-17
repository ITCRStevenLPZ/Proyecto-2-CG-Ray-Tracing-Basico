/*
 * Instituto Tecnologico de Costa Rica
 * Escuela de Ingenieria en Computacion
 * Computer Graphics
 *
 * Primer proyecto de CG
 * Ronald ESquivel Lopez
 * Ricardo Murillo Jimenez
 *
 * Programa: Ray Tracing Basico
 * Archivo:  rayTracingBasico.h
 */

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <math.h>  


#define H_SIZE 1008
#define V_SIZE 567

#define OBJ_C 8
#define LUZ_C 3

typedef struct{
  double X;
  double Y;
  double Z;
}VECTOR;

typedef struct {
  double r;
  double g;
  double b;
} COLOR;

typedef struct {
  COLOR *color_objeto;
  VECTOR* centro;
  double radio;
  double KD; //coeficiente relexion difusa
  double KA; //coeficiente iluminacion ambiente
  double KS; //coeficiente relexion especular
} ESFERA;

typedef struct {
  COLOR *color_objeto;
  VECTOR *objeto_tocado;
  double radio_objeto;
  double distancia_ojo_minima;
  double KD_objeto;
  double KA_objeto;
  double KS_objeto;
} INTERSECCION;

typedef struct{
  COLOR *color_objeto;
  double intensidad;
  VECTOR* posicion;
}LUZ;


//-----------------------FUNCIONES--------------------------

void crear_buffer();
void crear_objetos_escena();
void ray_tracing();
COLOR* averiguar_color(double xd,double yd, double zd);
INTERSECCION * primera_interseccion(double xd,double yd, double zd, double posE_X, double posE_Y, double posE_Z); //busca la primera interseccion que se encuentre el rayo 
void dibujar_escena();
double calcular_angulo(double AX, double AY, double AZ, double BX, double BY, double BZ);
void guardar_imagen();