/*
 * Instituto Tecnologico de Costa Rica
 * Escuela de Ingenieria en Computacion
 * Computer Graphics
 *
 * Segundo proyecto de CG
 * Ronald Esquivel Lopez
 * Ricardo Murillo Jimenez
 *
 * Programa: Ray Tracing Basico
 * Archivo:  rayTracingBasico.c
 */
#include "rayTracingBasico.h"
#include "malloc.h"

COLOR **buffer;
FILE *COORD;
LUZ *luces[LUZ_C];
ESFERA *esferas[OBJ_C];

#define INFINITE 999999999
#define X_MAX H_SIZE
#define X_MIN 0.0
#define Y_MAX V_SIZE
#define Y_MIN 0.0
//#define XW(x) ((x + 1/2)*(X_MAX - X_MIN) / H_SIZE + X_MIN)
//#define YW(x) ((x + 1/2)*(Y_MAX - Y_MIN) / V_SIZE + Y_MIN)

#define XE 504.0
#define YE 283.0
#define ZE -1000.0   //Coordenadas de los ojos

#define IA 0.5



#define BCKGR_R 0
#define BCKGR_G 0 //Colores del background [negro]
#define BCKGR_B 0

int cantidad_objetos = 0;


int main(int argc, char** argv){ //main del programa
    COORD=fopen ("objetosEscena.txt","r");
    if ( COORD == NULL )
    {
        //printf("No se puede abrir archivo") ;
    }
    crear_buffer();
    crear_objetos_escena();
    fclose(COORD);
    ray_tracing();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(H_SIZE,V_SIZE);
    glutCreateWindow("Draw CostaRica");
    glClear(GL_COLOR_BUFFER_BIT);
    gluOrtho2D(-0.5, H_SIZE +0.5, -0.5, V_SIZE + 0.5);
    glutDisplayFunc(dibujar_escena);
    //glutKeyboardFunc(teclado);
    glutMainLoop();
  	return 1;
}

void dibujar_escena(){
    int i , j;
    for (i = 0; i < H_SIZE; i++){
        for (j = 0; j < V_SIZE; j++){
            glBegin (GL_POINTS);
            double n_r = buffer[i][j].r * 1.0 / 255.0;
            double n_g = buffer[i][j].g * 1.0 / 255.0; 
            double n_b = buffer[i][j].b * 1.0 / 255.0;  
            ////printf("\n COLORES A PINTAR = R = %f, G = %f, B = %f",buffer[i][j].r,buffer[i][j].g,buffer[i][j].b);
            glColor3f (n_r,n_g,n_b);
            glVertex2i (i,j);
            glEnd();      
        }
  }
  glFlush();
}


void crear_buffer(){ //funcion encargada de liberar memoria y de crear el framebuffer (en ese orden)
    int i, j;
    if(buffer != NULL){
        for (i = 0; i < H_SIZE; i++){
            free(buffer[i]);
        }
         free(buffer);
    }
  	buffer = (COLOR **)malloc(H_SIZE * sizeof(COLOR*));
  	for (i = 0; i < H_SIZE; i++){
    	buffer[i] = (COLOR *)malloc(V_SIZE * sizeof(COLOR));
    }
  	for (i = 0; i < H_SIZE; i++){
    	for (j = 0; j < V_SIZE; j++){
            buffer[i][j].r = 255;
            buffer[i][j].g = 255;
            buffer[i][j].b = 255;
        }
   }
}

void crear_objetos_escena(){
    rewind(COORD);
    int i;
    for(i =0; i < LUZ_C; i++){ //almacena las luces en un sruct
        COLOR* color_luz = (COLOR *)malloc(sizeof(COLOR));
        VECTOR* pos_luz = (VECTOR *)malloc(sizeof(VECTOR));
        LUZ* nueva_luz = (LUZ *)malloc(sizeof(LUZ));
        double x,y,z,r,g,b,intensidad;
        fscanf(COORD, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &intensidad, &x, &y, &z, &r, &g, &b);
        //printf("\n LUZ Antes X = %f, Y = %f, Z = %f, r = %f, g = %f, b = %f, intensidad = %f",x,y,z,r,g,b,intensidad);
        nueva_luz->intensidad = intensidad;
        color_luz->r = r; color_luz->g = g; color_luz->b = b;
        pos_luz->X = x; pos_luz->Y = y; pos_luz->Z = z;
        nueva_luz->color_objeto = color_luz;
        nueva_luz->posicion = pos_luz;
        luces[i] = nueva_luz;
        ////printf("\n LUZ DESPUES X = %f, Y = %f, Z = %f, r = %f, g = %f, b = %f, intensidad = %f",nueva_luz->posicion->X,nueva_luz->posicion->Y,nueva_luz->posicion->Z,nueva_luz->color_objeto->r,nueva_luz->color_objeto->g,nueva_luz->color_objeto->b,nueva_luz->intensidad);
    }
    for(i = 0; i < OBJ_C; i++){ //almacena las esferas en un struct
        COLOR* color_objeto = (COLOR *)malloc(sizeof(COLOR));
        VECTOR* centro_objeto = (VECTOR *)malloc(sizeof(VECTOR));
        ESFERA* nueva_esfera = (ESFERA *)malloc(sizeof(ESFERA));
        double x,y,z,r,g,b,radio,kd,ka;
        fscanf(COORD, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &radio, &x, &y, &z, &r, &g, &b,&kd,&ka);
        //printf("\n ESFERA Antes X = %f, Y = %f, Z = %f, Radio = %f, r = %f, g = %f, b = %f, kd = %f, ka = %f",x,y,z,radio,r,g,b,kd,ka);
        nueva_esfera->radio = radio;
        centro_objeto->X = x; centro_objeto->Y = y; centro_objeto->Z = z;
        color_objeto->r = r; color_objeto->g = g; color_objeto->b = b;
        nueva_esfera->centro = centro_objeto;
        nueva_esfera->color_objeto = color_objeto;
        nueva_esfera->KD = kd;
        nueva_esfera->KA = ka;
        ////printf("\nDespues X = %f, Y = %f, Z = %f, Radio = %f, r = %f, g = %f, b = %f, kd = %f, ka = %f",nueva_esfera->centro->X,nueva_esfera->centro->Y,nueva_esfera->centro->Z,nueva_esfera->radio, nueva_esfera->color_objeto->r,nueva_esfera->color_objeto->g, nueva_esfera->color_objeto->b,nueva_esfera->KD,nueva_esfera->KA);
        esferas[i] = nueva_esfera;
    }
}

void ray_tracing(){
    int x,y;
    for(x = 0; x < H_SIZE; x++){
        for(y = 0; y < V_SIZE; y++){
            double dx = (double) x; double dy = (double) y;
            double xw = (double)((dx + 0.5)*(X_MAX - X_MIN) / H_SIZE + X_MIN);
            double yw = (double)((dy + 0.5)*(Y_MAX - Y_MIN) / V_SIZE + Y_MIN);
            double zw = 0.0; 
            double L = sqrt(pow((xw-XE),2)+pow((yw-YE),2)+pow((zw-ZE),2));
            double xd = (xw - XE) / L;
            double yd = (yw - YE) / L;
            double zd = (zw - ZE) / L;
            ////printf("\n x = %d, y = %d, xw = %f, yw = %f, zw = %f, L = %f, xd = %f, yd = %f, zd = %f",x,y,xw,yw,zw,L,xd,yd,zd);
            COLOR *color_resultado = averiguar_color(xd,yd,zd);
            buffer[x][y].r = color_resultado->r; buffer[x][y].g = color_resultado->g; buffer[x][y].b = color_resultado->b;
        }
    }
}

COLOR* averiguar_color(double xd,double yd, double zd){ //averigua el color de la interseccion del rayo, si es que hay. Si no hay, devuelve los colores definidos como background
    COLOR* color_buscado = (COLOR *)malloc(sizeof(COLOR));
    INTERSECCION * interseccion_rayo;
    interseccion_rayo = primera_interseccion(xd,yd,zd);
    if(interseccion_rayo == NULL){
        color_buscado->r = BCKGR_R;
        color_buscado->g = BCKGR_G;
        color_buscado->b = BCKGR_B;
    }else{
        double distancia_ojo_min = interseccion_rayo-> distancia_ojo_minima;
        double XC = interseccion_rayo->objeto_tocado->X; double YC = interseccion_rayo->objeto_tocado->Y; double ZC = interseccion_rayo->objeto_tocado->Z; double RADIO = interseccion_rayo->radio_objeto;
        double XI = XE + (distancia_ojo_min * xd); double YI = YE + (distancia_ojo_min * yd); double ZI = ZE + (distancia_ojo_min * zd);
        double KA = interseccion_rayo->KA_objeto; double KD = interseccion_rayo->KD_objeto;
        double intensidad_total = 0;
        int i;
        for (i = 0; i < LUZ_C; i++){
            //printf("\n\n------------------------------------------------");
            double xp = luces[i]->posicion->X; double yp = luces[i]->posicion->Y; double zp = luces[i]->posicion->Z;
            //printf("\nxp = %f,yp = %f,zp = %f",xp,yp,zp);
            double vx = xp - XI ;double vy = yp - YI;double vz = zp - ZI;
            //printf(" vx = %f,vy = %f,vz = %f",vx,vy,vz);
            double distancia = sqrt(pow(vx,2)+pow(vy,2)+pow(vz,2));
            //printf(" Distancia = %f",distancia);
            double LX = vx / distancia;  double LY = vy / distancia;  double LZ = vz / distancia;
            //printf(" LX = %f,LY = %f,LZ = %f",LX,LY,LZ);
            double NX = (XI - XC) / RADIO; double NY = (YI - YC) / RADIO; double NZ = (ZI - ZC) / RADIO;
            //printf(" NX = %f,NY = %f,NZ = %f",NX,NY,NZ);
            double producto_escalar = (LX*NX)+(LY*NY)+(LZ*NZ);
            //printf(" Producto Escalar = %f",producto_escalar);
            double modulo_L = sqrt(pow(LX,2) + pow(LY,2) + pow(LZ,2));
            double modulo_N = sqrt(pow(NX,2) + pow(NY,2) + pow(NZ,2));
            //printf(" Modulo L = %f , Modulo N = %f",modulo_L,modulo_N);
            double angulo_vectores = producto_escalar / (modulo_L*modulo_N);
            //printf(" Angulo Vectores = %f",angulo_vectores);
            double intensidad;
            if(angulo_vectores < 0){
                intensidad = 0;
            }else{
                intensidad = angulo_vectores;
            }
            //printf("\nIntensidad Antes = %f",intensidad);
            intensidad = (intensidad * KD * luces[i]->intensidad /** FATT*/);
            //printf(" Intensidad Despues = %f",intensidad);
            intensidad_total += intensidad;
            //printf("\n------------------------------------------------\n");
        }
        ////printf("\nIntensidad Total Antes = %f",intensidad_total);
        intensidad_total+= (IA * KA);
        //printf("\nIntensidad Total = %f",intensidad_total);
        if(intensidad_total > 1){
            intensidad_total = 1;
        }
        double r_editado = (interseccion_rayo->color_objeto->r * intensidad_total); double g_editado = (interseccion_rayo->color_objeto->g * intensidad_total); double b_editado = (interseccion_rayo->color_objeto->b * intensidad_total);
        //printf("\nR = %f,G = %f,B = %f",interseccion_rayo->color_objeto->r,interseccion_rayo->color_objeto->g,interseccion_rayo->color_objeto->b);
        //printf("\nR NUEVO = %f,G NUEVO = %f,B NUEVO = %f",r_editado,g_editado,b_editado);
        color_buscado->r = r_editado; color_buscado->g = g_editado; color_buscado->b = b_editado;
        free(interseccion_rayo);
    }
    return color_buscado;
}

INTERSECCION * primera_interseccion(double xd,double yd, double zd){ //busca la primera interseccion que se encuentre el rayo 
    INTERSECCION *nueva_interseccion = NULL;
    COLOR* color_buscado = (COLOR *)malloc(sizeof(COLOR));
    VECTOR* vector_esfera = (VECTOR *)malloc(sizeof(VECTOR));
    int r,g,b;
    double XC,YC,ZC,RADIO,KD,KA; //posiciones del centro y radio de la esfera con la minima distancia con respecto al ojo
    long double distancia_ojo_min = INFINITE; 
    int i,encontrado;
    encontrado = 0;
    for(i = 0; i < OBJ_C; i++){
        double distancia_ojo_1 = INFINITE;
        double distancia_ojo_2 = INFINITE;
        double xc = esferas[i]->centro->X;
        double yc = esferas[i]->centro->Y;
        double zc = esferas[i]->centro->Z;
        double R = esferas[i]->radio;
        //double alpha = pow(xd,2) + pow(yd,2) + pow(zd,2);
        double beta = 2*(xd * (XE - xc) + yd * (YE - yc) + zd * (ZE - zc));
        double Y = pow((XE - xc), 2) + pow((YE - yc),2) + pow((ZE - zc), 2) - pow(R, 2);
        double discriminante = pow(beta, 2) - (4.0 * Y); //se simplifica porque alpha es unitario
        if(discriminante > 0){
            if(discriminante > 0.00001){
                distancia_ojo_1 = ((-beta + sqrt(pow(beta,2) - 4 * Y))/2);
                distancia_ojo_2 = ((-beta - sqrt(pow(beta,2) - 4 * Y))/2);
                ////printf("\n xc = %f, yc = %f, zc = %f, x = %f, yd = %f, zd = %f, radio = %f, beta = %f, Y = %f, Discriminante = %lf, t1 = %f, t2 = %f",xc,yc,zc,xd,yd,zd,R,beta,Y,discriminante, distancia_ojo_1, distancia_ojo_2);
            }else{
                distancia_ojo_1 = ((-beta + sqrt(pow(beta,2) - 4 * Y))/2);
                ////printf("\n xc = %f, yc = %f, zc = %f, x = %f, yd = %f, zd = %f, radio = %f, beta = %f, Y = %f, Discriminante = %f, t1 = %f",xc,yc,zc,xd,yd,zd,R,beta,Y,discriminante, distancia_ojo_1);
            }
        }
        if(distancia_ojo_1 < distancia_ojo_min && distancia_ojo_1 >= 0){
            ////printf("\nSE HA ENCONTRADO INTERSECCION");
            encontrado = 1;
            distancia_ojo_min = distancia_ojo_1;
            XC = xc; YC = yc; ZC = zc; RADIO = R;
            r = esferas[i]->color_objeto->r; g = esferas[i]->color_objeto->g; b = esferas[i]->color_objeto->b;
            KD = esferas[i]->KD; KA = esferas[i]->KA;
            ////printf("\nMinimo = %Lf", distancia_ojo_min);
            //color_buscado = esferas[i]->color_objeto;
        }if(distancia_ojo_2 < distancia_ojo_min && distancia_ojo_2 >= 0){
            ////printf("\nSE HA ENCONTRADO INTERSECCION");
            encontrado = 1;
            distancia_ojo_min = distancia_ojo_2;
            XC = xc; YC = yc; ZC = zc; RADIO = R;
            r = esferas[i]->color_objeto->r; g = esferas[i]->color_objeto->g; b = esferas[i]->color_objeto->b;
            KD = esferas[i]->KD; KA = esferas[i]->KA;
            ////printf("\nMinimo = %Lf", distancia_ojo_min);
            //color_buscado = esferas[i]->color_objeto;
        }
    }
    if(encontrado){
        ////printf("\nMinimo FINAL= %Lf", distancia_ojo_min);
        nueva_interseccion = (INTERSECCION *)malloc(sizeof(INTERSECCION));
        vector_esfera->X = XC; vector_esfera->Y = YC; vector_esfera->Z = ZC;
        color_buscado->r = r; color_buscado->g = g; color_buscado->b = b;
        nueva_interseccion->distancia_ojo_minima = distancia_ojo_min;
        nueva_interseccion->color_objeto = color_buscado;
        nueva_interseccion->objeto_tocado = vector_esfera;
        nueva_interseccion->radio_objeto = RADIO;
        nueva_interseccion->KD_objeto = KD;
        nueva_interseccion->KA_objeto = KA; 
    }else{
        free(color_buscado);
        free(vector_esfera);
    }
    return nueva_interseccion;

}

void plot(int x, int y, int r, int g, int b){ //funcion encargada de colocar en x y coordenadas valores de RGB dentro del framebuffer

	////printf("Resultado = %d , %d\n", x, y);
    if(x >= 0 && x < V_SIZE && y >= 0 && y < H_SIZE){
        buffer[x][y].r = r;
        buffer[x][y].g = g;
        buffer[x][y].b = b;
    }
}


