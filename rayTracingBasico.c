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
    glutCreateWindow("Proyecto#2 Ronald Esquivel y Ricardo Murillo");
    glClear(GL_COLOR_BUFFER_BIT);
    gluOrtho2D(-0.5, H_SIZE +0.5, -0.5, V_SIZE + 0.5);
    glutDisplayFunc(dibujar_escena);
    //glutKeyboardFunc(teclado);
    glutMainLoop();
  	return 1;
}

void dibujar_escena(){ //dibuja el framebuffer, mediante el uso de glut
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
  guardar_imagen();
}

void guardar_imagen(){
  int i, j;
  FILE *fp = fopen("RayTrace-Esquivel-López-Murillo-Jiménez.ppm", "wb"); //funcion encargada de pasar el framebuffer a formato Netpbm en su tipo Portable Pixmap. Las imagenes output son Raster
  (void) fprintf(fp, "P6\n%d %d\n255\n", H_SIZE, V_SIZE);
  for (j = 0; j < V_SIZE; ++j)
  {
    for (i = 0; i < H_SIZE; ++i)
    {
      static unsigned char color[3];
      color[0] = buffer[i][j].r; 
      color[1] = buffer[i][j].g;  
      color[2] = buffer[i][j].b;  
      (void) fwrite(color, 1, 3, fp);
    }
  }
  (void) fclose(fp);
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

void crear_objetos_escena(){//funcion que crea los objetos, leyendo los parametros del archivo objetosEscena.txt y los almacena en structs
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
        double x,y,z,r,g,b,radio,kd,ka,ks;
        fscanf(COORD, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &radio, &x, &y, &z, &r, &g, &b,&kd,&ka,&ks);
        //printf("\n ESFERA Antes X = %f, Y = %f, Z = %f, Radio = %f, r = %f, g = %f, b = %f, kd = %f, ka = %f, ks = %f",x,y,z,radio,r,g,b,kd,ka,ks);
        nueva_esfera->radio = radio;
        centro_objeto->X = x; centro_objeto->Y = y; centro_objeto->Z = z;
        color_objeto->r = r; color_objeto->g = g; color_objeto->b = b;
        nueva_esfera->centro = centro_objeto;
        nueva_esfera->color_objeto = color_objeto;
        nueva_esfera->KD = kd;
        nueva_esfera->KA = ka;
        nueva_esfera->KS = ks;
        //printf("\nDespues X = %f, Y = %f, Z = %f, Radio = %f, r = %f, g = %f, b = %f, kd = %f, ka = %f, ks = %f",nueva_esfera->centro->X,nueva_esfera->centro->Y,nueva_esfera->centro->Z,nueva_esfera->radio, nueva_esfera->color_objeto->r,nueva_esfera->color_objeto->g, nueva_esfera->color_objeto->b,nueva_esfera->KD,nueva_esfera->KA,nueva_esfera->KS);
        esferas[i] = nueva_esfera;
    }
}

void ray_tracing(){//funcion main del ray tracing
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

double calcular_angulo(double AX, double AY, double AZ, double BX, double BY, double BZ){ //dado dos vectores, calcula el angulo de estos
        double producto_escalar = (AX*BX)+(AY*BY)+(AZ*BZ);
        //printf(" Producto Escalar = %f",producto_escalar);
        double modulo_A = sqrt(pow(AX,2) + pow(AY,2) + pow(AZ,2));
        double modulo_B = sqrt(pow(BX,2) + pow(BY,2) + pow(BZ,2));
        //printf(" Modulo L = %f , Modulo N = %f",modulo_L,modulo_N);
        double angulo_vectores = producto_escalar / (modulo_A*modulo_B);
        return angulo_vectores;
}

COLOR* averiguar_color(double xd,double yd, double zd){ //averigua el color de la interseccion del rayo, si es que hay. Si no hay, devuelve los colores definidos como background
    COLOR* color_buscado = (COLOR *)malloc(sizeof(COLOR));
    INTERSECCION * interseccion_rayo;
    interseccion_rayo = primera_interseccion(xd,yd,zd,XE,YE,ZE);
    if(interseccion_rayo == NULL){
        color_buscado->r = BCKGR_R;
        color_buscado->g = BCKGR_G;
        color_buscado->b = BCKGR_B;
    }else{
        double distancia_ojo_min = interseccion_rayo-> distancia_ojo_minima;
        double XC = interseccion_rayo->objeto_tocado->X; double YC = interseccion_rayo->objeto_tocado->Y; double ZC = interseccion_rayo->objeto_tocado->Z; double RADIO = interseccion_rayo->radio_objeto;
        double XI = XE + (distancia_ojo_min * xd); double YI = YE + (distancia_ojo_min * yd); double ZI = ZE + (distancia_ojo_min * zd);
        double KA = interseccion_rayo->KA_objeto; double KD = interseccion_rayo->KD_objeto; double KS = interseccion_rayo->KS_objeto;
        double intensidad_total = 0;
        double reflexion_especular_total = 0;
        double VX = XE - XI; double VY = YE - YI; double VZ = ZE - ZI;
        int i;
        for (i = 0; i < LUZ_C; i++){
            double xp = luces[i]->posicion->X; double yp = luces[i]->posicion->Y; double zp = luces[i]->posicion->Z;
            double vx = xp - XI ;double vy = yp - YI;double vz = zp - ZI;
            //printf(" vx = %f,vy = %f,vz = %f",vx,vy,vz);
            double distancia = sqrt(pow(vx,2)+pow(vy,2)+pow(vz,2));
            //printf(" Distancia = %f",distancia);
            double LX = vx / distancia;  double LY = vy / distancia;  double LZ = vz / distancia;
            //printf(" LX = %f,LY = %f,LZ = %f",LX,LY,LZ);
            double NX = (XI - XC) / RADIO; double NY = (YI - YC) / RADIO; double NZ = (ZI - ZC) / RADIO;
            //printf(" NX = %f,NY = %f,NZ = %f",NX,NY,NZ);
            double angulo_vectores = calcular_angulo(LX,LY,LZ,NX,NY,NZ);
            double FATT = (1/pow(distancia,2));
            INTERSECCION * interseccion_obstaculo = primera_interseccion(XI,YI,ZI,xp,yp,zp);
            int esSombra = 0;
            if(interseccion_obstaculo != NULL && angulo_vectores > 0){
                double XI_obs = xp + (interseccion_obstaculo->distancia_ojo_minima * XI); double YI_obs = yp + (interseccion_obstaculo->distancia_ojo_minima  * YI); double ZI_obs = zp + (interseccion_obstaculo->distancia_ojo_minima * ZI);
                double vx_obs = xp - XI_obs ;double vy_obs = yp - YI_obs;double vz_obs = zp - ZI_obs;
                double distancia_obs = sqrt(pow(vx_obs,2)+pow(vy_obs,2)+pow(vz_obs,2));
                //printf("\nCOLOR DE OBJETO INTERSECATO = R = %f, G = %f Y B = %f, distancia de objeto a luz = %f y distancia hacia el obstaculo = %f",interseccion_obstaculo->color_objeto->r,interseccion_obstaculo->color_objeto->g ,interseccion_obstaculo->color_objeto->b, distancia, distancia_obs);
                if(distancia_obs > distancia){
                    esSombra = 0;
                }else{
                    esSombra = 1;
                }
            }else if(interseccion_obstaculo == NULL){
                esSombra = 0;
            }
            if(esSombra == 0){
                double RX = (2 * NX * calcular_angulo(NX,NY,NZ,LX,LY,LZ)) - LX;
                double RY = (2 * NY * calcular_angulo(NX,NY,NZ,LX,LY,LZ)) - LY;
                double RZ = (2 * NZ * calcular_angulo(NX,NY,NZ,LX,LY,LZ)) - LZ;
                //printf("\nRX= %f, RY = %f, RZ = %f", RX,RY,RZ);
                double angulo_especular = calcular_angulo(RX,RY,RZ,VX,VY,VZ);
                if(angulo_especular > 0 && angulo_vectores > 0){
                    //printf("\nESPECULAR ANGULO = %f", calcular_angulo(NX,NY,NZ,LX,LY,LZ));
                    double reflexion_especular = (angulo_especular * KS * luces[i]->intensidad /** FATT*/);
                    //printf("\nreflexion_especular = %f, angulo_especular = %f, KS = %f, luces intensidad = %f", reflexion_especular,angulo_especular,KS,luces[i]->intensidad);
                    reflexion_especular_total += reflexion_especular;
                }
                double intensidad;
                if(angulo_vectores < 0){
                    intensidad = 0;
                }else{
                    intensidad = angulo_vectores;
                }
                //printf("\nIntensidad Antes = %f",intensidad);
                intensidad = (intensidad * KD * luces[i]->intensidad);
                intensidad_total += intensidad;
            }
            free(interseccion_obstaculo);
        }
        //printf("\nIntensidad Total Antes = %f",intensidad_total);
        intensidad_total += (IA * KA);
        //printf("\nIntensidad Total = %f",intensidad_total);
        if(intensidad_total > 1){
            intensidad_total = 1;
        }
        if(reflexion_especular_total > 1){
            reflexion_especular_total = 1;
        }
        //printf("\nESPECULAR TOTAL = %f", reflexion_especular_total);
        double r_difuso = (interseccion_rayo->color_objeto->r * intensidad_total); double g_difuso = (interseccion_rayo->color_objeto->g * intensidad_total); double b_difuso = (interseccion_rayo->color_objeto->b * intensidad_total);
        double r_especular = (r_difuso + (reflexion_especular_total*(255.0-r_difuso)));
        double g_especular = (g_difuso + (reflexion_especular_total*(255.0-g_difuso)));
        double b_especular = (b_difuso + (reflexion_especular_total*(255.0-b_difuso)));
        //printf("\nESPECULAR R = %f, G = %f, B = %f", r_especular, g_especular, b_especular);
        //printf("\nDIFUSO R = %f,G = %f,B = %f",r_difuso,g_difuso,b_difuso);
        //printf("\nESPECULAR R = %f,G = %f,B = %f",r_especular,g_especular,b_especular);
        color_buscado->r = r_especular; color_buscado->g = g_especular; color_buscado->b = b_especular;
        free(interseccion_rayo);
    }
    return color_buscado;
}

INTERSECCION * primera_interseccion(double xd,double yd, double zd, double posE_X, double posE_Y, double posE_Z){ //busca la primera interseccion que se encuentre el rayo 
    INTERSECCION *nueva_interseccion = NULL;
    COLOR* color_buscado = (COLOR *)malloc(sizeof(COLOR));
    VECTOR* vector_esfera = (VECTOR *)malloc(sizeof(VECTOR));
    int r,g,b;
    double XC,YC,ZC,RADIO,KD,KA,KS; //posiciones del centro y radio de la esfera con la minima distancia con respecto al ojo
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
        double beta = 2*(xd * (posE_X - xc) + yd * (posE_Y - yc) + zd * (posE_Z - zc));
        double Y = pow((posE_X - xc), 2) + pow((posE_Y - yc),2) + pow((posE_Z - zc), 2) - pow(R, 2);
        double discriminante = pow(beta, 2) - (4.0 * Y); //se simplifica porque alpha es unitario
        if(discriminante>0){
            if(discriminante > 0.0005){
                distancia_ojo_1 = ((-beta + sqrt(pow(beta,2) - 4 * Y))/2);
                distancia_ojo_2 = ((-beta - sqrt(pow(beta,2) - 4 * Y))/2);
                ////printf("\n xc = %f, yc = %f, zc = %f, x = %f, yd = %f, zd = %f, radio = %f, beta = %f, Y = %f, Discriminante = %lf, t1 = %f, t2 = %f",xc,yc,zc,xd,yd,zd,R,beta,Y,discriminante, distancia_ojo_1, distancia_ojo_2);
            }else{
            distancia_ojo_1 = ((-beta + sqrt(pow(beta,2) - 4 * Y))/2);
            ////printf("\n xc = %f, yc = %f, zc = %f, x = %f, yd = %f, zd = %f, radio = %f, beta = %f, Y = %f, Discriminante = %f, t1 = %f",xc,yc,zc,xd,yd,zd,R,beta,Y,discriminante, distancia_ojo_1);
            }
        }
        
        if(distancia_ojo_1 < distancia_ojo_min && distancia_ojo_1 > 0.0005){
            //printf("\nSE HA ENCONTRADO INTERSECCION");
            encontrado = 1;
            distancia_ojo_min = distancia_ojo_1;
            XC = xc; YC = yc; ZC = zc; RADIO = R;
            r = esferas[i]->color_objeto->r; g = esferas[i]->color_objeto->g; b = esferas[i]->color_objeto->b;
            KD = esferas[i]->KD; KA = esferas[i]->KA; KS = esferas[i]->KS;
            ////printf("\nMinimo = %Lf", distancia_ojo_min);
            //color_buscado = esferas[i]->color_objeto;
        }if(distancia_ojo_2 < distancia_ojo_min && distancia_ojo_2 > 0.0005){
            //printf("\nSE HA ENCONTRADO INTERSECCION");
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
        nueva_interseccion->KS_objeto = KS; 
    }else{
        free(color_buscado);
        free(vector_esfera);
    }
    return nueva_interseccion;

}


