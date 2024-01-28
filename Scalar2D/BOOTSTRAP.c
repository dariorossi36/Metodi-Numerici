#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include"/home/dario/Documents/UNI/Metodi/ran2.h"
#include"/home/dario/Documents/UNI/Metodi/usefulstuff.h"
#include"/home/dario/Documents/UNI/Metodi/listfunction.h"


/*NOTA BENE: i file _NO_RELAX sono identici a quelli normali eccetto per il fatto che sono stati
             generati solo da heathbath senza over-relaxation-----> commentare o scommentare tale 
             scritta a seconda del tipo di simulazione che si vuole fare, per non sovrascrivere 
             dati importanti
             */ 


int N_misure; //N = numero di misure----> lo leggo da file 

void bootstrap(FILE *f, FILE *BOOT, float m, int Nt, int Nx, float Ntm, float Nxm){
    unsigned long s; //serve per salvarsi il numero casuale enorme che produco
    int i=0,j,k,M=100,x,bin=1000,taglio=0;         /*M = numero di resampling, 
    x serve per scanfare i file, i e j indici ricorsivi, bin è il bin del bootstrap, 
    taglio è per levare le misure non termalizzate*/
    int V=N_misure-taglio;
    long int seed = 13;         // seed per il ran2()
    float e_mass, e_spat, e_temp, e_density; /*cinetica*/         // risultati del simulatore_boot.c e valore del beta corrente
    int iter;           // risultati del simulatore_boot.c
    float *dati_mass = malloc(sizeof(long double) * (N_misure-taglio));
    float *dati_dens = malloc(sizeof(long double) * (N_misure-taglio));
    float std_mass=0, std_dens=0, medie_mass[M], medie_dens[M];
    float media_definitiva_mass=0, media_definitiva_dens=0;

    x=fscanf(f, "%d %f  %f  %f  %f", &iter, &e_mass, &e_spat, &e_temp, &e_density);

    while (x!=EOF){
            dati_mass[i]=e_mass;
            dati_dens[i]=e_density;
            x=fscanf(f, "%d %f  %f  %f  %f", &iter, &e_mass, &e_spat, &e_temp, &e_density);
            i++;
        }


    for(j=0;j<M;j++){
        medie_mass[j]=0;
        medie_dens[j]=0;
        for (k=0;k<(N_misure-taglio)/bin;k++){

            s=ran2(&seed)*(N_misure-taglio);
            s=(int)s;

            s=s+taglio;
            //printf("%ld \n",s);
            for(int h=0;h<bin;h++){
              //  printf("terzo for \n");
                if(s+h>N_misure-1){
                    medie_mass[j]=medie_mass[j]+dati_mass[s+h+1-(N_misure-taglio)];
                    medie_dens[j]=medie_dens[j]+dati_dens[s+h+1-(N_misure-taglio)];
                }
                if(s+h<=N_misure-1){
                     medie_mass[j]=medie_mass[j]+dati_mass[s+h];
                     medie_dens[j]=medie_dens[j]+dati_dens[s+h];
                }

            }
            
        }
            medie_mass[j]=medie_mass[j]/(N_misure-taglio);
            media_definitiva_mass=media_definitiva_mass+medie_mass[j];
            medie_dens[j]=medie_dens[j]/(N_misure-taglio);
            media_definitiva_dens=media_definitiva_dens+medie_dens[j];
        }
        media_definitiva_dens=media_definitiva_dens/M;
        media_definitiva_mass=media_definitiva_mass/M;
    
        for(j=0;j<M;j++){
            std_dens=std_dens+pow((medie_dens[j]-media_definitiva_dens),2)/(M-1);
            std_mass=std_mass+pow((medie_mass[j]-media_definitiva_mass),2)/(M-1);
        }
        std_dens=sqrt(std_dens);
        std_mass=sqrt(std_mass);
        fprintf(BOOT,"%.9f %.9f %.9f %.9f %f %d %d %f %f \n", media_definitiva_dens, std_dens, media_definitiva_mass, std_mass, m, Nt, Nx, Ntm, Nxm);


    free(dati_mass);
    free(dati_dens);

    return;
}


/*==== NOTA: Nel file di input.txt classico non viene specificata la lunghezza del reticolo. Quindi è necessario
crearne uno che la contenga (io l'ho chiamato input2.txt)*/ 
void read_file_input(){ 

    FILE *in;
    in = fopen("/home/dario/Documents/UNI/Metodi/Modulo3/Scalar2D/input.txt", "r");
    control_file(in);
    int iflag, idec, extfield; 
    float enne, d_metro;
    int x;
    x=fscanf(in, "%d  %d  %d  %d", &iflag, &N_misure, &idec, &extfield);
    fclose(in);
    return ;
}


int main(){

    read_file_input();  // funzione per leggere file con variabili di input (numero di misure, lunghezza del lato del reticolo, iflag per la partenza e campo esterno e idecorrel)

    char path[500];

    FILE *f,*file_bootstrap;         // puntatore al file degli eta
    int p=0, L;         // p serve per scorrere il file degli eta e L è la lunghezza dell'array degli eta. 
    float Ntm;
    float Nxm;
    float singolom;
    float m;
    int Nx;
    char dir[100];


    
        printf("######################################################## \n");
        printf("    Per bootstrappare il limite al continuo (nt*m, ns*m fisso) digiti 1     \n");
        printf("######################################################## \n");
        printf(" Per bootstrappare la variazione di temperatura (Nt variabile, m fisso) digiti 0 \n");
        printf("######################################################## \n");
        printf("          Per bootstrappare entrambe digiti 2            \n");
        printf("######################################################## \n");

        int scelta;
        scanf("%d", &scelta);

     if(scelta==1 || scelta==2){
         printf("############################################################# \n");
         printf("Ha scelto di bootstrappare il limite al continuo, buona giornata\n");
         printf("############################################################# \n");

         p=0;

        f=fopen("/home/dario/Documents/UNI/Metodi/Modulo3/Scalar2D/valori_cont.txt","r");
        control_file(f);
        list *Nt=NULL;           // lista contenente i valori degli Nt, letti da file

        Nt = scan_file(f,Nt);         // lettura degli valori degli Nt, con funzione ("listfunction.h")
        L=count(Nt)-3;         // funzione che conta la lunghezza di una lista ("listfunction.h")
        Ntm=val_posizione(0,Nt);
        Nxm=val_posizione(1,Nt);
        singolom=val_posizione(2,Nt);
        Nt=move_to_position(3,Nt);


        char dir[50]="Continuo";

        char filebootstrap[900];
        sprintf(filebootstrap,"/home/dario/Documents/UNI/Metodi/Modulo3/Scalar2D/%s/Bootstrap/bootstrap_NO_RELAX_new.txt", dir);
        file_bootstrap=fopen(filebootstrap,"w");
        control_file(file_bootstrap);
        fprintf(file_bootstrap, "# media_definitiva_dens / std_dens / media_definitiva_mass / std_mass / m / Nt / Nx / Ntm / Nxm \n");

        FILE *misure[L];          // puntatore per i file che salvano gli obs, per ciascun eta
        while(Nt!=NULL){
            char filemisure[500];
            m=Ntm/(Nt->val);
            Nx=(int)(Nxm/m);

            sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo3/Scalar2D/Continuo/misure_NO_RELAX(Ntm=%.3f,Nxm=%.3f,Nt=%.0f,Nx=%d,m=%.4f).txt", Ntm, Nxm, Nt->val, Nx, m); //it modifies each time the name of the file to etae created
            misure[p]=fopen(filemisure, "r");
            control_file(misure[p]);
            

            bootstrap(misure[p], file_bootstrap, m, (int)Nt->val, Nx, Ntm, Nxm); 
             /*eseguo la simulazione per il valore di eta in causa. 
            Funzione di oscillazione.h Mi creerà automaticamente un file di misure contenente i valori delle osservabili che
            ci interessano (y^2medio, Dy^2medio, ymedio) e eta associato per ogni misura fatta. */
     
             printf("######################################################## \n");
             printf("              Ho terminato il valore %f                  \n", Nt->val);
             printf("######################################################## \n");      

            fclose(misure[p]);
            Nt=Nt->next;
            p++;
        }
        fclose(f);
        fclose(file_bootstrap);
      }




        if(scelta==0 || scelta==2) {
            printf("######################################################## \n");
            printf("Ha scelto di bootstrappare la temperatura, buona giornata\n");
            printf("######################################################## \n");


            p=0;

            f=fopen("/home/dario/Documents/UNI/Metodi/Modulo3/Scalar2D/valori_temp.txt","r");
            control_file(f);

            list *Nt=NULL;           // lista contenente i valori degli Nt, letti da file

            Nt = scan_file(f,Nt);         // lettura degli valori degli Nt, con funzione ("listfunction.h")
            L=count(Nt)-3;         // funzione che conta la lunghezza di una lista ("listfunction.h")
            Ntm=val_posizione(0,Nt);
            Nxm=val_posizione(1,Nt);
            singolom=val_posizione(2,Nt);
            Nt=move_to_position(3,Nt);


            char dir[200]="Temperatura/Simmetrico";

            char filebootstrap[900];
            sprintf(filebootstrap,"/home/dario/Documents/UNI/Metodi/Modulo3/Scalar2D/%s/Bootstrap/bootstrap.txt", dir);
            file_bootstrap=fopen(filebootstrap,"w");
            control_file(file_bootstrap);
            fprintf(file_bootstrap, "# media_definitiva_dens / std_dens / media_definitiva_mass / std_mass / m / Nt / Nx / Ntm / Nxm \n");

            FILE *misure[L];          // puntatore per i file che salvano gli obs, per ciascun eta
            while(Nt!=NULL){
                Ntm=Nt->val*singolom;
                Nx=(int)(Nxm/singolom);
                char filemisure[500];
                sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo3/Scalar2D/Temperatura/Simmetrico/misure(Ntm=%.3f,Nxm=%.3f,Nt=%.0f,Nx=%d,m=%.3f).txt", Ntm, Nxm,Nt->val, Nx, singolom); //it modifies each time the name of the file to etae created
                misure[p]=fopen(filemisure, "r");
                control_file(misure[p]);

                bootstrap(misure[p], file_bootstrap, singolom, (int)Nt->val, Nx, Ntm, Nxm);
                /*eseguo la simulazione per il valore di eta in causa. 
                Funzione di oscillazione.h Mi creerà automaticamente un file di misure contenente i valori delle osservabili che
                ci interessano (y^2medio, Dy^2medio, ymedio) e eta associato per ogni misura fatta. */
                printf("######################################################## \n");
                printf("              Ho terminato il valore %f                  \n", Nt->val);
                printf("######################################################## \n");

                fclose(misure[p]);
                Nt=Nt->next;
                p++;
            }
            fclose(f);
            fclose(file_bootstrap);
        }


    printf("############################################################# \n");
    printf("Ci auguriamo sia stata una piacevole esperienza, arrivederci!\n");
    printf("############################################################# \n");
    return 0;
}
