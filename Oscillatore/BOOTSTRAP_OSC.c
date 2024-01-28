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

void bootstrap(FILE *f, FILE *BOOT, float e, int enne){
    unsigned long s; //serve per salvarsi il numero casuale enorme che produco
    int i=0,j,k,M=100,x,bin=2000,taglio=0;         /*M = numero di resampling, 
    x serve per scanfare i file, i e j indici ricorsivi, bin è il bin del bootstrap, 
    taglio è per levare le misure non termalizzate*/
    int V=N_misure-taglio;
    long int seed = 13;         // seed per il ran2()
    float ics2, dics2, fild /*cinetica*/;         // risultati del simulatore_boot.c e valore del beta corrente
    int iter;           // risultati del simulatore_boot.c
    float *dati_x2 = malloc(sizeof(long double) * (N_misure));
    float *dati_dx2 = malloc(sizeof(long double) * (N_misure));
    float *dati_ene = malloc(sizeof(long double) * (N_misure));
    float *dati_y = malloc(sizeof(long double) * (N_misure));
    float std_x2=0, std_dx2=0, std_ene=0, std_y=0, medie_x2[M], medie_dx2[M], medie_ene[M], medie_y[M];
    float media_definitiva_x2=0, media_definitiva_dx2=0, media_definitiva_ene=0, media_definitiva_y=0;

    x=fscanf(f, "%f  %f  %f  %d", &ics2, &dics2, &fild, &iter);

    while (x!=EOF){
            dati_x2[i]=ics2;
            dati_dx2[i]=dics2;
            dati_y[i]=fild;
            dati_ene[i] = ics2/2 - dics2/(2*pow(e,2));
            x=fscanf(f, "%f  %f  %f  %d", &ics2, &dics2, &fild, /*&cinetica,*/ &iter);
            i++;
        }

    for(j=0;j<M;j++){
        medie_x2[j]=0;
        medie_dx2[j]=0;
        medie_ene[j]=0;
        medie_y[j]=0;
        for (k=0;k<(N_misure-taglio)/bin;k++){

            s=ran2(&seed)*(N_misure-taglio);
            s=(int)s;

            s=s+taglio;
            //printf("%ld \n",s);
            for(int h=0;h<bin;h++){
              //  printf("terzo for \n");
                if(s+h>N_misure-1){
                    medie_x2[j]=medie_x2[j]+dati_x2[s+h+1-(N_misure-taglio)];
                    medie_dx2[j]=medie_dx2[j]+dati_dx2[s+h+1-(N_misure-taglio)];
                    medie_ene[j]=medie_ene[j]+dati_ene[s+h+1-(N_misure-taglio)];
                    medie_y[j] = medie_y[j]+dati_y[s+h+1-(N_misure-taglio)];

                }
                if(s+h<=N_misure-1){
                    medie_x2[j]=medie_x2[j]+dati_x2[s+h];
                    medie_dx2[j]=medie_dx2[j]+dati_dx2[s+h];
                    medie_ene[j]=medie_ene[j]+dati_ene[s+h];
                    medie_y[j] = medie_y[j]+dati_y[s+h];
                }

            }
            
        }

            medie_x2[j]=medie_x2[j]/(N_misure-taglio);
            medie_dx2[j]=medie_dx2[j]/(N_misure-taglio);
            medie_ene[j]=medie_ene[j]/(N_misure-taglio);
            medie_y[j] = medie_y[j]/(N_misure-taglio);

            media_definitiva_x2=media_definitiva_x2+medie_x2[j];
            media_definitiva_dx2=media_definitiva_dx2+medie_dx2[j];
            media_definitiva_ene=media_definitiva_ene+medie_ene[j];
            media_definitiva_y=media_definitiva_y+medie_y[j];

        }
        media_definitiva_x2=media_definitiva_x2/M;
        media_definitiva_dx2=media_definitiva_dx2/M;
        media_definitiva_y=media_definitiva_y/M;
        media_definitiva_ene=media_definitiva_ene/M;
    
        for(j=0;j<M;j++){
            std_x2=std_x2+pow((medie_x2[j]-media_definitiva_x2),2)/(M-1);
            std_dx2=std_dx2+pow((medie_dx2[j]-media_definitiva_dx2),2)/(M-1);
            std_y=std_y+pow((medie_y[j]-media_definitiva_y),2)/(M-1);
            std_ene=std_ene+pow((medie_ene[j]-media_definitiva_ene),2)/(M-1);
            
        }
        std_x2=sqrt(std_x2);
        std_dx2=sqrt(std_dx2);
        std_y=sqrt(std_y);
        std_ene=sqrt(std_ene);
        fprintf(BOOT,"%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %f %d\n", media_definitiva_x2, std_x2, media_definitiva_dx2, std_dx2, media_definitiva_y, std_y, media_definitiva_ene, std_ene, e, enne);


    free(dati_x2);
    free(dati_dx2);
    free(dati_y);
    free(dati_ene);

    return;
}


/*==== NOTA: Nel file di input.txt classico non viene specificata la lunghezza del reticolo. Quindi è necessario
crearne uno che la contenga (io l'ho chiamato input2.txt)*/ 
void read_file_input(){ 

    FILE *in;
    in = fopen("/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/input.txt", "r");
    control_file(in);
    int x, idec, iflag, i_term; 
    float eta, d_metro;
    x=fscanf(in, "%f  %d  %d  %d  %d", &d_metro, &N_misure, &idec, &iflag, &i_term);
    fclose(in);
    return ;

}


int main(){

    read_file_input();  // funzione per leggere file con variabili di input (numero di misure, lunghezza del lato del reticolo, iflag per la partenza e campo esterno e idecorrel)

    char path[500];

    FILE *f, *file_bootstrap;         // puntatore al file degli eta
    int p=0, L;         // p serve per scorrere il file degli eta e L è la lunghezza dell'array degli Nlatt. 
    float ennepereta;   //è il prodotto di lunghezza path per valore di eta, è costante nel limite al continuo e sarà il primo elemento del file valori.txt nel caso 1
    float singoloeta;   //primo elemento del file valori.txt nel caso scelta = 0, eta fisso al variare degli N
    char dir[70];
    

    list *Nlatt=NULL;           // lista contenente i valori degli eta, letti da file, per caso scela=1
    list *N=NULL;              // lista contenente i valori degli N, letti da file, per caso scela=0
    
    
    printf("######################################################## \n");
    printf("       Per bootstrappare il limite al continuo digiti 1       \n");
    printf("######################################################## \n");
    printf("   Per bootstrappare la variazione di temperatura digiti 0    \n");
    printf("######################################################## \n");

    int scelta=0;
    //scanf("%d", &scelta);

    if(scelta==1){
        printf("############################################################# \n");
        printf("Ha scelto di bootstrappare il limite al continuo, buona giornata\n");
        printf("############################################################# \n");
        f=fopen("/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/valori_netafisso.txt","r");  //leggo il file valori.txt
        control_file(f);
        Nlatt = scan_file(f,Nlatt);         // lettura degli valori degli eta, con funzione ("listfunction.h")
        L=count(Nlatt)-1;         // funzione che conta la lunghezza di una lista ("listfunction.h")
        ennepereta = (int)val_posizione(0,Nlatt);
        Nlatt = move_to_position(1,Nlatt);

        //Cartella con i vari dati
        sprintf(dir,"N_ETA_fisso/Bootstrap");
    // Nome del file su cui salvare i dati bootstrappati
        char filemisure[700];
        sprintf(filemisure,"/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/%s/bootstrap_NEta=%.0f.txt", dir,ennepereta);
        file_bootstrap=fopen(filemisure,"w");
        control_file(file_bootstrap);
        fprintf(file_bootstrap, "# y2 / err_y2 / dy2 / err_dy2 / y / err_y / ene / err_ene / eta / Nlatt \n");

        FILE *misure[L];          // puntatore per i file che salvano gli obs, per ciascun eta
        while(Nlatt!=NULL){
            printf("Sono nel ciclo %d \n",p);
            char filemisure[700];
            float eta;
            eta=ennepereta/(Nlatt->val);

            // Lettura del file delle misure
            sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/N_ETA_fisso/N_Eta=%.0f/misure(N=%.0f)_eta=%.3f.txt", ennepereta,Nlatt->val, eta); //it modifies each time the name of the file to etae created
            misure[p]=fopen(filemisure, "r");
            control_file(misure[p]);


            printf("sono nella cartella %s \n", dir);
            int ennelatt;
            ennelatt=(int)Nlatt->val;


            bootstrap(misure[p], file_bootstrap, eta, ennelatt);

            printf("Ho letto il valore %f\n", Nlatt->val);     

            fclose(misure[p]);
            printf("chiuso misure\n");
            Nlatt=Nlatt->next;
            p++;
        }
    fclose(file_bootstrap);
    }

//CIao
// DOPO C'è DA METTERE SCELTA TRA LE OPZIONI DI HARMONIC_OSCILLATOR E FARGLI FARE COSE DIVERSE


    else {
        printf("######################################################## \n");
        printf("Ha scelto di variare la temperatura, buona giornata\n");
        printf("######################################################## \n");
        f=fopen("/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/valori_n.txt","r");  //leggo il file valori.txt
        control_file(f);
        N = scan_file(f,N);         // lettura degli valori degli eta, con funzione ("listfunction.h")
        L=count(N)-1;         // funzione che conta la lunghezza di una lista ("listfunction.h")
        singoloeta=val_posizione(0,N);
        N=move_to_position(1,N);

        //Cartella con i vari dati
        sprintf(dir,"N_variabile");
        // Nome del file su cui salvare i dati bootstrappati
        char filemisure[700];
        sprintf(filemisure,"/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/%s/Bootstrap/bootstrap_eta=%.3f.txt", dir, singoloeta);
        file_bootstrap=fopen(filemisure,"a");
        control_file(file_bootstrap);
        fprintf(file_bootstrap, "# y2 / err_y2 / dy2 / err_dy2 / y / err_y / ene / err_ene / eta / Nlatt \n");

        FILE *misure[L];          // puntatore per i file che salvano gli obs, per ciascun eta
        while(N!=NULL){
            char filemisure[500];
            sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/N_variabile/ETA=%.3f/misure(eta=%.3f)_N=%.3f.txt", singoloeta, singoloeta, N->val); //it modifies each time the name of the file to etae created
            misure[p]=fopen(filemisure, "r");
            control_file(misure[p]);

            sprintf(dir,"N_variabile/ETA=%.3f", singoloeta);

            bootstrap(misure[p], file_bootstrap, singoloeta, N->val);

            printf("Ho terminato il valore %f\n", N->val);

            fclose(misure[p]);
            N=N->next;
            p++;
        } 
        fclose(file_bootstrap);     
    }


    printf("############################################################# \n");
    printf("Ci auguriamo sia stata una piacevole esperienza, arrivederci!\n");
    printf("############################################################# \n");
    fclose(f);

    return 0;
}
