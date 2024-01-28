#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <errno.h>
#include"/home/dario/Documents/UNI/Metodi/usefulstuff.h"
#include"/home/dario/Documents/UNI/Metodi/ran2.h"
#include"/home/dario/Documents/UNI/Metodi/listfunction.h"

int n,idec; //n = numero di misure----> lo leggo da file 
int max_k=20;       // massimo valore di k (per il C(k))
float results[2]={0,0};

//# FUNZIONE PER FARE IL BOOTSTRAP
float * bootstrap(FILE *f, int N_misure){
    FILE *ckappa;
    unsigned long s; //serve per salvarsi il numero casuale enorme che produco
    int i=0,j,k,M=10,x,bin=1,taglio=0;         /*M = numero di resampling, 
    x serve per scanfare i file, i e j indici ricorsivi, bin è il bin del bootstrap, 
    taglio è per levare le misure non termalizzate*/
    int V=N_misure-taglio;
    long int seed = 13;         // seed per il ran2()
    float ck;         // risultati del simulatore_boot.c e valore del beta corrente
    int iter;           // risultati del simulatore_boot.c
    float *medie=malloc(sizeof(long double)*M);
    float media_definitiva=0, std=0;
    results[0]=0; results[1]=0;
    float *dati= malloc(sizeof(long double) * (N_misure));

    x=fscanf(f, "%f", &ck);

    while (x!=EOF){
            dati[i]=ck;
            x=fscanf(f, "%f", &ck);
            i++;
        }

    for(j=0;j<M;j++){
        medie[j]=0;
        for (k=0;k<(V)/bin;k++){

            s=ran2(&seed)*(V);
            s=(int)s;

            s=s+taglio;
            //printf("%ld \n",s);
            for(int h=0;h<bin;h++){
              //  printf("terzo for \n");
                if(s+h>N_misure-1){
                    medie[j]=medie[j]+dati[s+h+1-(V)];

                }
                if(s+h<=N_misure-1){
                    medie[j]=medie[j]+dati[s+h];
                }

            }
            
        }

        medie[j]=medie[j]/(V);
        media_definitiva=media_definitiva+medie[j];
    }

    media_definitiva=media_definitiva/M;
    
    for(j=0;j<M;j++){
        std=std+pow((medie[j]-media_definitiva),2)/(M-1);
            
    }
    std=sqrt(std);
    
    free(medie);
    free(dati);

    results[0]=media_definitiva;
    results[1]=std;

    return results;
    
}

// FUNZIONE PER LEGGERE FILE DI INPUT
/*==== NOTA: Nel file di input.txt classico non viene specificata la lunghezza del reticolo. Quindi è necessario
crearne uno che la contenga (io l'ho chiamato input2.txt)*/ 
void read_file_input(){ /*Mi serve solo per leggere il numero di misure N e la larghezza del reticolo L,
 gli altri non mi servono a niente, quindi non gli do nemmeno un nome significativo*/
    FILE *in;
    in = fopen("/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/input_split.txt", "r");
    control_file(in);
    int x, iflag, i_term; 
    float eta, d_metro;
    x=fscanf(in, "%f  %d  %d  %d  %d", &d_metro, &n, &idec, &iflag, &i_term);
    fclose(in);
    return ;
}


int main(){

    read_file_input();  // funzione per leggere file con variabili di input (numero di misure, lunghezza del lato del reticolo, iflag per la partenza e campo esterno e idecorrel)

    // NON SERVE char path[500];

    FILE *f, *ckappa;         // puntatore al file dei valori degli eta e al file dove stamperò il ck medio con errore
    int N_misure;
    int L;         // p serve per scorrere il file degli eta e L è la lunghezza dell'array degli Nlatt. 
    float ennepereta;   //è il prodotto di lunghezza path per valore di eta, è costante nel limite al continuo e sarà il primo elemento del file valori.txt nel caso 1
    float singoloeta;   //primo elemento del file valori.txt nel caso scelta = 0, eta fisso al variare degli N
    // NON SERVE char dir[70];
    int altro;

    list *Nlatt=NULL;           // lista contenente i valori degli eta, letti da file, per caso scela=1


    f=fopen("/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/valori_split.txt","r");  //leggo il file valori.txt
    control_file(f);
    Nlatt = scan_file(f,Nlatt);         // lettura dei valori degli eta, con funzione ("listfunction.h")
    L=count(Nlatt)-1;         // funzione che conta la lunghezza di una lista ("listfunction.h")
    ennepereta = val_posizione(0,Nlatt);       // il primo valore è il prodotto di Nlatt*eta
    Nlatt = move_to_position(1,Nlatt);      // dopodiché salvo i valori successivi (Nlatt)


    FILE *misure1, *misure2, *misure3, *misure4;          // puntatore per i file che salvano gli obs, per ciascun eta

    
    while(Nlatt!=NULL){
        char filemisure[700];
        float eta;
        eta=ennepereta/(Nlatt->val);
        printf("%.3f\n", eta);
        int ennelatt;
        ennelatt=(int)Nlatt->val;

        //apro il file in cui stampare i valori, che stampo con la funzione bootstrap direttamente
        char fileboot[700];
        sprintf(fileboot,"/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/Splitting/Bootstrap/bootstrap(eta=%f,N=%d).txt", eta, ennelatt); 
    
        ckappa=fopen(fileboot,"w");
        control_file(ckappa);

        float *risY=NULL, *risY2=NULL;  //puntatori che prendono i risultati del bootatrap per y e y2
        float Y=0, DY=0, Y2=0, DY2=0   //variabili dove mi salvo i risultati del bootstrap per y e y2

        //prima volta <Y>
        char string3[10]="<Y>";
        sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/Splitting/Campo/C(k)/%s(eta=%.3f,N_Eta=%.0f).txt",  string3, eta, ennepereta); //it modifies each time the name of the file to etae created
        misure3=fopen(filemisure, "r");
        control_file(misure3);
        altro=1;

        risY=bootstrap(misure3, n);
        Y=risY[0]; DY=risY[1];
        fclose(misure3); 
        //printf("%f\n", risCK[0]);

        //prima volta <Y2>
        char string4[10]="<Y2>";
        sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/Splitting/Campo/C(k)/%s(eta=%.3f,N_Eta=%.0f).txt", string4, eta, ennepereta); //it modifies each time the name of the file to etae created
        misure4=fopen(filemisure, "r");
        control_file(misure4);
        altro=1;

        risY2=bootstrap(misure4, n);
        Y2=risY2[0]; DY2=risY2[1];
        fclose(misure4);
        //printf("%f\n", risCK[0]);

        for(int k=0; k<max_k; k++){
            float *risCK=NULL, *risCK2=NULL;
            float ckY=0, DckY=0, ckY2=0, DckY2=0;
            altro = 0;      
            
            //ho 2 file quasi uguali di cui fare il bootstrap, cambia solo un 2 nel nome=> faccio 2 volte bootstrap a stesso eta

            // PRIMA VOLTA
            char string[10]="";     // stringa per inserire il 2 nel nome
            sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/Splitting/Campo/C(k)/Ck_Y%s(eta=%.3f,N_Eta=%.0f,k=%d).txt", string, eta, ennepereta,k); //it modifies each time the name of the file to etae created
            misure1=fopen(filemisure, "r");
            control_file(misure1);
            N_misure=n*idec;    //n è globale e anche idec
            
            risCK=bootstrap(misure1, N_misure);
            ckY=risCK[0]; DckY=risCK[1];
            fclose(misure1);
            //printf("%f\n", ckY);

            // SECONDA VOLTA
            char string2[10]="2";       // stringa per inserire il 2 nel nome
            sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/Splitting/Campo/C(k)/Ck_Y%s(eta=%.3f,N_Eta=%.0f,k=%d).txt", string2, eta, ennepereta,k); //it modifies each time the name of the file to etae created
            misure2=fopen(filemisure, "r");
            control_file(misure2);
            N_misure=n*idec;
            
            risCK2=bootstrap(misure2, N_misure);
            ckY2=risCK2[0]; DckY2=risCK2[1];
            fclose(misure2); 
            //printf("%f\n", ckY);

            printf("fine k=%d \n", k);
            fprintf(ckappa, "%d  %f  %f  %f  %f  %f  %f  %f  %f \n", k, ckY, DckY, ckY2, DckY2, Y, DY, Y2, DY2);
        }

        printf("Ho letto il valore %f\n", Nlatt->val);
        fclose(ckappa);
        
        Nlatt=Nlatt->next;

    }

    fclose(f);
    

    printf("############################################################# \n");
    printf("Ci auguriamo sia stata una piacevole esperienza, arrivederci!\n");
    printf("############################################################# \n");
    

    return 0;
}
