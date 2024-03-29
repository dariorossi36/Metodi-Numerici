#include "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/oscillazione.h"

/* IN QUESTO CODICE VIENE ESEGUITA LA SIMULAZIONE MONTECARLO PER VARI VALORI DI eta, RIPORTATI NEL FILE "eta.txt". 
PER OGNI VALORE NEL FILE "misure%.1f.txt" SONO RIPORTATE LE MEDIE DELLA MAGNETIZZAZIONE E DELL'ENERGIA PER CIASCUN 
VALORE DI eta */

/*La lunghezza del reticolo devve essere specificata nel file oscillazione.h,
        mentre il numero di misure nel file input.txt, che verrà letto in oscillazione.h*/


int main(void){
    FILE *f;         // puntatore al file degli eta
    int p=0, L;         // p serve per scorrere il file degli eta e L è la lunghezza dell'array degli eta . 
    float ennepereta;
    float singoloeta;
    float eta;

    list *Nlatt=NULL;           // lista contenente i valori degli eta, letti da file
    list *N=NULL;
    
    char valori[300];
    printf("######################################################## \n");
    printf("Per simulare il limite al continuo (sia eta che variabili, ma N*eta fisso) digiti 1\n");
    printf("######################################################## \n");
    printf("Per simulare la variazione di temperatura (N variabile, eta fisso) digiti 0\n");
    printf("######################################################## \n");

    int scelta=0;
    //scanf("%d", &scelta);

    if(scelta==1){
        printf("############################################################# \n");
        printf("Ha scelto di simulare il limite al continuo, buona giornata\n");
        printf("############################################################# \n");

        //apro il file valori_n.txt che contiene alla prima riga un valore di eta fissato, e le righe successive
        //valori di N (cioè passi di path) variabili
        sprintf(valori,"/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/valori_netafisso.txt");
        f=fopen(valori,"r");
        control_file(f);

        Nlatt = scan_file(f,Nlatt);         // lettura degli valori degli eta, con funzione ("listfunction.h")
        L=count(Nlatt)-1;         // funzione che conta la lunghezza di una lista ("listfunction.h")
        ennepereta=val_posizione(0,Nlatt);
        Nlatt=move_to_position(1,Nlatt);
        FILE *misure[L];          // puntatore per i file che salvano gli obs, per ciascun eta
        while(Nlatt!=NULL){
            char filemisure[500];
            eta=ennepereta/(Nlatt->val);
            int enne;
            enne=(int)Nlatt->val;
            sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/N_ETA_fisso/N_Eta=%.0f/misure(N=%d)_eta=%.3f.txt", ennepereta, enne, eta); //it modifies each time the name of the file to etae created
            misure[p]=fopen(filemisure, "w");
            control_file(misure[p]);

            Harmonic_metropolis(Nlatt->val, misure[p], scelta, ennepereta); 
             /*eseguo la simulazione per il valore di eta in causa. 
            Funzione di oscillazione.h Mi creerà automaticamente un file di misure contenente i valori delle osservabili che
            ci interessano (y^2medio, Dy^2medio, ymedio) e eta associato per ogni misura fatta. */
     
            printf("Ho terminato il valore %f\n", Nlatt->val);     

            fclose(misure[p]);
            Nlatt=Nlatt->next;
            p++;
        }
    }


// DOPO C'è DA METTERE SCELTA TRA LE OPZIONI DI HARMONIC_OSCILLATOR E FARGLI FARE COSE DIVERSE


    else {
        printf("######################################################## \n");
        printf("Ha scelto di variare la temperatura, buona giornata\n");
        printf("######################################################## \n");

        //apro il file valori_n.txt che contiene alla prima riga un valore di eta fissato, e le righe successive
        //valori di N (cioè passi di path) variabili
        sprintf(valori,"/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/valori_n.txt");
        f=fopen(valori,"r");
        control_file(f);
        printf("ho letto il file %s", valori);
        N = scan_file(f,N);         // lettura degli valori degli eta, con funzione ("listfunction.h")
        printf("ho letto il file %s", valori);
        L=count(N)-1;         // funzione che conta la lunghezza di una lista ("listfunction.h")
        singoloeta=val_posizione(0,N);
        N=move_to_position(1,N);
        FILE *misure[L];          // puntatore per i file che salvano gli obs, per ciascun eta
        while(N!=NULL){
            char filemisure[500];
            sprintf(filemisure, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/N_variabile/ETA=%.3f/misure(eta=%.3f)_N=%.3f.txt", singoloeta, singoloeta, N->val); //it modifies each time the name of the file to etae created
            misure[p]=fopen(filemisure, "w");
            control_file(misure[p]);

            Harmonic_metropolis(N->val, misure[p], scelta, singoloeta); 
             /*eseguo la simulazione per il valore di eta in causa. 
            Funzione di oscillazione.h Mi creerà automaticamente un file di misure contenente i valori delle osservabili che
            ci interessano (y^2medio, Dy^2medio, ymedio) e eta associato per ogni misura fatta. */
            printf("Ho terminato il valore %f\n", N->val);

            fclose(misure[p]);
            N=N->next;
            p++;
        }       
    }


    printf("############################################################# \n");
    printf("Ci auguriamo sia stata una piacevole esperienza, arrivederci!\n");
    printf("############################################################# \n");
    fclose(f);
    return 0;
}
