#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include"/home/dario/Documents/UNI/Metodi/ran2.h"
#include"/home/dario/Documents/UNI/Metodi/usefulstuff.h"
#include"/home/dario/Documents/UNI/Metodi/listfunction.h"
/* Programma per la simulazione dell'oscillatore armonico*/

#define Nmax 10000000
int Nlatt;
long int seed = 23456789;
float d_metro; //eta=a*omega = parametro reticolo * pulsazione;  d_metro = parametro del metropolis = 2*sqrt(eta) 
int iflag, measures, i_decorrel, i_term;  //vedi ising. i_term = passo di termalizzazione
int npp[Nmax], nmm[Nmax]; //array per definire le posizioni dei primi vicini del lattice
float field[Nmax];  //array in cui salviamo le Nlatt posizioni del percorso 
unsigned long acc=0;  //intero che mi calcola quante volte ho accettato il metropolis, per poi calcolare l'accettanza



/* per ogni coordinata definisco il passo in avanti o indietro con le opportune condizioni al bordo*/
void geometry(){
    
    for(int i=0; i<Nlatt; i++){
        npp[i] = i + 1;
        nmm[i] = i - 1;
    }
    
    npp[Nlatt-1] = 0; //Condizioni al
    nmm[0] = Nlatt-1; //bordo periodiche
    
    return;
        
}

/* --------------- REFUSO DAL MODELLO DI ISING ------------------*/

/* Assegno la configurazione di partenza della catena di Markov */


void initialize_lattice(int iflag){
    float x;
    
    // PARTENZA A FREDDO (tutti gli spin a 1, o -1, come se fosse T = 0) 
    if(iflag == 0){
        for(int i = 0; i<Nlatt; i++){
                
            field[i] = 0;
            
        }
    }
    // PARTENZA A CALDO (spin random, come se fosse T = infinito) 
    if(iflag == 1){
        for(int i = 0; i<Nlatt-1; i++){
            //float x;
            //seed = i;
            x = 1-2*ran2(&seed);
            field[i] = x;
            //printf("%f\n", field[i]);
        }
        field[Nlatt-1]=field[0];
    }
    
    // AGGIUNGERE IF PER RIPARTIRE DALLA CONFIGURAZIONE PRECEDENTE
    
    
    return;
}



/*funzione per avanzare col metropolis e modificare la configurazione*/

void update_metropolis(float eta){
    float c1, c2; //sono solo shortcut per scrivere 1/eta e altra funz di eta
    int ip, im; //coordinate dei 2 primi vicini
    float force, phi, phi_prova;  //force = forza del campo intorno; phi = valore attuale del campo; phi_prova = valore di prova del campo
    float p_rat, x; // p_rat = valore di energia con cui confrontare x che è estratto a caso e mi dirà se accettare il passo o no
    
    c1 = 1/eta;
    c2 = 1/eta + eta/2;

    /*loop su tutti i siti, qui il sito non è scelto a caso ma faccio una spazzata 
    iterativa su tutti i siti, si può dimostrare che va bene lo stesso per il bilancio dettagliato, 
    ma meno banale da provare*/
    
    for(int i = 0; i<Nlatt; i++){
        ip = npp[i];
        im = nmm[i];
        force = field[ip] + field[im];

        phi = field[i];
        phi_prova = phi + 2*d_metro*(0.5 - ran2(&seed));  // Il secondo termine serve per traslare nell'intervallo [-d_metro, d_metro]. Calcolo del phi di prova per trovare la probabilità e decidere se accettare il passaggio o no
        p_rat= c1 * phi_prova * force - c2*pow(phi_prova,2);
        p_rat= p_rat - c1*phi*force + c2*pow(phi,2);
        //printf("%f\n", p_rat);
        x = log(ran2(&seed));  //METRO-TEST! 
        
        if(x<p_rat){
          field[i]=phi_prova;  //test accettanza, se p_rat>1 accetto
          acc++;
        }
        
    }
    
    return ;
}


//funzione che prende le misure delle osservabili che ci interessano

void measure(float eta, FILE* campo1, FILE* campo2){
    for(int i=0; i<Nlatt; i++){

        fprintf(campo1,"%lf  ", field[i]);
        fprintf(campo2, "%lf  ", pow(field[i],2));
    }
        fprintf(campo1,"\n ");
        fprintf(campo2, "\n");
    return ;
}



/*=================================== SIMULAZIONE =============================================*/

// y è Nlatt o eta e viene letto come lista dopo primovalore dal file valori.txt

void Harmonic_metropolis(float y, float primovalore){
    FILE * lat, *input; // file in cui stampo il field (lat) e da cui prendo i valori iniziali ( init )
    int x; //per leggere l'init.txt; l è l'Nlatt che qui non serve a niente
    //OPERAZIONI PRELIMINARI
    float eta;
    input = fopen("/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/input_split.txt","r");
    //printf("1\n");

    control_file(input);
    //printf("4\n");
    x= fscanf(input, "%f  %d  %d  %d  %d", &d_metro, &measures, &i_decorrel, &iflag, &i_term);

    Nlatt=y;
    eta=primovalore/Nlatt;
    if(eta<0.05) d_metro=0.3;
    else if(eta<0.09 && eta>=0.05) d_metro=0.5;
    else if(eta>=0.09 && eta<0.15) d_metro=0.8;
    else if(eta>=0.15 && eta<0.45) d_metro=1;
    else if(eta>=0.45) d_metro=3;
    Nlatt=(int)Nlatt;
        

    FILE* campo1;
    FILE* campo2;
    char filename[700];
    sprintf(filename, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/Splitting/Campo/Valoricampo(Eta=%.3f)(NEta=%.0f).txt", eta, primovalore);
    campo1=fopen(filename,"w");
    control_file(campo1);
    char filenome[700];
    sprintf(filenome, "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/Splitting/Campo/Valoricampo2(Eta=%.3f)(NEta=%.0f).txt", eta, primovalore);
    campo2=fopen(filenome,"w");
    control_file(campo2);


    initialize_lattice(iflag);
    geometry();
    //printf("8\n");
    //SESSIONE ALL'EQUILIBRIO con MISURE

    for (int i = 1; i<i_term+1; i++){
        update_metropolis(eta);
    }

    acc=0;
    double obs[3]={0,0,0};
    for (int iter=0; iter<measures; iter++){
        // AGGIORNAMENTO CONFIGURAZIONE
        for(int idec=0; idec<i_decorrel; idec++ ){

            update_metropolis(eta);
    }
        measure(eta, campo1, campo2);
}

    float accettanza;
    unsigned long num_updates;
    num_updates=measures*i_decorrel*(Nlatt-2);
    accettanza=(float)acc*100/num_updates;

    printf("\n \n");

    printf("delta=%f, eta=%f, Nlatt=%d, accettanza=%.2f percento,\n numero di update = %lu, numero di accettati %lu\n", d_metro, eta, Nlatt, accettanza, num_updates, acc);
    //fclose(lat);
    fclose(input);
    fclose(campo1);
    fclose(campo2);
    return;

}