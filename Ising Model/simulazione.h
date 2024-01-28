#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include"/home/dario/Documents/UNI/Metodi/ran2.h"
#include "/home/dario/Documents/UNI/Metodi/usefulstuff.h"

/* Programma per la simulazione del modello di ising bidimensionale con 
possibilita` di inserire un campo magnetico esterno */


/* ------------------------------------- DICHIARAZIONE DELLE VARIABILI 
GLOBALI ------------------------------------------- */

/* !!!!!!!!!!!!!!!! CARTELLA ROOT, DA CAMBIARE SOLO QUESTA QUANDO CAMBIA COMPUTER !!!!!!!!!!!!!!!! */
char root_dir[500] = "/home/dario/Documents/UNI/Metodi/Modulo1";

/* PROVA PER LEGGERE N DA FILE */
#define Nmax 200       // dimensione del reticolo (eccessiva per allocare memoria sufficiente)
int field[Nmax][Nmax]; //campo
int npp[Nmax], nmm[Nmax]; //array per definire le posizioni dei primi vicini del lattice
/*FINE PROVA */
int N; //lunghezza lato reticolo
float extfield; //nomi dei parametri della simulazione (beta = 1/kT | extfield = campo esterno)
float xmagn = 0, xene = 0;

long int seed = 138;
float medie[2]={0,0};

unsigned long acc=0;  //intero che mi calcola quante volte ho accettato il metropolis, per poi calcolare l'accettanza

// ! RICORDARE DI CAMBIARE ANCHE LA LUNGHEZZA DELL'ARRAY DEI BETA, NEL MAIN !





/* ------------------------------------- INIZIO SUBROUTINES 
------------------------------------------- */


/* per ogni coordinata definisco il passo in avanti o indietro con le 
opportune condizioni al bordo*/
void geometry(){
    
    for(int i=0; i<N; i++){
        npp[i] = i + 1;
        nmm[i] = i - 1;
    }
    
    npp[N-1] = 0; //Condizioni al
    nmm[0] = N-1; //bordo periodiche
    
    return;
        
}

/* -------------------------*/

/* Assegno la configurazione di partenza della catena di Markov */
void initialize_lattice(int iflag){
    
    float x;
    
    /* PARTENZA A FREDDO (tutti gli spin a 1, o -1, come se fosse T = 0) */
    if(iflag == 0){
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                
                field[i][j] = 1;
            }
        }
    }
    /* PARTENZA A CALDO (spin random, come se fosse T = infinito) */
    if(iflag == 1){
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                
                //float x;
                //seed = i;
                x = ran2(&seed);
                field[i][j] = 1;
                if(x<0.5) field[i][j] = -1;
            }
        }
    }
    
    // AGGIUNGERE IF PER RIPARTIRE DALLA CONFIGURAZIONE PRECEDENTE
    
    
    return;
}


/* -------------------------*/

/* Calcolo magnetizzazione media del reticolo */
float magnetization(float xmagn){
    
    int nvol = N*N; //dimensioni del reticolo (volume)
    
    xmagn = 0; // inizializzazioine xmagn
    
    for(int i = 0; i<N; i++){            // faccio il loop su
        for(int j = 0; j<N; j++){        // su tutto il reticolo
            
            xmagn = xmagn + field[i][j];     // e sommo su tutti i valori del campo
        }
    }
    xmagn = xmagn/nvol;                      // normalizzo dividendo per il volume
    return xmagn;
}

/* -------------------------*/

/* Energia media (= 0 per configurazione ordinata e campo esterno 0) */
float energy(float xene){
    
    int ip, im, jp, jm;
    float force;
    int nvol = N*N; //dimensioni del reticolo (volume)
    
    xene = 0;                             // inizializzazioine xene
    
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            
            ip = npp[i];        // calcolo coordinate
            im = nmm[i];        // dei primi vicini
            jp = npp[j];        // del sito
            jm = nmm[j];
            
            force = field[i][jp] + field[i][jm] + field[ip][j] + 
field[im][j]; // somma dei 4 primi vicini
            xene = xene - 0.5*force*field[i][j]; // 1/2 per il conteggio giusto
            xene = xene - extfield*field[i][j]; // contributo campo esterno
        }
    }
    
    xene = xene/nvol;       // normalizzo dividendo per il volume
    return xene;
}

/* -------------------------*/

/* faccio aggiornamenti locali delle variabili di spin con metropolis, la 
variabile di spin di prova e` sempre quella opposta a quella attuale */
void update_metropolis(float beta){
    
    int i, j, ip, im, jp, jm;
    float force, phi, p_rat, x; //acc = 0;
    int nvol = N*N; //dimensioni del reticolo (volume)
    
    for(int ivol=0; ivol<nvol; ivol++){     // loop su tutti i siti del reticolo
        
        /* sclego a caso un sito del reticolo */
        
        //seed = ivol;
        
        i = ran2(&seed)*N;                       // scelgo a caso un sito del reticolo
        
    
        j = ran2(&seed)*N;
        
        /*
        i = rand()%nlatt;                   // genera direttamente un numero intero tra 0 e nlatt
        j = rand()%nlatt;
        */
        ip = npp[i];                        // calcolo le coordinate dei quattro primi vicini
        im = nmm[i];                        // del sito che ho scelto
        jp = npp[j];
        jm = nmm[j];
        
        force = field[i][jp] + field[i][jm] + field[ip][j] + field[im][j];
        force = beta*(force+extfield);
        
        phi = field[i][j];
        
        p_rat = exp(-2*phi*force);
        
    
        x = ran2(&seed);

        
        
        if(x<p_rat) {
            field[i][j] = - phi;
            acc++;
        }
        
    }
    
    return;
}


    /* ------------------------------------- FINE SUBROUTINES 
------------------------------------------- */

    /* ------------------------------------- INIZIO MAIN 
------------------------------------------- */

void simulazione(float beta, FILE *misure){
    
    FILE *f  ; //*l;
    int x; //variabili per leggere il file

    float m_magn = 0, m_ene = 0;
    
    // input
    int iflag, measures, i_decorrel; //nomi dei parametri della simulazione
                                    //(iflag per partenza caldo(1), freddo(0),...| measures = n° misure|
                                    // i_decorrel = numero di passi della catena Markov tra una misura e l'altra)


    
    /* PROVA PER LEGGERE N DA FILE
    int npp[N], nmm[N];
    int nlatt = N, nvol = N*N; //dimensioni del reticolo (lunghezza del lato e volume)
    int field[N][N]; //campo
     FINE PROVA */
    
    /* lettura file input */
    char fileinput[200];
    sprintf(fileinput, "%s/Ising/input.txt", root_dir);
    f=fopen(fileinput,"r");
    control_file(f);
    
    /* lettura dei parametri della simulazione */
    //int a; //è la larghezza di reticolo, non mi serve a niente qua
    x = fscanf(f, "%d  %d  %d  %f  %d", &iflag, &measures, &i_decorrel,
&extfield, &N);
    
    /* OPERAZIONI PRELIMINARI */
    geometry(); //inizializza condizioni al bordo
    initialize_lattice(iflag); //inizializza configurazione iniziale
    


    for(int iter=0; iter<measures; iter++){
        
        /* AGGIORNAMENTO CONFIGURAZIONE: i_decorrel spazza tutto il reticolo con l'algoritmo scelto */
        for(int idec=0; idec<i_decorrel; idec++){
            update_metropolis(beta);
        }
        fprintf(misure, "%f  %f  %d  %f\n", magnetization(xmagn), 
energy(xene), iter, beta );
        /* MISURA DELLE VARIABILI FISICHE */
        m_magn = m_magn + magnetization(xmagn);
        m_ene = m_ene + energy(xene);
        
    }
    
    m_magn = m_magn/measures;
    m_ene = m_ene/measures;
    medie[0]=m_magn;
    medie[1]=m_ene;
  //  medie = {m_magn, m_ene};

    /* --------------------------- TERMINE SIMULAZIONE MONTECARLO 
--------------------------- */

    /* ------------------ CALCOLO ACCETTANZA ----------------------------*/
    
    float accettanza = 0;
    unsigned long long num_updates = 0;
    num_updates=i_decorrel*(N*N);//(Nlatt-2);
    num_updates = measures*num_updates;
    accettanza=(float)acc*100/measures;
    accettanza=(float)accettanza/i_decorrel;
    accettanza=(float)accettanza/(N*N);
   
    printf("\n \n");

    printf("Nlatt=%d, accettanza=%.2f percento,\n numero di update = %llu, numero di accettati %lu\n", N, accettanza, num_updates, acc);
    
    acc = 0;
    

    fclose(f);

	
	return;
}
