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

/* !!!!!!!!!!!!!!!! CARTELLA ROOT, DA CAMBIARE SOLO QUESTA QUANDO CAMBIA COMPUTER !!!!!!!!!!!!!!!! */
char root_dir[200] = "/home/dario/Documents/UNI/Metodi/Modulo1";


int M=170,bin=1000;
int N, L; //N = numero di misure, L = lunghezza reticolo ----> li leggo da file

void bootstrap(FILE *f, FILE *boot){    //*f = file delle misure da bootsrappare
    unsigned long s;
    int i=0,j,k,x,taglio=0;         /*M = numero di resampling,
    x serve per scanfare i file, i e j indici ricorsivi, bin è il bin del bootstrap,taglio è per levare le misure non termalizzate*/
    long int seed = 17;         // seed per il ran2()
    float mag, ene, beta;         // risultati del simulatore_boot.c e valore del beta corrente
    int iter;           // risultati del simulatore_boot.c
    int vol=L*L;

    //allocazione dati da bootstrappare: i numeri in fondo indicano la potenza a cui è elevata la quantità in questione
    printf("Numero di misure: %d \n", N);
    
    float *dati_mag = malloc(sizeof(long double) * (N)); 
    float *dati_mag2 = malloc(sizeof(long double) * (N));
    float *dati_mag4 = malloc(sizeof(long double) * (N));
    float *dati_ene = malloc(sizeof(long double) * (N));
    float *dati_ene2 = malloc(sizeof(long double) * (N));
 
    //il singolo elem degli array medie_ è la media di un solo resampling.
    // Nelle medie_definitive_ ci vanno poi le medie delle medie dei resampling
    float media_def_ene=0, media_def_mag=0, std_ene=0, std_mag=0, media_calore=0, media_susce=0, media_binder=0, std_calore=0, std_susce=0, std_binder=0;
    float medie_ene[M],medie_ene2[M], medie_mag[M], medie_mag2[M], medie_mag4[M], susce[M], calore[M], binder[M];  
    
    // Lettura delle misure prodotte dalla simulazione
    x=fscanf(f, "%f  %f  %d  %f", &mag, &ene, &iter, &beta);
    while (x!=EOF){
            dati_mag[i] = fabs(mag);
            dati_mag2[i] = pow(mag,2);
            dati_mag4[i] = pow(mag,4);
            dati_ene[i] = fabs(ene);
            dati_ene2[i] = pow(ene,2);
            x=fscanf(f, "%f  %f  %d  %f", &mag, &ene, &iter, &beta);
            i++;
        }
    
    printf("Volume del reticolo: %d\n", L);

    //ciclo per fare media delle medie dei resampling
    for(j=0;j<M;j++){
        //inizializzazione a 0 delle medie di ogni resampling
        medie_ene[j]=0;
        medie_ene2[j]=0;
        medie_mag[j]=0;
        medie_mag2[j]=0;
        medie_mag4[j]=0;
        calore[j]=0;
        susce[j]=0;
        binder[j]=0;
        //creazione di un resampling e media
        for (k=0;k<(N-taglio)/bin;k++){

            s=ran2(&seed)*(N-taglio);
            s=(int)s;
            s=s+taglio;

            for(int h=0;h<bin;h++){ //condizioni al bordo per il binning
                if(s+h>N-1){
                    medie_ene[j]=medie_ene[j]+fabs(dati_ene[s+h+1-(N-taglio)]);
                    medie_mag[j]=medie_mag[j]+fabs(dati_mag[s+h+1-(N-taglio)]);
                    medie_ene2[j]=medie_ene2[j]+(dati_ene2[s+h+1-(N-taglio)]);
                    medie_mag2[j]=medie_mag2[j]+(dati_mag2[s+h+1-(N-taglio)]);
                    medie_mag4[j]=medie_mag4[j]+(dati_mag4[s+h+1-(N-taglio)]);
                }
                if(s+h<=N-1){
                    medie_ene[j]=medie_ene[j]+fabs(dati_ene[s+h]);
                    medie_mag[j]=medie_mag[j]+fabs(dati_mag[s+h]);
                    medie_ene2[j]=medie_ene2[j]+(dati_ene2[s+h]);
                    medie_mag2[j]=medie_mag2[j]+(dati_mag2[s+h]);
                    medie_mag4[j]=medie_mag4[j]+(dati_mag4[s+h]);
                }
            }
        }
        
        // media del resampling j
        medie_ene[j]=medie_ene[j]/(N-taglio);
        medie_mag[j]=medie_mag[j]/(N-taglio);
        medie_ene2[j]=medie_ene2[j]/(N-taglio);
        medie_mag2[j]=medie_mag2[j]/(N-taglio);
        medie_mag4[j]=medie_mag4[j]/(N-taglio);

        //definisco qui calore specifico, suscettività e cumulante di binder
        calore[j] = vol*(medie_ene2[j] - medie_ene[j]*medie_ene[j]);
        susce[j] = vol*(medie_mag2[j] - pow(medie_mag[j],2));
        binder[j] = medie_mag4[j]/pow(medie_mag2[j],2);

        media_def_ene=media_def_ene+medie_ene[j];
        media_def_mag=media_def_mag+medie_mag[j];
        media_calore=media_calore + calore[j];
        media_susce=media_susce + susce[j];
        media_binder=media_binder + binder[j];
    }

    media_def_ene=media_def_ene/M;
    media_def_mag=media_def_mag/M;
    media_calore=media_calore/M;
    media_susce=media_susce/M;
    media_binder=media_binder/M;

    //ciclo per deviazione standard dei resampling intorno a media_def_
    for(j=0;j<M;j++){
        std_ene=std_ene+pow((medie_ene[j]-media_def_ene),2)/(M-1);
        std_mag=std_mag+pow((medie_mag[j]-media_def_mag),2)/(M-1);
        std_calore=std_calore+pow((calore[j]-media_calore),2)/(M-1);
        std_susce=std_susce+pow((susce[j]-media_susce),2)/(M-1);
        std_binder=std_binder+pow((binder[j]-media_binder),2)/(M-1);
    }
    std_ene=sqrt(std_ene);
    std_mag=sqrt(std_mag);
    std_calore=sqrt(std_calore);
    std_susce=sqrt(std_susce);
    std_binder=sqrt(std_binder);

    //stampo su file medie_definitive_ e relative deviazioni standard e beta associato
    fprintf(boot,"%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f %.3f\n", media_def_ene, std_ene, media_def_mag, std_mag, media_calore, std_calore, media_susce, std_susce, media_binder, std_binder, beta);

    free(dati_mag);
    free(dati_ene);
    free(dati_mag2);
    free(dati_ene2);
    free(dati_mag4);

    return;
}


/*==== NOTA: Nel file di input.txt classico non viene specificata la 
lunghezza del reticolo. Quindi è necessario
crearne uno che la contenga (io l'ho chiamato input2.txt)*/ 
void read_file_input(){ /*Mi serve solo per leggere il numero di misure N 
e la larghezza del reticolo L,
 gli altri non mi servono a niente, quindi non gli do nemmeno un nome 
significativo*/
    FILE *in;
    char fileinput[200];
    sprintf(fileinput, "%s/Ising/input.txt", root_dir);
    in = 
fopen(fileinput, "r");
    control_file(in);
    int x, part, idec; // part = partenza (caldo/freddo)
    float c_ext; // valore del campo esterno
    x=fscanf(in, "%d  %d  %d  %f  %d", &part, &N, &idec, &c_ext, &L);
    fclose(in);
    return ;
}


int main()
{
    read_file_input();  // funzione per leggere file con variabili di input (numero di misure, lunghezza del lato del reticolo, iflag per la partenza e campo esterno e idecorrel)
    FILE *misure, *boot, *be;       // file dei dati da bootstrappare
    int  y; // variabile che serve per leggere da file medie.txt
    float beta;
    
    /* Apertura del file totale in cui salvare medie e std per ogni beta*/
    char fileboot[500];
    sprintf(fileboot,"%s/Ising/Nlatt=%d/Bootstrap/bootstrap(res=%d).txt",root_dir, L, M);
    boot=fopen(fileboot,"w");
    control_file(boot);
    fprintf(boot, "# ene / err_ene / mag / err_mag / calore / err_calore / susce / err_susce / binder / err_binder/ beta \n");
    
    //apertura e lettura file medie
    char filebeta[500];
    sprintf(filebeta, "%s/Ising/beta.txt", root_dir);
    be=fopen(filebeta,"r");
    printf("%s\n", filebeta);
    control_file(be);

    y=fscanf(be, "%f", &beta);
    printf("%f\n",beta);
    while(y!=EOF){
    /* Ciclo che scorre file nella cartella dei risultati */
        char path[500];
        sprintf(path, "%s/Ising/Nlatt=%d/Risultati/misure(%.3f).txt",root_dir,L, beta); // crea il nome del file dei risultati da analizzare
        printf("%s\n", path);
        misure = fopen(path, "r");
        control_file(misure);
        
        // Chiama funzione bootstrap, passando file dei dati da bootstrappare
        bootstrap(misure, boot);
        y=fscanf(be, "%f",&beta);
        // Chiusura file misure da bootstrappare
        fclose(misure);
    }

    fclose(boot);
    fclose(be);

    return(0);
}

