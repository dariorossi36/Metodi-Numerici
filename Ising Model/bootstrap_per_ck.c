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

int N, L; //N = numero di misure, L = lunghezza reticolo ----> li leggo da file


void bootstrap(FILE *f){
    FILE *m ,*e;
    unsigned long s;
    int i=0,j,k,M=100,x,bin=999,taglio=0;         /*M = numero di resampling,
    x serve per scanfare i file, i e j indici ricorsivi, bin è il bin del bootstrap,taglio è per levare le misure non termalizzate*/
    int V=N-taglio;
    long int seed = 13;         // seed per il ran2()
    float mag, ene, beta;         // risultati del simulatore_boot.c e valore del beta corrente
    int iter;           // risultati del simulatore_boot.c
    float *dati_mag = malloc(sizeof(long double) * (N));
    float *dati_ene = malloc(sizeof(long double) * (N));
    

    x=fscanf(f, "%f  %f  %d  %f", &mag, &ene, &iter, &beta);
    while (x!=EOF){
            dati_mag[i]=mag;
            dati_ene[i]=ene;
            x=fscanf(f, "%f  %f  %d  %f", &mag, &ene, &iter, &beta);
            i++;
        }
    /* Questo file va aperto quando sappiamo già il valore di beta*/
    char filemisure[500];
    sprintf(filemisure, 
"%s/Ising/Nlatt=%d/Bootstrap_noidec/Magnetizzazione/resample_binned%.3f.txt",
root_dir, L, beta);
    m=fopen(filemisure,"w");
    control_file(m);
    
    sprintf(filemisure, 
"%s/Ising/Nlatt=%d/Bootstrap_noidec/Energia/resample_binned%.3f.txt",
root_dir, L, beta);
    e=fopen(filemisure,"w");
    control_file(e);
    
    for(j=0;j<M;j++){
        for (k=taglio;k<N/(bin+1);k++){
            s = ran2(&seed) * RAND_MAX;
            s = (int)s;
            s <<= 15;
            s ^= (int)(ran2(&seed) * RAND_MAX);
            s %= V+1;
            s=s+taglio;
            //printf("%ld \n",s);
            for(int h=0;h<bin+1;h++){
                if(s+h>N-1){
                    fprintf(m,"%lf ",dati_mag[s+h+1-V]);
                    fprintf(e,"%lf ",dati_ene[s+h+1-V]);
                 //printf("sono dentro il primo if \n");
                }
                if(s+h<=N-1){
                    fprintf(m,"%lf ",dati_mag[s+h]);
                    fprintf(e,"%lf ",dati_ene[s+h]);
                 //printf("sono dentro il secondo if \n");
                }

            }
        }
        fprintf(m,"\n");
        fprintf(e,"\n");
    }

    free(dati_mag);
    free(dati_ene);
    fclose(f);
    fclose(m);
    fclose(e);
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
    DIR *folder;
    FILE *input;
    struct dirent *entry;
    int files = 0;
    
    char directory[500];
    sprintf(directory, "%s/Ising/Nlatt=%d/Risultati_noidec", root_dir, L);
    folder = opendir(directory);
    if(folder == NULL)
    {
        perror("Non riesco ad accedere alla cartella dei Risultati \n");
        return(1);
    }

    while( (entry=readdir(folder)) )
    {
        files++;
        
        printf("File %3d: %s\n", files, entry->d_name);

        if (!strcmp (entry->d_name, "."))
            continue;
        if (!strcmp (entry->d_name, ".."))
            continue;
        if (!strcmp (entry->d_name, ".DS_Store"))
            continue;
            
            
      
        char path[500];
        sprintf(path, 
"%s/Ising/Nlatt=%d/Risultati_noidec/%s",root_dir,L, entry->d_name); //it modifies each time the name of the file to be created
        
        input = fopen(path, "r");
        if (input == NULL){
            
            fprintf(stderr, "Error : Il file nella cartella Risultati è vuoto \n");
            exit(-1);
        }
        bootstrap(input);
    }

    closedir(folder);

    return(0);
}
