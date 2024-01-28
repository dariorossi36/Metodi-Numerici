#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"/home/dario/Documents/UNI/Metodi/listfunction.h"
#include"/home/dario/Documents/UNI/Metodi/Modulo1/Ising/simulazione.h"



/* IN QUESTO CODICE VIENE ESEGUITA LA SIMULAZIONE MONTECARLO PER VARI 
VALORI DI BETA, RIPORTATI NEL FILE "beta.txt". 
PER OGNI VALORE VIENE SALVATA LA CONFIGURAZIONE FINALE DEL RETICOLO NEI 
FILE "Lattice.txt" NELLA CARTELLA DELL'Nlatt RELATIVO E NEL FILE "medie.txt"
SONO RIPORTATE LE MEDIE DELLA MAGNETIZZAZIONE E DELL'ENERGIA PER CIASCUN 
VALORE DI BETA */

/*La lunghezza del reticolo devve essere specificata nel file 
simulazione_boot.h, mentre il numero di misure nel file input.txt, che verrà letto in simulazione_boot.h*/





int main(void){
	FILE *f,*g, *q;         // puntatori ai file dei beta e delle medie e dell'input
	int p=0, L, s=0;         // p serve per scorrere il file dei beta e L è la lunghezza dell'array dei beta. s serve a scorrere i file di beta di cui voglio stampare il lattice
    
    // variabili di input
    int iflag, measures, i_decorrel, N; //nomi dei parametri della simulazione
                                    //(iflag per partenza caldo(1), freddo(0),...| measures = n° misure|
                                    // i_decorrel = numero di passi della catena Markov tra una misura e l'altra)
    float extfield; //nomi dei parametri della simulazione (beta = 1/kT | extfield = campo esterno)

    
	float beta_interesting[]={0.43, 0.45, 0};   //array contenente i beta di cui voglio vedere il lattice. Assicurarsi che i valori inseriti sianopresenti nel file beta_init.txt
    printf("%f\n",beta_interesting[0]);
    // beta usati nella simulazione
    char filebeta[500];
    sprintf(filebeta, "%s/Ising/beta.txt", root_dir);
	f=fopen(filebeta,"r");
	control_file(f);
    
    // LETTURA VALORI DI INPUT + PROVA CARTELLA ROOT
    
  
    /* apertura dei file di input e di output */
    char fileinput[500];
    sprintf(fileinput, "%s/Ising/input.txt", root_dir);
    q=fopen(fileinput,"r");
    control_file(q);

    
    /* lettura dei parametri della simulazione */
    int x; // variabile per leggere il file
    x = fscanf(q, "%d  %d  %d  %f  %d", &iflag, &measures, &i_decorrel,
&extfield, &N);
    
    
    // file dove sono salvate medie su tutto il reticolo (media mag e ene per ogni beta). File unico con tanti beta diversi
    char filemedie[500];
    sprintf(filemedie,"%s/Ising/Nlatt=%d/medie.txt", root_dir,N);
    
	g=fopen(filemedie,"w");
	control_file(g);
    
    
  
    
    // LETTURA DEI BETA DA LISTA
	list *b=NULL;           // lista contenente i valori dei beta, letti da file
	b = scan_file(f,b);         // lettura dei valori dei beta, con funzione ("listfunction.h")
	L=count(b);         // funzione che conta la lunghezza di una lista ("listfunction.h")
    
    // FILE
	FILE *beta[3], *misure[L];          // puntatore per i file che salvano la configurazione, per ciascun beta
    
    // SCORRO BETA
	while(b!=NULL){

        char filemisure[500];
        sprintf(filemisure, "%s/Ising/Nlatt=%d/Risultati/misure(%.3f).txt", root_dir, N, b->val); // file diverso per ogni beta e per ogni N. Ci sono le misure per ogni update
        misure[p]=fopen(filemisure, "w");
        control_file(misure[p]);

        
		simulazione(b->val, misure[p]); /* eseguo la simulazione per il valore di beta in causa. Funzione di simulazione_boot.h Mi creerà automaticamente un file di misure contenente valore di magnetizz, energia e beta associato per ogni misura fatta. */
		fprintf(g," %f %f %f \n",medie[0],medie[1],b->val);  
//stampo medie di magn e energia con beta associato in file medie.txt
		
		if (b->val==beta_interesting[0]||b->val==beta_interesting[1]||b->val==beta_interesting[2]){ 
//faccio i file contenenti i lattice interessanti

			printf("%f\n", b->val);
	       	char filename[500];
			sprintf(filename, "%s/Ising/Nlatt=%d/Lattice/Lattice(beta%.3f).txt",root_dir, N, b->val); //it modifies each time the name of the file to be created
	      	beta[s]=fopen(filename, "w");
	      	control_file(beta[s]);

	       	for(int j=0; j<N; j++){
	       	    for(int k=0; k<N; k++){
	       	        fprintf(beta[s], "%d ", field[j][k]);
	        	}
	        	fprintf(beta[s],"\n");
	      	}
	     	fclose(beta[s]);
	     	s++;
      	}
     	
        fclose(misure[p]);
     	printf("ho fatto il beta numero %d \n", p); 
		b=b->next;
		p++;
	}
	printf("Tutto a posto \n");
	fclose(f);
	fclose(g);
	return 0;
}
