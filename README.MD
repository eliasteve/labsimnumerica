# Esercizi del Laboratorio di Simulazione Numerica

Questa cartella contiene le seguenti sottocartelle:

* Cartelle con nome `LSN_XX`, che contengono lo svolgimento dell'esercitazione `XX`. In particolare queste cartelle contengono:
    * File chiamati `lsnXX.y.cpp`, con lo svolgimento dell'esercizio `XX.y`. Questi programmi si compilano con `make lsnXX.y` e si eseguono (una volta compilati) con `./lsnXX.y.out`. L'eccezione sono le cartelle per le esercitazioni 4, 6 e 7 che non contengono sorgente, visto che usano il codice in `NSL_SIMULATOR`, e che contengono invece delle cartelle con gli output delle simulazioni in vari casi, e le esercitazioni 11 e 12, che sono svolte interamente nei notebook Jupyter.
    * Eventualmente file del tipo `my_lib.cpp` e `my_lib.h`, con librerie di supporto per l'esecuzione dei programmi.
    * File del tipo `my_data.dat`, che contengono i dati di output dei programmi.
    * Il Makefile e 
    * il notebook Jupyter con l'analisi dei risultati, con il nome di `LSN_XX.ipynb`, tranne che nell'esercitazione 11, dove i notebook per lo svolgimento degli esercizi sono separati e hanno il nome `Exercise 11.x.ipynb`.

* La cartella `lib_NSL`, che contiene librerie di funzioni in C++ che vengono usate in più esercitazioni.

* La cartella `lib_NSL_py`, che contiene librerie in Python per il controllo del simulatore in `NSL_SIMULATOR` e per l'analisi dei dati nel caso del problema del commesso viaggiatore.

* La cartella `NSL_SIMULATOR`, che contiene il codice modificato per lo svolgimento delle esercitazioni 4, 6 e 7.

* La cartella `random`, che contiene il codice modificato del generatore di numeri casuali.