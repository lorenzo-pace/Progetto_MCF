#####################################################
#                                                   #
#                   Lorenzo Pace                    #
#      Università degli Studi di Perugia            #
#   Corso di Metodi Computazionali per la Fisica    #
#                                                   #
#---------------------------------------------------#
#                                                   #
#              Progetto finale                      #
#  Deriva di particelle in campi elettromagnetici   #
#                                                   #
#        Script con il main dell'esecuzione         #                       
#                                                   #
#####################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
import scipy
import sys,os
import argparse
import mod 

#Campi di default
E_def = [0.01, 0.01, 0] #V/m
B_def = [0, 0, 0.0001] #T
E_def2 = [-0.02, -0.05, 0] #V/m
B_def2 = [0, 0, -0.0003] #T
E_def3 = [10, 4, 0] #V/m

def parse_arguments():
    """
    Funzione usata per gestire gli argomenti di argparse. 
    """
    parser = argparse.ArgumentParser(description="Simulazione del moto di particelle in campi elettrici e magnetici.")
    
    parser.add_argument("-c", "--configurazione", choices=['exb', 'grad', 'both'], required=True,
                        help="Specifica la configurazione della simulazione, con scelte possibili: 'exb', per regime con E e con B uniforme; 'grad' per regime senza E e con B variabile, 'both' per regime con E e con B variabile.")

    parser.add_argument("-s", "--singola", action="store_true",
                        help="Avvia la simulazione con particella singola permettendo l'inserimento dei parametri.")

    parser.add_argument("-m", "--multipla", action="store_true",
                        help="Avvia la simulazione con un numero N di particelle permettendo l'inserimento dei parametri.")
    
    parser.add_argument("-t", "--traiettorie", action="store_true",
                        help="Mostra la traiettoria di un numero fornito di particelle con velocità iniziale casuale nei campi E, B di default.")
    
    parser.add_argument("-stat", "--statistica", action="store_true",
                        help="Studia statisticamente la velocità di deriva per un campione adeguato di particelle con velocità iniziale casuale e per diverse configurazioni di E e B.")
 
    args = parser.parse_args()

    #controlli

    if not args.configurazione:
        print("\nErrore: l'argomento -c/--configurazione è obbligatorio.\n")
        print("Utilizzare -h per vedere l'elenco delle opzioni disponibili.\n")
        sys.exit(1)

    return args

def funzioni(args):
    """
    Associa agli argomenti di argparse le funzioni corrispondenti del mod.
    -------------------------------------------
    Parametri:
    args: argomenti di argparse
    """

    #Configurazione E X B : B uniforme, E diverso da 0.
    
    if args.configurazione == 'exb':

        if args.singola :

            print('Inserire i parametri della simulazione.')
            Ex = float(input('Inserire la componente x del campo elettrico E [V/m]: '))
            Ey = float(input('Inserire la componente y del campo elettrico E [V/m]: '))
            Bz = float(input('Inserire la componente z del campo magnetico B [T]: '))
            E0 = [Ex, Ey, 0]
            B0 = [0, 0, Bz]
            print('Inserire il vettore posizione iniziale x0 [m]: ')
            x0 = mod.inserimento_vettori()
            print('Inserire il vettore velocità iniziale v0 [m/s]: ')
            v0 = mod.inserimento_coppie()
            particella = mod.inserimento_particella()
            passi = int(input('Inserire il numero di passi N della simulazione: '))
            tr, dt = mod.avanzamento(E0, B0, x0, v0, particella, passi, 1, '0')
            mod.grafico2d (tr, E0, B0, particella, 0, dt)

        if args.multipla :

            print('Inserire i parametri della simulazione.')
            Ex = float(input('Inserire la componente x del campo elettrico E [V/m]: '))
            Ey = float(input('Inserire la componente y del campo elettrico E [V/m]: '))
            Bz = float(input('Inserire la componente z del campo magnetico B [T]: '))
            E0 = [Ex, Ey, 0]
            B0 = [0, 0, Bz]
            N = int(input ('Inserire il numero di particelle con cui avviare la simulazione: '))
            v0= mod.velocita_montecarlo(N)
            x0 = mod.scelta_posizione(N)
            particella = mod.inserimento_particella()
            passi = int(input('Inserire il numero di passi della simulazione: '))
            tr, dt = mod.avanzamento(E0, B0, x0, v0, particella, passi, N, '0')
            mod.grafico2d (tr, E0, B0, particella, 0, dt)
        
        if args.traiettorie :

            print('Traiettorie di un numero fornito di particelle nei campi di default (60000 passi):')
            print('E = [1, 1, 0] *10^(-2) V/m')
            print('B = [0, 0, 1] *10^(-4) T')
            N = int(input('Fornire il numero di particelle da visualizzare:'))
            particella = mod.inserimento_particella()
            x0 = mod.posizioni_montecarlo(N)
            v0 = mod.velocita_montecarlo(N)
            tr, dt = mod.avanzamento(E_def, B_def, x0, v0, particella, 60000, N, '0')
            mod.grafico2d (tr, E_def, B_def, particella, 0, dt)

        if args.statistica : 
            
            N = 300 #campione adeguato di particelle
            print(f'Studio statistico della velocità di drift su {N} particelle nei campi di default.')
            particella = mod.inserimento_particella()
            v0 = mod.velocita_montecarlo(N)
            x0 = mod.scelta_posizione(N)
            tr1, dt1 = mod.avanzamento(E_def, B_def, x0, v0, particella, 60000, N, '0')
            #mod.grafico2d (tr1, E_def, B_def, particella, 0, dt1)
            teorica = mod.teorica_exb(E_def, B_def)
            mod.istogramma_vdrift(tr1, dt1, '1', '1', teorica, particella)
            while True:
                risposta = input('Digitare 1 per continuare con la configurazione successiva, altro per uscire: ')
                if (risposta=='1'):
                    tr2, dt2 = mod.avanzamento(E_def2, B_def, x0, v0, particella, 60000, N, '0')
                    #mod.grafico2d (tr2, E_def2, B_def, particella, 0, dt2)
                    teorica2 = mod.teorica_exb(E_def2, B_def)
                    mod.istogramma_vdrift(tr2, dt2, '2', '1', teorica2, particella)
                    risposta2 = input('Digitare 1 per continuare con la configurazione successiva, altro per uscire: ')
                    if (risposta2=='1'):
                        tr3, dt3 = mod.avanzamento(E_def, B_def2, x0, v0, particella, 60000, N, '0')
                        #mod.grafico2d (tr3, E_def, B_def2, particella, 0, dt3)
                        teorica3 = mod.teorica_exb(E_def, B_def2)
                        mod.istogramma_vdrift(tr3, dt3, '1', '2', teorica3, particella)
                        break
                    else:
                        break 
                else:
                    break          
       
    #Configurazione grad : B dipendente dalla posizione, E assente.

    if args.configurazione == 'grad':

        E0 = [0,0,0]
        B0 = [0,0,0] #variabili non usate

        if args.singola :
        
            print('Inserire i parametri della simulazione.')
            var = mod.scelta_gradiente()
            print('Inserire il vettore posizione iniziale x0 [m]: ')
            x0 = mod.inserimento_vettori()
            print('Inserire il vettore velocità iniziale v0 [m/s]: ')
            v0 = mod.inserimento_coppie()
            particella = mod.inserimento_particella()
            passi = int(input('Inserire il numero di passi N della simulazione: '))
            tr, dt = mod.avanzamento(E0, B0, x0, v0, particella, passi, 1, var)
            mod.grafico2d (tr, E0, B0, particella, var, dt)

        if args.multipla :

            print('Inserire i parametri della simulazione.')
            var = mod.scelta_gradiente()
            N = int(input ('Inserire il numero di particelle con cui avviare la simulazione: '))
            v0= mod.velocita_montecarlo(N)
            x0 = mod.scelta_posizione(N)
            particella = mod.inserimento_particella()
            passi = int(input('Inserire il numero di passi della simulazione: '))
            tr, dt = mod.avanzamento(E0, B0, x0, v0, particella, passi, N, var)
            mod.grafico2d (tr, E0, B0, particella, var, dt)            
        
        if args.traiettorie :

            print('Traiettorie di un numero fornito di particelle nei campi di default (60000 passi):')
            N = int(input('Fornire il numero di particelle da visualizzare:'))
            x0 = mod.posizioni_montecarlo(N)
            v0 = mod.velocita_montecarlo(N)
            particella = mod.inserimento_particella()
            var = mod.scelta_gradiente()
            tr, dt = mod.avanzamento(E0, B0, x0, v0, particella, 60000, N, var)
            mod.grafico2d (tr, E0, B0, particella, var, dt)

        if args.statistica : 
            
            N = 300 #campione adeguato di particelle
            print(f'Studio statistico della velocità di drift su {N} particelle nei campi di default.')
            particella = mod.inserimento_particella()
            v0 = mod.velocita_montecarlo(N)
            x0 = mod.scelta_posizione(N)
            tr1, dt1 = mod.avanzamento(E0, B0, x0, v0, particella, 60000, N, '1')
            #mod.grafico2d (tr1, E0, B0, particella, 1, dt1)
            teorica = mod.teorica_grad(tr1, particella, '1')
            mod.istogramma_vdrift(tr1, dt1, '3', '31', teorica, particella)
            while True:
                risposta = input('Digitare 1 per continuare con la configurazione successiva, altro per uscire: ')
                if (risposta=='1'):
                    tr2, dt2 = mod.avanzamento(E0, B0, x0, v0, particella, 60000, N, '2')
                    #mod.grafico2d (tr2, E0, B0, particella, 2, dt2)
                    teorica2 = mod.teorica_grad(tr2, particella, '2')
                    mod.istogramma_vdrift(tr2, dt2, '3', '32', teorica2, particella)
                    risposta2 = input('Digitare 1 per continuare con la configurazione successiva, altro per uscire: ')
                    if (risposta2 =='1'):
                        tr3, dt3 = mod.avanzamento(E0, B0, x0, v0, particella, 60000, N, '3')
                        #mod.grafico2d (tr3, E0, B0, particella, 3, dt3)
                        teorica3 = mod.teorica_grad(tr3, particella, '3')
                        mod.istogramma_vdrift(tr3, dt3, '3', '33', teorica3, particella)
                        break
                    else:
                        break
                else:
                    break      

    #Configurazione both : B dipendente dalla posizione, E diverso da 0.

    if args.configurazione == 'both':

        B0 = [0,0,0]

        if args.singola :
        
            print('Inserire i parametri della simulazione.')
            var = mod.scelta_gradiente()
            Ex = float(input('Inserire la componente x del campo elettrico E [V/m]: '))
            Ey = float(input('Inserire la componente y del campo elettrico E [V/m]: '))
            E0 = [Ex, Ey, 0]
            print('Inserire il vettore posizione iniziale x0 [m]: ')
            x0 = mod.inserimento_vettori()
            print('Inserire il vettore velocità iniziale v0 [m/s]: ')
            v0 = mod.inserimento_coppie()
            particella = mod.inserimento_particella()
            passi = int(input('Inserire il numero di passi N della simulazione: '))
            tr, dt = mod.avanzamento(E0, B0, x0, v0, particella, passi, 1, var)
            mod.grafico2d (tr, E0, B0, particella, var, dt)

        if args.multipla :

            print('Inserire i parametri della simulazione.')
            var = mod.scelta_gradiente()
            Ex = float(input('Inserire la componente x del campo elettrico E [V/m]: '))
            Ey = float(input('Inserire la componente y del campo elettrico E [V/m]: '))
            E0 = [Ex, Ey, 0]
            N = int(input ('Inserire il numero di particelle con cui avviare la simulazione: '))
            v0= mod.velocita_montecarlo(N)
            x0 = mod.scelta_posizione(N)
            particella = mod.inserimento_particella()
            passi = int(input('Inserire il numero di passi N della simulazione: '))
            tr, dt = mod.avanzamento(E0, B0, x0, v0, particella, passi, N, var)
            mod.grafico2d (tr, E0, B0, particella, var, dt)       
        
        if args.traiettorie :

            print('Traiettorie di un numero fornito di particelle nei campi di default (60000 passi):')
            N = int(input('Fornire il numero di particelle da visualizzare:'))
            x0 = mod.posizioni_montecarlo(N)
            v0 = mod.velocita_montecarlo(N)
            particella = mod.inserimento_particella()
            var = mod.scelta_gradiente()
            tr, dt = mod.avanzamento(E_def, B0, x0, v0, particella, 60000, N, var)
            mod.grafico2d (tr, E_def, B0, particella, var, dt)

        if args.statistica : 
            
            N = 300 #campione adeguato di particelle
            print(f'Studio statistico della velocità di drift su {N} particelle nei campi di default.')
            particella = mod.inserimento_particella()
            v0 = mod.velocita_montecarlo(N)
            x0 = mod.scelta_posizione(N)
            tr1, dt1 = mod.avanzamento(E_def, B0, x0, v0, particella, 60000, N, '1')
            #mod.grafico2d (tr1, E_def, B0, particella, '1', dt1)
            teorica = mod.teorica_both(tr1, particella, '1', E_def)
            mod.istogramma_vdrift(tr1, dt1, '1', '31', teorica, particella)
            while True:
                risposta = input('Digitare 1 per continuare con la configurazione successiva, altro per uscire: ')
                if (risposta=='1'):
                    tr2, dt2 = mod.avanzamento(E_def2, B0, x0, v0, particella, 60000, N, '2')
                    #mod.grafico2d (tr2, E_def2, B0, particella, '2', dt2)
                    teorica2 = mod.teorica_both(tr1, particella, '2', E_def2)
                    mod.istogramma_vdrift(tr1, dt1, '2', '32', teorica2, particella)
                    risposta2 = input('Digitare 1 per continuare con la configurazione successiva, altro per uscire: ')
                    if (risposta2 =='1'):
                        tr3, dt3 = mod.avanzamento(E_def, B0, x0, v0, particella, 60000, N, '3')
                        #mod.grafico2d (tr3, E_def, B0, particella, '3', dt3)
                        teorica3 = mod.teorica_both(tr1, particella, '3', E_def)
                        mod.istogramma_vdrift(tr3, dt3, '1', '33', teorica3, particella)   
                        break
                    else:
                        break 
                else:
                    break

# Main del programma

def main():
    args = parse_arguments()
    funzioni(args)

if __name__ == "__main__":
    main()


