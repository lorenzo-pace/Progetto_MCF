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
#           Modulo con le funzioni usate            #                       
#                                                   #
#####################################################

import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy
import sys,os
import argparse
from scipy.stats import norm
import matplotlib.colors as mcolors

#Classi e funzioni varie utili

class Particella():
    """
    Classe che definisce una particella di cui simulare il moto.
    Massa e carica sono espresse in unità del sistema internazionale. 
    -------------------------------------------
    Parametri:
    singolare : nome singolare
    plurale : nome plurale
    m : massa [kg]
    q : carica [C] 
    """
    
    def __init__(self, nomesing, nomeplur, massa, carica):
        self.singolare = nomesing
        self.plurale = nomeplur
        self.m = massa
        self.q = carica
    
    def __repr__(self):
        return 'Particella: {:}, massa {:} [kg], carica {:} [C]'.format(self.singolare, self.m, self.q)

def periodo_larmor (particella, B):
    """
    Calcola il periodo associato al compimento dell'orbita ciclotronica.
    -------------------------------------------
    Parametri:
    particella: elemento della classe omonima
    B : vettore campo magnetico
    """
    if not isinstance(particella, Particella):
        raise TypeError("Il primo argomento deve essere un oggetto della classe Particella.")

    intB = np.linalg.norm(B) #modulo del vettore campo magnetico
    T = 2*np.pi*particella.m/(intB * np.abs(particella.q))

    return T

#Funzioni di inserimento e scelta

def inserimento_vettori ():
    """
    Permette di inserire in input un vettore speficiandone le tre componenti cartesiane.
    """
    vec = np.zeros(3)
    coordinate = ['x', 'y', 'z']
    for j in range(3):
        vj = float(input(f"Componente {coordinate[j]}: "))
        vec[j] = vj
    return vec

def inserimento_particella ():
    """
    Permette di inserire in input una particella fra quelle preimpostate.
    - 'e': elettrone
    - 'ae': positrone
    - 'p': protone
    - 'ap': antiprotone
    """
    lista = {
        'e': Particella('elettrone', 'elettroni', 9.10938356e-31, -1.60217663e-19),
        'ae': Particella('positrone', 'positroni', 9.10938356e-31, 1.60217663e-19),
        'p': Particella('protone', 'protoni', 1.6726219e-27, 1.60217663e-19),
        'ap': Particella('antiprotone', 'antiprotoni', 1.6726219e-27, -1.60217663e-19)
    }
    scelta = ""
    while scelta not in lista:
        scelta = input("Inserire la specie della particella (e: elettrone, ae: positrone, p: protone, ap: antiprotone): ").lower()
        if scelta not in lista:
            print("Scelta non valida. Inserire una delle seguenti: e, ae, p, ap.")

    return lista[scelta]

def scelta_gradiente():
    """
    Permette di scegliere la legge di dipendenza spaziale del campo magnetico fra 3 opzioni.
    """    
    print('Scegliere la legge di dipendenza spaziale del campo magnetico selezionando fra 1, 2 e 3.')
    print('1. Bz = 0,00001 + 0,00005*x (T)')
    print('2. Bz = 0,01 + 0,02*x (T)')
    print('3. Bz = -0,0005 - 0,0008*x (T)')
    while True:
        risposta = input('Scelta: ')
        if (risposta=='1' or risposta == '2' or risposta == '3'):
            break
        else:
            print("Input non valido, scegliere fra 1, 2 e 3.")    
    return risposta

#Funzione di simulazione

def avanzamento(E0, B0, x0, v0, particella, passi, N, grad):
    """
    Definisce l'avanzamento in senso vettoriale come ciclo for di iterazioni singole
    di moto rettilineo uniforme. Ad ogni iterazione aggiorna la forza di Lorentz e 
    di conseguenza i vettori velocità e posizione.
    Regime both con campo B variabile con la posizione secondo la funzione B_grad e campo E diverso da 0. 
    Calcola il dt come 1/10000 del periodo dell'orbita ciclotronica per rendere il numero di passi
    inseriti in input indipendente da particella e campo magnetico fissati. 
    -------------------------------------------
    Parametri:
    - E0 : campo elettrico (array 3d) [V/m]
    - x0 : vettore posizione iniziale [m]
    - v0 : vettore velocità (array 3d) [m/s]
    - particella: elemento della classe omonima definita
    - passi: numero di iterazioni
    - N : numero di particelle (con valore di default 1)
    - gradiente: char che identifica il tipo di gradiente (0 assente, 1(2,3) per le funzioni B_grad_1 (2,3))
    -------------------------------------------
    Restituisce:
    - x : array che rappresenta la traiettoria
    - dt : passo temporale della simulazione
    """

    #Caso 1: simulazione singola
    if N == 1: 
        x = np.zeros((passi, 3))
        v = np.zeros((passi, 3))
        x[0, :] = x0
        v[0, :] = v0
        for j in range(passi - 1):
            if grad=="0":
                B = B0
            else:
                B = gradiente(grad,x[j, 0]) #dipendenza di B dalla posizione
            dt = 0.0001 * periodo_larmor(particella, B)
            f_lorentz = particella.q * (E0 + np.cross(v[j, :], B))  #forza di Lorentz
            v[j + 1, :] = v[j, :] + f_lorentz / particella.m * dt  #aggiornamento del vettore velocità
            x[j + 1, :] = v[j, :] * dt + x[j, :]  #passi di moto rettilineo uniforme

     #Caso 2: simulazione multipla
    elif N > 1:  
        x = np.zeros((passi, N, 3))
        v = np.zeros((passi, N, 3))
        for i, v_init in enumerate(v0):
            x[0, i, :] = x0[i]
            v[0, i, :] = v_init
            for j in range(passi - 1):
                if grad=="0": #in assenza di gradiente
                    B = B0
                else: #presenza di gradiente
                    B = gradiente(grad, x[j, i, 0]) #dipendenza di B dalla posizione (particella i-esima)
                dt = 0.0001 * periodo_larmor(particella, B)
                f_lorentz = particella.q * (E0 + np.cross(v[j, i, :], B))  #forza di Lorentz (particella i-esima)
                v[j + 1, i, :] = v[j, i, :] + f_lorentz / particella.m * dt  #aggiornamento del vettore velocità (particella i-esima)
                x[j + 1, i, :] = v[j, i, :] * dt + x[j, i, :]  #passi di moto rettilineo uniforme (particella i-esima)

    return x, dt

#Rappresentazioni grafiche 

def grafico2d (tr, E, B0, particella, grad, dt):
    """
    Disegna l'orbita in uno spazio cartesiano 2d mostrando esplicitamente anche
    direzione e verso dei vettori campo elettrico e magnetico. 
    L'intensità del campo magnetico lungo la direzione di variazione per gradiente è resa da una scala di colori. 
    -------------------------------------------
    Parametri:
    - tr  :  traiettoria
    - E : campo elettrico che giace su xy
    - B0 : campo magnetico che giace su z
    - particella: elemento della classe omonima
    - grad : char che identifica il gradiente
    - dt : restituito da avanzamento per l'eventuale rappresentazione di v
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    #Caso 1: simulazione singola (1 traiettoria, l'array tr ha dimensione 2)
    if len(tr.shape) == 2:
            x0, y0 = tr[0, 0], tr[0, 1]
            ax.plot(tr[:, 0], tr[:, 1], label="Traiettoria", color='royalblue')
            ax.scatter(x0, y0, color='blue', label='Posizione iniziale')
            ax.set_title(f"Traiettoria di un {particella.singolare} nel piano xy")
            ax.legend()

    #Caso 2: simulazione multipla (N traiettorie, l'array tr ha dimensione 3)
    else:
        N = tr.shape[1] 
        for j in range(N):
            x0, y0 = tr[0, j, 0], tr[0, j, 1]
            ax.plot(tr[:, j, 0], tr[:, j, 1], color=plt.cm.viridis(j/N))  #per avere colori diversi per ogni particella
            ax.scatter(x0, y0, color='blue', marker='.') #label tolte per non appesantire il grafico
        ax.set_title(f"Traiettorie per {N} {particella.plurale} nel piano xy")

    xmin, xmax = ax.get_xlim() 
    ymin, ymax = ax.get_ylim()

    #rappresentazione del campo B
    a, b = np.linspace(xmin,xmax,7), np.linspace(ymin,ymax,7)
    a, b = np.meshgrid(a,b)
    Bx_pos = a[0, 1]  #posizionamento della lettera B accanto al secondo segno della fila alta
    By_pos = b[-1, 1]
    ax.text(Bx_pos -0.06 * (xmax - xmin), By_pos - 0.05 * (ymax - ymin), 'B', color='black', fontsize=12, fontweight='bold')

    if (grad == 0): #assenza di gradiente magnetico
            if B0[2]>0 :
                ax.scatter(a, b, marker='o', edgecolors='gray', s=150, c='white')
                ax.scatter(a, b, marker='.', c='gray')
            else:
                ax.scatter(a, b, marker='o', edgecolors='gray', s=150, c='white')
                ax.scatter(a, b, marker='x', c='gray')

    else : #caso con gradiente magnetico
        B_grad = gradiente(grad, a)
        B_abs = np.abs(B_grad[2])
        colori = [(0.9, 0.9, 0.9), (0, 0, 0)]  #da grigio chiaro a nero
        mappa = mcolors.LinearSegmentedColormap.from_list("custom_gray", colori) 
        #grigiochiaro= np.min(B_abs) + 0.1* (np.max(B_abs) - np.min(B_abs)) #partenza dal grigio invece che dal bianco
        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                #Br_value = B_grad(a[i, j])[2]
                Br_value = gradiente(grad, a[i, j])[2]
                intensita = (np.abs(Br_value) - np.min(B_abs)) / (np.max(B_abs) - np.min(B_abs))
                intensita = np.clip(intensita, 0, 1)  #normalizzazione intensità colore
                colore = mappa(intensita) #colormap in scala di grigi
                if Br_value >= 0:
                    ax.scatter(a[i, j], b[i, j], c=[colore], marker='o', s=100)
                    ax.scatter(a[i, j], b[i, j], c='white', marker='o', s=10)  #puntino cerchiato
                else: 
                    ax.scatter(a[i, j], b[i, j], c=[colore], marker='o', s=100)
                    ax.scatter(a[i, j], b[i, j], c='white', marker='x', s=60, linewidths=1.5)  #croce cerchiata

        #colorbar della scala di colori di B
        normalizzazione = mcolors.Normalize(vmin=np.min(B_abs), vmax=np.max(B_abs)) #per fare la colorbar normalizzata
        sm = plt.cm.ScalarMappable(cmap=mappa, norm=normalizzazione)  #mappatura dei valori di B
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label('Intensità del campo B [T]') 
    
    #rappresentazione del campo E
    if np.linalg.norm(E) != 0 :
        Ex_pos = xmin + 0.03 * (xmax - xmin)
        Ey_pos = ymin + 0.05 * (ymax - ymin)
        arrow_length_x = 0.2 * (xmax - xmin)  #lunghezza fissata della freccia
        arrow_length_y = 0.2 * (ymax - ymin) 
        E_norm = E / np.linalg.norm(E) #normalizzazione del campo elettrico cosicché la freccia sia lunga indipendentemente da E
        ax.quiver(Ex_pos, Ey_pos, E_norm[0] * arrow_length_x, E_norm[1] * arrow_length_y,
                angles='xy', scale=1, scale_units='xy', width=0.005, color='red', pivot='middle')
        ax.text(Ex_pos + 0.1 * (xmax - xmin), Ey_pos + 0.1 * (ymax - ymin), 'E', color='red', fontsize=12, fontweight='bold')

    #rappresentazione della velocità netta (solo particella singola)
    if len(tr.shape) == 2:
        vdrift = vel_drift(tr,dt)
        if np.linalg.norm(vdrift) != 0 :
            vx_pos = xmin + 0.05 * (xmax - xmin)
            vy_pos = ymin + 0.05 * (ymax - ymin)
            arrow_length_x = 0.2 * (xmax - xmin)  #lunghezza fissata della freccia
            arrow_length_y = 0.2 * (ymax - ymin) 
            v_norm = vdrift / np.linalg.norm(vdrift) #normalizzazione della velocità cosicché la freccia sia lunga indipendentemente
            ax.quiver(vx_pos, vy_pos, v_norm[0] * arrow_length_x, v_norm[1] * arrow_length_y,
                    angles='xy', scale=1, scale_units='xy', width=0.005, color='green', pivot='middle')
            ax.text(vx_pos + 0.1 * (xmax - xmin), vy_pos + 0.04 * (ymax - ymin), 'v', color='green', fontsize=12, fontweight='bold')

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    plt.show()

#Gradienti magnetici

def gradiente (tipo, x):
    """
    Funzione necessaria per distinguere i 3 tipi di gradiente magnetico tramite un dizionario. 
    -------------------------------------------
    Parametri:
    x : coordinata su cui varia l'intensità del vettore campo magnetico
    """
    opzioni = {
        "1": B_grad_1,
        "2": B_grad_2,
        "3": B_grad_3
    }
    if tipo in opzioni:
        return opzioni[tipo](x) 
    else:
        raise ValueError(f"Tipo di gradiente non valido: {tipo}. Inserire uno fra tra '1', '2', '3'.")

def B_grad_1 (x) : 
    """
    Legge di variazione spaziale del campo magnetico. 
    -------------------------------------------
    Parametri:
    x : coordinata su cui varia l'intensità del vettore campo magnetico
    """
    B0 = 0.00001 #T
    q = 0.00005
    Br = B0+q*(x**1)

    return [0,0,Br]

def B_grad_2 (x) :
    """
    Legge di variazione spaziale del campo magnetico. 
    -------------------------------------------
    Parametri:
    x : coordinata su cui varia l'intensità del vettore campo magnetico
    """
    B0 = 0.01 #T
    q = 0.02
    Br = B0+q*(x**1)

    return [0,0,Br]

def B_grad_3 (x) : 
    """
    Legge di variazione spaziale del campo magnetico. 
    -------------------------------------------
    Parametri:
    x : coordinata su cui varia l'intensità del vettore campo magnetico
    """
    B0 = -0.0005 #T
    q = -0.0008
    Br = B0+q*(x**1)

    return [0,0,Br]

#Generazione di parametri iniziali casuali

def velocita_montecarlo(N):
    """
    Genera N velocità iniziali casuali per particelle. Ciascuna componente è compresa in un range
    fra 10^3 m/s e +10^5 m/s con probabilità uniforme. La limitazione permette di mantenersi nel
    regime non relativistico ed è adatta ai valori scelti nel modulo per le intensità di B ed E.
    La simulazione genera metà velocità positive e metà negative ma con le precedenti limitazioni in modulo.
    -------------------------------------------
    Parametri:
    - N : numero dei vettori velocità inziale da generare
    """
    numero_positive = N // 2 + (N % 2)  #perché N potrebbe essere dispari
    numero_negative = N // 2 

    vxpos = np.random.uniform(1e3, 1e5, numero_positive)
    vxneg = np.random.uniform(-1e5, -1e3, numero_negative)
    vx = np.concatenate((vxpos, vxneg))
    vypos = np.random.uniform(1e3, 1e5, numero_positive)
    vyneg = np.random.uniform(-1e5, -1e3, numero_negative)
    vy = np.concatenate((vypos, vyneg))
    vzpos = np.random.uniform(1e3, 1e5, numero_positive)
    vzneg = np.random.uniform(-1e5, -1e3, numero_negative)
    vz = np.concatenate((vzpos, vzneg))
    array_v = np.column_stack((vx, vy, vz))

    return array_v

def posizioni_montecarlo(N):
    """
    Genera N posizioni iniziali casuali per particelle. Ciascuna componente è compresa in un range
    fra -1.5 m e +1.5 m con probabilità uniforme. La limitazione permette di mantenersi entro una
    risoluzione spaziale ragionevole ed è adatta ai valori scelti nel modulo per le intensità di B ed E.
    -------------------------------------------
    Parametri:
    - N : numero dei vettori velocità inziale da generare
    """
    px = np.random.uniform(-1.5, 1.5, N)
    py = np.random.uniform(-1.5, 1.5, N)
    pz = np.random.uniform(-1.5, 1.5, N)
    array_p = np.column_stack((px, py, pz))

    return array_p

def scelta_posizione(N):
    """
    Permette di scegliere nelle simulazioni multiple se generare casualmente le posizioni iniziali o
    usarne una comune, inserita, per tutte le particelle. 
    -------------------------------------------
    Parametri:
    - N : numero delle particelle presenti nella simulazione
    """
    print('Premere A per fissare una sola posizione iniziale;')
    print('Premere B per generare N posizioni iniziali casuali;')
    while True:
        risp2 = input('Risposta: ')
        if (risp2=='A'):
            print('Inserire il vettore posizione iniziale comune x0 [m]: ')            
            pos0 = inserimento_vettori()
            x0 = np.tile(pos0, (N, 1))
            break
        elif (risp2 == 'B'):
            x0= posizioni_montecarlo(N)
            break
        else:
            print("Input non valido, inserire A o B.")

    return x0

#Funzioni di statistica delle velocità di drift

def vel_drift (tr, dt):
    """
    Identifica la velocità netta di drift come differenza fra i vettori posizione
    iniziale e finale diviso il tempo netto dall'inizio alla fine del moto. 
    -------------------------------------------
    Parametri:
    - tr : array tridimensionale dei vettori posizione (traiettoria)
    -------------------------------------------
    Resistuisce:
    - v : vettore velocità di drift
    """
    passi = tr.shape[0]
    delta_t = np.sum(dt)
    delta_x = tr[-1]-tr[0]
    v = delta_x/delta_t

    return v

# Grafico di distribuzione delle velocità

def istogramma_vdrift(tr, dt, E, B, vel_teo):
    """
    Rappresenta un istogramma che mostra la distribuzione delle velocità di drift tenendo conto della loro
    media e deviazione standard.
    -------------------------------------------
    Parametri:
    - tr : array delle traiettorie restituite dalla simulazione con N particelle
    - dt : passo temporale della simulazione
    - E : char che identifica un tipo di campo elettrico per il sottotitolo
    - B : char che identifica un tipo di campo magnetico per il sottotitolo
    - vel_teo : array velocità teorica [m/s]
    """
    plt.figure(figsize=(10, 6))
    vel = vel_drift(tr, dt)
    velocita = []
    for vj in vel:
        velocita = np.append(velocita, np.linalg.norm(vj))
    n, bins, p = plt.hist(velocita, bins='auto', density=True, alpha=0.6, color='royalblue')
    mu = np.mean(velocita)
    plt.axvline(mu, color='red', linewidth=2, label=f'Velocità media: {mu:.2e} m/s')
    sigma = np.std(velocita)
    x_vals = np.linspace(min(bins), max(bins), 300)
    plt.plot(x_vals, norm.pdf(x_vals, mu, sigma), 'r-', lw=2, label='Distribuzione gaussiana')

    #Confronto grafico con la velocità teorica
    teorica = np.linalg.norm(vel_teo)
    plt.axvline(teorica, color='green', linewidth=2, label=f'Velocità teorica: {teorica:.2e} m/s')
    
    descrizioni_E = {
        "1": "E = (1; 1; 0) × 10⁻² V/m",
        "2": "E = (-2; -5; 0) × 10⁻² V/m",
        "3": "E = 0"
    }
    descrizioni_B = {
        "1": "B = (0; 0; 1) × 10⁻⁴ T",
        "2": "B = (0; 0; -3) × 10⁻⁴ T",
        "31": "B (x) = 10⁻⁵ + 5 × 10⁻⁵ x (T)",
        "32": "B (x) = 10⁻² + 2 × 10⁻² x (T)",
        "33": "B (x) = - 5 × 10⁻⁴ - 8 × 10⁻⁴ x (T)"
    }
    sottotitolo = f"{descrizioni_E[E]}; {descrizioni_B[B]}"
    plt.suptitle(sottotitolo, fontsize=12, color="gray")

    plt.xlabel('Velocità di deriva (m/s)')
    plt.ylabel('Densità di probabilità')
    plt.title('Distribuzione della velocità di deriva')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.show()

# Velocità di deriva teoriche

def teorica_exb(E0, B0):
    """
    Calcola la velocità del drift E x B.
    -------------------------------------------
    Parametri:
    - E0 : array campo elettrico [V/m]
    - B0 : array campo magnetico [T]
    """
    intB = np.linalg.norm(B0) ** 2
    velteorica = np.cross(E0, B0) / intB
    return velteorica

def teorica_grad (tr, particella, grad):
    """
    Calcola la velocità del drift di gradiente magnetico ortogonale.
    -------------------------------------------
    Parametri:
    - tr : array traiettoria
    - particella : elemento della classe omonima
    - grad : char che identifica il tipo di gradiente
    """

    posizione_mediana = tr.shape[0]//2
    #print(posizione_mediana)
    x_media = tr[posizione_mediana, :, :]
    x_precedente = tr[posizione_mediana-1, :, :]
    velocita = []
    for xi in x_media:
        B_medio = gradiente(grad, xi[0])[2]
        dt = 0.0001*periodo_larmor(particella, B_medio)
        velocita.append((x_media - x_precedente) / dt)
    velocita = np.array(velocita) #ritrasformo in array perché prima non funzionava veniva letto come lista

    #Derivate prime che compaiono nella formula della velocità
    if grad == "1":
        k = 0.00005 
    elif grad == "2":
        k = 0.02
    elif grad == "3":
        k = -0.0008
        
    # Stima della velocità perpendicolare media dalla simulazione
    v_ortogonale_media =  np.mean(np.linalg.norm(velocita[:, :2], axis = 1))

    y = (particella.m * v_ortogonale_media ** 2) / (2 * particella.q * B_medio ** 2)  

    return [0, y * k, 0]

def teorica_both (tr, particella, grad, E0):
    """
    Calcola la velocità del drift generale dove è presente sia campo elettrico uniforme
    che una legge di dipendenza spaziale ortogonale del campo magnetico. In particolare la
    velocità netta di deriva è la somma delle due velocità di drift. 
    -------------------------------------------
    Parametri:
    - tr : array traiettoria
    - dt : array dei tempi per il calcolo della v
    - particella : elemento della classe omonima
    - grad : char che identifica il tipo di gradiente
    - E0 : array campo elettrico [V/m]
    """
    posizione_mediana = tr.shape[0]//2
    x_media = tr[posizione_mediana, :, :]
    Bi = []
    for xi in x_media:
        Bi.append(gradiente(grad, xi[0]))
    B_medio = [0,0,np.mean(Bi[2])]
    v_soloexb = teorica_exb(E0, B_medio)
    v_solograd = teorica_grad(tr, particella, grad)

    return v_soloexb + v_solograd