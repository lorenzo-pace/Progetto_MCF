# Drift di particelle in campi elettromagnetici

La deriva di particelle cariche in regioni dello spazio permeate da campi magnetici è l'effetto dell'azione di campi elettrici o variazioni temporali o spaziali dei campi stessi, che si traduce nella presenza di una velocità netta di drift per il centro di guida dell'orbita ciclotronica che normalmente la particella compierebbe nel solo campo magnetico uniforme e costante. 

Nel presente programma si è tentato di simulare il moto compiuto da particelle in tre casi:
- sotto l'azione di un campo elettrico uniforme e costante e di un campo magnetico uniforme e costante (drift $E\times B$);
- sotto l'azione di un campo magnetico con una certa dipendenza spaziale (drift di gradiente magnetico $\nabla B$);
- nella situazione congiunta di campo elettrico e gradiente magnetico. 

L'obiettivo è stato da una parte l'implementazione della simulazione con gli adeguati strumenti e la sua rappresentazione grafica; in secondo luogo, lo studio statistico delle velocità di deriva nette appositamente calcolate. 

## Introduzione teorica

Nella fisica dei plasmi di particelle l'approccio cinetico alla descrizione del moto è di solito complicato perché necessita la risoluzione di un gran numero di equazioni del moto e di equazioni di Maxwell. Quando tuttavia si tratti una particella singola, sotto l'ipotesi che la particella puntiforme sia tale da non modificare col suo movimento i campi in cui è immersa, l'approccio cinetico è risolvibile e si parla di *teoria delle orbite*. Il moto di una particella di carica $q$ in una regione dello spazio permeata da campi elettromagnetici è regolato dalla forza di Lorentz:

$$\vec{f} = q\  (\vec{v}\times \vec{B} + \vec{E})$$

Una conseguenza immediata della forma della forza di Lorentz è il cosiddetto moto ciclotronico, ossia l'orbita circolare compiuta da particelle cariche immerse in un campo magnetico $\vec{B}$. Infatti, l'equazione del moto restituisce delle componenti della velocità che sono oscillanti con una certa frequenza detta frequenza di ciclotrone:

$$\omega = \frac{|q|B}{m}$$

Al moto ciclotronico è associato un raggio dell'orbita detto di Larmor e un periodo pari rispettivamente a:

$$r_L = \frac{mv_{\perp}}{|q|B} \qquad T = \frac{2\pi m}{|q|B}$$

Se $v_{\perp}$ è la componente della velocità trasversa alla direzione di $\vec{B}$. Nel moto ciclotronico semplice, il centro di guida, ossia il centro dell'orbita ciclotronica, ha velocità nulla. Ogni qualvolta il raggio di Larmor vada incontro a delle variazioni, il centro di guida può avere una velocità diversa da 0 che viene chiamata velocità di drift. Ciò può essere dovuto a campi non uniformi nello spazio (drift di gradiente magnetico parallelo o ortogonale, campi elettrici variabili nello spazio), a campi variabili nel tempo (drift di polarizzazione), alla presenza di forze che agiscono sulla particella. 

Nel caso di drift $\vec{E}\times \vec {B}$ la velocità di deriva è definita come la componente della velocità del centro di guida ortogonale alla direzione del campo magnetico. La velocità del drift $\vec{E} \times \vec{B}$ è costante purché i campi siano uniformi e costanti nel tempo. Sotto queste ipotesi, moltiplicando vettorialmente $\times \vec{B}$ l'espressione della forza di Lorentz si ottiene una velocità teorica:

$$\vec{v}_{E\times B} = \frac{\vec{E} \times \vec{B}}{B^2}$$

Nel caso di gradiente magnetico ortogonale, cioè quando esista un $\nabla B \perp \vec{B}$ ricavabile dalla legge di dipendenza spaziale del campo magnetico su una coordinata ortogonale alla sua direzione, la velocità si ottiene da:

$$\vec{v}_{\nabla B} = \frac{mv_{\perp}^2}{2qB^2} \cdot \frac{\vec{B} \times \nabla B}{B}$$

Se, come sarà trattato nel seguito, la coordinata su cui varia il campo magnetico è la $x$ e quella su cui giace è la $z$, allora la formula si semplifica a:

$$\vec{v}_{\nabla B} = \left( 0; \ \frac{mv_{\perp}^2}{2qB^2} \cdot \frac{dB_z}{dx}; \ 0 \right)$$

Nel caso, infine, in cui siano presenti sia $\vec{E} \neq 0$ che un gradiente magnetico la velocità netta di drift sarà semplicemente la somma delle due (vettorialmente). 

## File presenti

Di seguito il contenuto dei due file allegati, eccetto il presente README.

- ***mod.py*** è il modulo all'interno del quale sono presenti tutte le funzioni implementate, raggruppate per ambito. 
- ***progetto.py*** è lo script che contiene il main del programma, organizzato tramite argparse() per guidare l'esecuzione. 

## Esecuzione del programma

L'esecuzione del programma è gestita dal modulo argparse(). Per conoscere gli argomenti di ArgParse da utilizzare per il codice, si consiglia di eseguire il comando di `help`. In generale, l'avvio di una simulazione è subordinato alla scelta fra i tre regimi possibili di simulazione:

- -c exb : campo magnetico uniforme, campo elettrico diverso da 0;
- -c grad: gradiente magnetico, campo elettrico assente;
- -c both: gradiente magnetico, campo elettrico diverso da 0.

La scelta del regime è obbligatoria. Una volta digitato, si può scegliere di avviare 4 tipi di simulazione. 

- -s (--singola) : permette di simulare il moto di una singola particella inserendo i parametri relativi ai campi, il tipo di particella, il numero di passi, la velocità e la posizione iniziale;
- -m (--multipla) : permette di simulare il moto di un numero N, dato in input, di particelle inserendo i parametri relativi ai campi, il tipo di particella, il numero di passi, eventualmente la posizione iniziale con la possibilità di generarla casualmente come avviene per le velocità iniziali;
- -t (--traiettorie) : permette di simulare la traiettoria di un numero fornito di particelle con velocità iniziale casuale nei campi E, B di default;
- -stat (--statistica) : permette di studiare statisticamente la velocità di deriva producendo un istogramma delle velocità nette di drift e confrontando i risultati con la velocità teorica per quel tipo di deriva.

Esempio di avvio di una simulazione:

```bash
python3 progetto.py -c exb -stat
python3 progetto.py -c both -m
python3 progetto.py -c grad -s
```

### Simulazioni singole e multiple

Le istruzioni -s e -m permettono di avviare una simulazione completamente personalizzabile per il moto di una o N (numero in input) particelle nei campi scelti. Il tutto è subordinato alla struttura fisica del problema, con il campo elettrico che giace sul piano trasverso (xy) alla direzione (z) del campo magnetico. Interagendo con il terminale l'utente è in grado di inserire le componenti dei campi e, nel caso di gradiente magnetico, il tipo di dipendenza spaziale fra tre predefinite. Si può altresì scegliere fra 4 diverse particelle per dare avvio alla simulazione (protone, antiprotone, elettrone, positrone). Alla fine del calcolo delle traiettorie viene generato sullo schermo un grafico bidimensionale che rappresenta la traiettoria della particella o delle N particelle nel piano xy. 

### Visualizzazione di traiettorie

L'istruzione -t serve a  mostrare la traiettoria di alcune particelle con velocità iniziale casuale (in un intervallo appropriato) con un numero finito di passi all'interno dei campi di default. Può essere utile per verificare velocemente alcune proprietà senza inserire troppi parametri da terminale.

### Statistica delle velocità di drift

L'istruzione -stat permette di avviare tre simulazioni con un campione importante di particelle (N=100). In questo caso, nei campi di default, vengono calcolate le traiettorie delle N particelle per poi estrarne la velocità netta di drift; queste N velocità sono poi rappresentate in un istogramma dove si può apprezzare il confronto con una loro distribuzione normale e con la velocità teorica calcolata secondo le formule precedentemente descritte.

## Struttura del codice

In questa sezione si elencano brevemente le principali scelte di programmazione effettuate nel corso della costruzione del codice. Nel modulo ***mod.py*** sono riportate tutte le funzioni che poi verranno chiamate nel main del programma. La struttura di ***progetto.py*** è stata già più o meno delineata nel paragrafo precedente e comunque segue gli schemi convenzionali del modulo Argparse. 

### Classe particella

La classe `Particella()` definisce oggetti che contengono nome, massa e carica che identificano un protone, un positrone e così via. La sua utilità sta nella possibilità di richiamare proprietà diverse della particella nelle funzioni includendo un solo argomento. 

### Funzioni di inserimento

Al fine di facilitare l'avvio delle simulazioni sono state implementate alcune funzioni di inserimento controllate da cicli `while`. Nella fattispecie si tratta delle funzioni `inserimento_vettori()`, `inserimento_particella()` e `scelta_gradiente()`, che restituiscono rispettivamente un vettore tridimensionale (array di float), una particella fra le predefinite, con grandezze in unità S.I., a seconda della sigla (e, ae, p, ap) inserita, e infine un char che identifica uno fra i tre tipi di gradiente predefiniti. Per quanto riguarda i campi personalizzabili, si è preferito controllare l'inserimento direttamente da ***progetto.py*** per economia di funzioni da definire. 

### Simulazione

La simulazione è mediata da una funzione essenziale che è `avanzamento(E0, B0, x0, v0, particella, passi, N, grad)`. Questa funzione ha per argomenti, nell'ordine, il vettore campo elettrico e il vettore campo magnetico, le posizioni iniziali, le velocità iniziali, il tipo di particella, il numero di passi, il numero di particelle e il char che identifica il tipo di gradiente. Nelle configurazioni `-c grad` e `-c both` la variabile `B0` è inizializzata nel main al vettore nullo, così come `E0` in `-c grad`. La funzione è ottimizzata con un ciclo `if` per funzionare sia con 1 che con N particelle della stessa specie, con conseguente `shape` degli argomenti `x0` e `v0`. Il processo di base è il riempimento di una matrice vuota (la traiettoria) su un ciclo `for` dimensionato come i `passi`. Per ciascun ciclo, il codice esegue nell'ordine:

- l'eventuale aggiornamento del campo magnetico secondo la posizione (solo se `grad` è diverso da '0');
- il calcolo del dt (si veda il paragrafo successivo);
- l'aggiornamento della forza di Lorentz;
- l'aggiornamento del vettore velocità secondo la legge oraria del moto uniformemente accelerato;
- l'aggiornamento del vettore posizione secondo la legge oraria del moto rettilineo uniforme. 

In tal modo alla fine della simulazione, che può essere anche lunga in termini di tempo di esecuzione (soprattutto in presenza di gradiente magnetico), la funzione è in grado di restituire la matrice traiettoria e la matrice dei `dt` che servono poi per il calcolo delle velocità nette. Il tempo di esecuzione è collegato alla variabile `dt` che viene calcolata all'interno di `avanzamento(...)` tramite una funzione secondaria denominata `periodo_larmor(particella, B)` che calcola il periodo associato al compimento dell'orbita ciclotronica secondo la formula già vista. Il `dt` usato in `avanzamento(...)` è $1/10^4$ di questo periodo. Un valore così piccolo aumenta inevitabilmente il tempo di esecuzione della simulazione, ma per come è concepita la simulazione stessa un intervallo di tempo minore produrrebbe risultati non fisici. Infatti, soprattutto nella configurazione di gradiente magnetico, per avere una risoluzione corretta in termini di aggiornamento del campo è necessario un intervallo di tempo più piccolo possibile e un `dt` che è una parte su diecimila del periodo ciclotronico è sembrato, dopo vari tentativi, quello più efficace per conciliare tempo di esecuzione e risultati corretti. L'alternativa a una risoluzione temporale piccola sarebbe l'introduzione di metodi più elaborati: di solito nelle simulazioni di plasma si usa un algoritmo numerico chiamato Boris Push, un metodo di integrazione esplicito a tempo diviso che aggiorna la velocità e la posizione delle particelle in più fasi, garantendo una conservazione accurata dell'energia e del momento angolare in campi magnetici puri. Nei limiti dei tempi di esecuzione e dell'esperimento operato si può comunque ritenere valida la scelta del `dt` come frazione infinitesima dell'orbita ciclotronica. 

### Grafico

La funzione più lunga in ***mod.py*** è `grafico2d (tr, E, B0, particella, grad, dt)`. È programmata per funzionare per tutte le possibili configurazioni e riporta sul piano xy la traiettoria delle particelle simulate. Il campo elettrico è rappresentato da una freccia rossa in basso a sinistra (quando presente) nella sua direzionalità, e per particelle singole si riporta con una freccia verde anche il vettore velocità di deriva. Il campo magnetico giace sulla dimensione trasversa e il suo verso è identificato da una croce o un punto a seconda che sia entrante o uscente (cioè che sia negativa o positiva la coordinata z): nel caso di gradiente magnetico l'andamento sulla coordinata x è rappresentato da un'adeguata scala di colori. Le dimensioni delle frecce sono normalizzate e anche la scala di colori è normalizzata, rispetto a quelle di default, per partire da un grigio chiaro invece che dal bianco e rendere in ogni caso visibile il segno del campo magnetico. 

### Parametri casuali

Per generare posizioni iniziali e velocità iniziali con metodo Montecarlo si è fatto uso della `np.rand.uniform` che genera numeri reali (float) con probabilità uniforme entro dei limiti assegnati. Le due funzioni `posizioni_montecarlo(N)` e `velocita_montecarlo(N)` restituiscono matrici Nx3 di vettori da usare nelle simulazioni. Le limitazioni sono state scelte in accordo sia alla fisica del problema che alla natura del codice. In particolare a livello di studio statistico la dispersione delle posizioni è irrilevante, ma a livello grafico può essere fastidiosa e limitare fortemente la risoluzione soprattutto per elettroni. Un range fra $\pm$ 1,5 m permette di evitare uno sparpagliamento eccessivo. A livello di velocità, invece, i limiti assumono anche rilevanza statistica. Ciascuna componente ha modulo compreso in un range fra $10^3$ m/s e $10^5$ m/s con probabilità uniforme. La limitazione permette di mantenersi nel regime non relativistico ed è adatta ai valori scelti nel modulo per le intensità di B ed E. La funzione `scelta_posizione()` serve semplicemente, nelle simulazioni multiple, a far decidere all'utente se impostare un'unica posizione iniziale per tutte le particelle o a generarle casualmente con la funzione già menzionata. 

### Velocità di drift teoriche

Le tre funzioni `teorica_exb(E0, B0)`, `teorica_grad (tr, v0, particella, grad)` e `teorica_both (tr, v0, particella, grad, E0)` servono a calcolare la velocità di deriva teorica secondo le formule presentate. Nella funzione per la configurazione con gradiente il char `grad` corrisponde a tre valori di derivate prime che compaiono nella formula teorica. La velocità ortogonale media è calcolata con le giuste componenti di una media delle velocità dei singoli passi; il campo $B$ rappresentativo è stimato con una posizione media lungo la traiettoria. Questa stima però è abbastanza soggetta a fluttuazioni e può fallire in alcuni casi numericamente. Si è riscontrato infatti una maggiore aderenza fra velocità di drift media e teorica nella configurazione `-c exb`. 

### Funzioni di statistica delle velocità

La velocità di drift può stimarsi come differenza fra i vettori posizione iniziale e finale diviso il tempo netto dall'inizio alla fine del moto, e in questo senso è implementata nella funzione `vel_drift(tr, dt)`. A quel punto la `istogramma_vdrift(tr, dt, E, B, vel_teo)` riceve in input i parametri e visualizza la distribuzione delle N velocità della modalità statistica in un istogramma. All'interno del grafico sono presenti due barre verticali di colori diversi che identificano la posizione sull'asse delle velocità dei valori teorico e medio del drift. Sovrapposta all'istogramma c'è anche la curva gaussiana di distribuzione calcolata mediante media e deviazione standard della distribuzione. Il senso della rappresentazione statistica è valutare, quando sia apprezzabile, la normalità della distribuzione della velocità e il confronto del valor medio con quello teorico. Il valore del numero N di particelle simulate è stato scelto come miglior compromesso fra il tempo di esecuzione del programma e la rilevanza statistica della simulazione. 

## Conclusioni

Al netto di diverse prove, il codice è funzionante e replica in maniera soddisfacente il fenomeno fisico. Il problema principale del programma è il tempo di esecuzione della modalità statistica, che non permette di studiare e confrontare con rapidità i risultati. Come già accennato, questo è intrinsecamente legato al piccolo passo temporale della simulazione, che quando portata avanti con un numero alto di particelle aumenta notevolmente. Tanto più, nelle configurazioni con gradiente magnetico queste iterazioni sono aumentate per l'aggiornamento del vettore $\vec{B}$. La soluzione non è semplice, nel senso che andrebbero introdotti algoritmi più complessi usati nelle simulazioni di plasma, e con un approccio così classico l'unica alternativa possibile è aumentare la risoluzione temporale diminuendo il valore del `dt` come frazione del periodo ciclotronico. 

***Lorenzo Pace a.a. 2024-2025 Corso di Metodi Computazionali per la Fisica***

