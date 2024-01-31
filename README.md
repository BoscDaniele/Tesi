# Tesi
Definizione di Indicatori per la Caratterizzazione dello Stile di Guida di Veicoli Leggeri
## Target
1. Capire lo stato attuale del conducente che possono essere
   1. Fermo: accelerazione lungo l'asse x(asse direzione moto) pari a "zero" accelerazione gravitazionale distribuita lungo gli assi y e z (il peso viene spostato su una gamba, la bicicletta si inclina e così anche il sensore)
   2. Accelerazione:
      * Accelerometro
         - asse x: composta da una componente "costante" (nel senso che varia lentamente nel tempo, sembra essere molto inaffidabile) e da una componente variabile/oscillante (azione delle pedalate sull'accelerazione della bicicletta)
         - asse y: accelerazione oscillante dovuta all'azione del pedale
      * Giroscopio
         - rollio/roll (rotazione attorno a x): azione del pedale
   3. Azione Freno
      * Accelerometro
         - asse x: rapida diminuzione della parte costante accelerazione
      * Giroscopio
         - beccheggio/pitch (rotazione attorno a y): si dovrebbe notare l'azione del freno anteriore (sopratutto se la bici, come nel mio caso, è dotata di ammortizzatori)
   4. In Corsa
      * Accelerometro
         - asse x: l'accelerazione "costante" diminuisce in quanto la bicicletta è già a "regime" (ho mai detto che questa componente è altamente inaffidabile?), anche la parte oscillante dovuta alle pedalate dovrebbe ridursi in ampiezza (in quanto le pedalate riescono a trasferire una minore accelerazione) dovrebbe esserci anche un'aumento della frequenza delle pedalate in quanto la ruota non fornisce (quasi) più "resistenza"
      * Giroscopio
         - rollio/roll (rotazione attorno a x): azione del pedale
   5. Cambio Rapporti (marcia)
      * Accelerometro
         - asse x: l'accelerazione "costante" dovrebbe aumentare/diminuire, la componente oscillante, in caso di inserimento di una marcia superiore, aumenta in ampiezza e diminuisce in frequenza (l'azione del pedale trasferisce maggiore accelerazione ma diventa più "duro")
   6. Salita
      * Accelerometro
         - asse x: accelerazione "costante" diminuisce per effetto della forza di gravità (la componente gravitazionale viene registrata dal sensore anche se la bicicletta è ferma, quindi se sei fermo in salita/discesa per il sensore stai accelerando o in avanti(in discesa) o all'indietro (in salita)), la componente oscillante aumenta in ampiezza (a causa della forza di gravità c'è una maggiore differenza nella trasmissione dell'accelerazione a seconda del punto in cui si trova il pedale) ma diminuisce in frequenza (la bicicletta oppone più resistenza alle pedalate)
      * Giroscopio
         - beccheggio/pitch (rotazione attorno a y): rotazione positiva attorno all'angolo di pitch dovuta al cambiamento d'inclinazione della bicicletta
   7. Discesa (come salita ma al contrario)
   8. Buche/Dossi
      * Accelerometro
         - assi x e z: impulso
   9. Curva


## Procedure
   - db: ogni cartella è una sequenza di letture (di solito la prima o le prime due servono per ruotare il sistema di riferimento del sensore al fine di farlo coincidere con quello della bicicletta), per ogni rilievo il sensore raccoglie dati con una frequenza di 25Hz (una sequenza di dati ogni 0.04s) finchè il sensore non decide diversamente (nuovoApproccio2), nella prima colonna abbiamo il tempo (per renderlo leggibile sottrai la prima misura a ogni altra), nelle colonne 2, 3 e 4 abbiamo le accelerazioni in x, y e z espresse in mg (millesimi di g, quindi, in teoria ma non nella pratica, 1000 corrisponde a g, cioè 9.81m/s^2), nelle colonne 5, 6 e 7 abbiamo le velocità angolari attorno agli assi x, y e z espresse in mdegree al secondo(degree: 360° corrisponde all'angolo giro, dato che la misura è in milli degree l'angolo giro sono 360000) dovrai fare la conversione in radianti ovvero moltiplicare per pi/180. Nelle altre colonne abbiamo le misure del magnetometro (di cui possiamo fare a meno al momento), della pressione e della temperatura (che non penso ci serviranno)
   
   - Research: qui dentro metto gli articoli interessanti che trovo e il db degli articoli che potenzialmente possono essere interesanti ma che ancora non ho letto

   - GRot funzione matlab per calcolare la matrice di rotazione per portare l'accelerazione gravitazionale lungo l'asse z
   - GZRot come sopra ma aggiusta anche la rotazione attorno all'asse z
   
   - plotta3 esegue il plot di una matrice n*3 (per esempio gli passi la matrice delle accelerazioni e lui ti stampa 3 grafici, uno per la x, uno per la y e uno per la z)
   - multiPlotta3 come sopra ma stampa due matrici

   - RotMat funzione che riceve in ingresso un vettore di angoli e restituisce la matrice di rotazione, al momento non so se sarà ancora utiile

   - Gli altri file sono file Matlab (più o meno uno per ogni rilievo che ho fatto)
