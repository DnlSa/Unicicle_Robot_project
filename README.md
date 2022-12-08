# Unicicle_Robot_project
UnycleRobot

L'uniciclo è dotato di un sensore GPS caratterizzato da una certa risoluzione spaziale, il quale
fornisce una misura quantizzata e rumorosa della posizione dell'uniciclo. 
Più precisamente, le misure zx e zy della posizione (x,y) dell'uniciclo sono definite in questo modo: se qStep indica il passo di quantizzazione del GPS, la misura zx vale m*qStep se x (a cui si somma una perturbazione nx) è
compreso tra (m-1/2)qStep e (m+1/2)qStep, essendo m un intero. Lo stesso discorso vale per la y. In
formule questo può essere ottenuto definendo le due misure nel seguente modo:
zx = qStep*round((x+nx)/qStep),
zy = qStep*round((y+ny)/qStep),
essendo nx e ny errori gaussiani a media nulla con deviazione standard σp=10 pixel e round la funzione
che arrotonda all'intero più vicino.

Il robot possiede anche una bussola che gli fornisce una misura del suo orientamento. Cioè, nel
vettore z, oltre alle misure zx e zy, c'è una componente zb = θ + nb , che è una misura dell'orientamento
θ dell'uniciclo ed è caratterizzata da un errore nb da prendere gaussiano a media nulla con deviazione
standard σb=10 gradi.

La bussola può essere attivata o disattivata da tastiera, premendo per esempio il tasto 'B' per
attivarla e 'b' per disattivarla. Riportare a schermo l'indicazione sullo stato della bussola: ON o OFF.
Anche il passo qStep della misura del GPS deve essere modificabile da tastiera (per esempio con i
tasti + e -) da un minimo di 10 pixel senza limiti sul massimo.

Il reticolo che definisce la misura quantizzata del GPS va disegnato, insieme agli assi x=0 e y=0.
Lasciare attiva la funzionalità già presente nello sketch KalmanUniciclo.pde di modificare con le
frecce il tempo tStep tra una misura e l'altra. Come avviene in KalmanUniciclo.pde, quando in un
passo k non è disponibile nessuna misura, il filtro di Kalman usa solo la predizione basata sui passi
encoder delle ruote.
