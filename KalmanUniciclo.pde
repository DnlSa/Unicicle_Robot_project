/*

Pieraldo Santurro 0269843        Sabellico Danilo 0268692



Si disegna un robot mobile di tipo 'uniciclo'
e lo si controlla in velocità (senza limiti di saturazione)
verso il punto cliccato col mouse mediante la legge
di controllo proporzionale, tuttavia il calcolo del controllo
è fatto sulle variabili stimate anziché su quelle vere.
Sulla schermata viene anche disegnato (in giallo) il robot
nella posizione stimata (quello vero è disegnato in rosso).
In più, dei triangoli in bianco con bordo rosso rappresentano 
la posizione dei landmark presenti nell'ambiente. All'interno
di ciascun triangolo viene anche riportato un identificatore del
tipo Li per il landmark i-esimo.

Con le frecce è possibile modificare la frequenza delle
misure di distanza dai landmark o, più precisamente, il
tempo tStep che intercorre tra una misura e la successiva
(che comunque non può essere inferiore al ciclo di Processing,
per default dt = 1/60 secondo): con le frecce SU e GIU si 
incrementa e decrementa rispettivamente di 1 tale tempo mentre 
con le frecce DESTRA e SINISTRA lo si moltiplica e divide 
rispettivamente per 2 (per modificarlo più velocemente). Quando 
tale tempo è molto grande la ricostruzione diviene puramente 
odometrica.

Sulla schermata vengono riportati: le coordinate
(x,y,theta) dell'uniciclo, le coordinate (xDes,yDes)
del punto da raggiungere e il tempo impiegato per la
missione. Inoltre si riportano le velocità angolari 
(omegaR, omegaL) delle due ruote, la velocità longitudinale 
v1 e quella angolare v2.
Infine si indicano anche le grandezze stimate (xHat,yHat,thetaHat)
con la loro matrice di covarianza P e il tempo tStep (modificabile 
da tastiera) che intercorre tra una lettura e la successiva.
*/

// Dimensioni finestra
int sizeX = 1250;//1000;
int sizeY =  950;//550;

// Coordinate attuali uniciclo
float x = 0;
float y = 0;
float theta = 0;

// Coordinate desiderate uniciclo
float xDes = 0;
float yDes = 0;

// Caratteristiche fisiche uniciclo
float r = 8; // raggio ruote in pixel
float d = 25; // distanza tra le ruote in pixel
float w = 5; // spessore ruota in pixel
float R = 1.2*sqrt(pow(r,2)+pow(d/2+w/2,2)); // raggio robot

float dt = (float) 1/60; // tempo di campionamento

float e_p = 0; // errore posizionamento
float v1 = 0; // velocità lineare robot
float kv1 = 1; // costante legge proporzionale controllo v1
float v2 = 0; // velocità rotazionale robot
float kv2 = 2;  // costante legge proporzionale controllo v2
int nGiri = 0; // conta quanti giri su se stesso ha fatto il robot
float thetaDes; // orientamento desiderato (per raggiungere il target)
float t0,tempo; // tempo inizio e attuale missione
float tempoUltimaMisura; // tempo in cui si è effettuata l'ultima misura
long nStime = 0; // numero progressivo totale delle stime fatte

float omegaR = 0; // velocità angolare ruota destra
float omegaL = 0; // velocità angolare ruota sinistra
float uR, uL, uRe, uLe; // spostamenti ruote destra e sinistra veri e presunti (comandati)

// Variabili relative alla stima e al filtro di Kalman esteso

float sigmaX0=10;
float sigmaY0=10;
float sigmaGPS = 10; // deviazione standard dell'errore di stima del valore iniziale di y0
float sigmaThetaB = Gaussian(0,10); // deviazione standard (in rad) dell'errore di stima del valore iniziale di theta0
float sigmaTheta0 = 10*PI/180; // deviazione standard (in rad) dell'errore di stima del valore iniziale di theta0
float xHat = x + Gaussian(0,sigmaX0); // stima di x inizializzata al valore vero + una perturbazione gaussiana a media 0 e deviazione standard sigmaX0
float yHat = y + Gaussian(0,sigmaY0); // stima di y inizializzata al valore vero + una perturbazione gaussiana a media 0 e deviazione standard sigmaY0
float thetaHat = random(-PI,PI); // stima di theta inizializzata al valore vero + una perturbazione gaussiana a media 0 e deviazione standard sigmaTheta0 rad
float xHatMeno,yHatMeno,thetaHatMeno; // stime a priori di x,y,theta
float KR = 0.01; // coefficiente varianza errore odometrico ruota destra
float KL = KR; // coefficiente varianza errore odometrico ruota sinistra
float stDevR,stDevL; // deviazione standard errore odometrico ruota destra e sinistra



//// parte nuova ////
float n; 
int i;
int j ; 
float q;
float qStep=100;
int qStepLimit = 10 ; 
float xBase;
float yBase;
float nb = Gaussian(0,sigmaTheta0);
float Nx;
float Ny;
float mx;
float my;
int bussola=1;


// nuove matrici dichiarate 
float[][] RsOFF = idMat(2,pow(10,2)); // matrice di covarianza errore misura dai landmark
float[][] RsON = idMat(3,pow(10,2)); // matrice di covarianza errore misura dai landmark
float[][] innovazioneOFF = new float [2][1]; // innovazione EKF
float[][] innovazioneON = new float[3][1]; // innovazione EKF
float[][] Koff = new float[3][2]; // guadagno di Kalman
float[][] Kon = new float[3][3]; // guadagno di Kalman
float[][] Hoff = new float [2][3]; // matrice giacobiana H = dh/dx
float[][] Hon = new float[3][3] ;// matrice giacobiana H = dh/dx
float[][] F = {{1, 0, 0},{0, 1, 0},{0, 0, 1}}; // matrice giacobiana F=df/dx (alcuni elementi delle prime due righe vanno aggiustati durante l'esecuzione)
float[][] P = {{pow(sigmaX0,2), 0, 0},{0, pow(sigmaY0,2), 0},{0, 0, pow(sigmaTheta0,2)}}; // matrice di covarianza P inizializzata in base all'incertezza iniziale
float[][] Pmeno = new float[3][3]; // matrice di covarianza a priori P-
float[][] W = {{1, 0},{0,1},{1/d,-1/d}}; //  matrice giacobiana W = df/dw (gli elementi delle prime due righe vanno aggiustati durante l'esecuzione) 
float[][] Q = {{1, 0},{0,1}}; // matrice di covarianza del rumore di misura odometrico w (gli elementi sulla diagonale vanno aggiustati durante l'esecuzione)
float DeltaX,DeltaY,DeltaXY; // variabili di supporto
float[][] correzione = new float[3][1]; // termine correttivo stima
float tStep = 100000000; // tempo (in ms) tra una misura e la successiva (impostabile da tastiera)

//nx e ny errori gaussiani a media nulla con deviazione standard σp=10 pixel e round la funzione che arrotonda all'intero più vicino.
float nx; 
float ny;
float zx;
float zy;
float zb;



void setup() 
{
  size(1250,950);
  tempoUltimaMisura = 0; // Inizializzo a zero il tempo in cui è stata effettuata l'ultima misura
  xBase=width/2;
  yBase=height/2;
  
  bussola=0; //la bussola parte spenta
  
}

void draw() 
{
  
  background(0);
  pushMatrix();
  translate(sizeX/2,sizeY/2);

  if (keyPressed)
  {
    if (keyCode == UP) // aumento di 1 il tempo tra una misura e la successiva
    {
      tStep += 1;
    }
    if (keyCode == DOWN)  // decremento di 1 il tempo tra una misura e la successiva
    {
      tStep = max(0,tStep-1);
    }
    if (keyCode == RIGHT) // moltiplico per due il tempo tra una lettura e la successiva
    {
      tStep = tStep*2;
    }
    if (keyCode == LEFT) // divido per due il tempo tra una lettura e la successiva
    {
      tStep = tStep/2;
    }
    
    //mod 1
    if (key == '+') {
      qStep++; 
      
    }
    if (key == '-') {
      if(qStep > qStepLimit){
         qStep--;       
      }
    }
    if(key == 'b')
          bussola=0;
    if(key == 'B')
          bussola=1;    
  }
  
  if (mousePressed) // assegno target
  {
    xDes = mouseX - sizeX/2;
    yDes = sizeY/2 - mouseY;
    t0 = millis(); // inizio conteggio del tempo di missione
  }
  
////////////////LEGGE PROPORZIONALE DELLA VELOCITA  /////////////////////////////////////

// Calcolo errore e controllo basandomi sul valore stimato (xHat,yHat,thetaHat) delle variabili dell'uniciclo
  e_p = sqrt(pow(xDes-xHat,2)+pow(yDes-yHat,2));
  if (e_p > 1) // mi muovo solo se l'errore è maggiore di una certa quantità
  {
    tempo = (millis()-t0)/1000;  // tempo missione in secondi

    // assegno velocità secondo legge proporzionale (in termini delle quantità stimate!)
    v1 = -kv1*((xHat-xDes)*cos(thetaHat) + (yHat-yDes)*sin(thetaHat));

    // Calcolo l'angolo verso il target: scelgo il multiplo di 2PI 
    // più vicino all'orientamento corrente ***stimato*** del robot
    thetaDes = atan2(yDes-yHat,xDes-xHat) + nGiri*2*PI;
    if (abs(thetaDes+2*PI-thetaHat) < abs(thetaDes-thetaHat))
    {
      thetaDes = thetaDes+2*PI;
      nGiri += 1;
    }
    else
    {
      if (abs(thetaDes-2*PI-thetaHat) < abs(thetaDes-thetaHat))
      {
        thetaDes = thetaDes-2*PI;
        nGiri += -1;
      }
    }
    
   // assegno velocità angolare secondo legge proporzionale sempre in funzione delle quantità stimate   
    v2 = kv2*(thetaDes-thetaHat);
  }
  else // se penso di essere vicino al target mi fermo
  {
    v1 = 0;
    v2 = 0;
  }
  
  // Calcolo i movimenti da impartire alle ruote in base alle v1 e v2 trovate
  omegaR = (v1+v2*d/2)/r;
  omegaL = (v1-v2*d/2)/r;
  uRe = r*omegaR*dt; // spostamento comandato alla ruota destra (da considerarsi anche come informazione odometrica disponibile sul suo spostamento in (t,t+dt))
  uLe = r*omegaL*dt; // spostamento comandato alla ruota sinistra (da considerarsi anche come informazione odometrica disponibile sul suo spostamento in (t,t+dt))
  
  // Perturbo i due movimenti: gli spostamenti reali delle ruote non saranno mai esattamente uguali a quelli comandati 
  stDevR = sqrt(KR*abs(uRe));
  stDevL = sqrt(KL*abs(uLe));  
  uR = uRe + Gaussian(0,stDevR); // Spostamento vero ruota destra
  uL = uLe + Gaussian(0,stDevL); // Spostamento vero ruota sinistra
  
  
  
  // Dinamica effettiva dell'uniciclo
  x = x + ((uR+uL)/2)*cos(theta);
  y = y + ((uR+uL)/2)*sin(theta);
  theta = theta + (uR-uL)/d;

////////////////INIZIO FILTRO DI KALLMAN /////////////////////////////////////


  // STIMA FILTRO KALMAN ESTESO: PASSO di PREDIZIONE
  xHatMeno = xHat + ((uRe+uLe)/2)*cos(thetaHat);
  yHatMeno = yHat + ((uRe+uLe)/2)*sin(thetaHat);
  thetaHatMeno = thetaHat + (uRe-uLe)/d;

  //Aggiorno la giacobiana F (solo gli elementi variabili)
  F[0][2] = -(uRe+uLe)*sin(thetaHat)/2;
  F[1][2] = (uRe+uLe)*cos(thetaHat)/2;
  
  // Aggiorno W (solo gli elementi variabili)
  W[0][0] = .5*cos(thetaHat);
  W[0][1] = .5*cos(thetaHat);
  W[1][0] = .5*sin(thetaHat);
  W[1][1] = .5*sin(thetaHat);

  //Aggiorno Q (solo gli elementi variabili)
  Q[0][0] = KR*abs(uRe);
  Q[1][1] = KL*abs(uLe);
  
  // Calcolo matrice covarianza a priori
  Pmeno = mSum(mProd(mProd(F,P),trasposta(F)),mProd(mProd(W,Q),trasposta(W))); // Pmeno = F*P*F' + W*Q*W'

 // STIMA FILTRO DI KALMAN ESTESO: PASSO di CORREZIONE
  if (millis()-tempoUltimaMisura >= tStep) // attuo la correzione solo se ho le misure (che arrivano ogni tStep ms)
  {
    tempoUltimaMisura = millis(); // memorizzo il tempo in cui ho fatto l'ultima misura
    nStime++; // incremento il contatore delle stime fatte
    
      nx = Gaussian(0,10);
      ny = Gaussian (0,10);
      zx=(qStep*round((x+nx)/qStep));
      zy=(qStep*round((y+ny)/qStep));

   if(bussola == 0){          //bussola spenta
   
     if (xHatMeno < zx-qStep/2 || xHatMeno > zx + qStep/2 || yHatMeno < zy-qStep/2 || yHatMeno > zy+qStep/2){
     
      Hoff[0][0] =1;
      Hoff[0][1] =0;
      Hoff[0][2] =0; 
      Hoff[1][0] =0;
      Hoff[1][0] =1;
      Hoff[1][2] =0;
      
      innovazioneOFF[0][0]=  (zx - xHatMeno);  
      innovazioneOFF[1][0]= (zy - yHatMeno);
      
    
     }
     else
     {
       
       //robot dentro il quadratino della bussola
       
       innovazioneOFF[0][0]= 0;  
       innovazioneOFF[1][0]=0;
     }
     
     // Calcolo guadagno Kalman e aggiorno covarianza bussola OFF
     Koff = mProd(mProd(Pmeno,trasposta(Hoff)),invMat(mSum(mProd(mProd(Hoff,Pmeno),trasposta(Hoff)),RsOFF))); 
     P = mProd(mSum(idMat(3,1),mProd(idMat(3,-1),mProd(Koff,Hoff))),Pmeno);
      
     // correggo la stima
     correzione = mProd(Koff,innovazioneOFF);
   }
   else
   {
     
     //bussola attivata, cioè bussola=1
     
     nb=Gaussian(0,10);
     zb=theta + nb;
     
     if (xHatMeno < zx -qStep/2 || xHatMeno > zx+qStep/2 || yHatMeno < zy-qStep/2 || yHatMeno > zy+qStep/2) {
       
        // se il robot si trova fuori dal quadrato illuminato e la bussola è accesa
        
        Hon[0][0] = 1;
        Hon[0][1] = 0;
        Hon[0][2] = 0;
        Hon[1][0] = 0;
        Hon[1][1] = 1;
        Hon[1][2] = 0;
        Hon[2][0] = 0;
        Hon[2][1] = 0;
        Hon[2][2] = 1;
      
    
    
       innovazioneON[0][0]=(zx -xHatMeno); 
       innovazioneON[1][0]=(zy - yHatMeno);
       
   }
   else
   {
     
     //robot dentro quadratino e bussola ON
     
     innovazioneON[0][0]=0;
     innovazioneON[1][0]=0;
     
   }
   
   innovazioneON[2][0] = zb - thetaHatMeno;
   
   // Calcolo guadagno Kalman e aggiorno covarianza bussola ON
   Kon = mProd(mProd(Pmeno,trasposta(Hon)),invMat(mSum(mProd(mProd(Hon,Pmeno),trasposta(Hon)),RsON)));
   P = mProd(mSum(idMat(3,1),mProd(idMat(3,-1),mProd(Kon,Hon))),Pmeno);
   correzione=mProd(Kon,innovazioneON); 
   }
   
   //correggo i valori di xHat, yHat e thetaHat in ogni caso
   
   xHat = xHatMeno + correzione[0][0];
   yHat = yHatMeno + correzione[1][0];
   thetaHat = thetaHatMeno + correzione[2][0];
   
   }
   else
   {
    // se non ho misure non correggo nulla
    xHat = xHatMeno;
    yHat = yHatMeno;
    thetaHat = thetaHatMeno;
    P = Pmeno;
   }

/////////////FINE EKF/////////////////////////////////////////////////////

// disegno rettangolo che  si illumina 
  
  rect(zx-qStep/2,-zy-qStep/2,qStep,qStep);
  
    

// Disegno il robot vero e quello stimato
  robot(x,y,theta,1); // l'argomento 1 fa un robot rosso (robot reale)
  robot(xHat,yHat,thetaHat,0); // l'argomento 0 un robot giallo (robot nella posa stimata)

// FUNZIONE CHE DISEGNA IL RETICOLO  ////
  
   n = width/qStep ; // numero di linee verticali e orizzontali 
   round(n);
   drawline(n, qStep/2); // disegna la griglia
   popMatrix();

  ///// funzione di scrittura testi ///////////
  Testi(); // vedasi funzione nelle ultime righe 
  
}


void robot(float x, float y, float theta, int colore)
{
// funzione che disegna uniciclo in (x,y) con orientamento theta      
  pushMatrix();
  translate(x,-y);
  rotate(-theta);
  if (colore==1)
  {
    fill(255,0,0);
  }
  else
  {
    fill(255,255,0);
  }
  ellipse(0,0,2*R,2*R); // il robot
  fill(0,0,255);  
  rect(-r,-d/2-w/2,2*r,w);
  rect(-r,d/2-w/2,2*r,w);
  fill(255);
  ellipse(.8*R,0,.2*R,.2*R);
  ellipse(-.8*R,0,.2*R,.2*R);  
  fill(0,255,0);
  triangle(-.1*R,.3*R,-.1*R,-.3*R,.5*R,0);
  popMatrix();
}


void arco(float x,float y,float theta,float rMax, float betaMax) {                                        //funzione che realizza l'arco intorno al robot rosso (robot reale)
  fill(0,255,0);
  pushMatrix();
  translate(x,-y);
  rotate(-theta);
  arc(0, 0, 2*rMax, 2*rMax, -betaMax/2, betaMax/2);
  popMatrix();
}


/******************************************************
/******************************************************
  DA QUI IN POI CI SONO FUNZIONI DI SERVIZIO: 
  1) CALCOLO ERRORE GAUSSIANO
  2) ALGEBRA MATRICIALE
/******************************************************
/******************************************************/

float Gaussian(float media, float stDev) // Restituisce variabile N(media,stDev^2) approssimata sfruttando il teorema del limite centrale
{
  float somma = 0;
  for (int k=0; k<27; k+=1) // 27 in modo che sqrt(3/27)=1/3
  {
    somma = somma + random(-stDev/3,stDev/3);
  }
  return media + somma;
}

float[][] mProd(float[][] A,float[][] B) // Calcola prodotto di due matrici A e B (si assume, senza controllarlo, che numero colonne A = numero righe B!)
{
  int nA = A.length;
  int nAB = A[0].length;
  int nB = B[0].length;
  
  float[][] C = new float[nA][nB]; 

  for (int i=0; i < nA; i++) 
  {
    for (int j=0; j < nB; j++) 
    {  
      for (int k=0; k < nAB; k++) 
      {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return C;
}

float[][] mSum(float[][] A,float[][] B) // Calcola la somma di due matrici A e B (si assume, senza controllarlo, che A e B abbiano stesse dimensioni!)
{
  int nA = A.length;
  int nB = A[0].length;
  
  float[][] C = new float[nA][nB]; 

  for (int i=0; i < nA; i++) 
  {
    for (int j=0; j < nB; j++) 
    {  
      C[i][j] = A[i][j] + B[i][j];
    }
  }
  return C;
}

float[][] trasposta(float[][] A) // Calcola la trasposta di una matrice A
{
  int nR = A.length;
  int nC = A[0].length; 
  
  float[][] C = new float[nC][nR]; 

  for (int i=0; i < nC; i++) 
  {
    for (int j=0; j < nR; j++) 
    {  
      C[i][j] = A[j][i];
    }
  }
  return C;
}


float[][] minore(float[][] A, int i, int j) // Determina il minore (i,j) di una matrice A (si assume, senza controllarlo, che A sia quadrata!)
{
  int nA = A.length;
  float[][] C = new float[nA-1][nA-1];
  
  for (int iM = 0; iM < i; iM++)
  {
    for (int jM = 0; jM < j; jM++)
    {
      C[iM][jM] = A[iM][jM];
    } 
    for (int jM = j; jM < nA-1; jM++)
    {
      C[iM][jM] = A[iM][jM+1];
    } 
  }
  for (int iM = i; iM < nA-1; iM++)
  {
    for (int jM = 0; jM < j; jM++)
    {
      C[iM][jM] = A[iM+1][jM];
    } 
    for (int jM = j; jM < nA-1; jM++)
    {
      C[iM][jM] = A[iM+1][jM+1];
    } 
  }
  return C;
}


float det(float[][] A) // Calcola il determinante di A (si assume, senza controllarlo, che A sia quadrata!)
{
  int nA = A.length;
  float determinante = 0;
  
  if (nA == 1)
  {
    determinante = A[0][0];
  }
  else
  {
    for (int j=0; j < nA; j++) 
    {
      determinante = determinante + A[0][j]*pow(-1,j)*det(minore(A,0,j));
    }
  }
  return determinante;
}


float[][] invMat(float[][] A) // Calcola l'inversa di una matrice A (si assume, senza controllarlo, che A sia quadrata!)
{
  int nA = A.length;
  float[][] C = new float[nA][nA];
  float detA = det(A);

/*
  if (abs(detA)<0.001) // Per evitare casi singolari con determinanti vicini a 0
  {
    if (detA>0)
    {
      detA = 0.001;
    }
    else
    {
      detA = -0.001;
    }
  }
*/  
  
  if (nA == 1)
  {
    C[0][0] = 1/detA;
  }
  else
  {
    for (int i=0; i < nA; i++) 
    {
      for (int j=0; j < nA; j++) 
      {
        C[j][i] = pow(-1,i+j)*det(minore(A,i,j))/detA;
      }
    }
  }
  return C;
}

float[][] idMat(int nA, float sigma) // Assegna una matrice identità di ordine nA moltiplicata per una costante sigma
{
  float[][] I = new float[nA][nA]; 

  for (int i=0; i < nA; i++) 
  {
    for (int j=0; j < nA; j++) 
    {  
      I[i][j] = 0;
    }
    I[i][i] = sigma;
  }
  return I;
}

//tale funzione in realta disegna un reticolo scalato in quanto se predessi 10 pixel effettivi su schermo il reticolo si saturerebbe 
void drawline(float n,float  qStep)
 {
   float count = qStep*2;
   pushMatrix();
   stroke(#FF0505);
   line(0,-height/2,0,height); //linea verticale 
   line(-width/2,0,width,0); // linea orizzontale
   stroke(255);
   
   for(i=0  ; i<n ; i++){
      line(-width/2,-qStep,width,-qStep); // linea orizzontale 
      line(-width/2,+qStep,width,+qStep); // linea orizzontale 
   
   for(j=0; j<n ; j++){
      line(-qStep,-height/2,-qStep,height);
      line(+qStep,-height/2,+qStep,height);        
    }
       qStep += count;
  }
  popMatrix();
}


void Testi(){

  textSize(20);
  fill(0,0,255);
  text("v1 (pixel/s) = ",10,20); 
  text(v1,200,20);
  text("v2 (gradi/s) = ",10,50); 
  text(v2*180/PI,200,50);
  
  fill(255,0,0);  
  text("x = ",10,160); 
  text(x,80,160);
  text("y = ",10,190); 
  text(y,80,190);
  text("theta = ",10,220); 
  text(theta*180/PI,100,220);  

  fill(255,255,255);
  text("tempo = ",10,110); 
  text(tempo,120,110);  

  fill(0,0,255);
  text("omegaR (gradi/s) = ",650,20); 
  text(omegaL*180/PI,900,20);
  text("omegaL (gradi/s) = ",650,50); 
  text(omegaL*180/PI,900,50);
  
  fill(255,0,0);
  text("xDes = ",700,100); 
  text(xDes,800,100);
  text("yDes = ",700,130); 
  text(yDes,800,130);

  fill(255,255,0);  
  text("xHat = ",10,280); 
  text(xHat,120,280);
  text("yHat = ",10,310); 
  text(yHat,120,310);
  text("thetaHat = ",10,340); 
  text(thetaHat*180/PI,160,340);  

  fill(255,255,255);
  text("nStime = ",10,390); 
  text(nStime,120,390);  

  fill(255,255,255);
  text("tStep (ms) = ",10,420); 
  text(tStep,150,420);  

  fill(255,255,0);  
  text("P = ",10,460); 
  text(P[0][0],10,490); text(P[0][1],100,490); text(P[0][2],190,490);
  text(P[1][0],10,520); text(P[1][1],100,520); text(P[1][2],190,520); 
  text(P[2][0],10,550); text(P[2][1],100,550); text(P[2][2],190,550);   
    
  text("qStep = ",10,660); // SCALO LA MISURE DEL RETICOLO in quanto dovro riempirli 
  text(round(qStep),100,660); 
  
  text("Bussola = ",10,700); // SCALO LA MISURE DEL RETICOLO in quanto dovro riempirli 
  if(bussola == 1)
      text("OFF",130,700);
  else  
      text("ON",130,700);
  
  text("q=distanza= ", 10, 730);
  text(q,130,730);
    
 
 
}
