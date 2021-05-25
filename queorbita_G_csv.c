/*	====================================================
	
	Calcula todas as órbiras com valores dados de e, q, 
	m e M que satisfaçam as condições observacionais 
	fornecidas (posições e velocidades na linha de visada 
	de dois	objetos). 

	-> o plano do céu é o plano XY e o eixo Z aponta 
	para o observador.
	-> M esta no foco da órbita e é o centro do sist.
	de coordenadas (0,0,0);

	ATENÇÃO: Aqui modifiquei para G=43007.1 que serve para as galáxias
	do VOLKER

	Condições iniciais:
	-> Excentricidade da órbita (e)
	-> dist. de pericentro (q)
	-> Massas M e m

        SISTEMA DE UNIDADES: Os dados de entrada devem ser informados
        no sistema de coordenadas usado nas simulações com Gadget2
	, em que
	-> G=43007.1 
	-> Unid de distancia = 1 kpc
	-> unid de velocidade = 1 km/s
	-> Unid de massa = 1e10 M_sun
	-> Unid de tempo = 0.9779 Gyr


        ATENÇÃO: A velocidade na linha de visada (Vsys) deve ser
        informada no sistema de unidades acima, em que o sinal é
        invertido com respeito à à convenção observacional. Na
        convenção observacional o eixo z aponta do observador
        para o objeto (velocidade positiva para o objeto se afastando
        do observador). Aqui deve ser o contrário.

	Escrito em 1998.

	Modificado em fev/2004 

	Modificado em out/2005 - para inclui situações pré e
	pós-pericentro. Antes só levava em conta o
	pós-pericentro. Feitas verificações gerais, corrigidos
	pequenos erros e melhorada a escolha dos diversos angulos e
	transformações de coords para tornar mais
	inteligível. Adicionados muitos comentarios.

	====================================================
	compile instructions: 
	gcc queorbita_G_csv.c -o ../bin/queorbita_G_csv -lm 

*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void main()
{
  float	
    // G=1.,
    // Modifico para usar com os modelos do Volker:
    G=43007.1,
    qini,                       /* dist pericentro inicial */ 
    qfin,                       /* dist pericentro final */ 
    qstep,                      /* dist pericentro passo */ 
    eini,                       /* excentricidade inicial */ 
    efin,                       /* excentricidade final */ 
    estep,                      /* excentricidade passo */ 
    M,                          /* Massa galáxia central */
    m,                          /* Massa galáxia 2 */
    q,                          /* dist pericentro (índice do loop) */
    e,                          /* excentricidade (índice do loop) */
    rnow,
    mu,                         /* G*(M+m) */
    vq,                         /* velocidade no pericentro (módulo) */
    v,
    a,
    h,
    w,                          /* Ângulo entre rnow e q */
    w_dentro,                   
    theta,                      /* Ângulo entre rnow e a direção X do sistema órbita */
    phi,                        /* Ângulo entre rnow e V */
    phi_dentro,
    aux,                        
    aux1,
    gzini,                     /* coordenada Z inicial */
    gzfin,                     /* coordenada Z final */
    gzstep,                    /* coordenada Z passo */
    gx,                        /* coordenada X sistema céu (dado de entrada) */
    gy,                        /* coordenada Y sistema céu (dado de entrada) */
    gz,                        /* coordenada Z sistema céu (A ser calculado: índice do loop) */
    x,                         /* coordenada X no sistema órbita */
    y,                         /* coordenada Y no sistema órbita */
    z,                         /* coordenada Z=0 no sistema órbita */
    vx,                        /* velocidade X no sistema órbita */
    vy,                        /* velocidade Y no sistema órbita */
    thetav,                    /* Ângulo entre V e X do sistema órbita */
    vsys,                      /* Velocidade sistemática */
    vsyser,                    /* Erro na Velocidade sistemática */
    sphi,
    stheta,
    spsi,                      /* angulo para girar V em torno de RNOW */
    SA1,                       /* Elem. matriz transformação coord */
    SA2,                       /*                idem              */
    SA3,                       /*                idem              */
    SB1,                       /*                idem              */
    SB2,                       /*                idem              */
    SB3,                       /*                idem              */
    SC1,                       /*                idem              */
    SC2,                       /*                idem              */
    SC3,                       /*                idem              */
    svv[3],                    /* vetor velocidade no sist. rnow paral Z */
    vv[3],                     /* vetor velocidade sist ceu */
    spin[3],                   /* vetor spin da órbita */
    spinmod,                   /* módulo do vetor spin da órbita */
    qvec[3],                   /* vetor pericentro */
    sqvec[3],                  /* vetor pericentro no sist. rnow paral Z */
    qvecmod,                   /* módulo do vetor pericentro */
    spindsk[3],                /* vetor spin do disco M */
    apocenter;                 /* distância de apocentro */

  int orb,
    i,
    circfirstpass,
    pre_pos_peri;
  
  puts("	Excentricidade: (elipse: e < 1)");
  puts("			(paráb:  e = 1)");
  puts("			(hiperb: e > 1)");
  puts("	O programa vai varrer um intervalo de ");
  puts("	excentricidades. Indique os valores inicial");
  puts("	(eini), final (efin) e o passo (estep):");
  scanf("%f %f %f",&eini,&efin,&estep);

  puts("	O programa vai varrer um intervalo de ");
  puts("	pericentros. Indique os valores inicial");
  puts("	(qini), final (qfin) e o passo (qstep):");
  scanf("%f %f %f",&qini,&qfin,&qstep);

  puts("	Massas M e m:");
  scanf("%f %f",&M,&m);

  puts("	Particula de massa M em (0,0).");
  puts("	Indique a posição XY da partícula de massa m:");
  scanf("%f %f",&gx,&gy);

  /* Escolha dos gz */
  gzini=-10.*sqrt(gx*gx + gy*gy); 
  /*   gzini=-2.*sqrt(gx*gx + gy*gy) * 0.2 ;  */
  gzfin=-gzini; 
  gzstep=(gzfin - gzini)/1000.;  


  /* ATENÇÃO: Para indicar a velocidade sistemática, deve-se levar em
     conta que no sistema de referencias usado, o eixo z cresce para o
     observador, que é o contrário da convenção observacional. Isso
     significa que o sinal informado deve ser o contrário do 
     OBSERVACIONAL. */

  puts(" V_sys de m (em relação a M) e seu erro (+-V_sys_err):"); 
  scanf("%f %f",&vsys,&vsyser); 

  puts(" Indique as 3 coordenadas do vetor de spin do"); 
  puts(" disco principal (para ver se o movimento é"); 
  puts(" progrado ou retrógrado):"); 
  scanf("%f %f %f",&spindsk[0],&spindsk[1],&spindsk[2]);

  mu=G*(M+m);

  /*Inicia a saída da tabela de órbitas */
  puts("	");
  puts("	");  
  puts("	-----------------------------------------");  
  puts("	Tabela das possíveis órbitas encontradas:");
  puts("	-----------------------------------------");
  puts("	");
  puts("	Dados de entrada:");
  printf("	Galáxia A em (0,0,0) com massa %g\n",M);
  printf("	Galáxia B em (%g,%g,%g a %g) com massa %g\n",gx,gy,gzini,gzfin,m);
  printf("	Excentricidade: de %g a %g a passos de %g\n",eini,efin,estep);
  printf("	Pericentro: de %g a %g a passos de %g\n",qini,qfin,qstep);
  printf("	Velocidade sistemática de B em relação a A: V_sys=%g+-%g    (ATENÇÃO: No sistema de coordenadas da simulação, sinal inverso do observacional).\n",vsys,vsyser);
  puts("	-----------------------------------------");
  puts("	CULUNAS:");
  puts("	e = excentricidade");
  puts("	q = distancia de pericentro");
  puts("	rnow = distancia atual");
  puts("	vq = velocidade no pericentro");
  puts("	gx = coordenada X (sistema CÉU, dado de entrada)");
  puts("	gy = coordenada Y (sistema CÉU, dados de entrada)");
  puts("	gz = coordenada Z (sistema CÉU)");
  puts("	VX_ceu = Velocidade coordenada X (sistema CÉU)");
  puts("	VY_ceu = Velocidade coordenada Y (sistema CÉU)");
  puts("	VZ_ceu = Velocidade coordenada Z (sistema CÉU)");
  puts("	x = coordenada X no sistema da órbita (orb no plano XY, pericentro em Y, Z=0, antihorário)");
  puts("	y = coordenada Y no sistema da órbita ");
  puts("	vx = componente velocidade X no sistema da órbita ");
  puts("	vy = componente velocidade Y no sistema da órbita ");
  puts("	spx, spy, spz = componentes vetor spin da órbita no sistema CÉU ");
  puts("	qx, qy, qz = componentes vetor direção pericentro no sistema CÉU ");
  puts("	vsys = componente velocidade Z no sistema CÉU ");
  puts("	Dir = Direção da órbita: progrado ou retrógrado em rel ao disco M ");
  puts("	spin-orb = Ângulo entre vetor de spin do disco M e spin da órbita ");
  puts("	PERIC = Situação PRÉ ou PÓS PERICENTRO ");
  puts("	pos-peri = Ângulo entre vetor posição (g) e o vetor de pericentro (q) ");
  puts("	APOCENTRO = Distância de apocentro ");
  puts("	-----------------------------------------");
  puts("	");
  puts("e,q,rnow,vq,gx,gy,gz,VX_ceu,VY_ceu,VZ_ceu,x,y,vx,vy,spx,spy,spz,qx,qy,qz,vsys,Dir,spin-orb,PERIC,pos-peri,APOCENTRO");
  

  for (e=eini;e<=efin;e=e+estep) {        // loop sobre e 
    for (q=qini;q<=qfin;q=q+qstep) {        // loop sobre q
      if (e < 1.0) {
	/* eliptica */
	orb=(int)0;
	a=q/(1.-e);
	vq=sqrt(2*mu*( 1/q  - 1/(2*a) ));
      } else if (e == 1.0) { 
	/* parabolica */
	orb=(int)1;
	vq=sqrt(2.*mu/q);
      } else if (e > 1.0) {
	/* hiperbolica */
	orb=(int)2;
	a=q/(e-1.);
	vq=sqrt(2*mu*( 1/q  + 1/(2*a) ));
      }
  
      h=q*vq;

      /* Loop sobre possíveis valores da coordenada Z (gz) da galáxia de
	 massa m. Escolhi o step baseado na posição fornecida (gx,gy)
	 sobre o plano XY (atrás). */
      circfirstpass=1;
      if (e == 0.0) circfirstpass=0;     /* usado adiante no caso de orb circular */
      for (gz=gzini;gz<=gzfin;gz=gz+gzstep){         // loop sobre gz  
	rnow=sqrt(gx*gx + gy*gy + gz*gz);
	switch (orb){
	case 0:
	  /* elipse */
	  if ((e == 0.0) && (q >= sqrt(gx*gx + gy*gy)) ) {     
	    /* No caso de orbita circular, temos
	       somente dois valores possiveis de gz,
	       um positivo e um negativo, aque
	       satisfazem a condição de que
	       rnow=q=cte.

	       Neste caso, tenho que modificar o gz
	       para ir direto a esses valores,
	       passando pelo loop somente 2 vezes,
	       uma para o gz negativo e uma para o
	       positivo.
	    */
	    rnow = q;
	    if (circfirstpass == 0) { 
	      gz = -sqrt(rnow*rnow - gx*gx - gy*gy);
	      circfirstpass = 1;
	    } else {
	      if (circfirstpass == 1) { 
		gz = sqrt(rnow*rnow - gx*gx - gy*gy);
		circfirstpass = 2;
	      } else {
		if (circfirstpass == 2) { 
		  gz = gzfin+1.0;    /* para interromper o loop */
		}
	      }
	    }
	  }



	  apocenter=a*(1+e);

	  /* Tem que testar se o RNOW não está sendo maior que a
	     distância de apocentro. Se for, então faço v=0 só para
	     descartar a conta */
	  if (rnow < apocenter) {
	    v=sqrt(2.*(mu/rnow - mu/(2.*a)));
	    if (e == 0.0) {
	      aux = 1.0;    /* se a orbita for circular, o angulo vai
			       ser indefinido. Defino que seu cosseno
			       é 1 (angulo será zero) */
	    } else {
	      aux=(h*h/(rnow*mu) - 1)/e ;  
	      aux=aux-0.000000001;         /* Para evitar problemas 
					   já que às vezes no pericentro 
					   o aux que deveria ser 1.0 dá um
					   pouco maior que 1, por erros 
					   numéricos.  
					*/
	    }

	  } else {
	    v=0;
	  }

	  break;

	case 1:
	  /* parablola */
	  v=sqrt(2.*mu/rnow);
	  aux=2*q/rnow - 1;  
	  aux=aux-0.000000001;
	  break;
	case 2:
	  /* hiperbole */
	  v=sqrt(2.*(mu/rnow + mu/(2.*a)));
	  aux=(h*h/(rnow*mu) - 1)/e ;  
	  aux=aux-0.000000001;
	  break;
	} /* end switch */

	/* Calculo W: ângulo de rnow com a direção do pericentro (Y),
	   medido no sentido anti-horario de q para rnow. */
	w=acos(aux);     /* calculado assim o w sempre será um angulo
			    entre 0 e 180. Este angulo vale para a
			    situação PÓS-PERICENTRO. No PRÉ-PERI o
			    angulo seria -w */


	/* Calculo ângulo theta (ângulo de rnow com a direção X,
	   medido no sentido anti-horario). */
	theta=M_PI/2. - w ;  /* definiodo assim para a situação pré-pericentro */



	/* Calculo posicao - situação PRÉ-PERICENTRO*/
	x=rnow*cos(theta);
	y=rnow*sin(theta);
	/* Calculo angulo phi entre rnow e v */
	aux=h/(rnow*v);
	phi=asin(aux);   /* calculado assim o PHI sempre será um
			    angulo do 1º quadrante (0 a 90). Este
			    angulo vale para a situação
			    PÓS-PERICENTRO. No PRÉ-PERI o angulo seria
			    180-phi */

	/* Calculo o angulo do vetor veloc. com o eixo X medindo anti-horario	*/
	thetav=theta + M_PI - phi;  /* definiodo assim para a situação pré-pericentro */
	/* Calculo velocidade	*/
	vx=v*cos(thetav);
	vy=v*sin(thetav);
    
	/* Até aqui calculei as componentes X e Y da posição no plano
	   da órbita e as componentes VX e VY, do vetor velocidade no
	   plano da órbita, a partir dos angulos theta (posição da
	   particula rm rnow) e phi (angulo entre rnow e v). Os
	   valores obtidos aqui são calculados para o caso
	   PRÉ-PERICENTRO. 

	   A situação simétrica PÓS-PERICENTRO pode ser obtida
	   fazendo-se:
	   x = -x
	   y = y
	   vx = vx
	   vy = -vy

	   e para os angulos theta e phi:
	   w_pos = -w_pre
	   theta_pos = M_PI - theta_pre
	   phi_pos = M_PI - phi_pre
	   thetav_pos = theta_pos + M_PI - phi_pos    

	*/




	/* 	Até aqui calculei os vetores posição e velocidade num sistema
		de coordenadas em que M está no centro, o plano da órbita é o
		plano XY, o pricentro está na direção Y e a órbita é
		antihorária, ou seja, o vetor de spin da órbita aponta para o
		Z positivo. 
	
		O importante até aqui são os ângulos entre rnow e v
		(phi) e entre q e rnow (w). Agora devo girar o vetor v
		em torno da direção de rnow e medir sua projecao na
		linha de visada para ver se em alguma posição ele
		combina com o V_sys de entrada. Isso só será
		satisfeito, obviamente, se v for maior que V_sys.

		O jeito de fazer isso é assim: 

		- mudanca de coordenadas para um sistema em que o vetor RNOW
		esteja alinhado com um dos eixos. Para isso escrevo RNOW em
		coordenadas esfericas (s_r,stheta,sphi), giro um ângulo -sphi
		em torno de Z, e giro um ângulo stheta em torno de Y.

		- nesse sistema crio um vetor SVV com as coordenadas dos pontos
		que representam V, que estão sobre um circulo que é o corte
		de um cone formado pelo giro de V em redor de RNOW. 

		- mudanca de coordenada inversa sobre os pontos de V, para
		obter os possíveis vetores V no sistema do céu. A coordenada
		que interessa aqui é a alinhada com a linha de visada.
	*/

	/* testa pra ver se vale a pena fazer as contas */
	if ( (fabs(v) >= fabs(vsys)) && ((orb == 0 && rnow <= apocenter) || (orb > 0)) ) {    // if ( (fabs(v) >= fabs(vsys)) ...    
	  /* Escrevo o vetor RNOW em coordenadas esfericas para obter
	     suas coordenadas stheta e sphi */
	  sphi=atan2(gy,gx);
	  /* inverto o phi, por causa do sentido da transformação */
	  sphi=sphi-M_PI/2.0;

	  stheta=acos(gz/rnow);
	  /* troco o sinal do teta, por causa do sentido da transformação */
	  stheta=-stheta;

	  for (pre_pos_peri=1 ; pre_pos_peri <= 2 ; pre_pos_peri=pre_pos_peri+1){   // loop sobre situação pre e pos pericentro  
	    /* Este loop serve para levar em conta as 2 possiveis
	       situações: PRE-PERICENTRO e PÓS-PERICENTRO 
	       pre_pos_peri=1 para  PÓS-PERICENTRO
	       pre_pos_peri=2 para  PRE-PERICENTRO 
	    */
	    
	    if (pre_pos_peri==1) {  /* Se for PÓS-PERICENTRO, PHI e W
				       ficam como foranm definidos */
	      phi_dentro = phi;
	      w_dentro = w;
	    }
	    if (pre_pos_peri==2) {  /* Se for PRE-PERICENTRO tem que
				       acertar o PHI, e inverter o
				       W */
	      phi_dentro = M_PI - phi;
	      w_dentro = -w;
	    }


	    for (i=0;i<=359;i=i+1) {              //loop para rotação de V ao redor de RNOW 
	      /* Calculo o vetor velocidade no sistema em que RNOW
		 paralelo a z */
	      spsi=((float) i) * M_PI/180.;
	      /* v*sin(phi_dentro) é o raio do circulo descrito pelo
		 vetor V no sistema novo quando gira em torno de rnow*/
	      svv[0]=v*sin(phi_dentro) * cos(spsi); 
	      svv[1]=v*sin(phi_dentro) * sin(spsi);
	      svv[2]=v*cos(phi_dentro);

	      /* Calculo o vetor pericentro no sistema em que RNOW
		 paralelo a z 

		 q*sin(-w_dentro) é o raio do circulo descrito pelo
		 vetor PERI no sistema novo quando gira em torno de
		 rnow.

		 Tenho que trocar o sinal do W porque o W é definido
		 como o angulo entre q e rnow, medido de q para rnow,
		 mas para calcular o q quero o angulo de rnow para q.*/
	      sqvec[0]=q*sin(-w_dentro) * cos(spsi);
	      sqvec[1]=q*sin(-w_dentro) * sin(spsi);
	      sqvec[2]=q*cos(-w_dentro);



	      /* Calculo a transformação do vetor velocidade do sistema
		 rnow//z para o sistema do ceu. Uso os velhos angulos de
		 Euler SPHI, STHETA e SPSI, calculados acima.
	       
		 Elementos da matriz de transformação do sist. ceu para
		 para o sistema em que RNOW paralelo a Z */
	      SA1 =  cos(spsi)*cos(sphi) - cos(stheta)*sin(sphi)*sin(spsi);
	      SA2 = -sin(spsi)*cos(sphi) - cos(stheta)*sin(sphi)*cos(spsi);
	      SA3 =  sin(stheta)*sin(sphi);
	      SB1 =  cos(spsi)*sin(sphi) + cos(stheta)*cos(sphi)*sin(spsi);
	      SB2 = -sin(spsi)*sin(sphi) + cos(stheta)*cos(sphi)*cos(spsi);
	      SB3 = -sin(stheta)*cos(sphi);
	      SC1 =  sin(spsi)*sin(stheta);
	      SC2 =  cos(spsi)*sin(stheta);
	      SC3 =  cos(stheta);


	      /*  Calculo a transformação do vetor velocidade DO
		  sist. em que RNOW paralelo a Z para o sistema ceu*/
	      vv[0] = svv[0]*SA1 + svv[1]*SA2 + svv[2]*SA3 ;
	      vv[1] = svv[0]*SB1 + svv[1]*SB2 + svv[2]*SB3 ;
	      vv[2] = svv[0]*SC1 + svv[1]*SC2 + svv[2]*SC3 ;

	      /* testo pra ver se satisfaz a condicao imposta por vsys */
	      if (vv[2] >= vsys-vsyser && vv[2] <= vsys+vsyser) {        // if que testa se Vsys é respeitado  
		/* Calculo o "vetor de spin" da orbita (r X v) já no sistema ceu */
		spin[0] = gy*vv[2]-gz*vv[1];
		spin[1] = gz*vv[0]-gx*vv[2];
		spin[2] = gx*vv[1]-gy*vv[0];
		/* E seu módulo */
		spinmod=sqrt(spin[0]*spin[0] + spin[1]*spin[1] + spin[2]*spin[2]); 

		/*  Calculo o vetor de pericentro: transformação direta
		    do vetor de pericentro escrito no sistema em que
		    RNOW paralelo a Z para o sistema ceu*/
		qvec[0] = sqvec[0]*SA1 + sqvec[1]*SA2 + sqvec[2]*SA3 ;
		qvec[1] = sqvec[0]*SB1 + sqvec[1]*SB2 + sqvec[2]*SB3 ;
		qvec[2] = sqvec[0]*SC1 + sqvec[1]*SC2 + sqvec[2]*SC3 ;
		/* E seu módulo */
		qvecmod = sqrt(qvec[0]*qvec[0] + qvec[1]*qvec[1] + qvec[2]*qvec[2]); 


		/* TESTE com o vetor R */
		/* 	      printf("\ngx = %f\n",gx); */
		/* 	      printf("gy = %f\n",gy); */
		/* 	      printf("gz = %f\n\n",gz); */

		/* transf do sistema ceu para r//z - esta é a transf inversa*/
		/* 	      aaax=gx*SA1 + gy*SB1 + gz*SC1; */
		/* 	      aaay=gx*SA2 + gy*SB2 + gz*SC2; */
		/* 	      aaaz=gx*SA3 + gy*SB3 + gz*SC3; */

		/* 	      printf("x = %f\n",aaax); */
		/* 	      printf("y = %f\n",aaay); */
		/* 	      printf("z = %f\n\n",aaaz); */

		/* transf do r//z para o sistema ceu - esta é a transf direta*/
		/* 	      printf("gx = %f\n",aaax*SA1 + aaay*SA2 + aaaz*SA3); */
		/* 	      printf("gy = %f\n",aaax*SB1 + aaay*SB2 + aaaz*SB3); */
		/* 	      printf("gz = %f\n\n",aaax*SC1 + aaay*SC2 + aaaz*SC3); */




		/* 	      printf("\ne = %f  q = %f  apocenter = %f  rnow = %f  W = %f  THETA = %f  PHI = %f   VXceu = %f  VYceu = %f  VZceu = %f              spsi = %f\n",e,q,apocenter,rnow,w_dentro*180/M_PI,theta*180/M_PI,phi_dentro*180/M_PI, vv[0], vv[1], vv[2], spsi*180/M_PI); */


		printf("%g,",e);
		printf("%g,",q);
		printf("%g,",rnow);
		printf("%g,",vq);
		printf("%g,",gx);
		printf("%g,",gy);
		printf("%g,",gz);
		printf("%g,",vv[0]);
		printf("%g,",vv[1]);
		printf("%g,",vv[2]);
		if (pre_pos_peri==1) { 	 /* Se for PÓS-PERICENTRO inverto */
		  printf("%g,%g,",-x,y);
		  printf("%g,%g,",vx,-vy);
		}
		if (pre_pos_peri==2) { 	 /* Se for PRE-PERICENTRO fica como foi calculado lá encima */    
		  printf("%g,%g,",x,y);
		  printf("%g,%g,",vx,vy);
		}	
		printf("%g,%g,%g,",spin[0]/spinmod,spin[1]/spinmod,spin[2]/spinmod);
		printf("%g,%g,%g,",qvec[0]/qvecmod,qvec[1]/qvecmod,qvec[2]/qvecmod);
		printf("%g,",vv[2]);    /*  Vsys   */

		/* Testa se a órbita é direta ou retrógrada */
		aux1=spin[0]*spindsk[0]+spin[1]*spindsk[1]+spin[2]*spindsk[2];
		aux=aux1/(sqrt(spin[0]*spin[0] + 
			       spin[1]*spin[1] +
			       spin[2]*spin[2]) *
			  sqrt(spindsk[0]*spindsk[0] +
			       spindsk[1]*spindsk[1] + 
			       spindsk[2]*spindsk[2]));
		aux=acos(aux)*180./M_PI;
		if (aux1 >= 0){
		  printf("Pro,%g,",aux);
		} else {
		  printf("Ret,%g,",aux);
		}

		if (e>0){
		  if (pre_pos_peri==1) { 	 /* Se for PÓS-PERICENTRO */  
		    printf("POS,");
		  }	     
		  if (pre_pos_peri==2) { 	 /* Se for PRE-PERICENTRO */  
		    printf("PRE,");
		  }
		} else { /* orbita circular - sem pericentro... */
		  printf("---,");
		}



		/* Testa e imprime o ângulo entre a posição atual e o pericentro */
		aux1=qvec[0]*gx+qvec[1]*gy+qvec[2]*gz;
		aux=aux1/(sqrt(qvec[0]*qvec[0] + 
			       qvec[1]*qvec[1] +
			       qvec[2]*qvec[2]) *
			  sqrt(gx*gx +
			       gy*gy + 
			       gz*gz));
		aux=acos(aux)*180./M_PI;
		printf("%g,",aux);


		if (e<1){ /* orbita elíptica - escrevo o apocentro... */
		  printf("%g\n",apocenter);
		} else { /* as demais - sem apocentro... */
		  printf("---\n");
		}

	      
	      }   /* Fim do if que testa se Vsys é respeitado                    */
	    }   /* Fim do loop para rotação de V ao redor de RNOW                */
	  }   /* Fim do loop sobre situação pre e pos pericentro               */
	}   /* Fim do if ( (fabs(v) >= fabs(vsys)) ...                         */
      }   /* Fim do loop sobre gz                                              */
    }   /* Fim do loop sobre q                                                 */
  }   /* Fim do loop sobre e                                                   */
  
  
  exit(0);
} /* fim do main */
