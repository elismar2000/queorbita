/*
	====== ====== ====== ====== ====== ====== ======

		PROGRAM pot_hbd_din_fric

	====== ====== ====== ====== ====== ====== ======
	      
	 PROGRAM TO RECONSTRUCT THE ORBITS OF A DISK 
	 GALAXY POTENTIAL (Secondary) INTERACTING WITH 
	 A LARGER FIXED GALAXY POTENTIAL (Main) WITH 
         DYNAMICAL FRICTION ACTING BETWEEN THEM.

	 IT IS SUPPOSED THAT THE POTENTIALS HAVE DIFFERENT MASSES:
	 - THE MAIN GALAXY POTENTIAL (POT1) SHOULD BE MORE MASSIVE, 
	 POSITIONED AT THE ORIGIN WITH ZERO VELOCITY. 
	 - FOR THE SECOND (POT2), POSITION AND VELOCITY VECTORS 
	 SHOULD BE INPUT.

	 OUTPUT: THREE DIFFERENT ORBITS ARE PRESENTED:
	 - KEPLERIAN ORBIT
	 - POTENTIAL BASED ORBIT (NO DYN. FRIC.)
	 - ORBIT WITH DYNAMICAL FRICTION

	 IT WORKS EITHER FOR BACKWARD OR FORWARD 
	 INTEGRATION. THE SENSE OF INTEGRATION IS 
	 CONTROLED BY THE SIGNS OF tf AND dt

	====== ====== ====== ====== ====== ====== ======
	Irapuan - 27/08/2013 só halos.
	          02/10/2013 adicionei bojo e disco só parte radial. 
		  07/06/2014 separo a rotina do disco e incluo a parte em z
		             pois parece que o potencial do disco e sua 
			     inclinaćão são determinantes na forma da órbita.
		  21/10/2020 Arrumei a forma de escrever no arquivo de saída para 
     			     padronizar - facilita leitura por rotinas em python 
			     "calcula_orbitas_queorbita.py" e "plota_orbitas_queorbita.py"		     

	Para o mestrado do Fernando Silvério

	gcc -o pot_hbd_din_fric-2 pot_hbd_din_fric-2.c -O3 -lm

*/

/*
        Se comparar a órbita kepleriana com as demais (com e sem fricção), dá
        uma diferença enorme. Por exemplo, para uma órbita kepleriana
        parabólica, as massas das galáxias são consideradas massas
        puntuais. Nas demais, entra a massa acumulada da primária, que é muito
        menor, ou seja, a aceleração será muito menor.

	No programa pot_h_din_fric.c usei somente os potenciais dos
	halos. Comparando com uma órbita simulada, o resultado foi que a
	órbita da simulação, até o pericentro, segue melhor a órbita
	kepleriana do que as órbitas computadas com os potenciais de Hernquist
	(1990). 

	Uma hipótese é que ao considerar somente o halo, as componentes disco
	e bojo, cuja massa se concentra em r pequenos, faz falta. Considerando
	uma massa acumulada menor que a real a órbita calculada aqui fica mais
	aberta do que deveria.

	Por isso, nesta versão, incluo os potenciais do bojo e do disco. Os
	potenciais são:
	    Halo: Hernquist (1990)
	    Bojo: Hernquist (1990)
	    Disco: Hernquist (1993), mas ver a equação 3 de Springel et al. (2005)

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

float G, aHalo, MH, aBulge, aBfrac, MB, aDisk, MD, Z0, M200, R200, M_until_R200, dmfrac, bmfrac, MBH, bhfrac, m_p2, k_p2;   //global variables










/* This is the halo/bulge density function */
float dens_sphere(float r)
{
  extern float  G, aHalo, MH, aBulge, aBfrac, MB, aDisk, MD, Z0, M200, R200, M_until_R200, dmfrac, bmfrac, MBH, bhfrac, m_p2, k_p2;   //global variables
  float dens_s;        
  dens_s = MH*aHalo/(2*M_PI*r*pow((r+aHalo),3));              // Halo
  dens_s = dens_s + MB*aBulge/(2*M_PI*r*pow((r+aBulge),3)); // Bojo

  return dens_s;
}






/* This is the disk density function */
float dens_disk(float r[])
{
  extern float  G, aHalo, MH, aBulge, aBfrac, MB, aDisk, MD, Z0, M200, R200, M_until_R200, dmfrac, bmfrac, MBH, bhfrac, m_p2, k_p2;   //global variables
  float dens_d, R_cil; 

  R_cil = sqrt(r[0]*r[0] + r[2]*r[2]);
  dens_d = MD/(4.*M_PI*aDisk*aDisk*Z0)*exp(-R_cil/aDisk)/pow(cosh(r[1]/Z0),2);  

  // Uso essa:  sech(x)=1/cosh(x)

  return dens_d;
}







/* This is the galaxy cumulative mass */
float cumulative_mass(float r)
{
  extern float  G, aHalo, MH, aBulge, aBfrac, MB, aDisk, MD, Z0, M200, R200, M_until_R200, dmfrac, bmfrac, MBH, bhfrac, m_p2, k_p2;   //global variables
  float maccum;        
  maccum = MH*r*r/pow((r+aHalo),2);         /* Halo (Hernquist, 1990) */
  maccum = maccum * MH / (MH*R200*R200/pow((R200+aHalo),2));  
  /* Isto aqui acima é um artifício para fazer o halo ter a mesma massa
     do modelo gerado com o programa do Volker. O caso é que lá a massa
     do disco é calculada assim: 
     M_HALO = M200 - M_DISK - M_BULGE - M_BLACKHOLE
     e depois as massas das partículas do halo assim:
     mp_halo = M_HALO/N_HALO 
     Isso diz que a massa TOTAL do halo está toda contida até R200. 
     Mas isso não é verdade, pois se calculamos a
     massa do halo acumulada até R200, vai dar menor que M_HALO. Pra
     calcular a massa acumulada, de forma que dê M_HALO, eu tenho que
     escalar a massa, multiplicando por MH / cumulative_mass(R200), e
     essa é a mágica feita aqui... 
  */
  maccum = maccum + MB*r*r/pow((r+aBulge),2);      /* Bojo (Hernquist, 1990) */
  maccum = maccum + MD*(1.-exp(-r/aDisk)*(r/aDisk+1.)); /* Disco (Hernquist, 1993 e Springel, 2005) */
  return maccum;
}









/* Dynamical friction */
void din_fric(float r[], float v[], float a[])
{
  extern float  G, aHalo, MH, aBulge, aBfrac, MB, aDisk, MD, Z0, M200, R200, M_until_R200, dmfrac, bmfrac, MBH, bhfrac, m_p2, k_p2;   //global variables
  float rabs, vabs, Vrot, maccum, friction, dens, X;
  float COULOMB_S, COULOMB_D;

  //  float dens_sphere(float );
  //  float dens_disk(float );
  
  rabs = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  vabs = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

  if(rabs > aHalo*1.e-30) {                /* Para evitar problemas
					   quando um pot se aproximam
					   muito do centro do halo */

    maccum = cumulative_mass(rabs);

    Vrot = sqrt(G*maccum/rabs);
    //    Vrot = sqrt(G*M200/rabs); //testando
    X=vabs/Vrot;
   
    //    X=vabs/sqrt(2)/190.0/1.022;
    //    X=vabs/sqrt(2)/0.1;
    dens = dens_sphere(rabs) + dens_disk(r);
  


/*   COULOMB=log(rabs/aHalo); */
    COULOMB_S=2.4;             /* valor aconselhado por Taylor & Babul */
    COULOMB_D=0.5;             /* valor aconselhado por Taylor & Babul */

    // COULOMB_S=log(rabs/k_p2);

    friction = -4.*M_PI*G*G*m_p2/pow(vabs,3.);
    friction = friction*(erf(X) - 2.*X*exp(-X*X)/sqrt(M_PI));
    friction = friction*(dens_sphere(rabs)*COULOMB_S + dens_disk(r)*COULOMB_D);

    // friction=friction*0.3;
  } else {
    maccum = 0.0;
    Vrot = 0.0;
    X = 0.0;
    friction = 0.0;
  }

  printf("  HALO: ===>   COULOMB_S= %g rabs = %g  friction = %g\n", COULOMB_S, rabs, friction); 
     printf("        ===>  maccum = %g Vrot = %g  dens = %g  X = %g\n", maccum, Vrot, dens, X); 

  a[0] = friction*v[0];
  a[1] = friction*v[1];
  a[2] = friction*v[2];

     printf("        ===>  a[0] =  %g\n", a[0] ); 

}






/* Accelerations */
void accel(float r1[], float v1[], float a1[], int caso)
{
  /* 
     ENTRA: posição, velocidade, e massa da SECUNDÁRIA
     ÚLTIMO PARÂMETRO: caso - pode ser 1, 2 ou 3:
                       1=com fricção
		       2=sem fricção
		       3=kepleriano
     SAI: Aceleração na SECUNDÁRIA
  */
  float rabs, maccum;
  float ac_gal;        /* termo da Galaxia para potencial de Hernquist*/
  float ac_k;        /* termo da Galaxia para orbita kepleriana*/

  rabs = sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);

  /* Cálculo da aceleração entre a GALÁXIA PRINCIPAL e a SECUNDÁRIA. Aqui
     entra o cálculo da massa acumulada do potencial da PRIMÁRIA: maccum.
  */
  if(rabs > aHalo*1.e-30) {                /* Para evitar problemas
					   quando um pot se aproximam
					   muito do centro do halo */

    maccum = cumulative_mass(rabs);

    ac_gal = - G*maccum/(rabs*rabs)/rabs;  /* Calcula a aceleração por: 

					           GM    
					    a = - ---- * (vet.unit. da direção r)
					           r*r

					    Por isso o negocio é divivido por
					    r*r*r. as componentes do vetor r
					    serão multiplicadas depois.
					 */



  } else {
    maccum = 0.0;
    ac_gal = 0.0; 
  }


  a1[0] = 0.0;
  a1[1] = 0.0;
  a1[2] = 0.0;


  /* Sem fricção */
  if (caso==2){
    a1[0] = ac_gal*r1[0]  ;
    a1[1] = ac_gal*r1[1]  ;
    a1[2] = ac_gal*r1[2]  ;
  }
  /* Com fricção */
  if (caso==1){
    din_fric(r1, v1, a1);       /* termo da friccao dinamica:
				     sai o vetor a1 */
    printf("         ACCEL: ===>  maccum = %g ac_gal_x = %f  a_fric_x = %f\n", maccum, ac_gal*r1[0], a1[0]); 
    a1[0] = a1[0] + ac_gal*r1[0]  ;
    a1[1] = a1[1] + ac_gal*r1[1]  ;
    a1[2] = a1[2] + ac_gal*r1[2]  ;
  }

  /* Kepleriana */
  if (caso==3){
    ac_k = - G*M200/(rabs*rabs)/rabs;
    //ac_k = - G*M_until_R200/(rabs*rabs)/rabs; 
    a1[0] = ac_k*r1[0] ;
    a1[1] = ac_k*r1[1] ;
    a1[2] = ac_k*r1[2] ;
  }
}



main()
{
  extern float  G, aHalo, MH, aBulge, aBfrac, MB, aDisk, MD, Z0, M200, R200, M_until_R200, dmfrac, bmfrac, MBH, bhfrac, m_p2, k_p2;   //global variables
  int i,j;

  float r_nofric[3], r_kepler[3], r_p2[3],  
    v_nofric[3], v_kepler[3], v_p2[3], 
    a_nofric[3], a_kepler[3], a_p2[3];
  
  /* r_kepler, v_kepler, a_kepler : órbita kepleriana usando as massas totais 
                                    como massas puntuais. 
     r_nofric, v_nofric, a_nofric : órbita entre os potenciais de Hernquist, 
                                    mas sem fricção dinâmica
     r_p2, v_p2, a_p2             : órbita com fricção dinâmica
  */

  float rscale, vscale;                /* fatores de escala em r e v */
  float t, tf, dt;

  //  float dens_gal(float ); 

  float rabs,maccum;

  FILE *fp_out;


/*	=========================================
	READ INPUT PARAMETERS

	r_p2[t=0][x,y,z] = position of Pot2
	m_p2 = Pot2 mass
	k_p2 = Pot2 radial extension (aHalo)
*/
  puts("Program to reconstruct the orbit of a pair of galaxies (haloes");
  puts("only) with dynamical friction acting between them");
  puts("");
  puts("It is assumed that the potentials have different masses. ");
  puts("The MAIN GALAXY potential is the most massive, positioned ");
  puts("at the origin, with zero velocity. Both, MAIN and SECONDARY");
  puts("galaxies are represented by the same potentials used in");
  puts("Springel's MakeNewDisk code. ");
  puts("------------------------");
  puts("All the input data must be entered in the following physical system of units");
  puts("\tDistances in kpc;");
  puts("\tVelocities in km/s;");
  puts("\tMasses in Solar Masses;");
  puts("\tTime in Myr (t=0 reffers to the present time)");
  puts("------------------------");

  /* All the operations will be done internally in the following
     system of units:
           G=4.49828E-3 pc^3/Msol/Myr^2
	   distances in pc
	   velocities in pc/Myr
	   masses in solar masses
	   time in Myr
     So now we define the scale factors to convert distances and
     velocities from the input data system of units to the one 
     defined above. The conversion will be done as soon as the data
     are entered.
  */
  G=6.672e-8 * 1.989e33 / pow(3.08567802e18,3) * pow(3.15576e13,2);

  printf("G=%g\n",G);
  rscale=1000.;       // entrada em kpc, converte pra pc
  vscale=3.15576e+13*100000/3.08567802e18;  // entrada em km/s converte pra pc/Myr (vscale = 1.0227120197070982)

  puts("Input parameters at t=0");
  puts("-----------------------");
  puts("\tPotential 1 - MAIN GALAXY (Hernquist, 1990): ");
  puts("\t  MAIN Galaxy position at the origin: (0,0,0)");
  puts("\t  MAIN Galaxy velocity: (0,0,0). ");
  printf("\t  MAIN Galaxy Total mass (M200, MakeNewDisk output) [M_sun]: ");
  scanf("%f",&M200);
  printf("\t  MAIN Galaxy R200 (MakeNewDisk output) [kpc]: ");
  scanf("%f",&R200);
  R200 = R200 * rscale;

  printf("\t  HALO scalelength (aHalo, MakeNewDisk output - RH) [kpc]: ");
  scanf("%f",&aHalo);
  aHalo = aHalo * rscale;

  printf("\t  DISK mass fraction (dmfrac): ");
  scanf("%f",&dmfrac);
  MD = dmfrac*M200;
  printf("\t  DISK scalelength (aDisk, end of MakeNewDisk output) [kpc]: ");
  scanf("%f",&aDisk);
  aDisk = aDisk * rscale;
  printf("\t  DISK vertical scalelength (DiskHeight, from MakeNewDisk parameters file) [kpc]: ");
  scanf("%f",&Z0);
  Z0 = Z0 * aDisk;


  printf("\t  BULGE mass fraction (bmfrac, from MakeNewDisk parameters file - MD): ");
  scanf("%f",&bmfrac);
  MB = bmfrac*M200;
  printf("\t  BULGE scalelength fraction (aBulge, from MakeNewDisk  ");
  printf("\t                           parameters file - BulgeSize): ");
  scanf("%f",&aBfrac);
  aBulge = aBfrac*aDisk;

  printf("\t  BLACK HOLE mass fraction (bhfrac, from MakeNewDisk parameters file - MBH): ");
  scanf("%f",&bhfrac);
  MBH = bhfrac*M200;

  MH = M200 - MD - MB -MBH;

  M_until_R200 = cumulative_mass(R200); /* Calcula a massa REAL do
					   modelo até R200.  Isso dá
					   bem diferente do M200, por
					   causa da massa do halo: MH
					   é a massa total do halo,
					   que em R200 ainda não
					   terminou. Isso vai ser
					   usado no cálculo da
					   aceleração , caso
					   kepleriano, que pega a
					   massa total. Eu tava usando
					   M200 mas tava dando
					   discrepâncias grandes com
					   as outras órbitas. */


  printf("\n\n  MASSAS: M200=%g MH=%g MD=%g MB=%g MBH=%g\n\n", M200, MH, MD, MB, MBH);

  puts("-----------------------");
  puts("\tPotential 2 - SECONDARY GALAXY (Hernquist, 1990): ");

  printf("\t  Pot 2 Position (x,y,z) [kpc]: ");
  scanf("%f %f %f",&r_p2[0],&r_p2[1],&r_p2[2]);
  r_p2[0] = r_p2[0] * rscale;
  r_p2[1] = r_p2[1] * rscale;
  r_p2[2] = r_p2[2] * rscale;

  printf("\t  Pot 2 Velocity (vx,vy,vz) [km/s]: ");
  scanf("%f %f %f",&v_p2[0],&v_p2[1],&v_p2[2]);
  v_p2[0] = v_p2[0] * vscale;
  v_p2[1] = v_p2[1] * vscale;
  v_p2[2] = v_p2[2] * vscale;

  /* Initialize the other orbital quantities, for non-friction and kepler orbits */
  r_nofric[0] = r_p2[0]; r_nofric[1] = r_p2[1]; r_nofric[2] = r_p2[2];
  v_nofric[0] = v_p2[0]; v_nofric[1] = v_p2[1]; v_nofric[2] = v_p2[2];
  r_kepler[0] = r_p2[0]; r_kepler[1] = r_p2[1]; r_kepler[2] = r_p2[2];
  v_kepler[0] = v_p2[0]; v_kepler[1] = v_p2[1]; v_kepler[2] = v_p2[2];



  printf("\t  Pot 2 mass [M_sun]: ");
  scanf("%f",&m_p2);

  printf("\t  Pot 2 radial scalelength [kpc]: ");
  scanf("%f",&k_p2);
  k_p2 = k_p2 * rscale;

 
  puts("");
  puts("\tTime integration starts at t=0 [Myr]. ");
  printf("\tEnter the final time (tf [Myr]): ");
  scanf("%f",&tf);
  
  printf("\tEnter the timestep (dt [Myr]): ");
  scanf("%f", &dt);

  if ((tf < 0.0 && dt > 0) || (tf > 0.0 && dt < 0) || tf == 0.0 || dt == 0 ) { 
    puts("\n\n\n\nERROR:   Time step and final time must have the same sign. Please try again.\n\n");
    exit(1);
  }

  fp_out = fopen("orbits_pot_hbd_din_fric.dat","w");  

  puts("");
  puts("");
  puts("");
  puts("");
  puts("-----------------------");


	
  /* Now that we have the initial conditions, we proceed the
     time integration */
  

  fprintf(fp_out,"#\n");
  fprintf(fp_out,"#     	--------Com Fricção Dinâmica--------------	--------Sem Fricção Dinâmica--------------	------------ Kepleriana-------------------\n");
  fprintf(fp_out,"time\tx_f\ty_f\tz_f\tVx_f\tVy_f\tVz_f\tx_nf\ty_nf\tz_nf\tVx_nf\tVy_nf\tVz_nf\tx_kepl\ty_kepl\tz_kepl\tVx_kepl\tVy_kepl\tVz_kepl\n");
  fprintf(fp_out,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
	  t=0.0,
	  r_p2[0]/rscale, r_p2[1]/rscale, r_p2[2]/rscale, 
	  v_p2[0]/vscale, v_p2[1]/vscale, v_p2[2]/vscale,
	  r_nofric[0]/rscale, r_nofric[1]/rscale, r_nofric[2]/rscale,   // sem fricção
	  v_nofric[0]/vscale, v_nofric[1]/vscale, v_nofric[2]/vscale,
	  r_kepler[0]/rscale, r_kepler[1]/rscale, r_kepler[2]/rscale,   // Kepleriana
	  v_kepler[0]/vscale, v_kepler[1]/vscale, v_kepler[2]/vscale);

	  
  /* Here we start with the Leapfrog integrator */
  /* ATENÇÂO: uso fabs em t e tf abaixo para que o programa funcione forward and backward */
  
  for(t=0.0;fabs(t)<=fabs(tf);t+=dt){
    //    printf("   t=%f   dt=%f   tf=%f\n",t,dt,tf);
    for(i=0;i<=2;i++){         /* compute r(t+0.5*dt) from r(t) */
      r_p2[i]=r_p2[i] + 0.5*dt*v_p2[i];
      r_nofric[i]=r_nofric[i] + 0.5*dt*v_nofric[i];
      r_kepler[i]=r_kepler[i] + 0.5*dt*v_kepler[i];
    }
    
    /* the accelerations a(t+0.5*dt) on POT2 and POT3 from
       r(t+0.5*dt) and v(t) */
    accel(r_p2, v_p2, a_p2, 1);                       /* Último parâmetro = Caso 1 
							 Sai: a_p2 (com fricção)
						      */
    accel(r_nofric, v_nofric, a_nofric, 2);           /* Último parâmetro = Caso 2 
							 Sai: a_nofric (sem fricção)
						      */
    accel(r_kepler, v_kepler, a_kepler, 3);           /* Último parâmetro = Caso 3 
							 Sai: a_kepler (kepleriana)
						      */
    
    for(i=0;i<=2;i++){         /* compute v(t+dt) from a(t+0.5*dt) */
      v_p2[i]=v_p2[i] + dt*a_p2[i];
      v_nofric[i]=v_nofric[i] + dt*a_nofric[i];
      v_kepler[i]=v_kepler[i] + dt*a_kepler[i];
    }
    
    for(i=0;i<=2;i++){         /* compute r(t+dt) from
				  r(t+0.5*dt) and v(t+dt) */
      r_p2[i]=r_p2[i] + 0.5*dt*v_p2[i];
      r_nofric[i]=r_nofric[i] + 0.5*dt*v_nofric[i];
      r_kepler[i]=r_kepler[i] + 0.5*dt*v_kepler[i];
    }
    
    
    /* Send to the output file */  
    fprintf(fp_out,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
	    t+dt,
	    r_p2[0]/rscale, r_p2[1]/rscale, r_p2[2]/rscale, 
	    v_p2[0]/vscale, v_p2[1]/vscale, v_p2[2]/vscale,
	    r_nofric[0]/rscale, r_nofric[1]/rscale, r_nofric[2]/rscale,   // sem fricção
	    v_nofric[0]/vscale, v_nofric[1]/vscale, v_nofric[2]/vscale,
	    r_kepler[0]/rscale, r_kepler[1]/rscale, r_kepler[2]/rscale,   // Kepleriana
	    v_kepler[0]/vscale, v_kepler[1]/vscale, v_kepler[2]/vscale);
  } /* End of the loop over t */
		
  puts("-----------------------");
	
  puts("Massa TOTAL acumulada em função de r");
	

  for(rabs=0.;rabs<=160000.;rabs=rabs+1000.){
    maccum = cumulative_mass(rabs);

    printf("\t\t%g\t%g\n",rabs/rscale,maccum);
  }

  fclose(fp_out);


  exit(0);  
} /* End of main */

