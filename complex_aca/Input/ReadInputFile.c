/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                20.12.2008                                   *
\*****************************************************************************/

# include "Includes.h"

long ReadInputFile(void)
{
  char String[4096];
  FILE *InputFile;

  InputFile=fopen("Input.dat","r");
/*
  Surface file
*/
  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%s",SurfaceFile); 
  ReadDimension();
/*
  Minimal cluster size
*/
  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%ld",&NCluMin);
/*
  Admissibility parameters
*/
  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%lf",&ACAEta);

  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%lf",&ACAAlpha);

  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%lf",&ACABeta);

  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%lf",&ACAGamma);

  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%lf",&ACAEps);
/*
  Symmetry
*/
  fgets(String,sizeof(String),InputFile);
  if(strstr(String+40,"Non symmetric") > 0)
    Symmetry=NonSymmetric;
  else if (strstr(String+40,"Complex symmetric") > 0)
    Symmetry=ComplexSymmetric;
  else
    {
      printf("Wrong type of the symmetry !\n");
      exit(1);
    }
/*
  Maximal rank
*/
  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%ld",&MaxRank);
/*
  Example parameters
*/
  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%lf",&Kappa);

  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%lf",&TestEps);
/*
  Environment
*/
  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%s",Prefix); 
  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%s",Suffix); 
  fgets(String,sizeof(String),InputFile);   
  sscanf(String+40,"%s",User);

  fclose(InputFile);

  return 0;
}
