/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*       Deviation particle scheme for the Boltzmann equation in 0D            *
*                                                                             *
*                                08.12.2008                                   *
\*****************************************************************************/
/*
  Definition of cells
*/

typedef struct sCell
{
    Int Number;
    Int *List;
} Cell;
