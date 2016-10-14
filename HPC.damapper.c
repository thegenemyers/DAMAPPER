/*********************************************************************************************\
 *
 *  Produce a script to compute overlaps for all block pairs between two DBs, and then sort
 *    and merge *    them into as many .las files as their are blocks of the 1st DB.
 *
 *  Author:  Gene Myers
 *  Date  :  December 31, 2014
 *
 *********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "DB.h"

#undef  LSF  //  define if want a directly executable LSF script

static char *Usage[] =
  { "[-vbpCNz] [-k<int(20)>] [-t<int>] [-M<int>] [-e<double(.85)] [-s<int(100)]",
    "         [-n<double(1.)] [-m<track>]+ [-B<int( 4)>] [-T<int(4)>] [-f<name>]",
    "         <ref:db|dam> <reads:db|dam> [<first:int>[-<last:int>]]"
  };

#define LSF_ALIGN "bsub -q medium -n 4 -o DAMAPPER.out -e DAMAPPER.err -R span[hosts=1] -J align#%d"
#define LSF_SORT "bsub -q short -n 12 -o SORT.MAP.out -e SORT.MAP.err -R span[hosts=1] -J sort#%d"

int main(int argc, char *argv[])
{ int   nblocks2;
  int   useblock2;
  int   usepath1, usepath2;
  int   fblock, lblock;
#ifdef LSF
  int   jobid;
#endif

  FILE *out;
  char  name[100];
  char *pwd1, *root1;
  char *pwd2, *root2;

  int    FORWARD, REVERSE;
  int    BUNIT;
  int    VON, BON, PON, CON, NON;
  int    TINT, KINT, SINT, MINT;
  int    ZON;
  int    NTHREADS;
  double EREL, NEIGHBOR;
  int    MMAX, MTOP;
  char **MASK;
  char  *ONAME;

  { int    i, j, k;         //  Process options
    int    flags[128];
    char  *eptr;

    ARG_INIT("HPC.damapper")

    BUNIT    = 4;
    KINT     = 20;
    TINT     = 0;
    EREL     = 0.;
    SINT     = 100;
    MINT     = -1;
    NEIGHBOR = 1.;
    NTHREADS = 4;

    MTOP = 0;
    MMAX = 10;
    MASK = (char **) Malloc(MMAX*sizeof(char *),"Allocating mask track array");
    if (MASK == NULL)
      exit (1);
    ONAME = NULL;
    out   = stdout;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vbpCNz");
            break;
          case 'e':
            ARG_REAL(EREL)
            if (EREL < .7 || EREL >= 1.)
              { fprintf(stderr,"%s: Average correlation must be in [.7,1.) (%g)\n",Prog_Name,EREL);
                exit (1);
              }
            break;
          case 'k':
            ARG_POSITIVE(KINT,"K-mer length")
            break;
          case 'm':
            if (MTOP >= MMAX)
              { MMAX = 1.2*MTOP + 10;
                MASK = (char **) Realloc(MASK,MMAX*sizeof(char *),"Reallocating mask track array");
                if (MASK == NULL)
                  exit (1);
              }
            MASK[MTOP++] = argv[i]+2;
            break;
          case 'n':
            ARG_REAL(NEIGHBOR)
            if (NEIGHBOR < .7 || NEIGHBOR > 1.)
              { fprintf(stderr,"%s: Near optimal threshold must be in [.7,1.] (%g)\n",
                               Prog_Name,NEIGHBOR);
                exit (1);
              }
            break;
          case 's':
            ARG_POSITIVE(SINT,"Trace spacing")
            break;
          case 't':
            ARG_POSITIVE(TINT,"Tuple suppression frequency")
            break;
          case 'B':
            ARG_NON_NEGATIVE(BUNIT,"Blocks per command")
            break;
          case 'f':
            ONAME = argv[i]+2;
            break;
          case 'M':
            ARG_NON_NEGATIVE(MINT,"Memory allocation (in Gb)")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VON = flags['v'];
    BON = flags['b'];
    PON = flags['p'];
    CON = flags['C'];
    NON = flags['N'];
    ZON = flags['z'];

    if (argc < 3 || argc > 4)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        exit (1);
      }

    if (CON)
      if (NON)
        { FORWARD = 0;
          REVERSE = 1;
        }
      else
        { FORWARD = 1;
          REVERSE = 1;
        }
    else
      if (NON)
        { fprintf(stderr,"%s: Cannot specify both C and N flags together\n",Prog_Name);
          exit (1);
        }
      else
        { FORWARD = 1;
          REVERSE = 0;
        }
    
    for (j = 1; 2*j <= NTHREADS; j *= 2)
      ;
    NTHREADS = j;
  }

  //  Make sure DAM and DB exist and the DB is partitioned, get number of blocks in partition

  pwd1 = PathTo(argv[1]);
  if (strcmp(argv[1]+(strlen(argv[1])-4),".dam") == 0)
    root1 = Root(argv[1],".dam");
  else
    root1 = Root(argv[1],".db");

  { int i, nfiles;
    FILE *dbvis;

    dbvis = fopen(Catenate(pwd1,"/",root1,".dam"),"r");
    if (dbvis == NULL)
      { dbvis = Fopen(Catenate(pwd1,"/",root1,".db"),"r");
        if (dbvis == NULL)
          exit (1);
      }

    if (fscanf(dbvis,"files = %d\n",&nfiles) != 1)
      SYSTEM_ERROR
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];

        if (fgets(buffer,30000,dbvis) == NULL)
          SYSTEM_ERROR
      }

    usepath1 = (strcmp(pwd1,".") != 0);

    fclose(dbvis);
  }

  pwd2 = PathTo(argv[2]);
  if (strcmp(argv[2]+(strlen(argv[2])-4),".dam") == 0)
    root2 = Root(argv[2],".dam");
  else
    root2 = Root(argv[2],".db");

  if (strcmp(root2,root1) == 0 && strcmp(pwd1,pwd2) == 0)
    { fprintf(stderr,"%s: Comparing the same data base %s/%s against itself, use HPCdaligner\n",
                     Prog_Name,pwd1,root1);
      exit (1);
    }

  { int i, nfiles;
    FILE *dbvis;

    dbvis = fopen(Catenate(pwd2,"/",root2,".dam"),"r");
    if (dbvis == NULL)
      { dbvis = Fopen(Catenate(pwd2,"/",root2,".db"),"r");
        if (dbvis == NULL)
          exit (1);
      }

    if (fscanf(dbvis,"files = %d\n",&nfiles) != 1)
      SYSTEM_ERROR
    for (i = 0; i < nfiles; i++)
      { char buffer[30001];

        if (fgets(buffer,30000,dbvis) == NULL)
          SYSTEM_ERROR
      }

    useblock2 = 1;
    if (fscanf(dbvis,"blocks = %d\n",&nblocks2) != 1 || nblocks2 == 1)
      { useblock2 = 0;
        nblocks2  = 1;
      }

    usepath2 = (strcmp(pwd2,".") != 0);

    fclose(dbvis);
  }

  //  Set range fblock-lblock checking that DB.<fblock-1>.las exists & DB.<fblock>.las does not

  { char *eptr, *fptr, *src2;
    FILE *file;

    if (argc == 4)
      { fblock = strtol(argv[3],&eptr,10);
        if ((*eptr != '\0' && *eptr != '-') || eptr <= argv[3])
          { fprintf(stderr,"%s: final argument '%s' does not start with an integer\n",
                           Prog_Name,argv[3]);
            exit (1);
          }
        useblock2 = 1;
        if (*eptr == '-')
          { lblock = strtol(eptr+1,&fptr,10);
            if (*fptr != '\0' || fptr <= eptr+1)
              { fprintf(stderr,"%s: second part of range '%s' is not an integer\n",
                               Prog_Name,eptr+1);
                exit (1);
              }
          }
        else
          lblock = fblock;
        if (fblock < 1 || lblock > nblocks2 || fblock > lblock)
          { fprintf(stderr,"%s: range %d-%d is empty or out of bounds\n",Prog_Name,fblock,lblock);
            exit (1);
          }
      }
    else
      { fblock = 1;
        lblock = nblocks2;
      }

    if (usepath2)
      src2 = Strdup(Catenate(pwd2,"/",root2,""),"Allocating small string!");
    else
      src2 = Strdup(root2,"Allocating small string!");
    if (src2 == NULL)
      exit (1);

    if (fblock > 1)
      { file = fopen(Catenate(src2,Numbered_Suffix(".",fblock-1,"."),root1,".las"),"r");
        if (file == NULL)
          { fprintf(stderr,"%s: File %s.%d.%s.las should already be present!\n",
                           Prog_Name,src2,fblock-1,root1);
            exit (1);
          }
        else
          fclose(file);
      }
    if (useblock2)
      { file = fopen(Catenate(src2,Numbered_Suffix(".",fblock,"."),root1,".las"),"r");
        if (file != NULL)
          { fprintf(stderr,"%s: File %s.%d.%s.las should not yet exist!\n",
                           Prog_Name,src2,fblock,root1);
            exit (1);
          }
      }
    else
      { file = fopen(Catenate(src2,".",root1,".las"),"r");
        if (file != NULL)
          { fprintf(stderr,"%s: File %s.%s.las should not yet exist!\n",
                           Prog_Name,src2,root1);
            exit (1);
          }
      }

    free(src2);
  }

  { int njobs;
    int j, k;

    //  Produce all necessary daligner jobs ...

    if (ONAME != NULL)
      { sprintf(name,"%s.01.OVL",ONAME);
        out = fopen(name,"w");
      }

    njobs = (lblock-fblock)/BUNIT + 1;

    fprintf(out,"# Damapper jobs (%d)\n",njobs);

#ifdef LSF
    jobid = 1;
#endif
    { int bits;
      int low, hgh;

      bits = (lblock-fblock)/BUNIT+1;
      low  = fblock;
      for (j = 1; j <= bits; j++)
        {
#ifdef LSF
          fprintf(out,LSF_ALIGN,jobid++);
          fprintf(out," \"");
#endif
          fprintf(out,"damapper");
          if (VON)
            fprintf(out," -v");
          if (BON)
            fprintf(out," -b");
          if (PON)
            fprintf(out," -p");
          if (CON)
            fprintf(out," -C");
          if (NON)
            fprintf(out," -N");
          if (ZON)
            fprintf(out," -z");
          if (KINT != 20)
            fprintf(out," -k%d",KINT);
          if (TINT > 0)
            fprintf(out," -t%d",TINT);
          if (EREL > 0.)
            fprintf(out," -e%g",EREL);
          if (SINT != 100)
            fprintf(out," -s%d",SINT);
          if (NEIGHBOR > 0.)
            fprintf(out," -n%g",NEIGHBOR);
          if (MINT >= 0)
            fprintf(out," -M%d",MINT);
          if (NTHREADS != 4)
            fprintf(out," -T%d",NTHREADS);
          for (k = 0; k < MTOP; k++)
            fprintf(out," -m%s",MASK[k]);
          if (usepath1)
            fprintf(out," %s/%s",pwd1,root1);
          else
            fprintf(out," %s",root1);
	  hgh = fblock + (((lblock-fblock)+1)*j)/bits;
          for (k = low; k < hgh; k++)
            if (useblock2)
              if (usepath2)
                fprintf(out," %s/%s.%d",pwd2,root2,k);
              else
                fprintf(out," %s.%d",root2,k);
            else
              if (usepath2)
                fprintf(out," %s/%s",pwd2,root2);
              else
                fprintf(out," %s",root2);
#ifdef LSF
          fprintf(out,"\"");
#endif
          fprintf(out,"\n");
          low = hgh;
        }
    }

    //  Check .las (optional)

    if (ONAME != NULL)
      { fclose(out);
        sprintf(name,"%s.02.CHECK.OPT",ONAME);
        out = fopen(name,"w");
      }

    fprintf(out,"# Check all .las files (optional but recommended)\n");

    if (FORWARD)
      for (j = fblock; j <= lblock; j++)
        { fprintf(out,"LAcheck -vS");
          if (usepath2)
            fprintf(out," %s/%s",pwd2,root2);
          else
            fprintf(out," %s",root2);
          if (usepath1)
            fprintf(out," %s/%s",pwd1,root1);
          else
            fprintf(out," %s",root1);
          if (usepath2)
            if (useblock2)
              fprintf(out," %s/%s.%d.%s",pwd2,root2,j,root1);
            else
              fprintf(out," %s/%s.%s",pwd2,root2,root1);
          else
            if (useblock2)
              fprintf(out," %s.%d.%s",root2,j,root1);
            else
              fprintf(out," %s.%s",root2,root1);
          fprintf(out,"\n");
        }
    if (REVERSE)
      for (j = fblock; j <= lblock; j++)
        { fprintf(out,"LAcheck -vS");
          if (usepath1)
            fprintf(out," %s/%s",pwd1,root1);
          else
            fprintf(out," %s",root1);
          if (usepath2)
            fprintf(out," %s/%s",pwd2,root2);
          else
            fprintf(out," %s",root2);
          if (usepath2)
            if (useblock2)
              fprintf(out," %s/%s.%s.%d",pwd2,root1,root2,j);
            else
              fprintf(out," %s/%s.%s",pwd2,root1,root2);
          else
            if (useblock2)
              fprintf(out," %s.%s.%d",root1,root2,j);
            else
              fprintf(out," %s.%s",root1,root2);
          fprintf(out,"\n");
        }

    if (ONAME != NULL)
      fclose(out);
  }

  free(root2);
  free(pwd2);
  free(root1);
  free(pwd1);

  exit (0);
}
