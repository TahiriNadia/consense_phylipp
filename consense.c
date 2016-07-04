#include "phylip.h"
#include "cons.h"

/* version 3.6. (c) Copyright 1993-2008 by the University of Washington.
   Written by Joseph Felsenstein, Hisashi Horino,
   Akiko Fuseki, Dan Fineman, Sean Lamont, and Andrew Keeffe.
   Permission is granted
   to copy and use this program provided no fee is charged for it and
   provided that this copyright notice is not removed. */


/* The following extern's refer to things declared in cons.c */

extern int tree_pairing;

extern Char outfilename[FNMLNGTH], intreename[FNMLNGTH], intree2name[FNMLNGTH], outtreename[FNMLNGTH];
extern node *root;

extern long numopts, outgrno, col;
extern long maxgrp;               /* max. no. of groups in all trees found  */

extern boolean trout, firsttree, noroot, outgropt, didreroot, prntsets,
          progress, treeprint, goteof, strict, mr, mre, ml;
extern pointarray nodep;                 /* pointers to all nodes in tree */
extern group_type **grouping, **grping2, **group2;/* to store groups found  */
extern long **order, **order2, lasti;
extern group_type *fullset;
extern long tipy;

extern double trweight, ntrees, mlfrac;

#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   count_siblings(node **p);
void   treeout(node *, char [], int *);
int main_consense(char treeConsense [], char * treesNewick [],int trees_in,int tip_count);
/* function prototypes */
#endif


void getoptions()
{
  /* interactively set options */
  long loopcount, loopcount2;
  Char ch = 'Y';
  boolean done, done1;
  done = true;

  /* Initial settings */
  ibmpc          = IBMCRT;
  ansi           = ANSICRT;
  didreroot      = false;
  firsttree      = true;
  spp            = 0 ;
  col            = 0 ;

  /* This is needed so functions in cons.c work */
  tree_pairing   = NO_PAIRING ;

/*   fprintf(outfile, "\nConsensus tree");
  fprintf(outfile, " program, version %s\n\n", VERSION); */
  putchar('\n');
  strict = false;
  mr = true;
  mre = false;
  ml = false;
  mlfrac = 0.5;
  noroot = true;
  numopts = 0;
  outgrno = 1;
  outgropt = false;
  trout = true;
  prntsets = true;
  progress = true;
  treeprint = true;
  loopcount = 0;

}  /* getoptions */


void count_siblings(node **p)
{
  node *tmp_node;
  int i;

  if (!(*p)) {
    /* This is a leaf, */
    return;
  } else {
    tmp_node = (*p)->next;
  }

  for (i = 0 ; i < 1000; i++) {
    if (tmp_node == (*p)) {
      /* When we've gone through all the siblings, */
      break;
    } else if (tmp_node) {
      tmp_node = tmp_node->next;
    } else  {
      /* Should this be executed? */
      return ;
    }
  }
} /* count_siblings */

void treeout(node *p, char treeConsense[], int * compteur)
{
  /* write out file with representation of final tree */
  long i, n = 0;
  char c;
  node *q;
  double x;

  count_siblings (&p);

  if (p->tip) {
    /* If we're at a node which is a leaf, figure out how long the
       name is and print it out. */
    for (i = 1; i <= MAXNCH; i++) {
      if (p->nayme[i - 1] != '\0')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = p->nayme[i];
      if (c == ' ')
        c = '_';
      treeConsense[*compteur]=c;
	  (*compteur)++;
	  //putc(c, outtree);
    }
    col += n;
  } else {
    /* If we're at a furcation, print out the proper formatting, loop
       through all the children, calling the procedure recursively. */
    treeConsense[*compteur]='(';
	(*compteur)++;
	//putc('(', outtree);
    col++;
    q = p->next;
    while (q != p) {
      /* This should terminate when we've gone through all the
         siblings, */
      treeout(q->back,treeConsense,compteur);
      q = q->next;
      if (q == p)
        break;
      treeConsense[*compteur]=',';
	  (*compteur)++;
	  //putc(',', outtree);
      col++;
      /*if (col > 60) {
      treeConsense[*compteur]='\n';
	  (*compteur)++;
	  //putc('\n', outtree);
        col = 0;
      }*/
    }
    treeConsense[*compteur]=')';
	(*compteur)++;
	//putc(')', outtree);
    col++;
  }

  if (p->tip)
    x = ntrees;
  else
    x = (double)p->deltav;

  if (p == root) {
    /* When we're all done with this tree, */
    treeConsense[*compteur]=';';
	(*compteur)++;
	treeConsense[*compteur]='\n';
	(*compteur)++;
	//fprintf(outtree, ";\n");
    return;
  }

  /* Figure out how many characters the branch length requires: */
  else {
    if (!strict) {
      if (x >= 100.0) {
/* 		char arr[sizeof(x)];
		memcpy(&arr,&x,sizeof(x));
		printf("1) double = %d\n",x);
		printf("1) char * = %s\n",arr); */

		treeConsense[*compteur]=':';
		(*compteur)++;
		treeConsense[*compteur]='1';
		(*compteur)++;
		treeConsense[*compteur]='0';
		(*compteur)++;
		treeConsense[*compteur]='0';
		(*compteur)++;
		treeConsense[*compteur]='.';
		(*compteur)++;
		treeConsense[*compteur]='0';
		(*compteur)++;

        //fprintf(outtree, ":%5.1f", x);
        col += 4;
      } else if (x >= 10.0) {
/* 		  	char arr[sizeof(x)];
			memcpy(&arr,&x,sizeof(x));
			printf("2) double = %d\n",x);
			printf("2) char * = %s\n",arr); */

			treeConsense[*compteur]=':';
			(*compteur)++;
			treeConsense[*compteur]='1';
	        (*compteur)++;
			treeConsense[*compteur]='0';
	        (*compteur)++;
			treeConsense[*compteur]='.';
	        (*compteur)++;
			treeConsense[*compteur]='0';
	        (*compteur)++;

          //fprintf(outtree, ":%4.1f", x);
          col += 3;
        } else if (x >= 1.00) {
/* 			char arr[sizeof(x)];
			memcpy(&arr,&x,sizeof(x));
			printf("3) double = %d\n",x);
			printf("3) char * = %s\n",arr); */

			treeConsense[*compteur]=':';
			(*compteur)++;
			treeConsense[*compteur]='1';
	        (*compteur)++;
			treeConsense[*compteur]='.';
	        (*compteur)++;
			treeConsense[*compteur]='0';
	        (*compteur)++;
			treeConsense[*compteur]='0';
	        (*compteur)++;

            //fprintf(outtree, ":%4.2f", x);
            col += 3;
          }
    }
  }

}  /* treeout */

int main_consense(char treeConsense [], char * treesNewick [],int trees_in,int tip_count)
//int main(/* int argc, Char *argv[] */)
{
  /* Local variables added by Dan F. */
  pattern_elm  ***pattern_array;
  long i, j;
  node *p, *q;
  /*long trees_in = 0;
  long tip_count = 0;
  char treeConsense[100000];*/
  int compteur = 0;
  
  /*char * treesNewick [] ={"(4:0.6660,((1:1.2234,(10:0.8667,(13:0.7613,(2:0.6189,8:0.8176):1.4461):0.6259):0.5966):0.8599,((3:1.6857,7:0.6000):0.7439,(6:1.9101,15:1.5985):1.2384):1.3035):0.7627,(11:1.1495,((5:0.7647,9:0.6086):1.3284,(12:0.8095,(14:0.6788,16:0.8620):1.3349):0.8865):0.8509):1.4008);",
 "(4:0.6660,((2:1.2234,(10:0.8667,(13:0.7613,(1:0.6189,8:0.8176):1.4461):0.6259):0.5966):0.8599,((3:1.6857,7:0.6000):0.7439,(6:1.9101,15:1.5985):1.2384):1.3035):0.7627,(11:1.1495,((5:0.7647,9:0.6086):1.3284,(12:0.8095,(14:0.6788,16:0.8620):1.3349):0.8865):0.8509):1.4008);",
 "(13:0.8785,(4:1.4093,((10:0.9596,(1:1.3741,(6:1.0037,15:1.1345):1.6311):0.6104):0.7605,((7:0.7134,14:0.8217):0.8077,(9:0.7243,16:2.2580):0.8752):0.7479):0.7527):0.6755,(11:1.0347,(3:1.0933,(2:0.7588,(8:0.9821,(5:0.6977,12:0.6440):1.3816):1.3430):1.1131):0.8067):1.0069);",
 "(13:0.8785,(4:1.4093,((10:0.9596,(1:1.3741,(6:1.0037,15:1.1345):1.6311):0.6104):0.7605,((7:0.7134,14:0.8217):0.8077,(9:0.7243,16:2.2580):0.8752):0.7479):0.7527):0.6755,(11:1.0347,(12:1.0933,(2:0.7588,(8:0.9821,(5:0.6977,3:0.6440):1.3816):1.3430):1.1131):0.8067):1.0069);",
 "((10:0.8447,(6:3.1321,(1:0.9276,(16:1.2873,(11:0.5728,(2:0.9222,15:0.8928):1.8756):0.5969):0.7027):0.8222):0.5830):0.9111,(4:0.5815,((9:0.6427,13:0.9100):0.6754,(14:0.8118,(5:0.8935,7:1.0046):1.6148):1.1266):0.9732):1.1149,(12:1.1573,(3:1.0087,8:0.7018):0.7171):0.9952);",
 "((10:0.8447,(6:3.1321,(2:0.9276,(16:1.2873,(11:0.5728,(1:0.9222,15:0.8928):1.8756):0.5969):0.7027):0.8222):0.5830):0.9111,(4:0.5815,((9:0.6427,13:0.9100):0.6754,(14:0.8118,(5:0.8935,7:1.0046):1.6148):1.1266):0.9732):1.1149,(12:1.1573,(3:1.0087,8:0.7018):0.7171):0.9952);",
 "((1:0.8447,(6:3.1321,(10:0.9276,(16:1.2873,(11:0.5728,(2:0.9222,15:0.8928):1.8756):0.5969):0.7027):0.8222):0.5830):0.9111,(4:0.5815,((9:0.6427,13:0.9100):0.6754,(14:0.8118,(5:0.8935,7:1.0046):1.6148):1.1266):0.9732):1.1149,(12:1.1573,(3:1.0087,8:0.7018):0.7171):0.9952);",
 "(10:0.8316,((2:0.6145,(5:0.6205,6:0.9002):1.2925):1.0065,(16:1.6438,((7:1.2011,(15:0.8131,(1:0.6308,11:0.6147):0.9774):1.0724):0.6522,(12:1.7871,14:1.3660):1.6932):1.4504):1.0451):0.6352,(4:0.7807,(3:1.4190,(9:0.7055,(8:0.7526,13:0.6769):1.0528):0.7920):1.3507):0.6212);",
 "(10:0.8316,((2:0.6145,(5:0.6205,6:0.9002):1.2925):1.0065,(16:1.6438,((1:1.2011,(15:0.8131,(7:0.6308,11:0.6147):0.9774):1.0724):0.6522,(12:1.7871,14:1.3660):1.6932):1.4504):1.0451):0.6352,(4:0.7807,(3:1.4190,(9:0.7055,(8:0.7526,13:0.6769):1.0528):0.7920):1.3507):0.6212);",
 "(10:0.8316,((4:0.6145,(5:0.6205,6:0.9002):1.2925):1.0065,(16:1.6438,((7:1.2011,(15:0.8131,(1:0.6308,11:0.6147):0.9774):1.0724):0.6522,(12:1.7871,14:1.3660):1.6932):1.4504):1.0451):0.6352,(2:0.7807,(3:1.4190,(9:0.7055,(8:0.7526,13:0.6769):1.0528):0.7920):1.3507):0.6212);"};
*/
 #ifdef MAC
  argc = 1;                /* macsetup("Consense", "");        */
  /* argv[0] = "Consense"; */
#endif
  /* init(argc, argv); */

  /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
  //openfile(&intree, INTREE, "input tree file", "rb", "test", intreename);
  //openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  /* Initialize option-based variables, then ask for changes regarding
     their values. */
  getoptions();

  ntrees = 0.0;
  maxgrp = 32767;   /* initial size of set hash table */
  lasti  = -1;

/*   if (trout)
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename); */
  /* if (prntsets)
    fprintf(outfile, "Species in order: \n\n"); */

  //trees_in = countsemic(&intree);
  //trees_in is number of tree
  trees_in = 10;
  //countcomma(&intree,&tip_count);
  //tip_count++; /* countcomma does a raw comma count, tips is one greater */
  //tip_count is number of species
  tip_count = 16;

  //char ** treesNewick;
  char treesNewick1 [trees_in][maxgrp];
  strcpy(treesNewick1[0], "(4:0.6660,((1:1.2234,(10:0.8667,(13:0.7613,(2:0.6189,8:0.8176):1.4461):0.6259):0.5966):0.8599,((3:1.6857,7:0.6000):0.7439,(6:1.9101,15:1.5985):1.2384):1.3035):0.7627,(11:1.1495,((5:0.7647,9:0.6086):1.3284,(12:0.8095,(14:0.6788,16:0.8620):1.3349):0.8865):0.8509):1.4008);");
  strcpy(treesNewick1[1], "(4:0.6660,((2:1.2234,(10:0.8667,(13:0.7613,(1:0.6189,8:0.8176):1.4461):0.6259):0.5966):0.8599,((3:1.6857,7:0.6000):0.7439,(6:1.9101,15:1.5985):1.2384):1.3035):0.7627,(11:1.1495,((5:0.7647,9:0.6086):1.3284,(12:0.8095,(14:0.6788,16:0.8620):1.3349):0.8865):0.8509):1.4008);");
  strcpy(treesNewick1[2], "(13:0.8785,(4:1.4093,((10:0.9596,(1:1.3741,(6:1.0037,15:1.1345):1.6311):0.6104):0.7605,((7:0.7134,14:0.8217):0.8077,(9:0.7243,16:2.2580):0.8752):0.7479):0.7527):0.6755,(11:1.0347,(3:1.0933,(2:0.7588,(8:0.9821,(5:0.6977,12:0.6440):1.3816):1.3430):1.1131):0.8067):1.0069);");
  strcpy(treesNewick1[3], "(13:0.8785,(4:1.4093,((10:0.9596,(1:1.3741,(6:1.0037,15:1.1345):1.6311):0.6104):0.7605,((7:0.7134,14:0.8217):0.8077,(9:0.7243,16:2.2580):0.8752):0.7479):0.7527):0.6755,(11:1.0347,(12:1.0933,(2:0.7588,(8:0.9821,(5:0.6977,3:0.6440):1.3816):1.3430):1.1131):0.8067):1.0069);");
  strcpy(treesNewick1[4], "((10:0.8447,(6:3.1321,(1:0.9276,(16:1.2873,(11:0.5728,(2:0.9222,15:0.8928):1.8756):0.5969):0.7027):0.8222):0.5830):0.9111,(4:0.5815,((9:0.6427,13:0.9100):0.6754,(14:0.8118,(5:0.8935,7:1.0046):1.6148):1.1266):0.9732):1.1149,(12:1.1573,(3:1.0087,8:0.7018):0.7171):0.9952);");
  strcpy(treesNewick1[5], "((10:0.8447,(6:3.1321,(2:0.9276,(16:1.2873,(11:0.5728,(1:0.9222,15:0.8928):1.8756):0.5969):0.7027):0.8222):0.5830):0.9111,(4:0.5815,((9:0.6427,13:0.9100):0.6754,(14:0.8118,(5:0.8935,7:1.0046):1.6148):1.1266):0.9732):1.1149,(12:1.1573,(3:1.0087,8:0.7018):0.7171):0.9952);");
  strcpy(treesNewick1[6], "((1:0.8447,(6:3.1321,(10:0.9276,(16:1.2873,(11:0.5728,(2:0.9222,15:0.8928):1.8756):0.5969):0.7027):0.8222):0.5830):0.9111,(4:0.5815,((9:0.6427,13:0.9100):0.6754,(14:0.8118,(5:0.8935,7:1.0046):1.6148):1.1266):0.9732):1.1149,(12:1.1573,(3:1.0087,8:0.7018):0.7171):0.9952);");
  strcpy(treesNewick1[7], "(10:0.8316,((2:0.6145,(5:0.6205,6:0.9002):1.2925):1.0065,(16:1.6438,((7:1.2011,(15:0.8131,(1:0.6308,11:0.6147):0.9774):1.0724):0.6522,(12:1.7871,14:1.3660):1.6932):1.4504):1.0451):0.6352,(4:0.7807,(3:1.4190,(9:0.7055,(8:0.7526,13:0.6769):1.0528):0.7920):1.3507):0.6212);");
  strcpy(treesNewick1[8], "(10:0.8316,((2:0.6145,(5:0.6205,6:0.9002):1.2925):1.0065,(16:1.6438,((1:1.2011,(15:0.8131,(7:0.6308,11:0.6147):0.9774):1.0724):0.6522,(12:1.7871,14:1.3660):1.6932):1.4504):1.0451):0.6352,(4:0.7807,(3:1.4190,(9:0.7055,(8:0.7526,13:0.6769):1.0528):0.7920):1.3507):0.6212);");
  strcpy(treesNewick1[9], "(10:0.8316,((4:0.6145,(5:0.6205,6:0.9002):1.2925):1.0065,(16:1.6438,((7:1.2011,(15:0.8131,(1:0.6308,11:0.6147):0.9774):1.0724):0.6522,(12:1.7871,14:1.3660):1.6932):1.4504):1.0451):0.6352,(2:0.7807,(3:1.4190,(9:0.7055,(8:0.7526,13:0.6769):1.0528):0.7920):1.3507):0.6212);");


  /* Read the tree file and put together grouping, order, and timesseen */
  read_groups (&pattern_array,trees_in, tip_count,treesNewick);
  /* Compute the consensus tree. */
  //putc('\n', outfile);
  nodep      = (pointarray)Malloc(2*(1+spp)*sizeof(node *));
  for (i = 0; i < spp; i++) {
    nodep[i] = (node *)Malloc(sizeof(node));
    for (j = 0; j < MAXNCH; j++)
      nodep[i]->nayme[j] = '\0';
    strncpy(nodep[i]->nayme, nayme[i], MAXNCH);
  }
  for (i = spp; i < 2*(1+spp); i++)
    nodep[i] = NULL;
  consensus(pattern_array, trees_in);
  printf("\n");
  if (trout) {
    treeout(root,treeConsense,&compteur);
/*     if (progress)
      printf("Consensus tree written to file \"%s\"\n\n", outtreename); */
  }

/* 	treeConsense[compteur]=')';
	compteur++;
	treeConsense[compteur]=';';
	compteur++; */

	  printf("TREE CONSENUSE %s\n",treeConsense);
/*   if (progress)
    printf("Output written to file \"%s\"\n\n", outfilename); */
  for (i = 0; i < spp; i++)
    free(nodep[i]);
  for (i = spp; i < 2*(1 + spp); i++) {
    if (nodep[i] != NULL) {
      p = nodep[i]->next;
      do {
        q = p->next;
        free(p);
        p = q;
      } while (p != nodep[i]);
      free(p);
    }
  }
  free(nodep);
  FClose(outtree);
  FClose(intree);
  /* FClose(outfile); */

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
//printf("Done.\n\n");

#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif

return 0;
}  /* main */
