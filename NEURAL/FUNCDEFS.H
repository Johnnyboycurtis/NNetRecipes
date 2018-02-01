/******************************************************************************/
/*                                                                            */
/*  FUNCDEFS - Function definition headers                                    */
/*                                                                            */
/* Copyright (c) 1993 by Academic Press, Inc.                                 */
/*                                                                            */
/* All rights reserved.  Permission is hereby granted, until further notice,  */
/* to make copies of this diskette, which are not for resale, provided these  */
/* copies are made from this master diskette only, and provided that the      */
/* following copyright notice appears on the diskette label:                  */
/* (c) 1993 by Academic Press, Inc.                                           */
/*                                                                            */
/* Except as previously stated, no part of the computer program embodied in   */
/* this diskette may be reproduced or transmitted in any form or by any means,*/
/* electronic or mechanical, including input into storage in any information  */
/* system for resale, without permission in writing from the publisher.       */
/*                                                                            */
/* Produced in the United States of America.                                  */
/*                                                                            */
/* ISBN 0-12-479041-0                                                         */
/*                                                                            */
/******************************************************************************/

extern double actderiv ( double f ) ;
extern double activity ( double *input , double *coefs , int n ) ;
extern double act_func ( double x ) ;
extern void act_func_init () ;
extern double dotprod ( int n , double *vec1 , double *vec2 ) ;
extern void error_message ( char *msg ) ;
extern int getch () ;
extern void get_control_line ( char *control_line , int *control_file_number ,
                               FILE *control_files[] ) ;
extern double inverse_act ( double f ) ;
extern int kbhit () ;
extern long longrand () ;
extern long longrandmax () ;
extern void *memalloc ( unsigned int n ) ;
extern void memclose () ;
extern void memfree ( void *ptr ) ;
extern void *memrealloc ( void *ptr , unsigned int size ) ;
extern void memtext ( char *text ) ;
extern void memory_message ( char *msg ) ;
extern void normal_message ( char *msg ) ;
extern void nomemclose () ;
extern void notext ( char *text ) ;
extern double ParseDouble ( char **str ) ;
extern void progress_message ( char *msg ) ;
extern void shake ( int nvars , double *center , double *x , double temp ) ;
extern void slongrand ( long iseed ) ;
extern void split_control_line ( char *line , char **command , int *n_command ,
                                 char **rest , int *n_rest ) ;
extern void stack_control_file ( char *filename , int *control_file_number ,
                                 FILE *control_files[] ) ;
extern void str_to_upr ( char *str ) ;
extern double unifrand () ;
extern double veclen ( int n , double *v ) ;
extern void warning_message ( char *msg ) ;
extern Network *wt_restore ( char *name , int *net_model ) ;
extern void wt_save ( Network *net , int net_model , int updating , char *name);
