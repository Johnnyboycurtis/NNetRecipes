/******************************************************************************/
/*                                                                            */
/*  CLASSES - Headers for all classes                                         */
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

/*
--------------------------------------------------------------------------------

   TrainingSet - Collection of samples which will be used for training

   All training data is stored in 'data'.  If the output model is HETERO, the
   input for each sample is followed in 'data' by the number (1 through nout)
   of the output class to which that sample belongs.  If the output model is
   GENERAL, the output values follow the input values.  No output values are
   stored in AUTO mode, since the output is the input!

   The variable 'bufcnt' is the number of samples which could still fit in
   'data' without allocating more memory.  This allows us to allocate memory
   in large blocks, avoiding the overhead of often expensive operating system
   calls to malloc.

--------------------------------------------------------------------------------
*/

class TrainingSet {

public:

   TrainingSet ( int outmodel , int n_inputs , int n_outputs ) ;
   ~TrainingSet () ;
   void train ( char *file , int outclass ) ;

   unsigned ntrain ; // Number of samples in 'data'
   double *data ;    // Actual training data here
   int outmod ;      // Output model (see OUTMOD_? in CONST.H)
   int nin ;         // Number of input neurons
   int nout ;        // Number of output neurons

private:

   unsigned bufcnt ; // Sample areas remaining unused at end of 'data'
} ;

/*
--------------------------------------------------------------------------------

   SingularValueDecomp - Singular value decomposition of matrices

   Normally, the user would:
     1) Use 'new' to create a SingularValueDecomp object with all necessary
        memory (a, u?, w, v, work, b) allocated by the constructor.
     2) Fill in the public 'a' with the matrix to be decomposed.
     3) Call svdcmp to decompose a, replacing it with the U matrix if preserve
        is zero, else computing u.
        This will also compute w and v, which are normally not used but are
        public just in case the user wants to access them.
     4) Fill in the 'b' vector with the right hand side of the equations.
     5) Call backsub with a pointer to the cols vector which is where the
        solution will be placed.  This vector is NOT allocated by the
        constructor.  The outputs of svdcmp (a, u?, w, v) will not be disturbed.
     6) Repeat the above step as desired.
     7) Delete the SingularValueDecomp object, which frees all memory which
        was allocated by the constructor.

--------------------------------------------------------------------------------
*/

class SingularValueDecomp {

public:

   SingularValueDecomp ( int rows , int cols , int preserve ) ;
   ~SingularValueDecomp () ;
   void svdcmp () ;
   void backsub ( double thresh , double *x ) ;

   int ok ;         // Was memory allocation successful?

/*
   The following four input/output areas are allocated by the constructor
*/

   double *a ;      // Rows by cols input of 'A' matrix, output of U
   double *u ;      // unless preserve != 0, in which case U output here
   double *w ;      // Cols vector output of singular values, not sorted
   double *v ;      // Cols by cols output of 'V' matrix
   double *b ;      // Rows vector of RHS input to backsub


private:

   int rows ;       // Number of rows in 'A' matrix
   int cols ;       // and number of columns
   double *work ;   // Cols work vector (allocated by constructor)
} ;

/*
--------------------------------------------------------------------------------

   Network - General parent of things common to all networks

   The constructors for the child classes set 'ok' to one if all memory
   allocation succeeded.  Otherwise they set it to zero,

--------------------------------------------------------------------------------
*/

class Network {

public:
   void classify_from_file ( char *name , double thresh ) ;
   void execute_from_file ( char *inname , char *outname ) ;
   virtual void learn ( TrainingSet *tptr ,
                        struct LearnParams *lptr ) = 0 ;
   void reset_confusion () ;
   void save_confusion ( char *name ) ;
   void show_confusion () ;
   virtual void trial ( double *input ) = 0 ;
   virtual int wt_save ( FILE *fp ) = 0 ;
   virtual void wt_restore ( FILE *fp ) = 0 ;

   virtual ~Network() ;

   int ok ;         // Was all constructor memory allocation successful?
   double *out ;    // Output neuron activations here if exe nonzero
   double neterr ;  // Mean square error of the network if executed
   int outmod ;     // Output model (see OUTMOD_? in CONST.H)
   int nin ;        // Number of input neurons
   int nout ;       // Number of output neurons
   int *confusion ; // Handy work area avoids malloc/free (see CONFUSE.CPP)

protected:

   double trial_error ( TrainingSet *tptr ) ;  // Compute training set error

   int exe ;        // Are work areas allocated to allow executing it?
} ;

/*
--------------------------------------------------------------------------------

   LayerNet

   Nhid1 and nhid2 are the number of neurons in the first and second layers
   respectively.  Either or both may be zero.  If nhid1 is zero, it is assumed
   that nhid2 is also zero.

   Weights for a layer are stored as a two dimensional matrix strung out into
   a vector.  For example, the first element in hid1_coefs is the weight
   connecting the first input neuron to the first hidden neuron.  The second
   connects the second input to the first hidden.  The nin+1 element is the
   bias for the first hidden unit.  The nin+2 connects the first input to the
   second hidden, et cetera.

--------------------------------------------------------------------------------
*/

class LayerNet: public Network {

public:

   LayerNet ( int outmodel , int n_inputs , int n_hidden1 , int n_hidden2 ,
              int n_outputs , int exe , int zero ) ;
   ~LayerNet () ;
   void learn ( TrainingSet *tptr , struct LearnParams *lptr ) ;
   void trial ( double *input ) ;
   int wt_save ( FILE *fp ) ;
   void wt_restore ( FILE *fp ) ;

   int nhid1 ;              // Number of neurons in hidden layer 1
   int nhid2 ;              // Ditto for hidden layer 2

private:
    
   void anneal ( TrainingSet *tptr , struct LearnParams *lptr ,
                 LayerNet *bestnet , int init ) ;
   double conjgrad ( TrainingSet *tptr , int maxits ,
                     double reltol , double errtol ) ;
      void check_grad ( TrainingSet *tptr , double *grad ) ;
      double find_grad ( TrainingSet *tptr , double *hid2delta ,
                         double *outdelta , double *grad ) ;
      void find_new_dir ( double gam , double *g , double *h , double *grad ) ;
      double gamma ( double *g , double *grad ) ;
   void copy_weights ( LayerNet *dest , LayerNet *source ) ;
   double direcmin ( TrainingSet *tptr , double start_err ,
                     int itmax , double eps , double tol ,
                     double *base , double *direc ) ;
      void negate_dir ( double *direc ) ;
      void preserve ( double *base ) ;
      void step_out ( double step , double *direc , double *base ) ;
      void update_dir ( double step , double *direc ) ;
   void gen_init ( TrainingSet *tptr , struct LearnParams *lptr ) ;
   void perturb ( LayerNet *cent , LayerNet *perturbed , double temp ,
                  int regress) ;
   double regress ( TrainingSet *tptr , SingularValueDecomp *sptr ) ;
   void zero_weights () ;

   double *hid1_coefs ;     // nhid1 * (nin+1) weights (in changes fastest)
   double *hid2_coefs ;     // nhid2 * (nhid1+1) weights (hid1 changes fastest)
   double *out_coefs ;      // nout * (nhid?+1) weights (hid? changes fastest)
   double *hid1 ;           // Hid 1 neuron activations here if exe nonzero
   double *hid2 ;           // Ditto hidden layer 2
} ;

/*
--------------------------------------------------------------------------------

   KohNet

   There are no hidden units and no bias terms.
   The output weights are stored as a two dimensional matrix
   strung out into a vector.
   For example, the first element in out_coefs is the weight
   connecting the first input neuron to the first output neuron.
   The second connects the second input to the first output.
   Element nin+1 connects the first input to the second output, etc.

--------------------------------------------------------------------------------
*/

class KohNet: public Network {

public:

   KohNet ( int n_inputs , int n_outputs , KohParams *kp ,
            int exe , int zero ) ;
   ~KohNet () ;
   void learn ( TrainingSet *tptr , struct LearnParams *lptr ) ;
   void trial ( double *input ) ;
   int wt_save ( FILE *fp ) ;
   void wt_restore ( FILE *fp ) ;

   int normalization ;      // 0=multiplicative, 1=Z 

private:
    
   void copy_weights ( KohNet *dest , KohNet *source ) ;
   void epoch1 ( TrainingSet *tptr , double rate , int learn_method ,
                 int *won , double *bigerr , double *correc , double *work ) ;
   void epoch2 ( double rate , int learn_method , int *won ,
                 double *bigcorr , double *correc ) ;
   void force_win ( TrainingSet *tptr , int *won ) ;
   void initialize () ;
   void in_norm ( double *input , double *normfac , double *synth ) ;
   int winner ( double *input , double *normfac , double *synth ) ;
   void wt_norm ( double *w ) ;
   void zero_weights () ;

   double *out_coefs ;      // nout * (nin+1) weights
} ;
