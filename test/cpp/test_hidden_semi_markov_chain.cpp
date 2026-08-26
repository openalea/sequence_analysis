/****************************************************************
 *
 *  Test hidden semi-Markov chains
 */

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/discrete_mixture.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"

using namespace stat_tool;
using namespace sequence_analysis;

int main(void)
{

   register int var, t, u, nb_integral, cptm= 0; // j, i, nb_states,

   bool status, geometric_poisson=false, common_dispersion=false, counting_flag=true, state_sequence=true;
   stat_tool::rounding mode = ROUND;
   const int nb_sequence = 30 , length = 100;
   const stat_tool::process_type itype = stat_tool::ORDINARY;
   HiddenSemiMarkov *hsmc= NULL, *hsmc_ref= NULL, *hsmc_est_file= NULL;
   SemiMarkovData *hsmd= NULL;
   stat_tool::censoring_estimator estimator=stat_tool::COMPLETE_LIKELIHOOD;
   // Hidden_variable_order_markov *hmc= NULL, *hmc_init= NULL;
   MultiPlotSet *plotable=NULL;
   MarkovianSequences *mseq= NULL;
   Sequences *seq= NULL, *seq_cluster= NULL;
   StatError error;
   std::vector< int > select;
   const char * hsmcrefpath= "../../../share/data/pin_laricio_6.hsc";
   const char * datarefpath = "../../../share/data/pin_laricio_7x.seq";

 
   seq =  Sequences::ascii_read(error, datarefpath, false);
   seq_cluster = seq->cluster(error, 1, 10, mode);

   hsmc_ref = HiddenSemiMarkov::ascii_read(error, hsmcrefpath);
   mseq = new MarkovianSequences(*seq_cluster);
   hsmc_est_file = mseq->hidden_semi_markov_estimation(error, &cout, *hsmc_ref, geometric_poisson , common_dispersion, estimator, counting_flag, state_sequence);
   
   hsmc_est_file->ascii_write(cout);

   delete seq;
   delete seq_cluster;
   delete hsmc_est_file;
   delete hsmc_ref;
   delete mseq;

   return 0;
}
