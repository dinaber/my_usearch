# my_usearch
A clustering matlab algorithm for short DNA sequences
The main function: my_usearch is called with a cell array of sequences, threshold score and additional optional parameters.
my_usearch uses my_align_int and my_nt2int.

The outputs include: each bin centroner (cell array of sequences), bin number that was assigned to each sequence in 
the original cell array and a binary vector reporting which seuences had to be reversed in the process.

#############
Example:
[rep_grps,rep_inds,rep_rev] = my_usearch(seq_list ,score0 ,'swalign',-8,[],false,isrev) ;

