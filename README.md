Algorithms
===========

1. K Center problem in 1 dimension
    
        Choosing K Centers for input y (vector N) such that a metric d(..) defined
        on the input space summed over SUM [ d( center, i) ] over all i assigned
        to Centeri is minimized.
        d : is assumed to be fixed as Sum of Squares distance in the code.
        The problem is solved using Bellman's equation and backtracking.

        a. KCenter1D.jl  b. Kcenter1D.py
