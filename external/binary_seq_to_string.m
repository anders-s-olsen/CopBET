%BINARY_SEQ_TO_STRING String representation of a logical vector. 
%   This function takes a vector of logical values, and returns a string
%   representation of the vector.  e.g. [0 1 0 1 1] becomes '01011'.
%   
%   Usage: [s] = binary_seq_to_string(b)
%
%   INPUTS:
%   
%   b: 
%   A vector of logical values representing a binary sequence.
%   Numeric values will be converted to logical values depending on 
%   whether (0) or not (1) they are equal to 0.
%
%
%   OUTPUTS:
%
%   s:
%   A string representation of b.  The nth character in the string
%   corresponds with b(k).
%
%
%
%
%   Author: Quang Thai (qlthai@gmail.com)
%   Copyright (C) Quang Thai 2012


function [s] = binary_seq_to_string(b)


b = logical(b(:));
    

lookup_string = '01';

s = lookup_string(b + 1);
