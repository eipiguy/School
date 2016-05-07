function [ l ] = latex( M,Brac,NumStyle )
% latex: formats a matrix for LaTeX
%   Writes a string for a LaTeX matrix with a given bracket style:
%   Border = 
%       [blank], no brackets
%       b, box brackets
%       B, curly brackets
%       p, parenthesis
%       v, verticle line brackets
%       V, double verticle line brackets
%       small, small matrix for inline use
% Input:
%   M = matrix to be formatted
%   Brac = bracket style key as above
% Output:
%   latex = LaTeX formatted sting
%##########################################################################
% Pseudo Code: (flops per section)
%   Total Flops: 
%   ####
%   Input Format Check:
%   ====
%   Variable Declarations: 
%   ====
%   Main Algorithm: 
%   ####
%##########################################################################
% Input Format Check:

if nargin<1, error('Need input');end
if nargin<2 || isempty(Brack),b='';end

%==========================================================================
% Variable Declarations:

%==========================================================================
% Main Algorithm:

%##########################################################################

end

