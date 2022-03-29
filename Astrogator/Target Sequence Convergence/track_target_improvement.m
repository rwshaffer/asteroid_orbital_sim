function [cons_div] = track_target_improvement(percent_changes,div_threshold)
% Return the number of consecutive "divergences" for each result in a target sequence
% where divergence is defined as being more positive than a number
% -div_threshold.
% This function is used to determine whether corrective action is needed to get a 
% target sequence to converge, or whether it's going in the right direction by itself.
%
% INPUTS:
% percent_changes: n x m matrix, where n is the number of iterations a target sequence has
%   been run for, m is the number of results, and the matrix entries are the percent change
%   of that column's result from the last iteration to the current one.
% div_threshold: scalar, percentage that the result must improve by to be considered NOT
%   diverging. We don't want to consider a result that improved by 0.01% to be improving,
%   since that's basically no change from one iteration to the next.
%
% OUTPUTS:
% cons_div: 1 x m row vector giving how many consecutive divergences each result has recently had 
%   (i.e. since the last time it improved, how many times has it gotten worse?)


A = percent_changes >= -div_threshold;
% NOTES: A is a binary matrix, where 0 means improvement (change in the result at this iteration
% was more negative than the div_threshold)
%
% We want to count the 1's in each column, but reset the count every time there's a zero
%
% Clever way to do this: at each row in the matrix A, add 1 to the output vector in
% each column, then multiply by the element in that column in A. That way, every time we encounter
% a zero in a given column, that column in the output vector is multiplied by 0 (and thus reset).
%
% Example matrix: A = [0 0 1 1 0;
%                      0 0 0 1 1;
%                      1 1 1 0 0;
%                      1 0 1 1 1];
%
% Should return:      [2 0 2 1 1]

cons_div = zeros(1,size(A,2)); % Initialize output vector

% Iterate through rows of A
for n = 1:size(A,1)
    
    cons_div = A(n,:) .* (cons_div + 1); % Add one, then multiply by either 0 or 1 depending on A
    
end
    
    
    
end