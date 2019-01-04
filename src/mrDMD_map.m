function [ ptree, map, low_f_cutoff,Phi] = mrDMD_map(mrdmd)
% function [ ptree, map, low_f_cutoff,Phi] = mrDMD_map(mrdmd)

% Inputs:
% mrdmd     tree structure returned by mrDMD routine
% 
% Outputs:
% ptree        	matrix of identified mode amplitudes in order
% map         	matrix grid of identified mode amplitudes 
% low_f_cutoff  vector of low freq cutoff for each level of recursion
% Phi           matrix of mrDMD modes above the freq cutoff
%
% Code adapted from JN Kutz, SL Brunton, BW Brunton and JL Proctor,
% "Dynamic Mode Decomposition", SIAM
%
% Modified 2018/12/31

[levels, M] = size(mrdmd);

ptree = zeros(size(mrdmd));
map = zeros(levels, M);
Phi = [];

low_f_cutoff = zeros(levels+1, 1);
for i = 1:levels,
    chunks = 2^(i-1);
    K = M / chunks;
    for j = 1:chunks,
        f = imag(mrdmd{i, j}.omega);
        P = mrdmd{i, j}.P;

        P = P(abs(f) >= low_f_cutoff(i));
        Phi = [Phi, mrdmd{i,j}.Phi(:,f > low_f_cutoff(i))];
        
        if ~isempty(P),
            map(i, (1:K)+(j-1)*K) = mean(abs(P));
            ptree(i,j) = mean(abs(P));
        end;
    end;
    
    low_f_cutoff(i+1) = mrdmd{i,1}.rho;
end;