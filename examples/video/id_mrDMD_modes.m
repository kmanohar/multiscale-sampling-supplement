function [ modes ] = id_mrDMD_modes( mrdmd, omega_true, ptree )
%id_mrDMD_modes returns a matrix whose columns are the mrDMD modes 
%   associated with the given true frequencies of the data
%
%   INPUTS
%   mrdmd: cell array returned by mrDMD() routine
%   omega_true: true frequencies of the temporal dynamics
%   ptree: matrix storing the frequencies associated with mrdmd cell array
%
%   See script video_mrdmd_demo.m for example usage
% 
% Modified 2018/12/31


L = []; T = []; K = [];
Omega = [];
modes = [];

[I,J] = find(ptree>0);

for ii = 1:length(I)
    i = I(ii); j = J(ii);
    
    f = abs(imag(mrdmd{i, j}.omega));

    [k] = find(f>0);
    
    if ~isempty(k)
        Omega = [Omega; f(k(:))];
        L = [L; repmat(i,length(k),1)];
        T = [T; repmat(j,length(k),1)];
        K = [K; k(:)];
    end
end

for i= 1:length(omega_true)
    [~,ind] = min(abs(Omega -omega_true(i)));
    
    ind = ind(1);
    modes = [modes, real(mrdmd{L(ind),T(ind)}.Phi(:,K(ind)))];
    
end

