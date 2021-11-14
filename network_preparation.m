function [W,lb,ub] = network_preparation(A,normalization,diagonal)

%{
%=========================================================================%
PURPOSE          This function is meant to manipulate the I-O matrices, for
                 instance it can: a) make the matrix row-stochastic; b)
                 remove the main diagonal for preserving consistency. It
                 also calculates accordingly the support where the
                 associated spatial coefficient can live. 

%=========================================================================%

INPUT            A                 matrix (n x n)
                 normalization     string. It can be equal to 'yes' or
                                   'no'. If 'yes', the function carries out
                                   a row-normalization of the input matrix,
                                   otherwise it does not.
                 diagonal          string. It can be 'with' or 'without'.
                                   if 'without', the function removes the
                                   main diagonal, otherwise it does not.

OUTPUT           W                 is the resulting matrix after the
                                   aforementioned manipulations
                 lb                is the lower bound of the support of the
                                   associated spatial coefficient.
                 ub                is the upper bound of the support of the
                                   associated spatial coefficient.

%=========================================================================%

Favero, Karamysheva, Briganti 
"The network effects of fiscal adjustments" (2019)

%=========================================================================%
%}

% Size of the network:
[n1,n2] = size(A);
if n1~=n2
    error('Network is not symmetric')
end

% Remove the Main Diagonal:
switch diagonal
    case "without"
        Wd = A - diag(diag(A));
    case "with"
        Wd = A;
end

% Apply Row-Normalization:
switch normalization
    case "yes"
        W = Wd ./ sum(Wd,2) ;
    case "no"
        W = Wd;
end


% Calculate the bounds for the spacial coefficent: 1) calculate
% eigenvalues; 2) remove the complex ones; 3) sort them and invert
% them. The reason why we invert them is that we are looking for the
% solutions to: |eye(n) - rho*W|=0, that is, the value of rho which
% nullifies the determinant, which corresponds to the inverse of the
% eigenvalues of W. For further reference see LeSage&Pace (2009) -
% chapter 4 (intro. to spatial econometrics).

% Calculate Eigenvalues:
eigenvalues = eig(W);
% Retain only Real Eigenvalues:
eigenvalues_real = eigenvalues(imag(eigenvalues)==0); 
% Calculate positive eighenvalues:
eigenvalues_real_positive = eigenvalues_real(eigenvalues_real>0);
% Calculate upper bound:
ub = 1./max(eigenvalues_real_positive);
% Calculate negative eighenvalues:
eigenvalues_real_negative = eigenvalues_real(eigenvalues_real<0);
% Calculate lower bound:
if length(eigenvalues_real_negative) >= 1 
    lb = 1./min(eigenvalues_real_negative);
else
    lb = 0;
end



end