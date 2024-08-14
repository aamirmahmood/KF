clc;close all;clear all

C = [1 0.4; 0.4 1]; % define covariance matrix

[eigVec , eigVal ] = eig(C); % find eigenvalues and eigenvectors

if eigVal (1 ,1) > eigVal (2 ,2) % get the highest eigenvalue index

    a = sqrt ( eigVal (1 ,1)); % half - major axis length
    b = sqrt ( eigVal (2 ,2)); % half - minor axis length
    theta = atan ( eigVec (2 ,1) / eigVec (1 ,1)); % ellipse angle (radians )
else

    a = sqrt ( eigVal (2 ,2)); % half - major axis length
    b = sqrt ( eigVal (1 ,1)); % half - minor axis length
    theta = atan ( eigVec (2 ,2) / eigVec (2 ,1)); % ellipse angle (radians )
end