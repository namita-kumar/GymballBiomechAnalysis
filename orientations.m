function [ G ] = orientations( )
%This function reads data from the orientations of each contact point on the
%hand to the center of the device. It reads from orientation.xls file. 
%The first three columns of the spreadsheet corresponds to rotation matrices 
%and the last column corresponds to the position vector in (mm)
%There are 19 contact points and they are considered in the following order: 
%Thumb  (DP), Thumb proximal phalange (PP), Thumb metacarpal (MCP),
%Index DP, Index middle phalange (MP), Index PP, Mid finger DP, Mid finger
%MP, mid finger PP, mid finger MCP, Ring finger DP, Ring finger MP, Ring
%finger PP, Ring finger MCP, Little finger DP, Little finger MP, Little
%finger PP, Little finger MCP
%where DP is distal phalange, MP is middle phalange, PP is proximal phalange 
%MCP is metacarpal
data = xlsread('orientations',1,'A1:D94');      %Read data from spreadsheet

i = 1;  % initialize couter. Max = 19. There are 19 contact points

R = zeros(3, 3, 19);%Rotation matrix from contact point to object's center
p = zeros(3, 19);   %position vector from contact point to object's center
Ad = zeros(6, 6, 19);   %adjoint wrench transformation  
G = zeros(6, 4, 19);    %Grasp mapping between contact forces
B = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]; %wrench basis

while(i <= 19)
    
    %The rotation matrix data of xth contact point is given by the data in 
    %5x-4 to 5x-2 rows and the first 3 columns
    %The position vector of xth contact point is given by the data in the
    %4th column of 5x-4 to 5x-2 rows
    R(:, :, i) = data((5*i-4):(5*i-2), 1:3);
    p(:, i) = data((5*i-4):(5*i-2), 4);
    p(:, i) = p(:, i) / norm(p(:, i));      %normalize position vector
    
    %cross product of p and R where p is represented by a skew symmetric matrix
    XProd = [0 -p(3, i) p(2, i); p(3, i) 0 -p(1, i); -p(2, i) p(1, i) 0]...
        * R(:, :, i);   
    
    Ad(:, :, i) = [R(:, :, i), zeros(3); XProd, R(:, :, i)];
    G(:, :, i) = Ad(:, :, i)' * B;
    
    i = i+1;
    
end

