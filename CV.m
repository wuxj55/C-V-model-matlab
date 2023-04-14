%demo_CV.m
close all;
clear all;
clc;

% Read Image
c0 = 2;
Img = imread('skin_cancer_image1.bmp');
Temp = Img;
if ndims(Img) == 3
    Img = rgb2gray(Img);
end
Img = double(Img);

% Initial phi is the level set function
%figure;imagesc(Img,[0,255]);colormap(gray);hold on; axis off; axis equal;
figure;
imshow(Temp);colormap;
rNum = 1;% the cicle number in a row
cNum = 1;% the cicle number in a colcumn
[m,n] = size(Img);
phi= ones(m,n).*c0;
r = min(m/(2*rNum),n/(2*cNum))/2;
    for i = 1:rNum
        for j = 1:cNum
            px = (2*i-1)*(m/(2*rNum));
            py = (2*j-1)*(n/(2*cNum));%(px,py) is the centre of the initial zero level set cicle
            for x = 1:m
                for y = 1:n
                    d = (x-px)^2 + (y - py)^2;
                    if d < r^2
                       phi(x,y) = -c0;
                    end%if
                end%y
            end%x
        end%for j
     end%for i
hold on;
[c,h] = contour(phi,[0 0],'r');
hold off;
   
iterNum =600; % the total number of the iteration
lambda1 = 1;   % the weight parameter of the globe term which is inside the level set
lambda2 = 1;   % the weight parameter of the globe term which is ouside the level set
mu = 0.01*255*255; % the weight parameter of the length term
nu = 0; % the weight parameter of the area term
pu = 1.0; %the weight parameter of the penalizing term
timestep = 0.1; % the time step of the level set function evolution
epsilon = 1.0; % the parameter of the Heaviside function and the Dirac delta function
%tol=2*10^-5;

% all model's initial level set is same
phi_CV        = phi; 
phi_star      = phi;

figure;
imshow(uint8(Img)); colormap;

%start the level set evolution
%CV Model 
error=zeros(1,600);
tic
for iter = 1:iterNum
	numIter = 1; 
	% level set evolution. 
	phi_CV1 = EVOL_CV(phi_CV, Img, lambda1, lambda2, mu, nu,pu, timestep,epsilon, numIter); 
	if mod(iter, 10) == 0
	    contour(phi_CV1, [0,0], 'y'); 
	end 	
    error(iter)=norm(phi_CV1-phi_CV,2)/bwarea(Img);
    phi_CV = phi_CV1;
end 
toc
plot(error);

% Display Results
figure;
imshow(uint8(Img));
hold on;
contour(phi_star,[0,0],'r','linewidth',1);
title('Initial Level set');

figure; 
imshow(uint8(Img)); 
hold on; 
contour(phi_CV1, [0,0], 'r', 'linewidth', 1); 
title('Results of CV model'); 
