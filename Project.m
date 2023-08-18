clc ;
clear ;
 
%  Test data classification
X = [] ; 
Y = [] ;
  
imagefiles = dir('../../../Desktop/754 Project/photographic/*.jpg');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
     
     I = im2double((imread("../../../Desktop/754 Project/photographic/"+imagefiles(ii).name)));
     
%      Images of size less than 300x300 are rejected
     if length(I(:,1)) < 300 || length(I(1,:)) < 300
        continue ;
     end 
   
%      Take a 256x256 crop from centre
     r = centerCropWindow2d(size(I),[256 256]);
     I = imcrop(I,r);
     
     X = [X; find_stat(I)] ;
     Y = [Y; "photographic"] ;
end
  
imagefiles = dir('../../../Desktop/754 Project/photorealistic/*.jpg');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
     
     I = im2double((imread("../../../Desktop/754 Project/photorealistic/"+imagefiles(ii).name)));
     
%      Images of size less than 300x300 are rejected
     if length(I(:,1)) < 300 || length(I(1,:)) < 300
         continue ;
     end
    
%     Take a 256x256 crop from centre
     r = centerCropWindow2d(size(I),[256 256]);
     I = imcrop(I,r);
     
     X = [X; find_stat(I)] ;
     Y = [Y; "photorealistic"] ;
end
%  
% Fit SVM to Training dataset
SVM = fitcsvm(X,Y) ;
  
save('svm.mat', 'SVM') ;
  
load('svm.mat', 'SVM') ;

% Classification
Z1 = [] ;
Z2 = [] ;

fileID = fopen('wrong.txt','w');

imagefiles = dir('../../../Desktop/754 Project/test/photographic/*.jpg');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
    
   I = im2double((imread("../../../Desktop/754 Project/test/photographic/"+imagefiles(ii).name)));

%     Images of size less than 300x300 are rejected   
   if length(I(:,1)) < 300 || length(I(1,:)) < 300
       continue ;
   end 
   
%     Take a 256x256 crop from centre   
   r = centerCropWindow2d(size(I),[256 256]);
   I = imcrop(I,r);
   
%    Find the prediction
   Z1 = [Z1 ; predict(SVM, find_stat(I))] ;
   
   if(predict(SVM, find_stat(I))=="photorealistic") 
       fprintf(fileID,imagefiles(ii).name);
       fprintf(fileID,'\n');
   end
   
end

imagefiles = dir('../../../Desktop/754 Project/test/photorealistic/*.jpg');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
    
   I = im2double((imread("../../../Desktop/754 Project/test/photorealistic/"+imagefiles(ii).name)));

%     Images of size less than 300x300 are rejected   
   if length(I(:,1)) < 300 || length(I(1,:)) < 300
       continue ;
   end 
   
%     Take a 256x256 crop from centre   
   r = centerCropWindow2d(size(I),[256 256]);
   I = imcrop(I,r);
   
%    Find the prediction
   Z2 = [Z2 ; predict(SVM, find_stat(I))] ;
   
   if(predict(SVM, find_stat(I))=="photographic") 
       fprintf(fileID,imagefiles(ii).name);
       fprintf(fileID,'\n');
   end
   
end

fclose(fileID) ;

fprintf("correctly classified photographic images %d/%d \n", sum(Z1(:)=="photographic"), length(Z1)) ;
fprintf("correctly classified photorealistic images %d/%d \n", sum(Z2(:)=="photorealistic"), length(Z2)) ;


function stat = find_stat(I)

    [c,s] = wavedec2(I,3,'qd1');
    [H1,V1,D1] = detcoef2('all',c,s,1);

    [H2,V2,D2] = detcoef2('all',c,s,2);

    [H3,V3,D3] = detcoef2('all',c,s,3);
    
    stat = [] ;

%     Loop over colors RGB
    for i = 1:3
%         H1(x,y) base + error statistics
       [w,p] = func(H1(:,:,i), H2(:,:,i), D1(:,:,i), D2(:,:,i), H1(:,:,mod(i,3)+1), H1(:,:,mod(i+1,3)+1),1) ;
       stat = [stat mean(w) var(w) skewness(w) kurtosis(w)] ;
       stat = [stat mean(p) var(p) skewness(p) kurtosis(p)] ;

%         H2(x,y) base + error statistics
       [w,p] = func(H2(:,:,i), H3(:,:,i), D2(:,:,i), D3(:,:,i), H2(:,:,mod(i,3)+1), H2(:,:,mod(i+1,3)+1),1) ;
       stat = [stat mean(w) var(w) skewness(w) kurtosis(w)] ;
       stat = [stat mean(p) var(p) skewness(p) kurtosis(p)] ;

%         V1(x,y) base + error statistics
       [w,p] = func(V1(:,:,i), V2(:,:,i), D1(:,:,i), D2(:,:,i), V1(:,:,mod(i,3)+1), V1(:,:,mod(i+1,3)+1),1) ; 
       stat = [stat mean(w) var(w) skewness(w) kurtosis(w)] ;
       stat = [stat mean(p) var(p) skewness(p) kurtosis(p)] ;
       
%         V2(x,y) base + error statistics
       [w,p] = func(V2(:,:,i), V3(:,:,i), D2(:,:,i), D3(:,:,i), V2(:,:,mod(i,3)+1), V2(:,:,mod(i+1,3)+1),1) ;
       stat = [stat mean(w) var(w) skewness(w) kurtosis(w)] ;
       stat = [stat mean(p) var(p) skewness(p) kurtosis(p)] ;

%         D1(x,y) base + error statistics      
       [w,p] = func(D1(:,:,i), D2(:,:,i), H1(:,:,i), V1(:,:,i), D1(:,:,mod(i,3)+1), D1(:,:,mod(i+1,3)+1),0) ;
       stat = [stat mean(w) var(w) skewness(w) kurtosis(w)] ;
       stat = [stat mean(p) var(p) skewness(p) kurtosis(p)] ;
 
%         D2(x,y) base + error statistics       
       [w,p] = func(D2(:,:,i), D3(:,:,i), H2(:,:,i), V2(:,:,i), D2(:,:,mod(i,3)+1), D2(:,:,mod(i+1,3)+1),0) ;
       stat = [stat mean(w) var(w) skewness(w) kurtosis(w)] ;
       stat = [stat mean(p) var(p) skewness(p) kurtosis(p)] ;
       
    end
    
end

function [w, p] = func(A,a,D,d,A1,A2,flg)
   s = size(A) ;
   n = s(1) ;
   
   Q = zeros(n*n, 9) ;
   
%    Construct matrix Q
   for j = 1:n*n 
       x = floor((j-1)/n)+1 ;
       y = mod(j-1,n)+1 ;
       
       if x>1
           Q(j,1) = A(x-1,y) ;
       end
       
       if x<n
           Q(j,2) = A(x+1,y) ;
       end
       
       if y>1
           Q(j,3) = A(x,y-1) ;
       end
       
       if y<n
           Q(j,4) = A(x,y+1) ;
       end
       
       Q(j,5) = a(ceil(x/2),ceil(y/2)) ;
       Q(j,6) = D(x,y) ;
       if flg == 1
           Q(j,7) = d(ceil(x/2),ceil(y/2)) ;
       else 
           Q(j,7) = d(x,y) ;
       end
       Q(j,8) = A1(x,y) ;
       Q(j,9) = A2(x,y) ;
       
   end
      
%    w = inv(Q'Q)Q'v
   Q = abs(Q) ;
   v = abs(A(:)) ;
   w = inv(Q'*Q) * Q' * v ;
   
%    Log error
   p = log(v) - log(abs(Q*w)) ;
   
   
end
