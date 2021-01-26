function d = double(E,varargin)
% convert embedding to packed double
% 
% Syntax
%
%   d = double(E)
%   d = double(E,'full')
%
% Input
%  E - @embedding
%
% Options
%  full - return full embedding
%
% Output
%  d - double
%
% See also
% OrientationEmbeddings
%

if check_option(varargin,'full')
  
  for k = 1:length(E.u)
    
    d{k} = reshape(double(E.u{k}),[],size(E,1));
    
  end
  
  d = vertcat(d{:}).';
  return;
  
end

ind2 = [1 2 3 5 6];

% d = {1 1 2; 1 1 3; 1 2 2; 1 2 3; 1 3 3; 2 2 3; 2 3 3};
% ind3 = sub2ind([3 3 3],[d{:,1}],[d{:,2}],[d{:,3}])
ind3 = [10,19,13,22,25,23,26];

% d = {1 1 1 1; 1 1 1 2; 1 1 1 3; 1 1 2 2; 1 1 2 3; 1 1 3 3; 1 2 2 2; 1 2 2 3; ...
%      1 2 3 3; 1 3 3 3; 2 2 2 2; 2 2 2 3; 2 2 3 3; 2 3 3 3};
% ind4 = sub2ind([3 3 3 3],[d{:,1}],[d{:,2}],[d{:,3}],[d{:,4}]);
ind4 = [1 28 55 37 64 73 40 67 76 79 41 68 77 80];
  
ind43 = ind4([2 3 4 6 7 10 12 13 14]);

% d = unique(combnk([1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3],6),'rows');
% ind6 = sub2ind([3 3 3 3 3 3],[d(:,1)],[d(:,2)],[d(:,3)],[d(:,4)],[d(:,5)],[d(:,6)]);
ind6 = [1 244 487 325  568 649 352 595 676 703 361 604 685 712 721 364 607 688 715 724 727 365 608 689 716 725 728 729];

% d = [1 1 2; 1 1 3; 1 2 2; 1 2 3; 1 3 3; 2 2 3; 2 3 3 ];
% indT = sub2ind([3 3 3],[d(:,1)],[d(:,2)],[d(:,3)]);
indT = [10 19 13 22 25 23 26];
    
switch E.CS.Laue.id
    
  case 2 % C1
    
    d = cat(3,double(E.u{1}),double(E.u{2}),double(E.u{3}));
    d = reshape(d,[],9);
    
  case {5,8,11} % C2

    M2 = reshape(E.u{2}.M,9,[]).';
    M3 = reshape(E.u{3}.M,9,[]).';
    
    d = [reshape(double(E.u{1}),[],3), M2(:,ind2),M3(:,ind2)];
    d([5,6,8]) = sqrt(2) * d([5,6,8]);
    d([10,11,13]) = sqrt(2) *  d([10,11,13]);
    
  case 16 % D2
    
    M1 = reshape(E.u{1}.M,9,[]).';
    M2 = reshape(E.u{2}.M,9,[]).';
    M3 = reshape(E.u{3}.M,9,[]).';
    
    d = [M1(:,ind2),M2(:,ind2),M3(:,ind2)];
    d([2,3,5]) = sqrt(2) * d([2,3,5]);
    d([7,8,10]) = sqrt(2) *  d([7,8,10]);
    d([12,13,15]) = sqrt(2) *  d([12,13,15]);
    
       
  case 18 % C3

    M = reshape(E.u{1}.M,27,[]).';
    d = [M(:,ind3),reshape(double(E.u{2}),[],3)];
    d([1:7]) = sqrt(3) * d([1:7]);

  case {21,24} % D3
    
    M1 = reshape(E.u{1}.M,27,[]).';
    M2 = reshape(E.u{2}.M,9,[]).';
    
    d = [M1(:,ind3),M2(:,ind2)];
    d([1:7]) = sqrt(3) * d([1:7]);
    d([9,10,12]) = sqrt(2) * d([9,10,12]);
    
    
  case 27 % C4
    
    M1 = reshape(E.u{1}.M,3^4,[]).';   
    d = [M1(:,ind4),reshape(double(E.u{2}),[],3)];
    d([2,3,7,10,12,14]) =  sqrt(4) * d([2,3,7,10,12,14]);
    d([4,6,13]) = sqrt(6) * d([4,6,13]);
    d([5,8,9]) = sqrt(12) * d([5,8,9]);
    
  case 32 % D4
    
    M1 = reshape(E.u{1}.M,3^4,[]).';
    d = M1(:,ind4);
    d([2,3,7,10,12,14]) =  sqrt(4) * d([2,3,7,10,12,14]);
    d([4,6,13]) = sqrt(6) * d([4,6,13]);
    d([5,8,9]) = sqrt(12) * d([5,8,9]);
    
  case 35 % C6
    M1 = reshape(E.u{1}.M,3^6,[]).';
    d = [M1(:,ind6(1:end-1)),reshape(double(E.u{2}),[],3)];
    d([2,3,16,21,23,27]) = sqrt(6) * d([2,3,16,21,23,27]);
    d([4,6,11,15,24,26]) = sqrt(15) *  d([4,6,11,15,24,26]);
    d([5,17,20]) = sqrt(30) * d([5,17,20]);
    d([7,10,25]) = sqrt(20) * d([7,10,25]);
    d([8,9,12,14,18,19]) = sqrt(60) * d([8,9,12,14,18,19]);
    d([13]) = sqrt(15*6) * d([13]);
    
  case 40 % D6
    M1 = reshape(E.u{1}.M,3^6,[]).';
    d = M1(:,ind6(1:end-1));
    d([2,3,16,21,23,27]) = sqrt(6) * d([2,3,16,21,23,27]);
    d([4,6,11,15,24,26]) = sqrt(15) *  d([4,6,11,15,24,26]);
    d([5,17,20]) = sqrt(30) * d([5,17,20]);
    d([7,10,25]) = sqrt(20) * d([7,10,25]);
    d([8,9,12,14,18,19]) = sqrt(60) * d([8,9,12,14,18,19]);
    d([13]) = sqrt(15*6) * d([13]);
     
    
  case 42 % T
    M1 = reshape(E.u{1}.M,3^3,[]).';
    d = M1(:,indT);
    d = sqrt(3) * d;
    
  case 45 % O
    
    M = reshape(E.u{1}.M,3^4,[]).';
    d = M(:,ind43);
    d([1,2,5,6,7,9]) =  sqrt(4) * d([1,2,5,6,7,9]);
    d([3,4,8]) = sqrt(6) * d([3,4,8]);
            
end
end