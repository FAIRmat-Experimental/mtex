function fhat = convSO3(fhat1,fhat2)
%

% old sizes
s1 = size(fhat1);
s2 = size(fhat2);

% get bandwidth
L = min(dim2deg(s1(1)),dim2deg(s2(1)));

% new size
s = max(s1(2:end),s2(2:end));

% compute Fourier coefficients of the convolution
fhat = zeros([deg2dim(L+1),s]);
if prod(s) == 1 %simple SO3Fun
  for l = 0:L
    ind = deg2dim(l)+1:deg2dim(l+1);
    fhat(ind) = reshape(fhat2(ind),2*l+1,2*l+1) * ...
      reshape(fhat1(ind),2*l+1,2*l+1) ./ sqrt(2*l+1);     
  end
else % vector valued SO3Fun  
  for l = 0:L
    ind = deg2dim(l)+1:deg2dim(l+1);
    fhat_l = pagemtimes( reshape(SO3F2.fhat(ind,:),[2*l+1,2*l+1,s2]) , ...
      reshape(SO3F1.fhat(ind,:),[2*l+1,2*l+1,s1]) ) ./ sqrt(2*l+1);
    fhat(ind,:) = reshape(fhat_l,[],prod(s));
  end
end

end
