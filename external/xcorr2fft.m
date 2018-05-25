function [shiftx,shifty,c] = xcorr2fft(image1,image2)
    
    %   xcorr2fft computes offsets between images image1 and image2 based
    %   on Phase Correlation method. image1 & image2 are assumed to have
    %   the same dimensions.
    %   
    %   Written by: Sebastian J Streichan, EMBL, February 29, 2012
    % some comments and minor changes Idse Heemskerk
    
    % https://en.wikipedia.org/wiki/Cross-correlation
    % cross correlation c(y):
    % obeys F(f*g) = F(f)^* F(g) with ^* as complex conjugate
    % so if c(y) := f*g, then c = F^(-1)[F(f)^* F(g)] with y the
    % displacement vector
    
    F     = fftn(image1);
    Fc    = conj(fftn(image2));
    R     = F.*Fc; 
    
    % the cross correlation is real so its Fourier transform is hermitian
    % (conjugate symmetric), passing this option for speed
    c     = ifftshift(ifftn(R,'symmetric')); 
    [~,i] = max(c(:));
    [I,J] = ind2sub(size(c),i);
    shift = [I J] - size(c)/2 - 1;
    shifty = shift(1);
    shiftx = shift(2);
    
    c = c(I,J);
end