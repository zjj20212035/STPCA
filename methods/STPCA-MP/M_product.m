function [output] = M_product(A,B,M,trans)
%To conduct *M-product on third-order tensor.
%A is of size k*d*n, B is of size d*m*n
%trans = 0 or 1. 0 means C_M is the output, while 1 means C is the output.
n = size(A,3);
switch M
    case "FFT"
        A_M = fft(A,n,3);
        B_M = fft(B,n,3);
        C_M = pagemtimes(A_M,B_M);
        C = ifft(C_M,n,3);
    case "I"
        C = pagemtimes(A, B);
        C_M = C;
    otherwise
        A = tensor(A);
        B = tensor(B);

        A_M = double(ttm(A,M,3));
        B_M = double(ttm(B,M,3));
        C_M = pagemtimes(A_M,B_M);
        C = ttm(C_M,M^-1,3);

end
   
if trans == 0
    output = C_M;
else
    output = C;
end
end