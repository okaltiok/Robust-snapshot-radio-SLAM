classdef linearAlgebra
    methods (Static)
        function [L,flag] = choleskyFactorization(A)
            n = size(A,1);
            L = zeros(n);
            flag = true;
            
            for i = 1:n
                for j = 1:i
                    sum = 0;
                    for k = 1:j
                        sum = sum + L(i,k)*L(j,k);
                    end
                    
                    if i == j
                        tmp = A(i,i) - sum;
                        if tmp <= 0
                            flag = false;
                            return
                        end
                        L(i,j) = sqrt(tmp);
                    else
                        L(i,j) = (A(i,j)-sum)/L(j,j);
                    end
                end
            end
        end
        

        function [IA,detA] = matrixInverse(A)
            
            N = size(A,1);
            IA = zeros(N);
            P = zeros(N+1,1);
            
            for i = 1:N
                P(i) = i;
            end
            
            for i = 1:N
                maxA = 0.0;
                imax = i;
                
                for k = i:N
                    absA = abs(A(k,i));
                    
                    if (absA  > maxA)
                        maxA = absA;
                        imax = k;
                    end
                end
                
                %             if (maxA < Tol) return 0; //failure, matrix is degenerate
                
                if (imax ~= i)
                    % pivoting P
                    j = P(i);
                    P(i) = P(imax);
                    P(imax) = j;
                    
                    % pivoting rows of A
                    ptr = A(i,:);
                    A(i,:) = A(imax,:);
                    A(imax,:) = ptr;
                    
                    %counting pivots starting from N (for determinant)
                    P(N+1) = P(N+1) + 1;
                end
                
                
                for j = (i+1):N
                    A(j,i) = A(j,i)/A(i,i);
                    for k = (i+1):N
                        A(j,k) = A(j,k) - A(j,i) * A(i,k);
                    end
                end
            end
            
            
            
            detA = A(1,1);
            for j=2:N
                detA = detA*A(j,j);
            end
            
            if mod(P(end),2) ~= 0
                detA = -detA;
            end
            
            for j=1:N
                for i=1:N
                    if P(i) == j
                        IA(i,j) = 1;
                    else
                        IA(i,j) = 0;
                    end
                    
                    for k=1:i-1
                        IA(i,j) = IA(i,j) - A(i,k) * IA(k,j);
                    end
                end
                
                
                for i = N:-1:1
                    for k = (i+1):N
                        IA(i,j) = IA(i,j) - A(i,k) * IA(k,j);
                    end
                    IA(i,j) = IA(i,j)/A(i,i);
                end
            end
        end

    end
end

