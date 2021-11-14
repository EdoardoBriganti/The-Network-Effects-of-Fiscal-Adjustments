function [W_syn,lb_syn,ub_syn] = network_simulation(W,placebo_type)

% Calculate Network size:
[n1,n2] = size(W);
if n1~=n2
    error('Non-symmetric Network')
else
    n = n1;
end



switch placebo_type
       
        
    case "row_shuffle" % shuffle elements within a column
        
        W_syn = W(randperm(n),:);
        
    case "column_shuffle"    % shuffle elements within a row
        
        W_syn = W(:,randperm(n));
    
    case "match_first_and_second"   % Similar to Ozdagli&Weber --- it is WRONG!
        

        % Calculate first order indegree:
        w1 = sum(W,2);
        % Calculate second order indegree:
        w2 = sum(W^2,2);
        
        % Calculate number of zeros in the network and force values smaller
        % than 0.0001 to be zero:
        W_syn = W;
        mu = 0.0001;
        W_syn(W_syn<mu) = 0;
        nZ = sum(W_syn==0,2);
        
        % Estimate scale and shape parameters of Generalized Pareto curve:
        fit = W_syn(:);
        parmhat = gpfit( fit(fit>0) );
        % Store k: shape (or xi in wikipedia)
        k     = parmhat(1); 
        % Store sigma: scale (or sigma in wikipedia)
        sig = parmhat(2)+0.05;
        
        % since the choice of index H and K is arbitrary, set:
        K = 1;
        H = 2;
        except_K_and_H = setdiff(1:n,[H K]);
        
        % Construct synthetic matrix row by row:
        for i = 1:n
            
            a_iK = -10;
            a_iH = -10;
            count = 0;
            while a_iH<0 || a_iK<0 || a_iH>1 || a_iK>1 ||sum(a_ij>1)>0 
                
                % Count the number of iterations:
                count = 1 + count;
                
                % Draw n-2-nZ elements from the generalized pareto distribution
                % estimated above:
                a_ij = gprnd(k,sig,mu,[n-nZ(i)-2,1]);

                % Construct the last n-2 elements of row 1 of the synthetic
                % matrix:
                Ri = [a_ij'  zeros(1,nZ(i))];
                Ri = Ri(1,randperm(n-2));

                % calculate value of row column H:
                a_iH = (w2(i) - w1(i)*w1(K))/(w1(H)-w1(K)) + ...
                    sum(Ri .* ( (w1(K) - w1(except_K_and_H) )/ (w1(H)-w1(K)) )' );

                % Calculate the value in position K:
                a_iK = w1(i) - a_iH - sum(Ri);
                
                if count>1000
                    break;
                end
              
            end
            
            % Add element K and H to row i of synthetic matrix:
            W_syn(i,K) = a_iK;
            W_syn(i,H) = a_iH;
            W_syn(i,except_K_and_H) = Ri;
            
        end
        


        % Calculate synthetic first and second indegree:
        w1_syn = fix(sum(W_syn,2)   * 1e4) / 1e4;
        w2_syn = fix(sum(W_syn^2,2) * 1e4) / 1e4;
        w1 = fix(w1 * 1e4) / 1e4;
        w2 = fix(w2 * 1e4) / 1e4;
        
        if ~isequal(w1,w1_syn) || ~isequal(w2_syn,w2)
            disp('First or second order indegree of the network DONT match :-( ')
        else
            disp('First or second order indegree of the network match :-)')
        end
       
        
end



% Calculate bounds of the network:
[W_syn,lb_syn,ub_syn] = network_preparation(W_syn,"no","with");

end