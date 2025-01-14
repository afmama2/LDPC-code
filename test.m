matrix = genmatrix('rb_bravo.txt');
[G,H] = H2G(matrix);
[r,k] = size(G);
matrix=H;
num_v = size(H,2);
num_c = size(H,1);
num_iterations = 100;
% X = randi([0 1],1,r); 
% codeword = mod(X*G,2);
size(H)
count = 0;

num_points = 20;
num_runs = 10;

minSNR=1;
SNRstepsize=0.05;
maxSNR=minSNR+SNRstepsize*(num_points-1);

Bers = zeros(1,num_points);
Bers2 = Bers;
Bers3 = Bers;
Bers4 = Bers;

Rows = cell(num_c);
Columns = cell(num_v);

for n=1:num_c
    row = [];
    for m=1:num_v
        if matrix(n,m)==1
            row = [row m];
        end
    end
    Rows{n} = row;
end

for n=1:num_v
    column = [];
    for m=1:num_c
        if matrix(m,n)==1
            column = [column m];
        end
    end
    Columns{n} = column;
end

parfor Noise_lvls=1:num_points
    bers = zeros(1,num_runs);
    bers2 = bers;
    bers3 = bers;
    bers4 = bers;
    for l=1:num_runs
        counts = 0;
        countf = 0;
        countm = 0;
        countv = 0;
        X = randi([0 1],1,r); 
        codeword = mod(X*G,2);
    
        SNR = minSNR+((Noise_lvls-1)*SNRstepsize);
        N0 = 1/(10^(SNR/10));
        H = zeros(num_c,num_v);
        H2 = H;
        H3 = H;
        H4 = H;
        a_tx = -2*(codeword-0.5);
        a_rx= a_tx + sqrt(N0/2)*(randn(size(a_tx))+1i*randn(size(a_tx)));
        Channel = (abs(a_rx+1).^2-abs(a_rx-1).^2)/N0;
        for n = 1:num_iterations
            %% variable node calculations
        
            
            mat1 = [Channel;H];
            mat2 = repmat(sum(mat1),num_c,1);
            mat2(matrix==0) = 0;
            H = mat2 - H;
           
        
            %% check node calculations
            for c = 1:num_c
            H(c,Rows{c}) = check_node(H(c,Rows{c}));
            end
            countf = countf +1;
            HplusC = [Channel;H];
            sums = sum(HplusC);
            xhat = zeros(1,num_v);
            xhat(sums<0) = 1;
            %err = sum(abs(xhat-codeword));
            %bers(Noise_lvls,l) = err/k;
            
            s = (mod(xhat*matrix.',2));
            syndrome_output = sum(s);   
            if syndrome_output==0
              break
            end
        end
        for n = 1:num_iterations
            for c = 1:num_c
            for i=Rows{c}
                H2(c,i) = Channel(i) + sum(H2(Columns{i},i)) - H2(c,i);
            end
            H2(c,Rows{c}) = check_node(H2(c,Rows{c}));
            end
            
            counts = counts +1;
            H2plusC = [Channel;H2];
            sums2 = sum(H2plusC);
            xhat2 = zeros(1,num_v);
            xhat2(sums2<0) = 1;
            %err2 = sum(abs(xhat2-codeword));
            %bers2 = [bers2 err2/k];
        
            s = (mod(xhat2*matrix.',2));
            syndrome_output = sum(s);   
            if syndrome_output==0
              break
            end
        end
        for n = 1:num_iterations

            mat1 = [Channel;H3];
            mat2 = repmat(sum(mat1),num_c,1);
            mat2(matrix==0) = 0;
            H3 = mat2 - H3;
            %% check node calculations
            for c = 1:num_c
                H3(c,Rows{c}) = check_node1(H3(c,Rows{c}));
            end
            
            countm = countm + 1;
            H3plusC = [Channel;H3];
            sums3 = sum(H3plusC);
            xhat3 = zeros(1,num_v);
            xhat3(sums3<0) = 1;
            %err3 = sum(abs(xhat3-codeword));
            %bers3 = [bers3 err3/k];
        
            s = (mod(xhat3*matrix.',2));
            syndrome_output = sum(s);   
            if syndrome_output==0
              break
            end
        end
        
        for n = 1:num_iterations
%         tic
            for c = 1:num_c
            for i=Rows{c}
                H4(c,i) = Channel(i) + sum(H4(Columns{i},i)) - H4(c,i);
            end
            H4(c,Rows{c}) = check_node1(H4(c,Rows{c}));
            end
            countv = countv +1;
            HplusC = [Channel;H4];
            sums = sum(HplusC);
            xhat = zeros(1,num_v);
            xhat(sums<0) = 1;
            %err = sum(abs(xhat-codeword));
            %bers(Noise_lvls,l) = err/k;
            
            s = (mod(xhat*matrix.',2));
            syndrome_output = sum(s);   
            if syndrome_output==0
              break
            end
        end


        bers(l) = countf;
        bers2(l) = counts;
        bers3(l) = countm;
        bers4(l) = countv;
        
    end
    bers(bers==100) = [];
    bers2(bers2==100) = [];
    bers3(bers3==100) = [];
    bers4(bers4==100) = [];
    Bers(Noise_lvls) = sum(bers)/(length(bers));
    Bers2(Noise_lvls) = sum(bers2)/(length(bers2));
    Bers3(Noise_lvls) = sum(bers3)/(length(bers3));
    Bers4(Noise_lvls) = sum(bers4)/(length(bers4));
    disp(Noise_lvls)
end



figure

plot(minSNR:SNRstepsize:maxSNR,Bers2,"-square")
ylim([0 inf])
title("Average nuumber of iterations until convergence")
hold on
plot(minSNR:SNRstepsize:maxSNR,Bers4,"-.")
plot(minSNR:SNRstepsize:maxSNR,Bers,"-o")
plot(minSNR:SNRstepsize:maxSNR,Bers3,"-x")
legend("Serial_c SPA", "Serial_c MSA", "Flooding SPA", "Flooding MSA")
hold off

