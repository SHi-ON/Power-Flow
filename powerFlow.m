% Guass-Seidel Power Flow program.
% Developed by:
%                              Shayan Amani 
%                              shayanamani.com
% Contact:
%                              mail@shayanamani.com
clear all;
close all;

Vb = 10e3;       % Vbase (V)
Sb = 200e6;      % Sbase (VA)
Zb = (Vb.^2)./Sb;       % Zbase (Ohm)
Yb = 1./Zb;             % Ybase (Mho)

% the Line Impedances matrix given as follow:
%   Line    From    To  Z(ohm)
linedata = [ 1 1 2 0.04+0.0281*1i
             2 1 4 0.061+0.0304*1i
             3 1 5 0.0951+0.0064*1i
             4 2 3 0.002+0.0108*1i
             5 3 4 0.001+0.0297*1i
             6 4 5 0.005+0.0297*1i ];
         
% here Ybus matrix is formed :  
 ybus = zeros(5);
 for k = 1:5 
   for m = 1:5
         if linedata(m,2) == k || linedata(m,3) == k
             ybus(k,k)= ybus(k,k) + ((1/linedata(m,4))./Yb);
         else
             ybus(linedata(m,2),linedata(m,3)) = ((-1/linedata(m,4))./Yb);
             ybus(linedata(m,3),linedata(m,2)) = ybus(linedata(m,2),linedata(m,3));
         end
   end
 end

% p.u     Bus  Type    Pg      Qg      Pd      Qd       V        Qg_max        Qg_min
busdata = [ 1   1       0       0       0       0     1+0*1i         1.25        -0.5
            2   3       0       0      0.5    0.25    1+0*1i           0           0
            3   3     1.75     0.5     0.75   0.2     1+0*1i         0.5         -0.25
            4   2     1.75      0      2.5      1   1.015+0*1i       1.25        -0.5    
            5   3       1      0.35     2     0.5     1+0*1i         0.5         -0.3 
           ];  
 
P = busdata(:,3) - busdata(:,5);
Q = busdata(:,4) - busdata(:,6);
Vold = busdata(:,7);
Vnew = 2.*Vold;
err = 0.0001;
counter = 0;

ReCon = max(abs(real(Vnew-Vold)));
ImCon = max(abs(imag(Vnew-Vold)));

while ReCon > err || ImCon > err 
        counter = counter + 1;
        if counter == 1
            Vnew = Vold;
        end
        for b = 2:5
            if busdata(b,2) == 2
                Q(b) = -imag(conj(Vnew(b)) * ybus(b,:) * Vnew);
                Vnew(b) = (1/ybus(b,b))*(((P(b)-Q(b)*1i)/conj(Vnew(b))) - ( ybus(b,:)*Vnew - ybus(b,b)*Vnew(b)));
                Vnew(b) = (Vnew(b)/abs(Vnew(b))) * busdata(b,7);
                Q(b) = -imag(conj(Vnew(b)) * ybus(b,:) * Vnew);
                if Q(b) > busdata(b,9) - busdata(b,6)                   
                    Q(b) = busdata(b,9) - busdata(b,6);
                     busdata(b,2) = 3; 
                end
                if Q(b) > busdata(b,8) - busdata(b,6)
                    Q(b) = busdata(b,8) - busdata(b,6);
                    busdata(b,2) = 3;
                end
            end
            if busdata(b,2) == 3
                Vnew(b) = (1/ybus(b,b))*(((P(b)-Q(b)*1i)/conj(Vnew(b))) - ( ybus(b,:)*Vnew - ybus(b,b)*Vnew(b)));   
                
            end
        end
        ReCon = max(abs(real(Vnew-Vold)));
        ImCon = max(abs(imag(Vnew-Vold)));
        Vold = Vnew;
        dataVnew(:,counter)=Vold;       %#ok<*SAGROW>

end

 hold on
 title('Magnitudes of Bus Voltages Graph (V)');
 plot(abs(dataVnew(1,:)),'-xr');
 plot(abs(dataVnew(2,:)),'-sc');
 plot(abs(dataVnew(3,:)),'-om');
 plot(abs(dataVnew(4,:)),'-vb');
 plot(abs(dataVnew(5,:)),'-*g');
 xlabel('Iteration step');
 ylabel('Vbus (p.u.)');
 
 figure() ;
 hold on
 title('Angles of Bus Voltages Graph (deg)');
 plot(rad2deg(angle(dataVnew(1,:))),'-xr');
 plot(rad2deg(angle(dataVnew(2,:))),'-sc');
 plot(rad2deg(angle(dataVnew(3,:))),'-om');
 plot(rad2deg(angle(dataVnew(4,:))),'-vb');
 plot(rad2deg(angle(dataVnew(5,:))),'-*g');
 xlabel('Iteration step');
 ylabel('Vbus (p.u.)');
 
 
 Q(1) = -imag(conj(Vnew(1)) * ybus(1,:) * Vnew);
 P(1) = real(conj(Vnew(1)) * ybus(1,:) * Vnew);
 SlackP = [ 'The active power generation of Slack bus (PG Bus 1) is  :  ', num2str(P(1)), '  MW']; 
 SlackQ = [ '             The reactive power generation of Slack bus (QG Bus 1) is  :  ', num2str(Q(1)), '  MVar'];
 txt=[SlackP,SlackQ];
 msgbox(txt,'Results','warn')

 Qg4 = [ 'The reactive power generation of P-|V| bus (Qg Bus 4) is  :  ',  num2str(Q(4) - busdata(4,6)), '  MVar' ];
 msgbox(Qg4,'Results','warn')
 
 nIter = [ 'Number of iteration steps is : ', num2str(counter)];
 msgbox(nIter,' Iteration steps ','warn')
 
 VLastM = [ '|V Bus 1| = '
            '|V Bus 2| = '
            '|V Bus 3| = '
            '|V Bus 4| = '
            '|V Bus 5| = '
            ];     
 VLastMUnt = [ '  (V) ' 
               '  (V) '
               '  (V) '
               '  (V) '
               '  (V) '
               ];      
 txtVLM = [ VLastM, num2str(rad2deg(abs(Vnew))), VLastMUnt];
 msgbox(txtVLM,'|V|')
 
 VLastA = [ '<V Bus 1 = '
            '<V Bus 2 = '
            '<V Bus 3 = '
            '<V Bus 4 = '
            '<V Bus 5 = '
           ]; 
 VLastAUnt = [ '  deg ' 
               '  deg '
               '  deg '
               '  deg '
               '  deg '
             ]; 
 txtVLA = [ VLastA, num2str(rad2deg(angle(Vnew))), VLastAUnt];      
 msgbox(txtVLA,'<V')         
            
