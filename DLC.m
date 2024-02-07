clear
clc
close all

%PARAMETERS
xcoordinates=[];
ycoordinates=[];
lengthofthepath=[];
time=[];
Rforward=[];
delta_fforward=[];
Rbackward=[];
delta_fbackward=[];
Rcentral=[];
delta_fcentral=[];
Rinput=[];
delta_finput=[];

%INPUT DATA
V=25;
L=1.545;
stepsize=0.001; %While discretizing

%PATH DETAILS (DLC)
xi=5;
xlc1=5; ylc1=2.5;
xlc2=5; ylc2=2.5;
xr=20;
xlc3=5; ylc3=2.5;
xlc4=5; ylc4=2.5;
xe=5;

%POINTS ON THE PATH
x1=0; y1=0;
x2=xi; y2=0;
x3=x2+xlc1; y3=y2+ylc1;
x4=x3+xlc2; y4=y3+ylc2;
x5=x4+xr; y5=y4;
x6=x5+xlc3; y6=y5-ylc3;
x7=x6+xlc4; y7=y6-ylc4;
x8=x7+xe; y8=y7;

%RADIUS OF CURVATURE
R1=(xlc1^2+ylc1^2)/(2*ylc1);
R2=(xlc2^2+ylc2^2)/(2*ylc2);
R3=(xlc3^2+ylc3^2)/(2*ylc3);
R4=(xlc4^2+ylc4^2)/(2*ylc4);

%LENGTH OF THE PATH BETWEEN THE DEFINED POINTS
S12=xi;
S23=2*R1*asin((sqrt(ylc1^2+xlc1^2))/(2*R1));
S34=2*R2*asin((sqrt(ylc2^2+xlc2^2))/(2*R2));
S45=xr;
S56=2*R3*asin((sqrt(ylc3^2+xlc3^2))/(2*R3));
S67=2*R4*asin((sqrt(ylc4^2+xlc4^2))/(2*R4));
S78=xe;

S1=0;
S2=S1+S12;
S3=S2+S23;
S4=S3+S34;
S5=S4+S45;
S6=S5+S56;
S7=S6+S67;
S8=S7+S78;

%time
t1=0;
t2=S2/V;
t3=S3/V;
t4=S4/V;
t5=S5/V;
t6=S6/V;
t7=S7/V;
t8=S8/V;

%XY-COORDINATES W.R.T. TIME
for t=t1:stepsize:t8
    S=t*V;
    %   1 to 2
    S=t*V;
    x=S;
    y=y1;
    %   2 to 3
    if t>t2
        x=x2+R1*sin((t-t2)*V/R1);
        y=R1-sqrt(R1^2-(x-x2)^2);
        % 3 to 4
        if t>t3
            x=x3+xlc2-R2*sin((S34-(t-t3)*V)/R2);
            y=-R2+y4+sqrt(R2^2-(x-x4)^2);
            % 4 to 5
            if t>t4
                x=x4+(t-t4)*V;
                y=ylc1+ylc2;
                %   5 to 6
                if t>t5
                    x=x5+R3*sin((t-t5)*V/R3);
                    y=-R3+y4+sqrt(R3^2-(x-x5)^2);
                    %   6 to 7
                    if t>t6
                        x=x6+xlc4-R4*sin((S67-(t-t6)*V)/R4);
                        y=R4+ylc1+ylc2-ylc3-ylc4-sqrt(R4^2-(x-x7)^2);
                        %   7 to 8
                        if t>t7
                            x=x7+(t-t7)*V;
                            y=ylc1+ylc2-ylc3-ylc4;
                        end
                    end
                end
            end
        end
     end
    xcoordinates=[xcoordinates;x];
    ycoordinates=[ycoordinates;y];
    lengthofthepath=[lengthofthepath;S];
    time=[time;t];
    %xlswrite('DLCpath.xlsx',[time,xcoordinates,ycoordinates,lengthofthepath]);
end
 
%%RADIUS OF CURVATURE && FRONT STEERING ANGLE
    % Forward difference
for i=1:(((t8-t1)/stepsize)+1)
    if i>((t8-t1)/stepsize-1)
        Rfwd=1/0;
        delta_ffwd=0;
    else
        Rfwd=((1+((ycoordinates(i+1)-ycoordinates(i))/stepsize)^2)^1.5)/((ycoordinates(i+2)-2*ycoordinates(i+1)+ycoordinates(i))/stepsize^2);
        %disp(Rfwd);
        delta_ffwd=atan(L/Rfwd);
        %disp(delta_ffwd);  
    end
    Rforward=[Rforward;Rfwd];
    delta_fforward=[delta_fforward;delta_ffwd];
    %xlswrite('FwdDiffdelta_f.xlsx',[Rforward,delta_fforward]);
end

% Backward difference
for i=1:((t8-t1)/stepsize+1)
    if i<3
        Rbwd=1/0;
        delta_fbwd=0;
    else
        Rbwd=((1+((ycoordinates(i)-ycoordinates(i-1))/stepsize)^2)^1.5)/((ycoordinates(i)-2*ycoordinates(i-1)+ycoordinates(i-2))/stepsize^2);
        %disp(Rbwd);
        delta_fbwd=atan(L/Rbwd);
        %disp(delta_fbwd); 
    end
    Rbackward=[Rbackward;Rbwd];
    delta_fbackward=[delta_fbackward;delta_fbwd];
    %xlswrite('BwdDiffdelta_f.xlsx',[Rbackward,delta_fbackward]);
end

%Central difference
for i=1:((t8-t1)/stepsize+1)
    if i<2
        Rcnt=1/0;
        delta_fcnt=0;
    elseif i>((t8-t1)/stepsize)
       Rcnt=1/0;
       delta_fcnt=0;
    else
        Rcnt=((1+((ycoordinates(i+1)-ycoordinates(i-1))/(2*stepsize))^2)^1.5)/((ycoordinates(i+1)-2*ycoordinates(i)+ycoordinates(i-1))/stepsize^2);
        delta_fcnt=atan(L/Rcnt);
    end    
    Rcentral=[Rcentral;Rcnt];
    delta_fcentral=[delta_fcentral;delta_fcnt];
    disp([Rcnt,delta_fcnt]);
    xlswrite('CntDiffdelta_f.xlsx',[Rcentral,delta_fcentral]);
end

%Direct as INPUT _ Real Radius and delta_f
for t=t1:stepsize:t8
    Rip=1/0;
    if t>x2
        Rip=R1;
        if t==t3
            Rip=1/0;
        end    
        if t>t3
            Rip=-R2;
            if t>=t4
                Rip=1/0;
                if t>t5
                    Rip=-R3;
                    if t==t6
                        Rip=1/0;
                    end
                    if t>t6
                        Rip=R4;
                        if t>=t7
                            Rip=1/0;
                        end
                    end
                end
            end
        end
    end
    %disp([x,Rip]);
    delta_fip=atan(L/Rip);
    Rinput=[Rinput;Rip];
    delta_finput=[delta_finput;delta_fip];
    xlswrite('InputbBaseddelta_f.xlsx',[Rinput,delta_finput]);
end

%disp(["What next ??:";"1)repair central difference scheme";"2)compare with manual calculations at all 8 points and randomly between them"]);

%%PLOT
figure
plot(xcoordinates,time,'r-','LineWidth',3 );
xlabel('X(m)');
ylabel('t_r y_g S_b');
grid on;
title('DLC');
hold on
plot(xcoordinates,ycoordinates,'g-','LineWidth',3 );
plot(xcoordinates,lengthofthepath,'b-','LineWidth',3 );
hold off
