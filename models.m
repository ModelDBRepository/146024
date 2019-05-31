%MAPK Cascade Model
%Constructed by ZX.Mai.

%Species(20)
    %y1:MAP3K
    %y2:p-MAP3K
    %y3:MAPKK
    %y4:MAPKK:p-MAP3K
    %y5:p-MAPKK
    %y6:p-MAPKK:p-MAP3K
    %y7:pp-MAPKK
    %y8:MAPK
    %y9:MAPK:pp-MAPKK
    %y10:p-MAPK
    %y11:p-MAPK:pp-MAPKK
    %y12:pp-MAPK
    %y13:MKP2
    %y14:p-MAPKK:MKP2
    %y15:pp-MAPKK:MKP2
    %y16:MKP3
    %y17:p-MAPK:MKP3
    %y18:pp-MAPK:MKP3
    %y19:Signal_Ina
    %y20:Signal
    
%Parameters(28)
    %p1:k1/k-1   
    %p2:kb2      
    %p3:kd2     
    %p4:k2       
    %p5:kb3     
    %p6:kd3     
    %p7:k3      
    %p8:kb4     
    %p9:kd4     
    %p10:k4     
    %p11:kb5    
    %p12:kd5    
    %p13:k5     
    %p14:kb-2   
    %p15:kd-2   
    %p16:k-2    
    %p17:kb-3   
    %p18:kd-3   
    %p19:k-3    
    %p20:kb-4   
    %p21:kd-4   
    %p22:k-4    
    %p23:kb-5   
    %p24:kd-5   
    %p25:k-5    
    %p26:ks/k-s 
    %p27:wa     
    %p28:wb     
    
    
function rzt=sim()
  %generateInit();    %Parameters initialization
  options = [];
  y0 = importdata('iniCons.txt');
  p0 = importdata('iniP.txt'); 

  for k = 1:36
      tmpr1=[];
      tmpr2=[];
      tmpr3=[];
      ActMAPK = [];   
      for j = 1:2000
          y0r = y0(:,k)';  %initial cons
          t = [0];
          y0r(19) = 1; y0r(20) = 1;   %signal strength
          pr = p0(:,j);    %kinetic rates
          
          TA_Time = [0,0];
          TA_Amp = 0;
          OSC_Type = 0;
          OSC_Period = 0;
          OSC_Amp = [0,0];
          Section2nd = [];
          Section3rd = [];
          Section4th = [];
          IsSteady = 0;
          for l=1:100
              tspan = [(l-1)*50,l*50];  %simulation time
              [tr,y] = ode15s(@f,tspan,y0r(end,:)',options,pr);
              y0r = [y0r;y(2:end,:)];
              t = [t;tr(2:end)];
              if abs(1-y(1,12)/y(end,12)) < 0.002
                  IsSteady = 1;
                  break;
              end
          end
          if IsSteady == 0    %too slow
              tmpr3 = [tmpr3;[k,j]];
              continue;
          end
              %Activation strength and steady state
          ActMAPK = y0r(:,12);
          [MaxActMAPK,x1stMax] = max(ActMAPK);
          SteadyMAPK = ActMAPK(end);
          
          %Speed of activation
          xStartClimb = find(ActMAPK >= MaxActMAPK*0.1,1,'first');
          xEndClimb = find(ActMAPK >= MaxActMAPK*0.9,1,'first');
          StartClimbTime = t(xStartClimb);
          EndClimbTime = t(xEndClimb);
          ClimbGap = ActMAPK(xEndClimb)-ActMAPK(xStartClimb);
          
          %Drop back
          Section2nd = ActMAPK(x1stMax:end); %Points after 1st Max
          [FirstMinAfterMax,x1stMin] = min(Section2nd);
          x1stMin = x1stMin + x1stMax - 1;
          
          if (x1stMin > x1stMax) && (FirstMinAfterMax <= MaxActMAPK*0.75) %Trancient Activation
              TA_Amp = MaxActMAPK - FirstMinAfterMax;

              HalfUp = MaxActMAPK/2;
              HalfDrop = (MaxActMAPK+FirstMinAfterMax)/2;
              xHalfUp = find(ActMAPK >= HalfUp,1,'first');     
              xHalfDrop = find(Section2nd <= HalfDrop, 1, 'first');
              xHalfDrop = xHalfDrop + x1stMax - 1;

              TA_Time = [t(x1stMax) - t(xHalfUp), t(xHalfDrop) - t(x1stMax)];     %(Time of HalfMax->Max, Time of Max->HalfDropBack)

              Section3rd = ActMAPK(x1stMin:end); %Points after 1st min after 1st max
              [SecondMax,x2ndMax] = max(Section3rd);
              x2ndMax = x2ndMax + x1stMin - 1;

              if SecondMax > 1.1*SteadyMAPK     %Oscillation
                  if SecondMax > 0.9*MaxActMAPK
                      OSC_Type = 1;     %Sustained oscillation
                  else
                      OSC_Type = 2;     %Unsustained oscillation
                  end

                  Section4th = ActMAPK(x2ndMax:end);
                  [SecondMin,x2ndMin] = min(Section4th);
                  x2ndMin = x2ndMin + x2ndMax - 1;
                  OSC_Period = t(x2ndMax) - t(x1stMax);
                  OSC_Amp = [TA_Amp,(SecondMax - SecondMin)];
              end
          end
          if MaxActMAPK > (0.1*y0r(1,8))  %Activation          
              tmpr1 = [tmpr1;[k,j,y0r(1,8),MaxActMAPK/y0r(1,8),SteadyMAPK/y0r(1,8),StartClimbTime,EndClimbTime,ClimbGap/y0r(1,8),TA_Time,TA_Amp/y0r(1,8),OSC_Type,OSC_Period,OSC_Amp/y0r(1,8)]];
          else  %Inactivation
              tmpr2 = [tmpr2;[k,j,y0r(1,8),MaxActMAPK/y0r(1,8),SteadyMAPK/y0r(1,8),StartClimbTime,EndClimbTime,ClimbGap/y0r(1,8),TA_Time,TA_Amp/y0r(1,8),OSC_Type,OSC_Period,OSC_Amp/y0r(1,8)]];
          end
      end
      save('rzt_A_S1_2K.txt','tmpr1','-ascii','-append');
      save('rzt_I_S1_2K.txt','tmpr2','-ascii','-append');
      save('rzt_nonSteady_S1_2K.txt','tmpr3','-ascii','-append');
  end
  %save -ascii 'rzt_A_S10_1W.txt' tmpr1;
  %save -ascii 'rzt_I_S10_1W.txt' tmpr2;
  %save -ascii 'rzt_nonSteady_S10_1W.txt' tmpr3;

%   nAct = importdata('rzt_A_S10_2K.txt');
%   nCon = nAct(:,1);
%   nK = nAct(:,2);
%   nP = length(nCon);
%   tmpMAP3K = [];
%   tmpMAPKK = [];
%   tmpMAPK = [];
%   tmpVector3K = [];
%   tmpVector2K = [];
%   tmpVector = [];
%   
%   tmpSave = 0;
%   
%   for i=1:nP
%       y0n = nCon(i); kn = nK(i);
%       pr = p0(:,kn);
% 
%       s = [0,y0(:,y0n)'];
%       gradient = 0;
%       gradient3K = 0; 
%       gradient2K = 0;
% 
%       j = 0.001;
%       step = 0.001;
%       it = 1;
%       while j < 12
%           y0r = y0(:,y0n);
%           y0r(19) = j; y0r(20) = j;
%           IsSteady = 0;
%           for l=1:100
%               tspan = [(l-1)*50,l*50];  %simulation time
%               [tr,y] = ode15s(@f,tspan,y0r,options,pr);
%               y0r = y(end,:)';
%               if abs(1-y(1,12)/y(end,12)) < 0.002
%                   IsSteady = 1;
%                   break;
%               end
%           end
%           if IsSteady == 0
%             continue;
%           else
%             s=[s;[j,y0r']];
%           end
%           if s(end-1,12) == 0
%               j = j + step;
%               it = it + 1;
%               continue;
%           end
%           delta = abs(1-s(end,12)/s(end-1,12));
%           if (delta < 0.002) || (it>10)
%               step = step*5;
%               it = 1;
%           end
%           if (delta > 0.2) && (step > 0.05)
%               step = step/2;
%               it = 1;
%           end
%           if step > 1 
%               break;
%           end
%           j = j + step;
%           it = it + 1;
%       end
%       
%       
%       [maxMAP3K,max3K100] = max(s(:,3));
%       max3KSignal = s(max3K100,1);
%       
%       max3K10 = find(s(:,3)>=maxMAP3K*0.1,1,'first');
%       max3K90 = find(s(:,3)>=maxMAP3K*0.9,1,'first');
%       max3K10Signal = s(max3K10,1);
%       max3K90Signal = s(max3K90,1);
%       max10MAP3K = s(max3K10,3)/y0(1,y0n);
%       max90MAP3K = s(max3K90,3)/y0(1,y0n);
%       
%       [maxMAPKK,max2K100] = max(s(:,8));
%       max2KSignal = s(max2K100,1);
%       
%       max2K10 = find(s(:,8)>=maxMAPKK*0.1,1,'first');
%       max2K90 = find(s(:,8)>=maxMAPKK*0.9,1,'first');
%       max2K10Signal = s(max2K10,1);
%       max2K90Signal = s(max2K90,1);
%       max10MAPKK = s(max2K10,8)/y0(3,y0n);
%       max90MAPKK = s(max2K90,8)/y0(3,y0n);
%       
%       [maxMAPK,max100] = max(s(:,13));
%       maxSignal = s(max100,1);
%       
%       max10 = find(s(:,13)>=maxMAPK*0.1,1,'first');
%       max90 = find(s(:,13)>=maxMAPK*0.9,1,'first');
%       max10Signal = s(max10,1);
%       max90Signal = s(max90,1);
%       max10MAPK = s(max10,13)/y0(8,y0n);
%       max90MAPK = s(max90,13)/y0(8,y0n);
%       
%       y2 = s(end,2:end)';
%       y2(19)=max10Signal;y2(20)=max10Signal;
%       IsSteady = 0;
%       for l=1:100
%           tspan = [(l-1)*50,l*50];  %simulation time
%           [tr,y] = ode15s(@f,tspan,y2,options,pr);
%           y2 = y(end,:)';
%           if abs(1-y(1,12)/y(end,12)) < 0.002
%               IsSteady = 1;
%               break;
%           end
%       end
%       bistable = IsSteady * y2(12)/y0(8,y0n)/max10MAPK;
% 
%       y3 = s(end,2:end)';
%       y3(19)=max3K10Signal;y0b(20)=max3K10Signal;
%       IsSteady = 0;
%       for l=1:100
%           tspan = [(l-1)*50,l*50];  %simulation time
%           [tr,y] = ode15s(@f,tspan,y3,options,pr);
%           y3 = y(end,:)';
%           if abs(1-y(1,12)/y(end,12)) < 0.002
%               IsSteady = 1;
%               break;
%           end
%       end
%       bistable3K = IsSteady * y3(2)/y0(1,y0n)/max10MAP3K;
% 
%       y4 = s(end,2:end)';
%       y4(19)=max2K10Signal;y4(20)=max2K10Signal;
%       IsSteady = 0;
%       for l=1:100
%           tspan = [(l-1)*50,l*50];  %simulation time
%           [tr,y] = ode15s(@f,tspan,y4,options,pr);
%           y4 = y(end,:)';
%           if abs(1-y(1,12)/y(end,12)) < 0.002
%               IsSteady = 1;
%               break;
%           end
%       end
%       bistable2K = IsSteady * y4(7)/y0(3,y0n)/max10MAPKK;
% 
%       if max3K10Signal~=max3K90Signal
%           gradient3K = (max90MAP3K-max10MAP3K)/(max3K90Signal-max3K10Signal);
%       else
%           gradient3K = 0;
%       end
%       tmpMAP3K = [tmpMAP3K;[y0n,kn,maxMAP3K/y0(1,y0n),s(end,2)/y0(1,y0n),max3KSignal,max3K90Signal,gradient3K,bistable3K,max3K10Signal]];
%       
%       if max2K10Signal~=max2K90Signal
%           gradient2K = (max90MAPKK-max10MAPKK)/(max2K90Signal-max2K10Signal);
%       else
%           gradient2K = 0;
%       end
%       tmpMAPKK = [tmpMAPKK;[y0n,kn,maxMAPKK/y0(3,y0n),s(end,8)/y0(3,y0n),max2KSignal,max2K90Signal,gradient2K,bistable2K,max2K10Signal]];
% 
%       if max10Signal~=max90Signal
%           gradient = (max90MAPK-max10MAPK)/(max90Signal-max10Signal);
%       else
%           gradient = 0;
%       end
%       tmpMAPK = [tmpMAPK;[y0n,kn,maxMAPK/y0(8,y0n),s(end,13)/y0(8,y0n),maxSignal,max90Signal,gradient,bistable,max10Signal]];
%       
%       tmpVector3K = [tmpVector3K;[y0n,kn,s(max3K10,2:end)];[y0n,kn,y3']];
%       tmpVector2K = [tmpVector2K;[y0n,kn,s(max2K10,2:end)];[y0n,kn,y4']];
%       tmpVector = [tmpVector;[y0n,kn,s(max10,2:end)];[y0n,kn,y2']];
%       
%       tmpSave = tmpSave + 1;
%       if tmpSave >= 500
%           save('rzt_US_S10_2K_MAP3K.txt','tmpMAP3K','-ascii','-append');
%           save('rzt_US_S10_2K_MAPKK.txt','tmpMAPKK','-ascii','-append');
%           save('rzt_US_S10_2K_MAPK.txt','tmpMAPK','-ascii','-append');
%           save('rzt_Bistabale_Vector3K.txt','tmpVector3K','-ascii','-append');
%           save('rzt_Bistabale_Vector2K.txt','tmpVector2K','-ascii','-append');
%           save('rzt_Bistabale_Vector.txt','tmpVector','-ascii','-append');

%           tmpSave = 0;
%           tmpMAP3K = [];
%           tmpMAPKK = [];
%           tmpMAPK = [];
%           tmpVector3K = [];
%           tmpVector2K = [];
%           tmpVector = [];
%       end
%   end
%   save('rzt_US_S10_2K_MAP3K.txt','tmpMAP3K','-ascii','-append');
%   save('rzt_US_S10_2K_MAPKK.txt','tmpMAPKK','-ascii','-append');
%   save('rzt_US_S10_2K_MAPK.txt','tmpMAPK','-ascii','-append');
%   save('rzt_Bistabale_Vector3K.txt','tmpVector3K','-ascii','-append');
%   save('rzt_Bistabale_Vector2K.txt','tmpVector2K','-ascii','-append');
%   save('rzt_Bistabale_Vector.txt','tmpVector','-ascii','-append');


%             for i=0.1:0.1:5.0
%                y0r = y0(:,25);
%                kr = p0(:,159);
%                s=[];
%                y0r(8) = i;  %MAPK
%                gradient = 0;
%               for j = 0.001:0.005:0.05
%                   y0r(19) = j; y0r(20) = j;
%                   [t,y]=ode15s(@f,tspan,y0r,options,kr);
%                   s=[s;[j,y(end,12)/y0r(8)]];
%               end
% 
%               for j = 0.1:0.05:1.0
%                   y0r(19) = j; y0r(20) = j;
%                   [t,y]=ode15s(@f,tspan,y0r,options,kr);
%                   s=[s;[j,y(end,12)/y0r(8)]];
%               end
% 
%               for j = 1.5:0.5:5.0
%                   y0r(19) = j; y0r(20) = j;
%                   [t,y]=ode15s(@f,tspan,y0r,options,kr);
%                   s=[s;[j,y(end,12)/y0r(8)]];
%               end
% 
%               for j = 6.0:12.0
%                   y0r(19) = j; y0r(20) = j;
%                   [t,y]=ode15s(@f,tspan,y0r,options,kr);
%                   s=[s;[j,y(end,12)/y0r(8)]];
%               end
% 
%               [maxMAPK,max100] = max(s(:,2));
%               maxSignal = s(max100,1);
% 
% 
%               max10 = find(s(:,2)>=maxMAPK*0.1,1,'first');
%               max90 = find(s(:,2)>=maxMAPK*0.9,1,'first');
%               max10Signal = s(max10,1);
%               max90Signal = s(max90,1);
%               max10MAPK = s(max10,2);
%               max90MAPK = s(max90,2);
% 
%               y0b = y(end,:)';
%               y0b(19)=max10Signal;y0b(20)=max10Signal;
%               [t2,y2] = ode15s(@f,tspan,y0b,options,kr);    %Stable
%               if y2(end,12)/y0r(8) > max10MAPK*1.1
%                   IsBistable = 1;
%               else
%                   IsBistable = 0;
%               end
%               if max10Signal~=max90Signal
%                   gradient = (max90MAPK-max10MAPK)/maxMAPK/(max90Signal-max10Signal);
%               end
%                   tmpr3 = [tmpr3;[i,gradient,IsBistable,maxMAPK,s(end,2),maxSignal,max90Signal,max10Signal]];
%            end
%            save -ascii 'US_MAPK.txt' tmpr3;



function generateInit()
    i=0;
    iniY=[];
    AmpSubstance=[0.2;1;5.0];
    AmpEnzyme=[0.1;1.0];
    for i1=1:3
        tmpY = zeros(20,1);
        tmpY(1) = 1;    %MAP3K
        tmpY(3)=AmpSubstance(i1);   %MAPKK
        
        for i2=1:3
            tmpY(8)=AmpSubstance(i2);   %MAPK
            for j1=1:2
                tmpY(13)=tmpY(3)*AmpEnzyme(j1); %MKP2
                for j2=1:2
                    tmpY(16)=tmpY(8)*AmpEnzyme(j2); %MKP3
                    iniY = [iniY,tmpY];
                end
            end
        end
    end
    save -ascii 'iniConsNew.txt' iniY;

    i=0;
    iniP=[];
    while i<10000
        tmpP = generateP(floor(rand*3));
        IsDiff = 1;
        for j=1:i
            if tmpP==iniP(:,j)
                IsDiff = 0;
                break;
            end
        end
        if IsDiff == 1
            iniP=[iniP,tmpP];
            i=i+1;
        end
    end
    save -ascii 'iniPNew.txt' iniP;
    
    
function P=generateP(fbType)
    tmpP=zeros(28,1);
    %k1=[0.01;0.1;0.5;1.0;2.0;5.0];
    kOri=1;
    %k24=[0.01;0.05;0.1;0.5;1.0;2.0];
    %kb=[0.1;0.5;1.0;2.0;4.0;10.0];
    
    AmpKcat=[0.01;0.1;0.5;1.0;2.0;5.0];
    AmpBD=[0.1;1.0;10;100];
    AmpKs=[0.01;0.1;1.0;5.0];
    AmpWab=[0.2;1.0;5.0];

    tmpP(1)=kOri;
    
    tmpP(2)=kOri*AmpBD(1+floor(rand*4));
    tmpP(3)=kOri*AmpBD(1+floor(rand*4));
    %tmpP(3)=tmpP(2)*AmpBD(1+floor(rand*3));
    tmpP(4)=kOri*AmpKcat(1+floor(rand*6));
    
    tmpP(5)=kOri*AmpBD(1+floor(rand*4));
    tmpP(6)=kOri*AmpBD(1+floor(rand*4));
%    tmpP(6)=tmpP(5)*AmpBD(1+floor(rand*3));
    tmpP(7)=kOri*AmpKcat(1+floor(rand*6));
    
    tmpP(8)=kOri*AmpBD(1+floor(rand*4));
    tmpP(9)=kOri*AmpBD(1+floor(rand*4));
    %tmpP(9)=tmpP(8)*AmpBD(1+floor(rand*3));
    tmpP(10)=kOri*AmpKcat(1+floor(rand*6));
    
    tmpP(11)=kOri*AmpBD(1+floor(rand*4));
    tmpP(12)=kOri*AmpBD(1+floor(rand*4));
    %tmpP(12)=tmpP(11)*AmpBD(1+floor(rand*3));
    tmpP(13)=kOri*AmpKcat(1+floor(rand*6));
    
    tmpP(14)=kOri*AmpBD(1+floor(rand*4));
    tmpP(15)=kOri*AmpBD(1+floor(rand*4));
    %tmpP(15)=tmpP(14)*AmpBD(1+floor(rand*3));
    tmpP(16)=kOri*AmpKcat(1+floor(rand*6));
    %tmpP(16)=tmpP(4)*AmpKcat(3+floor(rand*3));
    
    tmpP(17)=kOri*AmpBD(1+floor(rand*4));
    tmpP(18)=kOri*AmpBD(1+floor(rand*4));
%    tmpP(18)=tmpP(17)*AmpBD(1+floor(rand*3));
    tmpP(19)=kOri*AmpKcat(1+floor(rand*6));
%    tmpP(19)=tmpP(7)*AmpKcat(3+floor(rand*3));
    
    tmpP(20)=kOri*AmpBD(1+floor(rand*4));
    tmpP(21)=kOri*AmpBD(1+floor(rand*4));
    %tmpP(21)=tmpP(20)*AmpBD(1+floor(rand*3));
    tmpP(22)=kOri*AmpKcat(1+floor(rand*6));
%    tmpP(22)=tmpP(10)*AmpKcat(3+floor(rand*3));
    
    tmpP(23)=kOri*AmpBD(1+floor(rand*4));
    tmpP(24)=kOri*AmpBD(1+floor(rand*4));
    %tmpP(24)=tmpP(23)*AmpBD(1+floor(rand*3));
    tmpP(25)=kOri*AmpKcat(1+floor(rand*6));
%   tmpP(25)=tmpP(13)*AmpKcat(3+floor(rand*3));
    
    if fbType>0 
        tmpP(26)=kOri*AmpKs(1+floor(rand*4));
        if fbType==1
            tmpP(27)=kOri*AmpWab(1+floor(rand*3));
        elseif fbType==2
            tmpP(28)=kOri*AmpWab(1+floor(rand*3));
        end
    end
    
    P=tmpP;



function dydt=f(t,y,p)
    v1=p(1)*y(20)*y(1);
    v21=p(2)*y(2)*y(3)-p(3)*y(4);
    v22=p(4)*y(4);
    v31=p(5)*y(2)*y(5)-p(6)*y(6);
    v32=p(7)*y(6);
    v41=p(8)*y(7)*y(8)-p(9)*y(9);
    v42=p(10)*y(9);
    v51=p(11)*y(7)*y(10)-p(12)*y(11);
    v52=p(13)*y(11);
    
    vm1=p(1)*y(2);
    vm21=p(14)*y(13)*y(5)-p(15)*y(14);
    vm22=p(16)*y(14);
    vm31=p(17)*y(13)*y(7)-p(18)*y(15);
    vm32=p(19)*y(15);
    vm41=p(20)*y(16)*y(10)-p(21)*y(17);
    vm42=p(22)*y(17);
    vm51=p(23)*y(16)*y(12)-p(24)*y(18);
    vm52=p(25)*y(18);
    
    va=(p(26)+p(27)*y(12))*y(19);
    vb=(p(26)+p(28)*y(12))*y(20);
    
    %-------------------------------------
    dy1=-v1+vm1;
    dy2=v1-vm1-v21+v22-v31+v32;
    dy3=-v21+vm22;
    dy4=v21-v22;
    dy5=v22-v31-vm21+vm32;
    dy6=v31-v32;
    dy7=v32-v41+v42-v51+v52-vm31;
    dy8=-v41+vm42;
    dy9=v41-v42;
    dy10=v42-v51-vm41+vm52;
    dy11=v51-v52;
    dy12=v52-vm51;
    dy13=-vm21+vm22-vm31+vm32;
    dy14=vm21-vm22;
    dy15=vm31-vm32;
    dy16=-vm41+vm42-vm51+vm52;
    dy17=vm41-vm42;
    dy18=vm51-vm52;
    dy19=vb-va;
    dy20=va-vb;
    
    dydt=[dy1;dy2;dy3;dy4;dy5;dy6;dy7;dy8;dy9;dy10;dy11;dy12;dy13;dy14;dy15;dy16;dy17;dy18;dy19;dy20];
    