%%%%%%%%%%%%%%%%             Cycle tests     wind and solar and Salt  UQS
lplot = 1;
nplot = 6;
ptitle = 'Sensitivities to wind, SWmax and Salt ';
%                                           a (wind   )  and  b(solar  )     parameter ranges
amin =  1.;
amax =  6.0;
bmin =  100.;
bmax =  320.;
na   = 26;
nb   = 12;
a = linspace(0,1,na) * (amax-amin) / (na-1)   + amin; % various U10 value
b = linspace(0,1,nb) * (bmax-bmin) / (nb-1)   + bmin;  % 10 times U10
%                                                constants and psuedo constants
SST = 30.0;
molvisc  =  1.8e-6 * exp(-SST/35.);
molPr    =  11.8   * exp(-SST/39.);
rhocpw = 4.1e6; % rho * heat capacity of water
grav   = 9.81;
alpha  = 3.0e-4; % coefficient of thermal expansion of water sw_alpha(salinity,temperature,pressure,'temp')
betaS  = 7.5e-4; %saline contraction coefficient  sw_beta(salinity,temperature,pressure,'temp')
sal    = 34.7; % salinity
vonk   = 0.4; % von Karman...
tiny =  1.0e-9;
%                                                 parameters
lambdaC  =  6.5;
d  = 3.5;
p  = 0.2;
pv = p;
LambdaL =   0.;
doLmax   = 1.0; % maximum ratio of d/L

kappa0  =  0.2e-4;
NUzero  = 40.0e-4;
Prandtl =   1.;
Rizero  =  1.0;
F0   =  0.60;
R1      =  10.;
F1   =  0.1;
%                                                       Equations (18) (19)
AMP   =  1./F0 - 1.;
pexp  =  log( ( 1./F1 - F0) / (1.-F0) )  /  log(R1);
%                                                         time stepping
niter = 2;
ndays = 1;
nstps = 24. * 60.;
dt    = 86400./(nstps);
tsec  = linspace(0,1,nstps) * dt  + dt;
thr   = tsec/3600.;
hour  = 0.5 + linspace(0,1,nstps+1) * 24. / nstps;
%                                                              forcing
C_D = 0.0012; % drag coefficient of air
rhoair = 1.22; % density of air
rhowater = 1030; % density of water
Qsmax  = 875.;
Pmmpd  =  100.0; % precipitation (milimeter per day)
precip = Pmmpd/0.084; % precipitation m/s
evap   =  -000/ 2.5; % evaporation

% output arrays
Tpeak = NaN( na, nb );
T1400 = NaN( na, nb );
Vpeak = NaN( na, nb );
Tave  = NaN( na, nb );
Speak = NaN( na, nb );
hpeak = NaN( na, nb );
hfirst = NaN( na, nb );
hlast  = NaN( na, nb );
hplot  = NaN( na, nb );
R1f    = NaN( na, nb );
R2f    = NaN( na, nb );
R3f    = NaN( na, nb );
R4f    = NaN( na, nb );
Smix   = NaN( na, nb );
%                                               wind   (a)  and  solar(b)  loops
for id = 1:na
    zd = -d;
    %                                               Appendix  EQ(20)
    fabsd  = 0.67 * exp(zd/1.) + 0.33 * exp(zd/17.); % Fraction of solar radiation absorbed.
    U10 =  a(id);
    ustar = U10 * sqrt((rhoair * C_D)/ rhowater);
    for ip = 1:nb
        Qsmax  =  4.*  b(ip);
        
%         print , ' '
%         print , ' U10   = ',  U10    , '    Swmax  = ' ,   Qsmax , ' d, p,  pv = ',d, p,  pv
        %                                     1 day integration loop
        Reg = 1.0;
        Wn  = 0.0;
        Sn  = 0.0;
        Vn  = 0.0;
        Tpeak( id , ip ) = 0.0;
        Tave (id  , ip ) = 0.0;
        hpeak( id , ip ) = 0.0;
        hfirst( id , ip ) = 0.0;
        hlast( id , ip ) = 24.0;
        R1f( id , ip ) = 0.0;
        R2f( id , ip ) = 0.0;
        R3f( id , ip ) = 0.0;
        R4f( id , ip ) = 0.0;
        Smix(id , ip ) = 0.0;
        FofRi = 0.0;
        Qbar = 0.0;
        
        
        for ns = 1:nstps
            if mod(id,10)==0&&mod(ns,500)==0&&mod(ip,10)==0
                id
                ns
            end
            Qsol = 0.0;
            %    Qnsol = -75.  - 10. * U10
            Qnsol = -135;
            precip = 0.0;
            
            if ( ( hour(ns) > 6. ) && ( hour(ns) < 18. ) )
                fsine = sin ( 2. * pi * (hour(ns) - 6.) / 24. );
                Qsol = Qsmax * fsine * fsine;
            end
            Qbar = Qbar + Qsol / nstps;
            %                                                                   initialize i = 1
            Wnp1  = Wn;
            Snp1  = Sn;
            Vnp1  = Vn;
            %                                                                     Iteration
            for it = 1:niter
                
                %                                       Cool Skin Equations (2) (3)
                Dcool = lambdaC *  molvisc / ustar;
                Qcool = Qnsol + Qsol * (0.137 + 11.*Dcool - 6.6e-5 / Dcool * ( 1. - exp(-Dcool/8.e-4) ) );
                DTcool = min( [0.0 ,  molPr * Qcool * (lambdaC / ustar) /  rhocpw ] );
                %                                            Warming Equations (7), (A22*), (A23), (9), (15)
                Hd = (Qnsol + Qsol * (1.-fabsd) ) / rhocpw;
                Fd = (precip+evap) * (sal+Snp1)  /1028. / 10e6;
                doL   = d*vonk*grav*(alpha*Hd + betaS*Fd ) / ustar.^3;
                Rid = grav*(alpha*Wnp1-betaS*Snp1)*p*d / (pv*(Vnp1+tiny)).^2;
                Rid = max( [ 0.0 ,  Rid ] );
                %                                                REGIMES
                if ( doL <= 0.0)
                    %                                                                           dol <= 0.0
                    if ( Wn <= 0.0 )
                        REG   = 1;
                        R1f( id , ip )   = R1f( id , ip ) + 100./nstps/niter;
                        Smult = 0.0;
                        Sfact = 0.0;
                        difs  = 0.0;
                        visc  = 0.0;
                    else
                        REG   = 4;
                        R4f( id , ip )   = R4f( id , ip ) + 100./nstps/niter;
                        Smult = ( p + 1 ) / ( d * p );
                        Sfact = ( p + 1 ) /  d^2;
                        %                                             Equation (16)
                        ghat = 0.0;
                        difs  =  vonk * ustar * d * (1. - 50. * doL)^(1/3);
                        visc  =  difs;
                        Smult = Smult * ( 1. + difs * ghat );
                    end
                else
                    %                                                                       0.0 < doL
                    phid  =   min ( [   1. +  5. * doL , 5. + doL ] );
                    % Equation (11),  (18)
                    FofRi  = 1. / (1. + Amp * (Rid/Rizero)^pexp  );
                    dif3 =   (kappa0 + NUzero * FofRi );
                    
                    if  ( (doL < LambdaL ) &&  (REG <= 2) )
                        REG =2;
                        R2f( id , ip )   = R2f( id , ip ) + 100./nstps/niter;
                        Smult = ( p + 1 ) / ( d * p );
                        Sfact = ( p + 1 ) /  d^2;
                        % Equation (10)
                        difs  =  vonk * ustar * d /  phid;
                        visc  = difs * (1. - doL/LambdaL )^2 + dif3 * (doL/LambdaL)^2 * ( 3. - 2. * doL / LambdaL);
                        difs  = visc;
                    else
                        % LambdaL <= doL
                        REG   = 3;
                        R3f( id , ip )   = R3f( id , ip ) + 100./nstps/niter;
                        if ( FofRi  > 0.0 )
                            Smix (id, ip) = Smix (id, ip) + 100./nstps/niter;
                        end
                        Smult = ( p + 1 ) / ( d * p );
                        Sfact = ( p + 1 ) /  d^2;
                        difs =             kappa0 + NUzero * FofRi;
                        visc  =  Prandtl * kappa0 + NUzero * FofRi;
                    end
                end
                %  Equations (24, 25, 26% Appendix)
                Wnp1 = (Wn + ( Smult * dt *  Hd) )       / (1.+  dt * Sfact  * difs);
                Snp1 = (Sn - ( Smult * dt *  Fd) )       / (1.+  dt * Sfact  * difs);
                Vnp1 = (Vn + ( Smult * dt * ustar^2) )   / (1.+  dt * Sfact  * visc);
                % iterate
            end
            
            
            %   print, hour(ns),Reg,Wnp1,Vnp1,Snp1,Tpeak(id,ip),Speak(id,ip)
            %   print, hour(ns),Reg,Wnp1,Vnp1,  Rid, FofRI, difs , Tpeak(id,ip)
            %   print, hour(ns),Reg,Wnp1,Vnp1, doL, LambdaL ,  Tpeak(id,ip)
            % print, hour(ns),Reg,Vnp1,Wnp1, Hd*rhocpw,  doL,  Tpeak(id,ip)
            % print, hour(ns),Reg, Wnp1, DTcool
            %                                                                 plotting results for a timestep
            Tave ( id, ip )  =  Tave ( id, ip )  +  Wnp1 / nstps;
            if ( Wnp1  > Tpeak(id,ip) )
                Tpeak(id, ip) = Wnp1;
                hpeak(id, ip) = hour(ns);
            end
            if ( ns == (14 * 60)  )
                T1400(id,ip) = Wnp1;
            end
            signmFd = 1. -  (abs(Fd) + Fd ) / Fd;
            if ( signmFd*Snp1  > signmFd*Speak(id,ip) )
                Speak(id, ip) = Snp1;
            end
            if ( Vnp1  > Vpeak(id,ip) )
                Vpeak(id, ip) = Vnp1;
            end
            
            Wn = max( [ Wnp1 , 0.0 ] );
            Sn =        Snp1;
            Vn = max( [ Vnp1 , 0.0 ] );
            %                                                                   end timestep after saving updates
        end
%         print, 'Qsmax =', Qsmax, '    Qbar= ', Qbar , '   U10 = ', U10
        %print,  Tpeak
    end
end




