%A magnetic dipole system-Ising model.
%Consider a two-dimensional system of dipoles arrayed on a square lattice.
%At each lattice point, the dipole can assume only one of two spin states, 
%and the dipoles interact with an externally applied constant magnetic field B and also with their nearest neighbors. 
%For a given temperature T, we shall observe the magnetization of the system.
%We use Monte Carlo M(RT)^2 to come up with samples.
%Y.S 2017


clear all;
N=100; %%% Number of samples
sample=80; %%% every k pass to get a sample
J=1; %%% The interaction constant
H=0; %%% The magematic constant-no external field
L=25;%%% The size of lattice
T=[1.0:0.1:4.0]; %%%  tempreature
%%% Generate the lattice of spins 1 or -1
Matrix=round(rand(L,L));
Matrix(Matrix==0)=-1;
%%% burn -in to eliminate effects of initial condition 
burn=300;
KbT=1;

for i=1:25*25*burn
            position_x=unidrnd(25);
            position_y=unidrnd(25);
            Update=Matrix(position_x,position_y);
            %%% Talk about the boundary condition
            position_x_up=position_x+1;
            position_x_down=position_x-1;
            position_y_up=position_y+1;
            position_y_down=position_y-1;

            if position_x_up==26
                position_x_up=1;
            end
            if position_x_down==0
                position_x_down=25;
            end
            if position_y_up==26
                position_y_up=25;
            end
            if position_y_down==0
                position_y_down=25;
            end
            %evaluate the change in energy caused by one fip 
            delta_u=2*J*Update*(Matrix(position_x,position_y_up)+Matrix(position_x,position_y_down)+Matrix(position_x_up,position_y)+Matrix(position_x_down,position_y));
            
            %M(RT)^2
            if delta_u<=0
                Matrix(position_x,position_y)=-Update;
            else
                ratio=min(exp(-delta_u/KbT),1);
                if rand()<ratio
                    Matrix(position_x,position_y)=-Update;
                end
            end
end


for t=1:31
    for count=1:sample*N
        Kb=1;
        KbT=Kb*T(t);    

        for i=1:25*25
            position_x=unidrnd(25);
            position_y=unidrnd(25);
            Update=Matrix(position_x,position_y);
            position_x_up=position_x+1;
            position_x_down=position_x-1;
            position_y_up=position_y+1;
            position_y_down=position_y-1;

            if position_x_up==26
                position_x_up=1;
            end
            if position_x_down==0
                position_x_down=25;
            end
            if position_y_up==26
                position_y_up=25;
            end
            if position_S_down==0
                position_S_down=25;
            end
            
            %Metropolis
            delta_u=2*J*Update*(Matrix(position_x,position_y_up)+Matrix(position_x,position_y_down)+Matrix(position_x_up,position_y)+Matrix(position_x_down,position_y));
            if delta_u<=0
                Matrix(position_x,position_y)=-Update;
            else
                ratio=min(exp(-delta_u/KbT),1);
                if rand()<ratio
                    Matrix(position_x,position_y)=-Update;
                end
            end
        end
        
        if mod(count,sample)==0
            %%%compute energy
                  u=0;
                  for i=1:L
                        i_up=i+1;
                        i_down=i-1;
                        if i_up==26
                            i_up=1;
                        end
                         if i_down==0
                                i_down=25;
                         end
                        for k=1:L
                            k_up=k+1;
                            k_down=k-1;
                                if k_up==26
                                    k_up=1;
                                end
                                if k_down==0
                                    k_down=25;
                                end
                            u=-J*0.5*Matrix(i,k)*(Matrix(i_up,k)+Matrix(i_down,k)+Matrix(i,k_up)+Matrix(i,k_down))+u;
                        end
                  end
              
                    Engergy(count/sample)=u;
                    Mag(count/sample)=sum(sum(Matrix))/(L*L);
        end
    end
    Fin_Engergy(t)=mean(Engergy);
    Specific(t)=var(Engergy)/(KbT*KbT);
    New_M=abs(Mag);
    Fin_Mag(t)=mean(New_M);
    Suscep(t)=var(New_M)/(KbT*KbT);
end

figure (1);
plot (T,Fin_Engergy,'-o');
figure (2);
plot (T,Specific,'-o');
figure (3);
plot (T,Fin_Mag,'-o');
figure (4);
plot (T,Suscep,'-o');
