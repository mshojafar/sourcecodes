for i=1:1:n
        if(S(i).E>0)
            temp_rand=rand;
            if ( (S(i).G)<=0)
                
                %Election of Cluster Heads for normal nodes
                if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )
                    
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                    
                    S(i).type='C';
                    S(i).G=100;
                    C(cluster).xd=S(i).xd;
                    C(cluster).yd=S(i).yd;
                    plot(S(i).xd,S(i).yd,'k*','MarkerEdgeColor','b',...
                    'MarkerFaceColor','b',...
                    'MarkerSize',25);
                    
                    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                    C(cluster).distance=distance;
                    C(cluster).id=i;
                    X(cluster)=S(i).xd;
                    Y(cluster)=S(i).yd;
                    cluster=cluster+1;
                    
                    %Calculation of Energy dissipated
                    distance;
                    if (distance>do)
                        S(i).E=S(i).E-((ETX+EDA)*(4000)+ Emp*4000*( distance^4));
                    end
                    if (distance<=do)
                        S(i).E=S(i).E-((ETX+EDA)*(4000)+ Efs*4000*( distance^2 ));
                    end
                end
                
                
                
                %Election of Cluster Heads for Advanced nodes
                if( ( S(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
                    
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                    
                    S(i).type='C';
                    S(i).G=100;
                    C(cluster).xd=S(i).xd;
                    C(cluster).yd=S(i).yd;
                    plot(S(i).xd,S(i).yd,'k*','MarkerEdgeColor','r',...
                    'MarkerFaceColor','r',...
                    'MarkerSize',25);
                    
                    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                    C(cluster).distance=distance;
                    C(cluster).id=i;
                    X(cluster)=S(i).xd;
                    Y(cluster)=S(i).yd;
                    cluster=cluster+1;
                    
                    %Calculation of Energy dissipated
                    distance;
                    if (distance>do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                    end
                    if (distance<=do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                    end
                end
                
            end
        end
    end