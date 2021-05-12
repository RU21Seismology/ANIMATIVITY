function Lambda=BuildGeneral(para)
%% Documentation
% Function to build up the general elastic tensor using different
% parameterizations
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
    eye=diag(ones(1,3));
    I=kron(eye,eye);
    
    %allocate the memory for Lambda
    Lambda=zeros(6,6,size(para,2));

for kk=1:size(para,2)
%read in the parameters
    A=para(1,kk);
    B=para(2,kk);
    C=para(3,kk);
    D=para(4,kk);
    E=para(5,kk);
    xi_deg=para(6,kk);
%     rho=para(7,i);
    tilt_deg=para(8,kk);
    baz_deg=para(9,kk);
  
    % add tilting and baz here
    phi_deg=xi_deg-baz_deg;
    phi=phi_deg/180*pi;
    theta=tilt_deg/180*pi;

    %calculate the new w;
    w=zeros(3,1);
    w(1)=sin(theta)*cos(phi);
    w(2)=sin(theta)*sin(phi);
    w(3)=cos(theta);

    %calculate the W and I
    W=kron(w,w.')-0.5*eye;
    Lambda0=zeros(6,6);

 for i=1:3
      for j=1:3
          for k=1:3
              for l=1:3
                 LA=I(i,j)*I(k,l);
                 LB=W(i,j)*I(k,l)+I(i,j)*W(k,l);
                 LC=8*W(i,j)*W(k,l)-I(i,j)*W(k,l);
                 LD=I(k,j)*I(i,l)+I(l,j)*I(k,i)-2*I(i,j)*I(k,l);
                 LB13=W(k,j)*I(i,l)+I(k,j)*W(i,l);
                 LB14=W(l,j)*I(k,i)+I(l,j)*W(k,i);
                 LE=2*(LB13+LB14-2*LB)+LD;
                 
                 L=A*LA+B*LB+C*LC+D*LD+E*LE;
                 
                 m=i*eq(i,j)+(1-eq(i,j))*(9-i-j);
                 n=k*eq(k,l)+(1-eq(k,l))*(9-k-l);
                 Lambda0(m,n)=L;
              end
          end
      end
  end

Lambda(:,:,kk)=Lambda0;

end