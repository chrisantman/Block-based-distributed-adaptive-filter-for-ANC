
IN=buff_entrada;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Buffer de tamaño bloque de la señal de entrada   
    buffx=[buffx(B+1:tam_fft);IN];
    
    E=fft([zeros(B,1);IN_e(:,nodo)],tam_fft); 
    Et=repmat(E,1,P);  
    
    % Normalización mu
    Pot(:,:,nodo)=gama*Pot(:,:,nodo)+(1-gama)*(V(:,:,nodo).*conj(V(:,:,nodo))); 
    D=1./Pot(:,:,nodo); 
    phi=D.*Et.*conj(V(:,:,nodo)); 
    step=ifft(phi); 
    O=fft([step(1:B,:);zeros(B,P)],tam_fft);
    
    % Actualización de coeficientes
    W(:,:,nodo)=W(:,:,nodo)-2*mu.*O; 
%     aux=ifft(W(:,:,nodo));
%     W(:,:,nodo)=fft([aux(1:B,:);zeros(B,P)],tam_fft);
       
    % Filtrado estima
    x_f(:,2:end,nodo)=x_f(:,1:F-1,nodo);
    x_f(:,1,nodo)=buffx;  
    X_f=fft(x_f,tam_fft);
    V_aux=X_f(:,:,nodo).*CSEC_nodo(:,:,nodo);
    Vaux=sum(V_aux,2);   
    V(:,2:end,nodo)=V(:,1:P-1,nodo);
    V(:,1,nodo)=Vaux;
      
    % Filtrado adaptativo 
    x_p(:,2:end,nodo)=x_p(:,1:P-1,nodo);
    x_p(:,1,nodo)=buffx;  
    X_p=fft(x_p,tam_fft);
    Y=X_p(:,:,nodo).*W(:,:,nodo);                  
    y=ifft(Y); 
    yb=sum(y,2);     
    OUT_y(:,nodo)=real(yb(B+1:tam_fft,:));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buff_salida=OUT_Y;
    
    