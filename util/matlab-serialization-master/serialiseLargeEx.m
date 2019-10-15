%Serialize data - only work for even number of frame

RcvData{1}=randi(2^16,[1024*100,128,46]);
S=whos('RcvData');
if S.bytes>2^32
    Nsize=size(RcvData{1});
    Nhalf=floor(Nsize(3)/2);
    RcvData{2,1}=RcvData{1}(:,:,Nhalf+1:end);
    RcvData{1,1}=RcvData{1}(:,:,1:Nhalf);
    RcvData{1,1}=serialize(RcvData(1,1));
    RcvData{2,1}=serialize(RcvData(2,1));
    RcvData=cat(3,RcvData{1},RcvData{2});
end

    savefast('test4','RcvData');
    
    %deserialize data
    temp=cell(size(RcvData,3),1);
    for i=1:size(RcvData,3)
        temp(i)=deserialize(RcvData(:,:,i));
    end
    
    RcvData=cat(3,temp{1},temp{2});
    
    %% Complex number serialization
    
  A=randi(2^16,[1024*100,128,25]);
  B=randi(2^16,[1024*100,128,25]);
  Data=complex(A,B);
  
  if ~isreal(Data)
      complexDataReal=real(Data);
      complexDataImag=imag(Data);
      clear Data;
      complexDataReal=serialize(complexDataReal);
      complexDataImag=serialize(complexDataImag);
  end
  
   savefast('test4','complexDataReal;','complexDataImag');
   
   complexDataReal=deserialize(complexDataReal);
   complexDataImag=deserialize(complexDataImag);
   
   data=complex(complexDataReal,complexDataImag);
   clear complexDataReal complexDataImag; 
   
  
  