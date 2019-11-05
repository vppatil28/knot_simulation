%Colormap "fiber" obtained from optomechanical fiber experiments

function [cm_data]=fiber(m)

%RGB values for strains from 0 to to 0.7 in steps of 0.05
e = [
255     0    65
255     0    40
255     0    17
255     0     5
255    23     0
255    96     0
255   135     0
238   174     0
189   207     0
111   231     0
 0   245     0
 0   251    84
 0   238   138
 0   211   178
 0   172   206];

Rave = e(:,1)'./256; Gave = e(:,2)'./256; Bave = e(:,3)'./256;
mymap=[Rave; Gave; Bave]'; %colormap specified by strain RGB values above

%interpolate colormap to get smooth gradient
cm_data=zeros(71,3); 
k=5;
for i=1:14
    v1 = mymap(i,:);
    v2 = mymap(i+1,:);
    cm_data(k*(i-1)+1,:) = v1;
    for j=1:k-1
        cm_data(k*(i-1)+j+1,:) = ((k-j)/k)*v1+(j/k)*v2;
    end
end
cm_data(end,:) = mymap(end,:);

end

    
    
    





