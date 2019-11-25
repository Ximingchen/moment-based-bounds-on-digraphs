%{
function A = RdSocialNetwork(filename)
    fileID = fopen(filename,'r');
    formatSpec = '%f %f';
    Data = fscanf(fileID,formatSpec,[2 inf]);
    n = max(max(Data(1:2,:)));
    A = zeros(n,n);
    for i = 1:size(Data,2)
        if Data(1,i)*Data(2,i) > 0
            A(Data(1,i),Data(2,i)) =  1;
        end
    end
end

%}
function A = RdSocialNetwork(filename)
    fileID = fopen(filename,'r');
    formatSpec = '%f %f';
    Data = fscanf(fileID,formatSpec,[2 inf]);
    
    indices = unique(Data);
    n = numel(indices);
    
    A = zeros(n,n);
    for i = 1:size(Data,2)
        if Data(1,i)*Data(2,i) > 0
            A(find(indices == Data(1,i)), find(indices == Data(2,i))) =  1;
        end
    end
end
%}