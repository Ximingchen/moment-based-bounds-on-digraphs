function A = RdNetwork(filename)
    fileID = fopen(filename,'r');
    formatSpec = '%f %f';
    Data = fscanf(fileID,formatSpec,[3 inf]);
    n = max(max(Data(1:2,:)));
    A = zeros(n,n);
    for i = 1:size(Data,2)
        A(Data(1,i),Data(2,i)) = Data(3,i);
    end
end