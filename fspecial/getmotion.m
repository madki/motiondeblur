% This function takes in name of the data set adn gives [dt, R, x]
function [dt, v, x, Rmat] = getmotion(name)
    fid = fopen(strcat(name,'acc.csv'));
    C = textscan(fid, '%f %f %f %f', 'Delimiter',',','EndOfLine','\n');
    fclose(fid);

    Xacc = cell2mat(C(:,1:4));

    fid = fopen(strcat(name,'rot.csv'));
    C = textscan(fid, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter',',','EndOfLine','\n');
    fclose(fid);

    Xrot = cell2mat(C(:,1:10));
    Rmat = Xrot(:,2:10);
    apt = Xacc(:,2:4) ;
    api = 0*apt;
    dt  = zeros(size(Xacc(:,1),1),1);
    for n = 1:size(Xacc(:,1),1)-1,
        dt(n)  = (Xacc(n+1,1)-Xacc(n,1))*1e-9;
    end
    x   = zeros(size(dt,1)+1,3);
    v   = zeros(size(dt,1)+1,3);

    for n = 1:size(dt,1),
        R = [Xrot(n,2:4)', Xrot(n,5:7)', Xrot(n,8:10)'];
        api(n,:) = (R*apt(n,:)')';
    end

    g = mean(api);

    for n = 1:size(dt,1),
        v(n+1,:) = v(n,:)+(api(n,:)-g)*dt(n);
        x(n+1,:) = x(n,:)+v(n,:)*dt(n)+0.5*(api(n,:)-g)*dt(n)*dt(n); 
    end
end