

mu=398600.4418; %standard gravitational parameter in km^3/s^2
Re=6371; %Earth Radius in km
alt=300; %orbital altitude in km

vel=sqrt(mu/(alt+Re))*1000; %magnitude of the orbital velocity in m/s

N=100;

for i=1:N

velocity=[vel; 0; 0;]; 

system('./ptetra>>/dev/null');
printf('%d%% complete...\n',(i/N*100))

fid=fopen('momout.txt' ,'r') ;

drag(i,1)=fscanf(fid,'Drag Force:\nx=   %f\n');
drag(i,2)=fscanf(fid,'y=   %f\n');
drag(i,3)=fscanf(fid,'z=   %f\n');

torq(i,1)=fscanf(fid,'Torques:\nx=   %f\n');
torq(i,2)=fscanf(fid,'y=   %f\n');
torq(i,3)=fscanf(fid,'z=   %f\n');
fflush(stdout);
end

fprintf('mean = %e, std = %e',mean(drag(:,1)),std(drag(:,1)));

hist(drag(:,1))

save('stat.mat','drag','torq')

